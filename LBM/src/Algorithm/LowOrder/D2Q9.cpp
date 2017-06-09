/*
 * D2Q9.cpp
 *
 *  Created on: 9 Jun 2015
 *      Author: thomas
 */

#include "D2Q9.h"
#include <iomanip>
D2Q9::D2Q9() {
	MultiBlock_=0;
	parallel=0;
	Writer=0;
	PtrParameters=0;
	f=0;
	nbvelo=9;
	Nb_VelocityCollide=nbvelo;
	nbnode=0;
	ftmp=0;
	tmp=0;
	tmpreturn=0;
	intTmpReturn=0;
	doubleTmpReturn=0;
	buf_recv=0;
	buf_send=0;
	size_buf=0;
	DiagConnect=0;
	Nd_variables_sync=0;
	buf_MacroSend=0;
	buf_MacroRecv=0;
	size_MacroBuf=0;
	Nd_MacroVariables_sync=0;
}
D2Q9::D2Q9(MultiBlock* MultiBlock__,ParallelManager* parallel__,WriterManager* Writer__, Parameters* Parameters_ ,InitLBM& ini) {


	f=new DistriFunct(MultiBlock__->Get_nnodes(),Parameters_->Get_NbVelocities());
	Set_Solver(MultiBlock__,parallel__,Writer__,Parameters_);
	ftmp=new double[nbnode];
	PtrFiStream=f;
	PtrFiCollide=PtrFiStream;
	InvTau=1.0/PtrParameters->Get_Tau();
	intTmpReturn=0;
	doubleTmpReturn=0;

	EiCollide=Ei;
	Nb_VelocityCollide=nbvelo;
	omegaCollide=omega;
	//if(PtrParameters->Get_WallType()==BounceBack)
		Dic->AddSync("Density",Rho);
	if(PtrParameters->Get_UserForceType()== ModelEnum::BodyForce ||PtrParameters->Get_HeleShawBodyForce()!= PorousMediaEnum::no)
	{
		F=new double*[2];
		Dic->AddVar(Vector,"BodyForce",true,false,false,F[0],F[1]);
		for(int i=0;i<nbnodes_total;i++)
			{F[0][i]=0;	F[1][i]=0;}
	}

	//Set Pointers On Functions for selecting the right model dynamically
		Set_PointersOnFunctions();
	// Initialise domain
		init(ini);
	//Initialise boundary conditions.
		InitD2Q9Bc(Dic, Parameters_,Ei);
	//Set_Convergence
		Set_Convergence();
	//Initialise communication between processors
		IniComVariables();
	//Set the variables names and the variable pointers for output in solution
		Solution2D::Set_output();
	//Set the variables names and the variable pointers for breakpoints in solution
		Solution2D::Set_breakpoint();

}
void D2Q9::Set_PointersOnFunctions(){
// Select the model for two-phase operator in the collision step
	Set_Collide();
// Select the macrocospic function for the colour fluid model depending of the output variables and models
	Set_Macro();
}
void D2Q9::Set_Collide(){

	PtrDicCollide=Dic;
	if(PtrParameters->Get_FluidType()==ModelEnum::Newtonian )
	{
		if(PtrParameters->Get_HeleShawBodyForce()!= PorousMediaEnum::no)
		{
			Select_Collide_2D(Std2DBody,PtrParameters->Get_cs2(),PtrParameters->Get_ReferenceDensity(),PtrParameters->Get_ModelOfFluid());
			Select_Collide_2D_V2(*PtrParameters,Std2DBody,PtrParameters->Get_cs2(),PtrParameters->Get_ReferenceDensity(),PtrParameters->Get_ModelOfFluid());
		}
		else if(PtrParameters->Get_UserForceType()== ModelEnum::LocalForce)
			{
			Select_Collide_2D(Std2DLocal,PtrParameters->Get_cs2(),PtrParameters->Get_ReferenceDensity(),PtrParameters->Get_ModelOfFluid());
			Select_Collide_2D_V2(*PtrParameters,Std2DLocal,PtrParameters->Get_cs2(),PtrParameters->Get_ReferenceDensity(),PtrParameters->Get_ModelOfFluid());
			}
		else if(PtrParameters->Get_UserForceType()== ModelEnum::BodyForce)
			{
			Select_Collide_2D(Std2DBody,PtrParameters->Get_cs2(),PtrParameters->Get_ReferenceDensity(),PtrParameters->Get_ModelOfFluid());
			Select_Collide_2D_V2(*PtrParameters,Std2DBody,PtrParameters->Get_cs2(),PtrParameters->Get_ReferenceDensity(),PtrParameters->Get_ModelOfFluid());
			}
		else
			{
			Select_Collide_2D(Std2D,PtrParameters->Get_cs2(),PtrParameters->Get_ReferenceDensity(),PtrParameters->Get_ModelOfFluid());
			Select_Collide_2D_V2(*PtrParameters,Std2D,PtrParameters->Get_cs2(),PtrParameters->Get_ReferenceDensity(),PtrParameters->Get_ModelOfFluid());
			}
	}
	else
	{
		if(PtrParameters->Get_HeleShawBodyForce()!= PorousMediaEnum::no)
		{
			Select_Collide_2D(Std2DNonCstTauBody,PtrParameters->Get_cs2(),PtrParameters->Get_ReferenceDensity(),PtrParameters->Get_ModelOfFluid());
			Select_Collide_2D_V2(*PtrParameters,Std2DNonCstTauBody,PtrParameters->Get_cs2(),PtrParameters->Get_ReferenceDensity(),PtrParameters->Get_ModelOfFluid());
		}
		if(PtrParameters->Get_UserForceType()== ModelEnum::LocalForce)
			{
			Select_Collide_2D(Std2DNonCstTauLocal,PtrParameters->Get_cs2(),PtrParameters->Get_ReferenceDensity(),PtrParameters->Get_ModelOfFluid());
			Select_Collide_2D_V2(*PtrParameters,Std2DNonCstTauLocal,PtrParameters->Get_cs2(),PtrParameters->Get_ReferenceDensity(),PtrParameters->Get_ModelOfFluid());
			}
		else if(PtrParameters->Get_UserForceType()== ModelEnum::BodyForce)
			{
			Select_Collide_2D(Std2DNonCstTauBody,PtrParameters->Get_cs2(),PtrParameters->Get_ReferenceDensity(),PtrParameters->Get_ModelOfFluid());
			Select_Collide_2D_V2(*PtrParameters,Std2DNonCstTauBody,PtrParameters->Get_cs2(),PtrParameters->Get_ReferenceDensity(),PtrParameters->Get_ModelOfFluid());
			}
		else
			{
			Select_Collide_2D(Std2DNonCstTau,PtrParameters->Get_cs2(),PtrParameters->Get_ReferenceDensity(),PtrParameters->Get_ModelOfFluid());
			Select_Collide_2D_V2(*PtrParameters,Std2DNonCstTau,PtrParameters->Get_cs2(),PtrParameters->Get_ReferenceDensity(),PtrParameters->Get_ModelOfFluid());
			}
	}
	if(PtrParameters->Get_UserForceType()== ModelEnum::BodyForce||PtrParameters->Get_HeleShawBodyForce()!= PorousMediaEnum::no)
		PtrCollision=&D2Q9::CollideD2Q9_WithBodyForce;
	else
		PtrCollision=&D2Q9::CollideD2Q9_NoBodyForce;

}
void D2Q9::Set_Macro(){
/*
		if(PtrParameters->Get_UserForceType()== ModelEnum::BodyForce||PtrParameters->Get_HeleShawBodyForce()!= PorousMediaEnum::no)
			PtrMacro=&D2Q9::UpdateMacroVariables_WithBodyForce;
		else*/
			PtrMacro=&D2Q9::UpdateMacroVariables_NoBodyForce;

}
D2Q9::~D2Q9() {
	delete [] size_buf;

	for (int i=0;i<Nd_variables_sync;i++)
	{
		delete [] buf_send[i][0];
		delete [] buf_send[i][1];
		delete [] buf_send[i][2];
		delete [] buf_send[i][3];
		delete [] buf_recv[i][0];
		delete [] buf_recv[i][1];
		delete [] buf_recv[i][2];
		delete [] buf_recv[i][3];
	}
	for (int i=0;i<Nd_variables_sync;i++)
	{
		delete [] buf_send[i];
		delete [] buf_recv[i];
	}
	if(PtrParameters->Get_UserForceType()== ModelEnum::BodyForce ||PtrParameters->Get_HeleShawBodyForce()!= PorousMediaEnum::no)
		delete [] F;
}
void D2Q9::init(InitLBM& ini){
	double* pos =new double[2];
	double* U_=new double[2];
	for (int j=0;j<NodeArrays->NodeInterior.size();j++)
	{
		pos[0]=NodeArrays->NodeInterior[j].get_x();
		pos[1]=NodeArrays->NodeInterior[j].get_y();
		ini.IniDomainSinglePhase(parallel->getRank(),NodeArrays->NodeInterior[j],0, NodeArrays->NodeInterior[j].Get_index(),pos,Rho[NodeArrays->NodeInterior[j].Get_index()],U_);
		U[0][NodeArrays->NodeInterior[j].Get_index()]=U_[0];
		U[1][NodeArrays->NodeInterior[j].Get_index()]=U_[1];
		if(Rho[j]==0)
			std::cerr<<" Density set to 0"<< pos[NodeArrays->NodeInterior[j].Get_index()]<<" "<<pos[1]<<std::endl;
	}
	for (int j=0;j<NodeArrays->NodeCorner.size();j++)
	{
		pos[0]=NodeArrays->NodeCorner[j].get_x();
		pos[1]=NodeArrays->NodeCorner[j].get_y();
		ini.IniDomainSinglePhase(parallel->getRank(),NodeArrays->NodeCorner[j],0, NodeArrays->NodeCorner[j].Get_index(),pos,Rho[NodeArrays->NodeCorner[j].Get_index()],U_);
		U[0][NodeArrays->NodeCorner[j].Get_index()]=U_[0];
		U[1][NodeArrays->NodeCorner[j].Get_index()]=U_[1];
	}
	for (int j=0;j<NodeArrays->NodeGlobalCorner.size();j++)
	{
		pos[0]=NodeArrays->NodeGlobalCorner[j].get_x();
		pos[1]=NodeArrays->NodeGlobalCorner[j].get_y();
		ini.IniDomainSinglePhase(parallel->getRank(),NodeArrays->NodeGlobalCorner[j],0, NodeArrays->NodeGlobalCorner[j].Get_index(),pos,Rho[NodeArrays->NodeGlobalCorner[j].Get_index()],U_);
		U[0][NodeArrays->NodeGlobalCorner[j].Get_index()]=U_[0];
		U[1][NodeArrays->NodeGlobalCorner[j].Get_index()]=U_[1];
		NodeArrays->NodeGlobalCorner[j].Set_UDef(U_[0],U_[1]);
		NodeArrays->NodeGlobalCorner[j].Set_RhoDef(Rho[NodeArrays->NodeGlobalCorner[j].Get_index()]);
	}
	for (int j=0;j<NodeArrays->NodeVelocity.size();j++)
	{
		pos[0]=NodeArrays->NodeVelocity[j].get_x();
		pos[1]=NodeArrays->NodeVelocity[j].get_y();
		ini.IniDomainSinglePhase(parallel->getRank(),NodeArrays->NodeVelocity[j],0, NodeArrays->NodeVelocity[j].Get_index(),pos,Rho[NodeArrays->NodeVelocity[j].Get_index()],U_);
		U[0][NodeArrays->NodeVelocity[j].Get_index()]=U_[0];
		U[1][NodeArrays->NodeVelocity[j].Get_index()]=U_[1];
	}

	for (int j=0;j<NodeArrays->NodePressure.size();j++)
	{
		pos[0]=NodeArrays->NodePressure[j].get_x();
		pos[1]=NodeArrays->NodePressure[j].get_y();
		ini.IniDomainSinglePhase(parallel->getRank(),NodeArrays->NodePressure[j],0, NodeArrays->NodePressure[j].Get_index(),pos,Rho[NodeArrays->NodePressure[j].Get_index()],U_);
		U[0][NodeArrays->NodePressure[j].Get_index()]=U_[0];
		U[1][NodeArrays->NodePressure[j].Get_index()]=U_[1];
	}
	for (int j=0;j<NodeArrays->NodeWall.size();j++)
	{
		pos[0]=NodeArrays->NodeWall[j].get_x();
		pos[1]=NodeArrays->NodeWall[j].get_y();
		ini.IniDomainSinglePhase(parallel->getRank(),NodeArrays->NodeWall[j],0, NodeArrays->NodeWall[j].Get_index(),pos,Rho[NodeArrays->NodeWall[j].Get_index()],U_);
		U[0][NodeArrays->NodeWall[j].Get_index()]=U_[0];
		U[1][NodeArrays->NodeWall[j].Get_index()]=U_[1];
	}
	for (int j=0;j<NodeArrays->NodeSpecialWall.size();j++)
	{
		pos[0]=NodeArrays->NodeSpecialWall[j].get_x();
		pos[1]=NodeArrays->NodeSpecialWall[j].get_y();
		ini.IniDomainSinglePhase(parallel->getRank(),NodeArrays->NodeSpecialWall[j],0, NodeArrays->NodeSpecialWall[j].Get_index(),pos,Rho[NodeArrays->NodeSpecialWall[j].Get_index()],U_);
		U[0][NodeArrays->NodeSpecialWall[j].Get_index()]=U_[0];
		U[1][NodeArrays->NodeSpecialWall[j].Get_index()]=U_[1];
	}
	for (int j=0;j<NodeArrays->NodeSymmetry.size();j++)
	{
		pos[0]=NodeArrays->NodeSymmetry[j].get_x();
		pos[1]=NodeArrays->NodeSymmetry[j].get_y();
		ini.IniDomainSinglePhase(parallel->getRank(),NodeArrays->NodeSymmetry[j],0, NodeArrays->NodeSymmetry[j].Get_index(),pos,Rho[NodeArrays->NodeSymmetry[j].Get_index()],U_);
		U[0][NodeArrays->NodeSymmetry[j].Get_index()]=U_[0];
		U[1][NodeArrays->NodeSymmetry[j].Get_index()]=U_[1];
	}
	for (int j=0;j<NodeArrays->NodePeriodic.size();j++)
	{
		pos[0]=NodeArrays->NodePeriodic[j].get_x();
		pos[1]=NodeArrays->NodePeriodic[j].get_y();
		ini.IniDomainSinglePhase(parallel->getRank(),NodeArrays->NodePeriodic[j],0, NodeArrays->NodePeriodic[j].Get_index(),pos,Rho[NodeArrays->NodePeriodic[j].Get_index()],U_);
		U[0][NodeArrays->NodePeriodic[j].Get_index()]=U_[0];
		U[1][NodeArrays->NodePeriodic[j].Get_index()]=U_[1];
	}
	for (int j=0;j<NodeArrays->NodeGhost.size();j++)
	{
		pos[0]=NodeArrays->NodeGhost[j].get_x();
		pos[1]=NodeArrays->NodeGhost[j].get_y();
		ini.IniDomainSinglePhase(parallel->getRank(),NodeArrays->NodeGhost[j],0, NodeArrays->NodeGhost[j].Get_index(),pos,Rho[NodeArrays->NodeGhost[j].Get_index()],U_);
		U[0][NodeArrays->NodeGhost[j].Get_index()]=U_[0];
		U[1][NodeArrays->NodeGhost[j].Get_index()]=U_[1];
	}
	for (int j=0;j<NodeArrays->NodeSolid.size();j++)
	{
		pos[0]=NodeArrays->NodeSolid[j].get_x();
		pos[1]=NodeArrays->NodeSolid[j].get_y();
		ini.IniDomainSinglePhase(parallel->getRank(),NodeArrays->NodeSolid[j],0, NodeArrays->NodeSolid[j].Get_index(),pos,Rho[NodeArrays->NodeSolid[j].Get_index()],U_);
		U[0][NodeArrays->NodeSolid[j].Get_index()]=U_[0];
		U[1][NodeArrays->NodeSolid[j].Get_index()]=U_[1];
	}

	if(PtrParameters->IsInitFromFile())
	{
		for(int i=0;i<PtrParameters->Get_NumberVariableToInit();i++)
		{
			Read_Variable(PtrParameters->Get_VariableNameToInit(i),PtrParameters->Get_FileNameToInit(i));
		}
	}
	for (int i=0;i<nbvelo;i++)
	{
		for (int j=0;j<nbnode;j++)
		{
			f->f[i][j]=CollideLowOrder::CollideEquillibrium(Rho[j], U[0][j], U[1][j], &Ei[i][0], omega[i]);
		}
	}

	Set_BcType();
	delete [] pos;
	delete [] U_;
	if(CalPressure)
		UpdatePressure();

}
void D2Q9::InitAllDomain(InitLBM& ini){

	InitDomainBc(ini);
	InitWall(ini);
	InitInterior(ini);

	double* pos =new double[2];
	double* U_=new double[2];
	int idx=0;
	for (int j=0;j<NodeArrays->NodeSolid.size();j++)
	{
		idx=NodeArrays->NodeSolid[j].Get_index();
		pos[0]=NodeArrays->NodeSolid[j].get_x();
		pos[1]=NodeArrays->NodeSolid[j].get_y();
		ini.IniDomainSinglePhase(parallel->getRank(),NodeArrays->NodeSolid[j],0, idx,pos,Rho[idx],U_);
		U[0][idx]=U_[0];
		U[1][idx]=U_[1];
	}
	delete [] pos;
	delete [] U_;
}
void
void D2Q9::InitDomainBc(InitLBM& ini){
	double* pos =new double[2];
	double* U_=new double[2];
	int idx=0;
	for (int j=0;j<NodeArrays->NodeGlobalCorner.size();j++)
	{
		idx=NodeArrays->NodeGlobalCorner[j].Get_index();
		pos[0]=NodeArrays->NodeGlobalCorner[j].get_x();
		pos[1]=NodeArrays->NodeGlobalCorner[j].get_y();
		ini.IniDomainSinglePhase(parallel->getRank(),NodeArrays->NodeGlobalCorner[j],0, idx,pos,Rho[idx],U_);
		U[0][idx]=U_[0];
		U[1][idx]=U_[1];
		NodeArrays->NodeGlobalCorner[j].Set_UDef(U_[0],U_[1]);
		NodeArrays->NodeGlobalCorner[j].Set_RhoDef(Rho[idx]);
	}
	for (int j=0;j<NodeArrays->NodeVelocity.size();j++)
	{
		idx=NodeArrays->NodeVelocity[j].Get_index();
		pos[0]=NodeArrays->NodeVelocity[j].get_x();
		pos[1]=NodeArrays->NodeVelocity[j].get_y();
		ini.IniDomainSinglePhase(parallel->getRank(),NodeArrays->NodeVelocity[j],0, idx,pos,Rho[idx],U_);
		U[0][idx]=U_[0];
		U[1][idx]=U_[1];
		NodeArrays->NodeVelocity[j].Set_UDef(U_[0],U_[1]);
	}

	for (int j=0;j<NodeArrays->NodePressure.size();j++)
	{
		idx=NodeArrays->NodePressure[j].Get_index();
		pos[0]=NodeArrays->NodePressure[j].get_x();
		pos[1]=NodeArrays->NodePressure[j].get_y();
		ini.IniDomainSinglePhase(parallel->getRank(),NodeArrays->NodePressure[j],0, idx,pos,Rho[idx],U_);
		U[0][idx]=U_[0];
		U[1][idx]=U_[1];
		NodeArrays->NodePressure[j].Set_RhoDef(Rho[idx]);
	}
	for (int j=0;j<NodeArrays->NodeSymmetry.size();j++)
	{
		idx=NodeArrays->NodeSymmetry[j].Get_index();
		pos[0]=NodeArrays->NodeSymmetry[j].get_x();
		pos[1]=NodeArrays->NodeSymmetry[j].get_y();
		ini.IniDomainSinglePhase(parallel->getRank(),NodeArrays->NodeSymmetry[j],0, idx,pos,Rho[idx],U_);
		U[0][idx]=U_[0];
		U[1][idx]=U_[1];
	}
	for (int j=0;j<NodeArrays->NodePeriodic.size();j++)
	{
		idx=NodeArrays->NodePeriodic[j].Get_index();
		pos[0]=NodeArrays->NodePeriodic[j].get_x();
		pos[1]=NodeArrays->NodePeriodic[j].get_y();
		ini.IniDomainSinglePhase(parallel->getRank(),NodeArrays->NodePeriodic[j],0, idx,pos,Rho[idx],U_);
		U[0][idx]=U_[0];
		U[1][idx]=U_[1];
	}
	delete [] pos;
	delete [] U_;
}
void D2Q9::InitWall(InitLBM& ini){
	double* pos =new double[2];
	double* U_=new double[2];
	int idx=0;
	for (int j=0;j<NodeArrays->NodeCorner.size();j++)
	{
		idx=NodeArrays->NodeCorner[j].Get_index();
		pos[0]=NodeArrays->NodeCorner[j].get_x();
		pos[1]=NodeArrays->NodeCorner[j].get_y();
		ini.IniDomainSinglePhase(parallel->getRank(),NodeArrays->NodeCorner[j],0, idx,pos,Rho[idx],U_);
		U[0][idx]=U_[0];
		U[1][idx]=U_[1];
		NodeArrays->NodeCorner[j].Set_UDef(U_[0],U_[1]);
		NodeArrays->NodeCorner[j].Set_RhoDef(Rho[idx]);
	}
	for (int j=0;j<NodeArrays->NodeWall.size();j++)
	{
		idx=NodeArrays->NodeWall[j].Get_index();
		pos[0]=NodeArrays->NodeWall[j].get_x();
		pos[1]=NodeArrays->NodeWall[j].get_y();
		ini.IniDomainSinglePhase(parallel->getRank(),NodeArrays->NodeWall[j],0, idx,pos,Rho[idx],U_);
		U[0][idx]=U_[0];
		U[1][idx]=U_[1];
	}
	for (int j=0;j<NodeArrays->NodeSpecialWall.size();j++)
	{
		idx=NodeArrays->NodeSpecialWall[j].Get_index();
		pos[0]=NodeArrays->NodeSpecialWall[j].get_x();
		pos[1]=NodeArrays->NodeSpecialWall[j].get_y();
		ini.IniDomainSinglePhase(parallel->getRank(),NodeArrays->NodeSpecialWall[j],0, idx,pos,Rho[idx],U_);
		U[0][idx]=U_[0];
		U[1][idx]=U_[1];
		NodeArrays->NodeSpecialWall[j].Set_UDef(U_[0],U_[1]);
		NodeArrays->NodeSpecialWall[j].Set_RhoDef(Rho[idx]);
	}

	delete [] pos;
	delete [] U_;
}
void D2Q9::InitInterior(InitLBM& ini){
	double* pos =new double[2];
	double* U_=new double[2];
	int idx=0;

	for (int j=0;j<NodeArrays->NodeInterior.size();j++)
	{
		idx=NodeArrays->NodeInterior[j].Get_index();
		pos[0]=NodeArrays->NodeInterior[j].get_x();
		pos[1]=NodeArrays->NodeInterior[j].get_y();
		ini.IniDomainSinglePhase(parallel->getRank(),NodeArrays->NodeInterior[j],0, idx,pos,Rho[idx],U_);
		U[0][idx]=U_[0];
		U[1][idx]=U_[1];
	}
	for (int j=0;j<NodeArrays->NodeGhost.size();j++)
	{
		idx=NodeArrays->NodeGhost[j].Get_index();
		pos[0]=NodeArrays->NodeGhost[j].get_x();
		pos[1]=NodeArrays->NodeGhost[j].get_y();
		ini.IniDomainSinglePhase(parallel->getRank(),NodeArrays->NodeGhost[j],0, idx,pos,Rho[idx],U_);
		U[0][idx]=U_[0];
		U[1][idx]=U_[1];
	}
	delete [] pos;
	delete [] U_;
}
//! Initialise the distributions.
void D2Q9::InitDistAllDomain(){
	InitDistDomainBc();
	InitDistWall();
	InitDistInterior();

// Init Solid
	double* pos =new double[2];
	double* U_=new double[2];
	int idx=0;

	for (int j=0;j<NodeArrays->NodeSolid.size();j++)
	{
// Set Index
		idx=NodeArrays->NodeSolid[j].Get_index();
// Get initialise values from the user

		for (int i=0;i<nbvelo;i++)
		{
			f->f[i][idx]=0;
		}
	}
	delete [] pos;
	delete [] U_;
}
void D2Q9::InitDistDomainBc(){
	double* pos =new double[2];
	double* U_=new double[2];
	int idx=0;
	for (int j=0;j<NodeArrays->NodeGlobalCorner.size();j++)
	{
// Set Index
		idx=NodeArrays->NodeGlobalCorner[j].Get_index();
		for (int i=0;i<nbvelo;i++)
		{
			f->f[i][idx]=CollideLowOrder::CollideEquillibrium(Rho[idx], U[0][idx], U[1][idx], &Ei[i][0], omega[i]);
		}
	}
	for (int j=0;j<NodeArrays->NodeSpecialWall.size();j++)
	{
// Set Index
		idx=NodeArrays->NodeSpecialWall[j].Get_index();
		for (int i=0;i<nbvelo;i++)
		{
			f->f[i][idx]=CollideLowOrder::CollideEquillibrium(Rho[idx], U[0][idx], U[1][idx], &Ei[i][0], omega[i]);
		}
	}
	for (int j=0;j<NodeArrays->NodeVelocity.size();j++)
	{
// Set Index
		idx=NodeArrays->NodeVelocity[j].Get_index();
		for (int i=0;i<nbvelo;i++)
		{
			f->f[i][idx]=CollideLowOrder::CollideEquillibrium(Rho[idx], U[0][idx], U[1][idx], &Ei[i][0], omega[i]);
		}
	}

	for (int j=0;j<NodeArrays->NodePressure.size();j++)
	{
// Set Index
		idx=NodeArrays->NodePressure[j].Get_index();
		for (int i=0;i<nbvelo;i++)
		{
			f->f[i][idx]=CollideLowOrder::CollideEquillibrium(Rho[idx], U[0][idx], U[1][idx], &Ei[i][0], omega[i]);
		}
	}
	for (int j=0;j<NodeArrays->NodeSymmetry.size();j++)
	{
// Set Index
		idx=NodeArrays->NodeSymmetry[j].Get_index();
		for (int i=0;i<nbvelo;i++)
		{
			f->f[i][idx]=CollideLowOrder::CollideEquillibrium(Rho[idx], U[0][idx], U[1][idx], &Ei[i][0], omega[i]);
		}
	}
	for (int j=0;j<NodeArrays->NodePeriodic.size();j++)
	{
// Set Index
		idx=NodeArrays->NodePeriodic[j].Get_index();
		for (int i=0;i<nbvelo;i++)
		{
			f->f[i][idx]=CollideLowOrder::CollideEquillibrium(Rho[idx], U[0][idx], U[1][idx], &Ei[i][0], omega[i]);
		}
	}
	delete [] pos;
	delete [] U_;
}
void D2Q9::InitDistWall(){
	double* pos =new double[2];
	double* U_=new double[2];
	int idx=0;

	for (int j=0;j<NodeArrays->CornerConcave.size();j++)
	{
// Set Index
		idx=NodeArrays->NodeCorner[NodeArrays->CornerConcave[j]].Get_index();
		for (int i=0;i<nbvelo;i++)
		{
			f->f[i][idx]=CollideLowOrder::CollideEquillibrium(Rho[idx], U[0][idx], U[1][idx], &Ei[i][0], omega[i]);
		}
	}

	for (int j=0;j<NodeArrays->NodeWall.size();j++)
	{
// Set Index
		idx=NodeArrays->NodeWall[j].Get_index();
		for (int i=0;i<nbvelo;i++)
		{
			f->f[i][idx]=CollideLowOrder::CollideEquillibrium(Rho[idx], U[0][idx], U[1][idx], &Ei[i][0], omega[i]);
		}
	}

	for (int j=0;j<NodeArrays->CornerConvex.size();j++)
	{
// Set Index
		idx=NodeArrays->NodeCorner[NodeArrays->CornerConvex[j]].Get_index();
		for (int i=0;i<nbvelo;i++)
		{
			f->f[i][idx]=CollideLowOrder::CollideEquillibrium(Rho[idx], U[0][idx], U[1][idx], &Ei[i][0], omega[i]);
		}
	}
	delete [] pos;
	delete [] U_;
}
void D2Q9::InitDistInterior(){
	double* pos =new double[2];
	double* U_=new double[2];
	int idx=0;

	for (int j=0;j<NodeArrays->NodeInterior.size();j++)
	{
// Set Index
		idx=NodeArrays->NodeInterior[j].Get_index();
		for (int i=0;i<nbvelo;i++)
		{
			f->f[i][idx]=CollideLowOrder::CollideEquillibrium(Rho[idx], U[0][idx], U[1][idx], &Ei[i][0], omega[i]);
		}
	}

	for (int j=0;j<NodeArrays->NodeGhost.size();j++)
	{
// Set Index
		idx=NodeArrays->NodeGhost[j].Get_index();
		for (int i=0;i<nbvelo;i++)
		{
			f->f[i][idx]=CollideLowOrder::CollideEquillibrium(Rho[idx], U[0][idx], U[1][idx], &Ei[i][0], omega[i]);
		}
	}
	delete [] pos;
	delete [] U_;
}
void D2Q9::UpdateAllDomain(Parameters* UpdatedParam,InitLBM& ini){
	//init field
	UpdateDomainBc(UpdatedParam,ini);
	UpdateWall(UpdatedParam,ini);
	UpdateInterior(UpdatedParam,ini);
	//init distri
	InitDistDomainBc();
	InitDistWall();
	InitDistInterior();
}
void D2Q9::UpdateDomainBc(Parameters* UpdatedParam,InitLBM& ini){
	//init field
	InitDomainBc(ini);
	//init distri
	InitDistDomainBc();
}
void D2Q9::UpdateWall(Parameters* UpdatedParam,InitLBM& ini){
	//init field
	InitWall(ini);
	//init distri
	InitDistWall();
}
void D2Q9::UpdateInterior(Parameters* UpdatedParam,InitLBM& ini){
	//init field
	InitInterior(ini);
	//init distri
	InitDistInterior();
}
void D2Q9::run(){

	int NbStep=PtrParameters->Get_NbStep();
	int OutPutNStep=PtrParameters->Get_OutPutNSteps();
	int listing=PtrParameters->Get_listing();
	double time_inirun=parallel->getTime();
	double time_run=0;
	int it=0;
	double max_error=PtrParameters->Get_ErrorMax()*listing;

	if(parallel->isMainProcessor())
		std::cout<<"Error max (listing*"<<PtrParameters->Get_ErrorMax()<<")= "<<max_error<<std::endl<< "Convergence will be checked every: "<<listing<<" iterations"<<std::endl;

	if(parallel->getSize()>1)
		SyncToGhost();
	UpdateMacroVariables();
	if(CalPressure)
			UpdatePressure();
	if(parallel->getSize()>1)
		SyncMacroVarToGhost();
	Convergence::Calcul_Error(it);
	Writer->Write_Output(it);
	it++;
 /*	char buffer[50]; // make sure it's big enough
 	std::ofstream myFlux;
 	snprintf(buffer, sizeof(buffer), "after_streaming_%d.txt", parallel->getRank());
 	myFlux.open(buffer);myFlux.close();
 	snprintf(buffer, sizeof(buffer), "before_streaming_%d.txt", parallel->getRank());
 	myFlux.open(buffer);myFlux.close();*/
	if(parallel->getSize()>1)
	{
		while(it<NbStep+1)
		{
			CollideD2Q9();
			SyncToGhost();
			StreamD2Q9();
			ApplyBc();
			UpdateMacroVariables();
			if(CalGradP)
				UpdatePressure();
			SyncMacroVarToGhost();

			if(it%OutPutNStep==0)
			{
				if(CalPressure&&!CalGradP)
						UpdatePressure();
				Writer->Write_Output(it);
			}
			if(it%listing==0  )
			{
				if(CalPressure&&!CalGradP)
						UpdatePressure();
				Convergence::Calcul_Error(it);
				if(parallel->isMainProcessor())
				{
					time_run=parallel->getTime()-time_inirun;

					std::cout<<"Iteration number: "<<it<< " Running time: ";
					if(time_run>60.0)
						std::cout<<trunc(time_run/60)<<"Minutes ";
					else //less than 1 minute
					std::cout<<trunc(time_run)<<"Seconds ";
					std::cout<< " Time per iteration: "<<time_run/it<<"s"<<std::endl;
					std::cout<<"Error is: "<<Get_Error()<<std::endl;
				}
				if(Get_Error()<max_error)
				{
					if(it<NbStep)
					{
						if(CalPressure&&!CalGradP)
							UpdatePressure();
						Writer->Write_Output(it);
					}
					it=NbStep;
				}
			}

			it++;
		}
	}
	else
	{
		while(it<NbStep+1)
		{
			CollideD2Q9();
			StreamD2Q9();
			ApplyBc();
			UpdateMacroVariables();
			if(CalGradP)
					UpdatePressure();
			if(it%OutPutNStep==0)
			{
				if(CalPressure&&!CalGradP)
					UpdatePressure();
				Writer->Write_Output(it);
			}
			if(it%listing==0  )
			{
				if(CalPressure&&!CalGradP)
						UpdatePressure();
				Convergence::Calcul_Error(it);
				time_run=parallel->getTime()-time_inirun;
				std::cout<<"Iteration number: "<<it<< " Running time: ";
				if(time_run>60.0)
					std::cout<<trunc(time_run/60)<<"Minutes ";
				else //less than 1 minute
				std::cout<<trunc(time_run)<<"Seconds ";
				std::cout<< " Time per iteration: "<<time_run/it<<"s"<<std::endl;
				std::cout<<"Error is: "<<Get_Error()<<std::endl;

				if(Get_Error()<max_error)
				{
					if(CalPressure&&!CalGradP)
							UpdatePressure();
					Writer->Write_Output(it);
					it=NbStep;
				}
			}
			it++;
		}
	}
	if((Get_Error()>=max_error)&&(it%OutPutNStep!=0))
	{
		if(CalPressure&&!CalGradP)
			UpdatePressure();
		it--;
		Writer->Write_Output(it);
	}
	Write_Breakpoint(PtrParameters);
	//Writer->Write_breakpoint(*PtrParameters);
}
void D2Q9::run(Parameters* UpdatedParam){

	PtrParameters=UpdatedParam;
	IniTau(PtrParameters);
	InvTau=Get_InvTau();
	D2Q9::run();
}

void D2Q9::StreamD2Q9() {

	for (unsigned int i=1;i<(unsigned int)nbvelo;i++) //No need to stream direction 0
	{
		for (int j=0;j<NodeArrays->NodeInterior.size();j++)
		{
			if (NodeArrays->NodeInterior[j].Get_connect()[i]!=NodeArrays->NodeInterior[j].Get_index())
			{
				ftmp[NodeArrays->NodeInterior[j].Get_connect()[i]]=f->f[i][NodeArrays->NodeInterior[j].Get_index()];
			}
		}

		for (int j=0;j<NodeArrays->NodeCorner.size();j++)
		{
			if (NodeArrays->NodeCorner[j].stream()[i])
			{
				ftmp[NodeArrays->NodeCorner[j].Get_connect()[i]]=f->f[i][NodeArrays->NodeCorner[j].Get_index()];
			}		}
		for (int j=0;j<NodeArrays->NodeGlobalCorner.size();j++)
		{
			if (NodeArrays->NodeGlobalCorner[j].stream()[i])
			{
				ftmp[NodeArrays->NodeGlobalCorner[j].Get_connect()[i]]=f->f[i][NodeArrays->NodeGlobalCorner[j].Get_index()];
			}

		}
		for (int j=0;j<NodeArrays->NodeVelocity.size();j++)
		{
			if (NodeArrays->NodeVelocity[j].stream()[i])
					{
						ftmp[NodeArrays->NodeVelocity[j].Get_connect()[i]]=f->f[i][NodeArrays->NodeVelocity[j].Get_index()];
					}
		}

		for (int j=0;j<NodeArrays->NodePressure.size();j++)
		{
			if (NodeArrays->NodePressure[j].stream()[i])
					{
						ftmp[NodeArrays->NodePressure[j].Get_connect()[i]]=f->f[i][NodeArrays->NodePressure[j].Get_index()];
					}
		}
		for (int j=0;j<NodeArrays->NodeWall.size();j++)
		{
			if (NodeArrays->NodeWall[j].stream()[i])
			{
				ftmp[NodeArrays->NodeWall[j].Get_connect()[i]]=f->f[i][NodeArrays->NodeWall[j].Get_index()];
			}
		}
		for (int j=0;j<NodeArrays->NodeSpecialWall.size();j++)
		{
			if (NodeArrays->NodeSpecialWall[j].stream()[i])
			{
				ftmp[NodeArrays->NodeSpecialWall[j].Get_connect()[i]]=f->f[i][NodeArrays->NodeSpecialWall[j].Get_index()];
			}
		}
		for (int j=0;j<NodeArrays->NodeSymmetry.size();j++)
		{
			if (NodeArrays->NodeSymmetry[j].stream()[i])
			{
				ftmp[NodeArrays->NodeSymmetry[j].Get_connect()[i]]=f->f[i][NodeArrays->NodeSymmetry[j].Get_index()];
			}
		}
		for (int j=0;j<NodeArrays->NodePeriodic.size();j++)
		{
			if (NodeArrays->NodePeriodic[j].stream()[i])
			{
				ftmp[NodeArrays->NodePeriodic[j].Get_connect()[i]]=f->f[i][NodeArrays->NodePeriodic[j].Get_index()];
			}
		}
		for (int j=0;j<NodeArrays->NodeGhost.size();j++)
		{
			if (NodeArrays->NodeGhost[j].stream()[i])
			{
				ftmp[NodeArrays->NodeGhost[j].Get_connect()[i]]=f->f[i][NodeArrays->NodeGhost[j].Get_index()];
			}
		}
		D2Q9::TmptoDistri(i);
	}

}
void D2Q9::CollideD2Q9(){
	(this->*PtrCollision)();
}
void D2Q9::CollideD2Q9_NoBodyForce(){
	double wtmp=0;
	double Fx, Fy, InvTau_tmp;
	double *fi_tmp; fi_tmp=new double [nbvelo];
	double *localforce;localforce=new double [nbvelo];
	 		for (int j=0;j<NodeArrays->NodeInterior.size();j++)
		{
			Fx=0;Fy=0;
			for (int i=0;i<9;i++)
			{
				fi_tmp[i]=f->f[i][NodeArrays->NodeInterior[j].Get_index()];localforce[i]=0;
			}
			Collide_2D_V2(fi_tmp,Rho[NodeArrays->NodeInterior[j].Get_index()], U[0][NodeArrays->NodeInterior[j].Get_index()], U[1][NodeArrays->NodeInterior[j].Get_index()],localforce, Fx, Fy, Get_InvTau(),Get_Mu());
			for (int i=0;i<9;i++)
			{
				f->f[i][NodeArrays->NodeInterior[j].Get_index()]=fi_tmp[i];
			}
		}
		for (int j=0;j<NodeArrays->NodeCorner.size();j++)
		{
			Fx=0;Fy=0;
			for (int i=0;i<9;i++)
			{
				fi_tmp[i]=f->f[i][NodeArrays->NodeCorner[j].Get_index()];localforce[i]=0;
			}
			Collide_2D_V2(fi_tmp,Rho[NodeArrays->NodeCorner[j].Get_index()], U[0][NodeArrays->NodeCorner[j].Get_index()], U[1][NodeArrays->NodeCorner[j].Get_index()],localforce, Fx, Fy, Get_InvTau(),Get_Mu());
			for (int i=0;i<9;i++)
			{
				f->f[i][NodeArrays->NodeCorner[j].Get_index()]=fi_tmp[i];
			}
		}
		for (int j=0;j<NodeArrays->NodeGlobalCorner.size();j++)
		{
			Fx=0;Fy=0;
			for (int i=0;i<9;i++)
			{
				fi_tmp[i]=f->f[i][NodeArrays->NodeGlobalCorner[j].Get_index()];localforce[i]=0;
			}
			Collide_2D_V2(fi_tmp,Rho[NodeArrays->NodeGlobalCorner[j].Get_index()], U[0][NodeArrays->NodeGlobalCorner[j].Get_index()], U[1][NodeArrays->NodeGlobalCorner[j].Get_index()],localforce, Fx, Fy, Get_InvTau(),Get_Mu());
			for (int i=0;i<9;i++)
			{
				f->f[i][NodeArrays->NodeGlobalCorner[j].Get_index()]=fi_tmp[i];
			}
		}
		for (int j=0;j<NodeArrays->NodeVelocity.size();j++)
		{
			Fx=0;Fy=0;
			for (int i=0;i<9;i++)
			{
				fi_tmp[i]=f->f[i][NodeArrays->NodeVelocity[j].Get_index()];localforce[i]=0;
			}
			Collide_2D_V2(fi_tmp,Rho[NodeArrays->NodeVelocity[j].Get_index()], U[0][NodeArrays->NodeVelocity[j].Get_index()], U[1][NodeArrays->NodeVelocity[j].Get_index()],localforce, Fx, Fy, Get_InvTau(),Get_Mu());
			for (int i=0;i<9;i++)
			{
				f->f[i][NodeArrays->NodeVelocity[j].Get_index()]=fi_tmp[i];
			}
		}
		for (int j=0;j<NodeArrays->NodePressure.size();j++)
		{
			Fx=0;Fy=0;
			for (int i=0;i<9;i++)
			{
				fi_tmp[i]=f->f[i][NodeArrays->NodePressure[j].Get_index()];localforce[i]=0;
			}
			Collide_2D_V2(fi_tmp,Rho[NodeArrays->NodePressure[j].Get_index()], U[0][NodeArrays->NodePressure[j].Get_index()], U[1][NodeArrays->NodePressure[j].Get_index()],localforce, Fx, Fy, Get_InvTau(),Get_Mu());
			for (int i=0;i<9;i++)
			{
				f->f[i][NodeArrays->NodePressure[j].Get_index()]=fi_tmp[i];
			}
		}
		for (int j=0;j<NodeArrays->NodeWall.size();j++)
		{
			Fx=0;Fy=0;
			for (int i=0;i<9;i++)
			{
				fi_tmp[i]=f->f[i][NodeArrays->NodeWall[j].Get_index()];localforce[i]=0;
			}
			Collide_2D_V2(fi_tmp,Rho[NodeArrays->NodeWall[j].Get_index()], U[0][NodeArrays->NodeWall[j].Get_index()], U[1][NodeArrays->NodeWall[j].Get_index()],localforce, Fx, Fy, Get_InvTau(),Get_Mu());
			for (int i=0;i<9;i++)
			{
				f->f[i][NodeArrays->NodeWall[j].Get_index()]=fi_tmp[i];
			}
		}
		for (int j=0;j<NodeArrays->NodeSpecialWall.size();j++)
		{
			Fx=0;Fy=0;
			for (int i=0;i<9;i++)
			{
				fi_tmp[i]=f->f[i][NodeArrays->NodeSpecialWall[j].Get_index()];localforce[i]=0;
			}
			Collide_2D_V2(fi_tmp,Rho[NodeArrays->NodeSpecialWall[j].Get_index()], U[0][NodeArrays->NodeSpecialWall[j].Get_index()], U[1][NodeArrays->NodeSpecialWall[j].Get_index()],localforce, Fx, Fy, Get_InvTau(),Get_Mu());
			for (int i=0;i<9;i++)
			{
				f->f[i][NodeArrays->NodeSpecialWall[j].Get_index()]=fi_tmp[i];
			}
		}
		for (int j=0;j<NodeArrays->NodeSymmetry.size();j++)
		{
			Fx=0;Fy=0;
			for (int i=0;i<9;i++)
			{
				fi_tmp[i]=f->f[i][NodeArrays->NodeSymmetry[j].Get_index()];localforce[i]=0;
			}
			Collide_2D_V2(fi_tmp,Rho[NodeArrays->NodeSymmetry[j].Get_index()], U[0][NodeArrays->NodeSymmetry[j].Get_index()], U[1][NodeArrays->NodeSymmetry[j].Get_index()],localforce, Fx, Fy, Get_InvTau(),Get_Mu());
			for (int i=0;i<9;i++)
			{
				f->f[i][NodeArrays->NodeSymmetry[j].Get_index()]=fi_tmp[i];
			}
		}
		for (int j=0;j<NodeArrays->NodePeriodic.size();j++)
		{
			Fx=0;Fy=0;
			for (int i=0;i<9;i++)
			{
				fi_tmp[i]=f->f[i][NodeArrays->NodePeriodic[j].Get_index()];localforce[i]=0;
			}
			Collide_2D_V2(fi_tmp,Rho[NodeArrays->NodePeriodic[j].Get_index()], U[0][NodeArrays->NodePeriodic[j].Get_index()], U[1][NodeArrays->NodePeriodic[j].Get_index()],localforce, Fx, Fy, Get_InvTau(),Get_Mu());
			for (int i=0;i<9;i++)
			{
				f->f[i][NodeArrays->NodePeriodic[j].Get_index()]=fi_tmp[i];
			}
		}

}
void D2Q9::CollideD2Q9_WithBodyForce(){
	double wtmp=0;
	//double InvTau_tmp;
	double *fi_tmp; fi_tmp=new double [nbvelo];
	double *localforce;localforce=new double [nbvelo];

	 		for (int j=0;j<NodeArrays->NodeInterior.size();j++)
		{
	 			F[0][NodeArrays->NodeInterior[j].Get_index()]=0;
	 			F[1][NodeArrays->NodeInterior[j].Get_index()]=0;
			for (int i=0;i<9;i++)
			{
				fi_tmp[i]=f->f[i][NodeArrays->NodeInterior[j].Get_index()];localforce[i]=0;
			}
			Collide_2D_V2(fi_tmp,Rho[NodeArrays->NodeInterior[j].Get_index()], U[0][NodeArrays->NodeInterior[j].Get_index()], U[1][NodeArrays->NodeInterior[j].Get_index()],localforce, F[0][NodeArrays->NodeInterior[j].Get_index()], F[1][NodeArrays->NodeInterior[j].Get_index()], Get_InvTau(),Get_Mu());
			for (int i=0;i<9;i++)
			{
				f->f[i][NodeArrays->NodeInterior[j].Get_index()]=fi_tmp[i];
			}
		}
		for (int j=0;j<NodeArrays->NodeCorner.size();j++)
		{
 			F[0][NodeArrays->NodeCorner[j].Get_index()]=0;
 			F[1][NodeArrays->NodeCorner[j].Get_index()]=0;
			for (int i=0;i<9;i++)
			{
				fi_tmp[i]=f->f[i][NodeArrays->NodeCorner[j].Get_index()];localforce[i]=0;
			}
			Collide_2D_V2(fi_tmp,Rho[NodeArrays->NodeCorner[j].Get_index()], U[0][NodeArrays->NodeCorner[j].Get_index()], U[1][NodeArrays->NodeCorner[j].Get_index()],localforce, F[0][NodeArrays->NodeCorner[j].Get_index()], F[1][NodeArrays->NodeCorner[j].Get_index()], Get_InvTau(),Get_Mu());
			for (int i=0;i<9;i++)
			{
				f->f[i][NodeArrays->NodeCorner[j].Get_index()]=fi_tmp[i];
			}
		}
		for (int j=0;j<NodeArrays->NodeGlobalCorner.size();j++)
		{
 			F[0][NodeArrays->NodeGlobalCorner[j].Get_index()]=0;
 			F[1][NodeArrays->NodeGlobalCorner[j].Get_index()]=0;
			for (int i=0;i<9;i++)
			{
				fi_tmp[i]=f->f[i][NodeArrays->NodeGlobalCorner[j].Get_index()];localforce[i]=0;
			}
			Collide_2D_V2(fi_tmp,Rho[NodeArrays->NodeGlobalCorner[j].Get_index()], U[0][NodeArrays->NodeGlobalCorner[j].Get_index()], U[1][NodeArrays->NodeGlobalCorner[j].Get_index()],localforce, F[0][NodeArrays->NodeGlobalCorner[j].Get_index()], F[1][NodeArrays->NodeGlobalCorner[j].Get_index()], Get_InvTau(),Get_Mu());
			for (int i=0;i<9;i++)
			{
				f->f[i][NodeArrays->NodeGlobalCorner[j].Get_index()]=fi_tmp[i];
			}
		}
		for (int j=0;j<NodeArrays->NodeVelocity.size();j++)
		{
 			F[0][NodeArrays->NodeVelocity[j].Get_index()]=0;
 			F[1][NodeArrays->NodeVelocity[j].Get_index()]=0;
			for (int i=0;i<9;i++)
			{
				fi_tmp[i]=f->f[i][NodeArrays->NodeVelocity[j].Get_index()];localforce[i]=0;
			}
			Collide_2D_V2(fi_tmp,Rho[NodeArrays->NodeVelocity[j].Get_index()], U[0][NodeArrays->NodeVelocity[j].Get_index()], U[1][NodeArrays->NodeVelocity[j].Get_index()],localforce, F[0][NodeArrays->NodeVelocity[j].Get_index()], F[1][NodeArrays->NodeVelocity[j].Get_index()], Get_InvTau(),Get_Mu());
			for (int i=0;i<9;i++)
			{
				f->f[i][NodeArrays->NodeVelocity[j].Get_index()]=fi_tmp[i];
			}
		}
		for (int j=0;j<NodeArrays->NodePressure.size();j++)
		{
 			F[0][NodeArrays->NodePressure[j].Get_index()]=0;
 			F[1][NodeArrays->NodePressure[j].Get_index()]=0;
			for (int i=0;i<9;i++)
			{
				fi_tmp[i]=f->f[i][NodeArrays->NodePressure[j].Get_index()];localforce[i]=0;
			}
			Collide_2D_V2(fi_tmp,Rho[NodeArrays->NodePressure[j].Get_index()], U[0][NodeArrays->NodePressure[j].Get_index()], U[1][NodeArrays->NodePressure[j].Get_index()],localforce, F[0][NodeArrays->NodePressure[j].Get_index()], F[1][NodeArrays->NodePressure[j].Get_index()], Get_InvTau(),Get_Mu());
			for (int i=0;i<9;i++)
			{
				f->f[i][NodeArrays->NodePressure[j].Get_index()]=fi_tmp[i];
			}
		}
		for (int j=0;j<NodeArrays->NodeWall.size();j++)
		{
 			F[0][NodeArrays->NodeWall[j].Get_index()]=0;
 			F[1][NodeArrays->NodeWall[j].Get_index()]=0;
			for (int i=0;i<9;i++)
			{
				fi_tmp[i]=f->f[i][NodeArrays->NodeWall[j].Get_index()];localforce[i]=0;
			}
			Collide_2D_V2(fi_tmp,Rho[NodeArrays->NodeWall[j].Get_index()], U[0][NodeArrays->NodeWall[j].Get_index()], U[1][NodeArrays->NodeWall[j].Get_index()],localforce, F[0][NodeArrays->NodeWall[j].Get_index()], F[1][NodeArrays->NodeWall[j].Get_index()], Get_InvTau(),Get_Mu());
			for (int i=0;i<9;i++)
			{
				f->f[i][NodeArrays->NodeWall[j].Get_index()]=fi_tmp[i];
			}
		}
		for (int j=0;j<NodeArrays->NodeSpecialWall.size();j++)
		{
 			F[0][NodeArrays->NodeSpecialWall[j].Get_index()]=0;
 			F[1][NodeArrays->NodeSpecialWall[j].Get_index()]=0;
			for (int i=0;i<9;i++)
			{
				fi_tmp[i]=f->f[i][NodeArrays->NodeSpecialWall[j].Get_index()];localforce[i]=0;
			}
			Collide_2D_V2(fi_tmp,Rho[NodeArrays->NodeSpecialWall[j].Get_index()], U[0][NodeArrays->NodeSpecialWall[j].Get_index()], U[1][NodeArrays->NodeSpecialWall[j].Get_index()],localforce, F[0][NodeArrays->NodeSpecialWall[j].Get_index()], F[1][NodeArrays->NodeSpecialWall[j].Get_index()], Get_InvTau(),Get_Mu());
			for (int i=0;i<9;i++)
			{
				f->f[i][NodeArrays->NodeSpecialWall[j].Get_index()]=fi_tmp[i];
			}
		}
		for (int j=0;j<NodeArrays->NodeSymmetry.size();j++)
		{
 			F[0][NodeArrays->NodeSymmetry[j].Get_index()]=0;
 			F[1][NodeArrays->NodeSymmetry[j].Get_index()]=0;
			for (int i=0;i<9;i++)
			{
				fi_tmp[i]=f->f[i][NodeArrays->NodeSymmetry[j].Get_index()];localforce[i]=0;
			}
			Collide_2D_V2(fi_tmp,Rho[NodeArrays->NodeSymmetry[j].Get_index()], U[0][NodeArrays->NodeSymmetry[j].Get_index()], U[1][NodeArrays->NodeSymmetry[j].Get_index()],localforce, F[0][NodeArrays->NodeSymmetry[j].Get_index()], F[1][NodeArrays->NodeSymmetry[j].Get_index()], Get_InvTau(),Get_Mu());
			for (int i=0;i<9;i++)
			{
				f->f[i][NodeArrays->NodeSymmetry[j].Get_index()]=fi_tmp[i];
			}
		}
		for (int j=0;j<NodeArrays->NodePeriodic.size();j++)
		{
 			F[0][NodeArrays->NodePeriodic[j].Get_index()]=0;
 			F[1][NodeArrays->NodePeriodic[j].Get_index()]=0;
			for (int i=0;i<9;i++)
			{
				fi_tmp[i]=f->f[i][NodeArrays->NodePeriodic[j].Get_index()];localforce[i]=0;
			}
			Collide_2D_V2(fi_tmp,Rho[NodeArrays->NodePeriodic[j].Get_index()], U[0][NodeArrays->NodePeriodic[j].Get_index()], U[1][NodeArrays->NodePeriodic[j].Get_index()],localforce, F[0][NodeArrays->NodePeriodic[j].Get_index()], F[1][NodeArrays->NodePeriodic[j].Get_index()], Get_InvTau(),Get_Mu());
			for (int i=0;i<9;i++)
			{
				f->f[i][NodeArrays->NodePeriodic[j].Get_index()]=fi_tmp[i];
			}
		}
		delete [] fi_tmp;delete [] localforce;
}
void D2Q9::UpdateMacroVariables(){
	(this->*PtrMacro)();
}
void D2Q9::UpdateMacroVariables_NoBodyForce(){
		for (int j=0;j<NodeArrays->NodeInterior.size();j++)
		{
			MacroVariables(NodeArrays->NodeInterior[j].Get_index());
		}
		for (int j=0;j<NodeArrays->NodeCorner.size();j++)
		{
			MacroVariables(NodeArrays->NodeCorner[j].Get_index());
		}
		for (int j=0;j<NodeArrays->NodeGlobalCorner.size();j++)
		{
			MacroVariables(NodeArrays->NodeGlobalCorner[j].Get_index());
		}
		for (int j=0;j<NodeArrays->NodeVelocity.size();j++)
		{
			MacroVariables(NodeArrays->NodeVelocity[j].Get_index());
		}
		for (int j=0;j<NodeArrays->NodePressure.size();j++)
		{
			MacroVariables(NodeArrays->NodePressure[j].Get_index());
		}
		for (int j=0;j<NodeArrays->NodeWall.size();j++)
		{
			MacroVariables(NodeArrays->NodeWall[j].Get_index());
		}
		for (int j=0;j<NodeArrays->NodeSpecialWall.size();j++)
		{
			MacroVariables(NodeArrays->NodeSpecialWall[j].Get_index());
		}
		for (int j=0;j<NodeArrays->NodeSymmetry.size();j++)
		{
			MacroVariables(NodeArrays->NodeSymmetry[j].Get_index());
		}
		for (int j=0;j<NodeArrays->NodePeriodic.size();j++)
		{
			MacroVariables(NodeArrays->NodePeriodic[j].Get_index());
		}
}

void D2Q9::UpdateMacroVariables_WithBodyForce(){
		for (int j=0;j<NodeArrays->NodeInterior.size();j++)
		{
			MacroVariablesWithForce(NodeArrays->NodeInterior[j].Get_index());
		}
		for (int j=0;j<NodeArrays->NodeCorner.size();j++)
		{
			MacroVariablesWithForce(NodeArrays->NodeCorner[j].Get_index());
		}
		for (int j=0;j<NodeArrays->NodeGlobalCorner.size();j++)
		{
			MacroVariablesWithForce(NodeArrays->NodeGlobalCorner[j].Get_index());
		}
		for (int j=0;j<NodeArrays->NodeVelocity.size();j++)
		{
			MacroVariablesWithForce(NodeArrays->NodeVelocity[j].Get_index());
		}
		for (int j=0;j<NodeArrays->NodePressure.size();j++)
		{
			MacroVariablesWithForce(NodeArrays->NodePressure[j].Get_index());
		}
		for (int j=0;j<NodeArrays->NodeWall.size();j++)
		{
			MacroVariablesWithForce(NodeArrays->NodeWall[j].Get_index());
		}
		for (int j=0;j<NodeArrays->NodeSpecialWall.size();j++)
		{
			MacroVariablesWithForce(NodeArrays->NodeSpecialWall[j].Get_index());
		}
		for (int j=0;j<NodeArrays->NodeSymmetry.size();j++)
		{
			MacroVariablesWithForce(NodeArrays->NodeSymmetry[j].Get_index());
		}
		for (int j=0;j<NodeArrays->NodePeriodic.size();j++)
		{
			MacroVariablesWithForce(NodeArrays->NodePeriodic[j].Get_index());
		}
}
void D2Q9::MacroVariables(int& idx){

		U[0][idx]=0;
		U[1][idx]=0;
		Rho[idx]=0;
		for (int k=0; k<nbvelo;k++)
		{
			Rho[idx]+=f->f[k][idx];
			for (int j=0;j<2;j++)
			{
				U[j][idx]+=f->f[k][idx]*Ei[k][j];
			}
		}
		U[0][idx]=U[0][idx]/Rho[idx];
		U[1][idx]=U[1][idx]/Rho[idx];

}
void D2Q9::MacroVariablesWithForce(int& idx){

		U[0][idx]=0;
		U[1][idx]=0;
		Rho[idx]=0;
		for (int k=0; k<nbvelo;k++)
		{
			Rho[idx]+=f->f[k][idx];
			for (int j=0;j<2;j++)
			{
				U[j][idx]+=f->f[k][idx]*Ei[k][j];
			}
		}
		U[0][idx]=(U[0][idx]+0.5*F[0][idx])/Rho[idx];
		U[1][idx]=(U[1][idx]+0.5*F[1][idx])/Rho[idx];

}
void D2Q9::TmptoDistri(unsigned int& direction){
	tmp= f->f[direction];
	f->f[direction]=ftmp;
	ftmp=tmp;
}


void D2Q9::Set_BcType(){
	for (int j=0;j<NodeArrays->NodeCorner.size();j++)
	{
		Set_CornerType(NodeArrays->NodeCorner[j]);
	}
	for (int j=0;j<NodeArrays->NodeGlobalCorner.size();j++)
	{
		Set_CornerType(NodeArrays->NodeGlobalCorner[j]);
	}
	for (int j=0;j<NodeArrays->NodeVelocity.size();j++)
	{
		Set_VelocityType(NodeArrays->NodeVelocity[j]);
	}

	for (int j=0;j<NodeArrays->NodePressure.size();j++)
	{
		Set_PressureType(NodeArrays->NodePressure[j]);
	}
	for (int j=0;j<NodeArrays->NodeWall.size();j++)
	{
		Set_WallType(NodeArrays->NodeWall[j]);
	}
	for (int j=0;j<NodeArrays->NodeSpecialWall.size();j++)
	{
		Set_WallType(NodeArrays->NodeSpecialWall[j]);
	}
	for (int j=0;j<NodeArrays->NodeSymmetry.size();j++)
	{
		Set_SymmetryType(NodeArrays->NodeSymmetry[j]);
	}
	for (int j=0;j<NodeArrays->NodePeriodic.size();j++)
	{
		Set_PeriodicType(NodeArrays->NodePeriodic[j]);
	}
	for (int j=0;j<NodeArrays->NodeGhost.size();j++)
	{
		Set_GhostType(NodeArrays->NodeGhost[j]);
	}
}
void D2Q9::Set_GhostType(NodeGhost2D& NodeIn){
	bool GhostStreaming[9];
	StreamingOrientation(NodeIn,GhostStreaming);
	NodeIn.Set_stream(GhostStreaming,nbvelo);

}
void D2Q9::Set_WallType(NodeWall2D& NodeIn){
	bool WallStreaming[9];
	StreamingOrientation(NodeIn,WallStreaming);
	NodeIn.Set_stream(WallStreaming,nbvelo);
	NodeIn.Set_RhoDef(Rho[NodeIn.Get_index()]);
	NodeIn.Set_UDef(U[0][NodeIn.Get_index()],U[1][NodeIn.Get_index()]);
}
void D2Q9::Set_SymmetryType(NodeSymmetry2D& NodeIn){
	bool SymmetryStreaming[9];
	StreamingOrientation(NodeIn,SymmetryStreaming);
	NodeIn.Set_stream(SymmetryStreaming,nbvelo);
}

void D2Q9::Set_PeriodicType(NodePeriodic2D& NodeIn){
	bool PeriodicStreaming[9];
	StreamingOrientation(NodeIn,PeriodicStreaming);
	NodeIn.Set_stream(PeriodicStreaming,nbvelo);
}
void D2Q9::Set_CornerType(NodeCorner2D& NodeIn){
	bool CornerStreaming[9];
	StreamingOrientation(NodeIn,CornerStreaming);
	NodeIn.Set_stream(CornerStreaming,nbvelo);
	NodeIn.Set_UDef(U[0][NodeIn.Get_index()],U[1][NodeIn.Get_index()]);
	NodeIn.Set_RhoDef(Rho[NodeIn.Get_index()]);
}
void D2Q9::Set_VelocityType(NodeVelocity2D& NodeIn){
	bool VelocityStreaming[9];
	StreamingOrientation(NodeIn,VelocityStreaming);
	NodeIn.Set_stream(VelocityStreaming,nbvelo);
	NodeIn.Set_UDef(U[0][NodeIn.Get_index()],U[1][NodeIn.Get_index()]);
}
void D2Q9::Set_PressureType(NodePressure2D& NodeIn){
	bool PressureStreaming[9];
	StreamingOrientation(NodeIn,PressureStreaming);
	NodeIn.Set_stream(PressureStreaming,nbvelo);
	NodeIn.Set_RhoDef(Rho[NodeIn.Get_index()]);
}
void D2Q9::ApplyPatchVelocity(VelocityPatchBc& VelPatchBc){
	SetVelocity(VelPatchBc.Get_VelocityModel(),VelPatchBc.Get_VelocityType());
	std::vector<int> NodeIdx=VelPatchBc.Get_NodeIndexByType();
	std::vector<int> NodeIdxSpecialWalls=VelPatchBc.Get_NodeIndexByTypeSpecialWalls();
	std::vector<int> NodeIdxGlobalCorner=VelPatchBc.Get_NodeIndexByTypeGlobalCorner();

	for (int j=0;j<NodeIdx.size();j++)
	{
		ApplyVelocity(NodeArrays->NodeVelocity[NodeIdx[j]].Get_BcNormal(),NodeArrays->NodeVelocity[NodeIdx[j]].Get_connect(),NodeArrays->NodeVelocity[NodeIdx[j]].Get_UDef(), f,Rho,U[0],U[1]);
	}
	for (int j=0;j<NodeIdxSpecialWalls.size();j++)
	{
		ExtrapolationOnCornerConcave(Rho,NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]].Get_connect(),NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]].Get_BcNormal());
		ApplyVelocityWall(NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]],U[0][NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]].Get_index()],U[1][NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]].Get_index()],NodeArrays->TypeOfNode,f,Rho,U[0],U[1]);
	}
	for (int j=0;j<NodeIdxGlobalCorner.size();j++)
	{
		ApplyGlobalCorner(NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]],NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]].Get_RhoDef(),NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]].Get_UDef(),NodeArrays->TypeOfNode,f,Rho,U[0],U[1]);
	}

}
void D2Q9::ApplyPatchPressure(PressurePatchBc& PresPatchBc){
	SetPressure(PresPatchBc.Get_PressureModel(),PresPatchBc.Get_PressureType());
	std::vector<int> NodeIdx=PresPatchBc.Get_NodeIndexByType();
	std::vector<int> NodeIdxSpecialWalls=PresPatchBc.Get_NodeIndexByTypeSpecialWalls();
	std::vector<int> NodeIdxGlobalCorner=PresPatchBc.Get_NodeIndexByTypeGlobalCorner();

	for (int j=0;j<NodeIdx.size();j++)
	{
		ApplyPressure(NodeArrays->NodePressure[NodeIdx[j]].Get_BcNormal(),NodeArrays->NodePressure[NodeIdx[j]].Get_connect(),NodeArrays->NodePressure[NodeIdx[j]].Get_RhoDef(), f,Rho,U[0],U[1]);
	}
	for (int j=0;j<NodeIdxSpecialWalls.size();j++)
	{
		ExtrapolationOnCornerConcave(Rho,NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]].Get_connect(),NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]].Get_BcNormal());
		ApplyPressureWall(NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]],Rho[NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]].Get_index()],NodeArrays->TypeOfNode,f,Rho,U[0],U[1]);
	}
	for (int j=0;j<NodeIdxGlobalCorner.size();j++)
	{
		ApplyGlobalCorner(NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]],NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]].Get_RhoDef(),NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]].Get_UDef(),NodeArrays->TypeOfNode,f,Rho,U[0],U[1]);
	}

}
void D2Q9::ApplyPatchSymmetry(SymmetryPatchBc& SymPatchBc){
	SetSymmetry(SymPatchBc.Get_SymmetryType());
	std::vector<int> NodeIdx=SymPatchBc.Get_NodeIndexByType();
	std::vector<int> NodeIdxSpecialWalls=SymPatchBc.Get_NodeIndexByTypeSpecialWalls();
	std::vector<int> NodeIdxGlobalCorner=SymPatchBc.Get_NodeIndexByTypeGlobalCorner();

	for (int j=0;j<NodeIdx.size();j++)
	{
			ApplySymmetry(NodeArrays->NodeSymmetry[NodeIdx[j]].Get_BcNormal(),NodeArrays->NodeSymmetry[NodeIdx[j]].Get_connect(),NodeArrays->NodeSymmetry[NodeIdx[j]].Get_RhoDef(),NodeArrays->NodeSymmetry[NodeIdx[j]].Get_UDef(),f,Rho,U[0],U[1]);
	}
	for (int j=0;j<NodeIdxSpecialWalls.size();j++)
	{
		ApplySymmetryWall(NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]],Rho[NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]].Get_index()],U[0][NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]].Get_index()],U[1][NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]].Get_index()],NodeArrays->TypeOfNode,f,Rho,U[0],U[1]);
	}
	for (int j=0;j<NodeIdxGlobalCorner.size();j++)
	{
		ApplyGlobalCorner(NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]],NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]].Get_RhoDef(),NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]].Get_UDef(),NodeArrays->TypeOfNode,f,Rho,U[0],U[1]);
	}

}
void D2Q9::ApplyPatchPeriodic(PeriodicPatchBc& PerPatchBc){
	SetPeriodic(PerPatchBc.Get_PeriodicType());
	std::vector<int> NodeIdx=PerPatchBc.Get_NodeIndexByType();
	std::vector<int> NodeIdxSpecialWalls=PerPatchBc.Get_NodeIndexByTypeSpecialWalls();
	std::vector<int> NodeIdxGlobalCorner=PerPatchBc.Get_NodeIndexByTypeGlobalCorner();
	for (int j=0;j<NodeIdx.size();j++)
	{
			ApplyPeriodic(NodeArrays->NodePeriodic[NodeIdx[j]].Get_BcNormal(),NodeArrays->NodePeriodic[NodeIdx[j]].Get_connect(),NodeArrays->NodePeriodic[NodeIdx[j]].Get_RhoDef(),NodeArrays->NodePeriodic[NodeIdx[j]].Get_UDef(),f,Rho,U[0],U[1]);
	}
	for (int j=0;j<NodeIdxSpecialWalls.size();j++)
	{
		ApplyPeriodicWall(NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]],Rho[NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]].Get_index()],U[0][NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]].Get_index()],U[1][NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]].Get_index()],NodeArrays->TypeOfNode,f,Rho,U[0],U[1]);
	}
	for (int j=0;j<NodeIdxGlobalCorner.size();j++)
	{
		ApplyGlobalCorner(NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]],NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]].Get_RhoDef(),NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]].Get_UDef(),NodeArrays->TypeOfNode,f,Rho,U[0],U[1]);
	}
}
void D2Q9::ApplyBc(){
	std::vector<SolverEnum::PatchType> PatchType=PatchsBc->Get_PatchTypeInType();
	std::vector<int> PatchIdInType=PatchsBc->Get_PatchIdInType();

	for (int i=0;i<PatchsBc->Get_NumberOfPatchBc();i++)
	{
		switch(PatchType[i])
		{
		case SolverEnum::Periodic:
			ApplyPatchPeriodic(PatchsBc->Get_PeriodicPatch()[PatchIdInType[i]]);
			break;
		case SolverEnum::Symmetry:
			ApplyPatchSymmetry(PatchsBc->Get_SymmetryPatch()[PatchIdInType[i]]);
			break;
		case SolverEnum::Pressure:
			ApplyPatchPressure(PatchsBc->Get_PressurePatch()[PatchIdInType[i]]);
			break;
		case SolverEnum::Velocity:
			ApplyPatchVelocity(PatchsBc->Get_VelocityPatch()[PatchIdInType[i]]);
			break;
		case SolverEnum::Wall:
			break;
		}
	}
/*
	for (int j=0;j<NodeArrays->NodeVelocity.size();j++)
	{
		ApplyVelocity(NodeArrays->NodeVelocity[j].Get_BcNormal(),NodeArrays->NodeVelocity[j].Get_connect(),NodeArrays->NodeVelocity[j].Get_UDef(), f);

	}

	for (int j=0;j<NodeArrays->NodePressure.size();j++)
	{
		ApplyPressure(NodeArrays->NodePressure[j].Get_BcNormal(),NodeArrays->NodePressure[j].Get_connect(),NodeArrays->NodePressure[j].Get_RhoDef(), f);
	}
	*/
	for (int j=0;j<NodeArrays->NodeWall.size();j++)
	{
		ApplyWall(NodeArrays->NodeWall[j].Get_BcNormal(),NodeArrays->NodeWall[j].Get_connect(),f);
	}
	/*
	for (int j=0;j<NodeArrays->NodeSpecialWall.size();j++)
	{
		ExtrapolationOnCornerConcave(Rho,NodeArrays->NodeSpecialWall[j].Get_connect(),NodeArrays->NodeSpecialWall[j].Get_BcNormal());
		ApplySpecialWall(NodeArrays->NodeSpecialWall[j],Rho[NodeArrays->NodeSpecialWall[j].Get_index()],U[0][NodeArrays->NodeSpecialWall[j].Get_index()],U[1][NodeArrays->NodeSpecialWall[j].Get_index()],NodeArrays->TypeOfNode,f,Rho,U[0],U[1]);

	}
	*/
	for (int j=0;j<NodeArrays->NodeCorner.size();j++)
	{
		ApplyCornerWall(NodeArrays->NodeCorner[j], f);
	}
/*
	for (int j=0;j<NodeArrays->NodeSymmetry.size();j++)
	{
			ApplySymmetry(NodeArrays->NodeSymmetry[j].Get_BcNormal(),NodeArrays->NodeSymmetry[j].Get_connect(),NodeArrays->NodeSymmetry[j].Get_RhoDef(),NodeArrays->NodeSymmetry[j].Get_UDef(),f);
	}
	for (int j=0;j<NodeArrays->NodePeriodic.size();j++)
	{
			ApplyPeriodic(NodeArrays->NodePeriodic[j].Get_BcNormal(),NodeArrays->NodePeriodic[j].Get_connect(),NodeArrays->NodePeriodic[j].Get_RhoDef(),NodeArrays->NodePeriodic[j].Get_UDef(),f);
	}
	*/
	/*
	for (int j=0;j<NodeArrays->NodeGlobalCorner.size();j++)
	{
		ApplyGlobalCorner(NodeArrays->NodeGlobalCorner[j],NodeArrays->TypeOfNode,f);
	}
*/

   parallel->barrier();

}
double D2Q9::Cal_RhoCorner(NodeCorner2D& nodeIn){

	unsigned int direction1,direction2;
	switch(nodeIn.Get_BcNormal())
	{
	case 5:
		direction1=1;
		direction2=2;
		doubleTmpReturn=(Rho[nodeIn.Get_connect()[direction1]]+Rho[nodeIn.Get_connect()[direction2]])*0.5;//Connect(nodenumber,direction1)]+Rho[Connect(nodenumber,direction2)])*0.5;
		break;
	case 6:
		direction1=2;
		direction2=3;
		doubleTmpReturn=(Rho[nodeIn.Get_connect()[direction1]]+Rho[nodeIn.Get_connect()[direction2]])*0.5;//(Rho[Connect(nodenumber,direction1)]+Rho[Connect(nodenumber,direction2)])*0.5;
		break;
	case 7:
		direction1=3;
		direction2=4;
		doubleTmpReturn=(Rho[nodeIn.Get_connect()[direction1]]+Rho[nodeIn.Get_connect()[direction2]])*0.5;//(Rho[Connect(nodenumber,direction1)]+Rho[Connect(nodenumber,direction2)])*0.5;
		break;
	case 8:
		direction1=1;
		direction2=4;
		doubleTmpReturn=(Rho[nodeIn.Get_connect()[direction1]]+Rho[nodeIn.Get_connect()[direction2]])*0.5;//(Rho[Connect(nodenumber,direction1)]+Rho[Connect(nodenumber,direction2)])*0.5;
		break;
	}

	return doubleTmpReturn;
}
void D2Q9::StreamingOrientation(NodeGhost2D& nodeIn, bool Streaming[9]){
	Streaming[0]=false;
	for (unsigned int i=1;i<(unsigned int)nbvelo;i++)
		if((nodeIn.Get_connect()[i]==nodeIn.Get_index())||(nodeIn.get_NodeType()==Solid)||(NodeArrays->TypeOfNode[nodeIn.Get_connect()[i]]==Solid)||(NodeArrays->TypeOfNode[nodeIn.Get_connect()[i]]==Ghost))
			Streaming[i]=false;
		else
			Streaming[i]=true;
}
void D2Q9::StreamingOrientation(NodeCorner2D& nodeIn, bool Streaming[9]){
	if(		NodeArrays->TypeOfNode[nodeIn.Get_connect()[1]]!=Periodic &&
			NodeArrays->TypeOfNode[nodeIn.Get_connect()[2]]!=Periodic &&
			NodeArrays->TypeOfNode[nodeIn.Get_connect()[3]]!=Periodic &&
			NodeArrays->TypeOfNode[nodeIn.Get_connect()[4]]!=Periodic	)
	{

	if(nodeIn.Get_CornerType()==Concave)
	{
		switch(nodeIn.Get_BcNormal())
			{
			case 5:
				Streaming[0]=false;
				Streaming[1]=true;
				Streaming[2]=true;
				Streaming[3]=false;
				Streaming[4]=false;
				Streaming[5]=true;
				Streaming[6]=false;
				Streaming[7]=false;
				Streaming[8]=false;
				break;
			case 6:
				Streaming[0]=false;
				Streaming[1]=false;
				Streaming[2]=true;
				Streaming[3]=true;
				Streaming[4]=false;
				Streaming[5]=false;
				Streaming[6]=true;
				Streaming[7]=false;
				Streaming[8]=false;
				break;
			case 7:
				Streaming[0]=false;
				Streaming[1]=false;
				Streaming[2]=false;
				Streaming[3]=true;
				Streaming[4]=true;
				Streaming[5]=false;
				Streaming[6]=false;
				Streaming[7]=true;
				Streaming[8]=false;
				break;
			case 8:
				Streaming[0]=false;
				Streaming[1]=true;
				Streaming[2]=false;
				Streaming[3]=false;
				Streaming[4]=true;
				Streaming[5]=false;
				Streaming[6]=false;
				Streaming[7]=false;
				Streaming[8]=true;
				break;
			}
	}
	else
	{
		switch(nodeIn.Get_BcNormal())
			{
			case 5:
				Streaming[0]=false;
				Streaming[1]=true;
				Streaming[2]=true;
				Streaming[3]=true;
				Streaming[4]=true;
				Streaming[5]=true;
				Streaming[6]=true;
				Streaming[7]=false;
				Streaming[8]=true;
				break;
			case 6:
				Streaming[0]=false;
				Streaming[1]=true;
				Streaming[2]=true;
				Streaming[3]=true;
				Streaming[4]=true;
				Streaming[5]=true;
				Streaming[6]=true;
				Streaming[7]=true;
				Streaming[8]=false;
				break;
			case 7:
				Streaming[0]=false;
				Streaming[1]=true;
				Streaming[2]=true;
				Streaming[3]=true;
				Streaming[4]=true;
				Streaming[5]=false;
				Streaming[6]=true;
				Streaming[7]=true;
				Streaming[8]=true;
				break;
			case 8:
				Streaming[0]=false;
				Streaming[1]=true;
				Streaming[2]=true;
				Streaming[3]=true;
				Streaming[4]=true;
				Streaming[5]=true;
				Streaming[6]=false;
				Streaming[7]=true;
				Streaming[8]=true;
				break;
			}
		}
	}
	else
	{
		for (unsigned int i=1;i<(unsigned int)nbvelo;i++)
				if((nodeIn.Get_connect()[i]==nodeIn.Get_index())||(nodeIn.get_NodeType()==Solid)||(NodeArrays->TypeOfNode[nodeIn.Get_connect()[i]]==Solid))
					Streaming[i]=false;
				else
					Streaming[i]=true;
	}

/*	for (unsigned int i=1;i<(unsigned int)nbvelo;i++)
		if((nodeIn.Get_connect()[i]==nodeIn.Get_index())||(nodeIn.get_NodeType()==Solid)||(NodeArrays->TypeOfNode[nodeIn.Get_connect()[i]]==Solid))
			Streaming[i]=false;
		else
			Streaming[i]=true;*/

}
void D2Q9::StreamingOrientation(NodeWall2D& nodeIn, bool Streaming[9]){
	switch(nodeIn.Get_BcNormal())
		{
		case 1:
			Streaming[0]=false;
			Streaming[1]=true;
			Streaming[2]=true;
			Streaming[3]=false;
			Streaming[4]=true;
			Streaming[5]=true;
			Streaming[6]=false;
			Streaming[7]=false;
			Streaming[8]=true;
			break;
		case 2:
			Streaming[0]=false;
			Streaming[1]=true;
			Streaming[2]=true;
			Streaming[3]=true;
			Streaming[4]=false;
			Streaming[5]=true;
			Streaming[6]=true;
			Streaming[7]=false;
			Streaming[8]=false;
			break;
		case 3:
			Streaming[0]=false;
			Streaming[1]=false;
			Streaming[2]=true;
			Streaming[3]=true;
			Streaming[4]=true;
			Streaming[5]=false;
			Streaming[6]=true;
			Streaming[7]=true;
			Streaming[8]=false;
			break;
		case 4:
			Streaming[0]=false;
			Streaming[1]=true;
			Streaming[2]=false;
			Streaming[3]=true;
			Streaming[4]=true;
			Streaming[5]=false;
			Streaming[6]=false;
			Streaming[7]=true;
			Streaming[8]=true;
			break;
// For special Walls
		case 5:
			Streaming[0]=false;
			Streaming[1]=true;
			Streaming[2]=true;
			Streaming[3]=false;
			Streaming[4]=false;
			Streaming[5]=true;
			Streaming[6]=false;
			Streaming[7]=false;
			Streaming[8]=false;
			break;
		case 6:
			Streaming[0]=false;
			Streaming[1]=false;
			Streaming[2]=true;
			Streaming[3]=true;
			Streaming[4]=false;
			Streaming[5]=false;
			Streaming[6]=true;
			Streaming[7]=false;
			Streaming[8]=false;
			break;
		case 7:
			Streaming[0]=false;
			Streaming[1]=false;
			Streaming[2]=false;
			Streaming[3]=true;
			Streaming[4]=true;
			Streaming[5]=false;
			Streaming[6]=false;
			Streaming[7]=true;
			Streaming[8]=false;
			break;
		case 8:
			Streaming[0]=false;
			Streaming[1]=true;
			Streaming[2]=false;
			Streaming[3]=false;
			Streaming[4]=true;
			Streaming[5]=false;
			Streaming[6]=false;
			Streaming[7]=false;
			Streaming[8]=true;
			break;
		default:
			std::cerr<<" Problem in setup streaming; Normal not found."<<std::endl;
		}
}
void D2Q9::StreamingOrientation(NodeSymmetry2D& nodeIn, bool Streaming[9]){
	switch(nodeIn.Get_BcNormal())
		{
		case 1:
			Streaming[0]=false;
			Streaming[1]=true;
			Streaming[2]=true;
			Streaming[3]=false;
			Streaming[4]=true;
			Streaming[5]=true;
			Streaming[6]=false;
			Streaming[7]=false;
			Streaming[8]=true;
			break;
		case 2:
			Streaming[0]=false;
			Streaming[1]=true;
			Streaming[2]=true;
			Streaming[3]=true;
			Streaming[4]=false;
			Streaming[5]=true;
			Streaming[6]=true;
			Streaming[7]=false;
			Streaming[8]=false;
			break;
		case 3:
			Streaming[0]=false;
			Streaming[1]=false;
			Streaming[2]=true;
			Streaming[3]=true;
			Streaming[4]=true;
			Streaming[5]=false;
			Streaming[6]=true;
			Streaming[7]=true;
			Streaming[8]=false;
			break;
		case 4:
			Streaming[0]=false;
			Streaming[1]=true;
			Streaming[2]=false;
			Streaming[3]=true;
			Streaming[4]=true;
			Streaming[5]=false;
			Streaming[6]=false;
			Streaming[7]=true;
			Streaming[8]=true;
			break;
		}
}
void D2Q9::StreamingOrientation(NodeVelocity2D& nodeIn, bool Streaming[9]){
	Streaming[0]=false;
	for (unsigned int i=1;i<(unsigned int)nbvelo;i++)
		if((nodeIn.Get_connect()[i]==nodeIn.Get_index())||(nodeIn.get_NodeType()==Solid)||(NodeArrays->TypeOfNode[nodeIn.Get_connect()[i]]==Solid))
			Streaming[i]=false;
		else
			Streaming[i]=true;
}
void D2Q9::StreamingOrientation(NodePressure2D& nodeIn, bool Streaming[9]){
	Streaming[0]=false;
	for (unsigned int i=1;i<(unsigned int)nbvelo;i++)
		if((nodeIn.Get_connect()[i]==nodeIn.Get_index())||(nodeIn.get_NodeType()==Solid)||(NodeArrays->TypeOfNode[nodeIn.Get_connect()[i]]==Solid))
			Streaming[i]=false;
		else
			Streaming[i]=true;
}
void D2Q9::StreamingOrientation(NodePeriodic2D& nodeIn, bool Streaming[9]){
	Streaming[0]=false;
	for (unsigned int i=1;i<(unsigned int)nbvelo;i++)
		if((nodeIn.Get_connect()[i]==nodeIn.Get_index())||(nodeIn.get_NodeType()==Solid)||(NodeArrays->TypeOfNode[nodeIn.Get_connect()[i]]==Solid))
			Streaming[i]=false;
		else
			Streaming[i]=true;
}



//Communication Functions
void D2Q9::IniComVariables(){
	MultiBlock_->Get_Connect_Node(IdNodeN,IdNodeE,IdNodeS,IdNodeW,IdNodeSW,IdNodeSE,IdNodeNW,IdNodeNE);
	MultiBlock_->Get_Connect_Node(IdRNodeN,IdRNodeE,IdRNodeS,IdRNodeW,IdGNodeN,IdGNodeE,IdGNodeS,IdGNodeW,
			IdRNodeSW,IdRNodeSE,IdRNodeNW,IdRNodeNE,IdGNodeSW,IdGNodeSE,IdGNodeNW,IdGNodeNE);
	MultiBlock_->Get_Connect_SolidNode(SolidIdRNodeN,SolidIdRNodeE,SolidIdRNodeS,SolidIdRNodeW,SolidIdGNodeN,SolidIdGNodeE,SolidIdGNodeS,SolidIdGNodeW,
			SolidIdRNodeSW,SolidIdRNodeSE,SolidIdRNodeNW,SolidIdRNodeNE,SolidIdGNodeSW,SolidIdGNodeSE,SolidIdGNodeNW,SolidIdGNodeNE);

	Nd_variables_sync=3;
	buf_send=new double** [Nd_variables_sync];
	buf_recv=new double** [Nd_variables_sync];
	for (int i=0;i<Nd_variables_sync;i++)
	{
		buf_send[i]=new double* [4];
		buf_recv[i]=new double* [4];
	}

	size_buf=new int [4];

	for (int i=0;i<Nd_variables_sync;i++)
	{
		buf_send[i][0]=new double[IdRNodeE.size()];
		buf_send[i][1]=new double[IdRNodeW.size()];
		buf_send[i][2]=new double[IdRNodeS.size()];
		buf_send[i][3]=new double[IdRNodeN.size()];
		buf_recv[i][0]=new double[IdGNodeW.size()];
		buf_recv[i][1]=new double[IdGNodeE.size()];
		buf_recv[i][2]=new double[IdGNodeN.size()];
		buf_recv[i][3]=new double[IdGNodeS.size()];
	}
	size_buf[0]=IdRNodeE.size();
	size_buf[1]=IdRNodeW.size();
	size_buf[2]=IdRNodeS.size();
	size_buf[3]=IdRNodeN.size();
// Macro sync
	Nd_MacroVariables_sync=Dic->Get_NbSyncVar();//6;
	std::cout<<"Synchromisation ; number of variable: "<<Nd_MacroVariables_sync<<std::endl;
	SyncVar=Dic->Get_SyncVar();
	buf_MacroSend=new double** [Nd_MacroVariables_sync];
	buf_MacroRecv=new double** [Nd_MacroVariables_sync];
	for (int i=0;i<Nd_MacroVariables_sync;i++)
	{
		buf_MacroSend[i]=new double* [4];
		buf_MacroRecv[i]=new double* [4];
	}

	size_MacroBuf=new int [8];

	for (int i=0;i<Nd_MacroVariables_sync;i++)
	{
		buf_MacroSend[i][0]=new double[IdRNodeE.size()];
		buf_MacroSend[i][1]=new double[IdRNodeW.size()];
		buf_MacroSend[i][2]=new double[IdRNodeS.size()];
		buf_MacroSend[i][3]=new double[IdRNodeN.size()];
		buf_MacroRecv[i][0]=new double[IdGNodeW.size()];
		buf_MacroRecv[i][1]=new double[IdGNodeE.size()];
		buf_MacroRecv[i][2]=new double[IdGNodeN.size()];
		buf_MacroRecv[i][3]=new double[IdGNodeS.size()];
	}
	size_MacroBuf[0]=IdRNodeE.size();
	size_MacroBuf[1]=IdRNodeW.size();
	size_MacroBuf[2]=IdRNodeS.size();
	size_MacroBuf[3]=IdRNodeN.size();
	size_MacroBuf[4]=IdGNodeW.size();
	size_MacroBuf[5]=IdGNodeE.size();
	size_MacroBuf[6]=IdGNodeN.size();
	size_MacroBuf[7]=IdGNodeS.size();
}
void D2Q9::GhostNodesSyncFromGhost(){

	for (unsigned int i=0;i<IdNodeE.size();i++)
	{
		buf_send[0][0][i]=f->f[1][IdNodeE[i]];
		buf_send[1][0][i]=f->f[5][IdNodeE[i]];
		buf_send[2][0][i]=f->f[8][IdNodeE[i]];
	}

	for (unsigned int i=0;i<IdNodeN.size();i++)
	{
		buf_send[0][3][i]=f->f[2][IdNodeN[i]];
		buf_send[1][3][i]=f->f[5][IdNodeN[i]];
		buf_send[2][3][i]=f->f[6][IdNodeN[i]];
	}

	int itmp=0;
	for (int i=0;i<Nd_variables_sync;i++)
	{
		MultiBlock_->CommunicationFromGhost(buf_send[i],buf_recv[i],size_buf);

	}
	for (unsigned int i=0;i<IdNodeW.size();i++)
	{
		f->f[1][IdNodeW[i]]=buf_recv[0][0][i];
		f->f[5][IdNodeW[i]]=buf_recv[1][0][i];
		f->f[8][IdNodeW[i]]=buf_recv[2][0][i];
	}

	for (unsigned int i=0;i<IdNodeS.size();i++)
	{
		f->f[2][IdNodeS[i]]=buf_recv[0][3][i];
		f->f[5][IdNodeS[i]]=buf_recv[1][3][i];
		f->f[6][IdNodeS[i]]=buf_recv[2][3][i];
	}

}
void D2Q9::GhostNodesSyncToGhost(){

	for (unsigned int i=0;i<IdRNodeW.size();i++)
	{
		buf_send[0][1][i]=f->f[3][IdRNodeW[i]];
		buf_send[1][1][i]=f->f[6][IdRNodeW[i]];
		buf_send[2][1][i]=f->f[7][IdRNodeW[i]];
	}
	for (unsigned int i=0;i<IdRNodeS.size();i++)
	{
		buf_send[0][2][i]=f->f[4][IdRNodeS[i]];
		buf_send[1][2][i]=f->f[7][IdRNodeS[i]];
		buf_send[2][2][i]=f->f[8][IdRNodeS[i]];
	}
	for (unsigned int i=0;i<IdRNodeE.size();i++)
	{
		buf_send[0][0][i]=f->f[1][IdRNodeE[i]];
		buf_send[1][0][i]=f->f[5][IdRNodeE[i]];
		buf_send[2][0][i]=f->f[8][IdRNodeE[i]];
	}
	for (unsigned int i=0;i<IdRNodeN.size();i++)
	{
		buf_send[0][3][i]=f->f[2][IdRNodeN[i]];
		buf_send[1][3][i]=f->f[5][IdRNodeN[i]];
		buf_send[2][3][i]=f->f[6][IdRNodeN[i]];
	}

	int tag_x_r=1;
	int tag_x_l=2;
	int tag_y_t=3;
	int tag_y_b=4;
	MPI_Status status;
	if(IdRNodeE.size()>=1)
	{
		MultiBlock_->Send(&buf_send[0][0][0],IdRNodeE.size(),1,tag_x_r);
		MultiBlock_->Send(&buf_send[1][0][0],IdRNodeE.size(),1,tag_x_r);
		MultiBlock_->Send(&buf_send[2][0][0],IdRNodeE.size(),1,tag_x_r);
	}
	if(IdGNodeW.size()>=1)
	{
		MultiBlock_->Recv(&buf_recv[0][0][0],IdGNodeW.size(),3,tag_x_r,status);
		MultiBlock_->Recv(&buf_recv[1][0][0],IdGNodeW.size(),3,tag_x_r,status);
		MultiBlock_->Recv(&buf_recv[2][0][0],IdGNodeW.size(),3,tag_x_r,status);
	}

	if(IdRNodeW.size()>=1)
	{
		MultiBlock_->Send(&buf_send[0][1][0],IdRNodeW.size(),3,tag_x_l);
		MultiBlock_->Send(&buf_send[1][1][0],IdRNodeW.size(),3,tag_x_l);
		MultiBlock_->Send(&buf_send[2][1][0],IdRNodeW.size(),3,tag_x_l);
	}
	if(IdGNodeE.size()>=1)
	{
		MultiBlock_->Recv(&buf_recv[0][1][0],IdGNodeE.size(),1,tag_x_l,status);
		MultiBlock_->Recv(&buf_recv[1][1][0],IdGNodeE.size(),1,tag_x_l,status);
		MultiBlock_->Recv(&buf_recv[2][1][0],IdGNodeE.size(),1,tag_x_l,status);
	}
		if(IdRNodeN.size()>=1)
	{
		MultiBlock_->Send(&buf_send[0][3][0],IdRNodeN.size(),0,tag_y_t);
		MultiBlock_->Send(&buf_send[1][3][0],IdRNodeN.size(),0,tag_y_t);
		MultiBlock_->Send(&buf_send[2][3][0],IdRNodeN.size(),0,tag_y_t);
	}
	if(IdGNodeS.size()>=1)
	{
		MultiBlock_->Recv(&buf_recv[0][3][0],IdGNodeS.size(),2,tag_y_t,status);
		MultiBlock_->Recv(&buf_recv[1][3][0],IdGNodeS.size(),2,tag_y_t,status);
		MultiBlock_->Recv(&buf_recv[2][3][0],IdGNodeS.size(),2,tag_y_t,status);
	}
	if(IdRNodeS.size()>=1)
	{
		MultiBlock_->Send(&buf_send[0][2][0],IdRNodeS.size(),2,tag_y_b);
		MultiBlock_->Send(&buf_send[1][2][0],IdRNodeS.size(),2,tag_y_b);
		MultiBlock_->Send(&buf_send[2][2][0],IdRNodeS.size(),2,tag_y_b);
	}
	if(IdGNodeN.size()>=1)
	{
		MultiBlock_->Recv(&buf_recv[0][2][0],IdGNodeN.size(),0,tag_y_b,status);
		MultiBlock_->Recv(&buf_recv[1][2][0],IdGNodeN.size(),0,tag_y_b,status);
		MultiBlock_->Recv(&buf_recv[2][2][0],IdGNodeN.size(),0,tag_y_b,status);
	}

	for (unsigned int i=0;i<IdGNodeE.size();i++)
	{
		f->f[3][IdGNodeE[i]]=buf_recv[0][1][i];
		f->f[6][IdGNodeE[i]]=buf_recv[1][1][i];
		f->f[7][IdGNodeE[i]]=buf_recv[2][1][i];
	}
	//if(block[0]>=0)
	for (unsigned int i=0;i<IdGNodeN.size();i++)
	{
		f->f[4][IdGNodeN[i]]=buf_recv[0][2][i];
		f->f[7][IdGNodeN[i]]=buf_recv[1][2][i];
		f->f[8][IdGNodeN[i]]=buf_recv[2][2][i];
	}
	//if(block[3]>=0)
	for (unsigned int i=0;i<IdGNodeW.size();i++)
	{
		f->f[1][IdGNodeW[i]]=buf_recv[0][0][i];
		f->f[5][IdGNodeW[i]]=buf_recv[1][0][i];
		f->f[8][IdGNodeW[i]]=buf_recv[2][0][i];
	}
	//if(block[2]>=0)
	for (unsigned int i=0;i<IdGNodeS.size();i++)
	{
		f->f[2][IdGNodeS[i]]=buf_recv[0][3][i];
		f->f[5][IdGNodeS[i]]=buf_recv[1][3][i];
		f->f[6][IdGNodeS[i]]=buf_recv[2][3][i];
	}
}
void D2Q9::CornerNodesSyncFromGhost(){
	 int tag_x_r=1;
	 int tag_y_t=2;
	 int tag_d_tr=3;

	//N=0,E=1,S=2,W=3,NE=4, SE=5, NW=6, SW=7;
	MPI_Status status;
	if(IdNodeNE.size()>=1)
		MultiBlock_->Send(&f->f[5][IdNodeNE[0]],1,4,tag_d_tr);
	if(IdNodeSW.size()>=1)
		MultiBlock_->Recv(&f->f[5][IdNodeSW[0]],1,7,tag_d_tr,status);
	if(IdNodeNW.size()>=1)
		MultiBlock_->Send(&f->f[2][IdNodeNW[0]],1,0,tag_y_t);
	if(IdNodeSW.size()>=1)
		MultiBlock_->Recv(&f->f[2][IdNodeSW[0]],1,2,tag_y_t,status);
	if(IdNodeNW.size()>=1)
		MultiBlock_->Send(&f->f[6][IdNodeNW[0]],1,0,tag_y_t);
	if(IdNodeSW.size()>=1)
		MultiBlock_->Recv(&f->f[6][IdNodeSW[0]],1,2,tag_y_t,status);
	if(IdNodeSE.size()>=1)
		MultiBlock_->Send(&f->f[1][IdNodeSE[0]],1,1,tag_x_r);
	if(IdNodeSW.size()>=1)
		MultiBlock_->Recv(&f->f[1][IdNodeSW[0]],1,3,tag_x_r,status);
	if(IdNodeSE.size()>=1)
		MultiBlock_->Send(&f->f[8][IdNodeSE[0]],1,1,tag_x_r);
	if(IdNodeSW.size()>=1)
		MultiBlock_->Recv(&f->f[8][IdNodeSW[0]],1,3,tag_x_r,status);
}
void D2Q9::CornerNodesSyncToGhost(){
	 int tag_x_l=1;
	 int tag_y_b=2;
	 int tag_d_bl=3;
	//N=0,E=1,S=2,W=3,NE=4, SE=5, NW=6, SW=7;
	MPI_Status status;

	if(IdRNodeSE.size()>=1)
	{
		MultiBlock_->Send(&f->f[8][IdRNodeSE[0]],1,5,tag_x_l);
	}
	if(IdGNodeNW.size()>=1)
	{
		MultiBlock_->Recv(&f->f[8][IdGNodeNW[0]],1,6,tag_x_l,status);
	}
	if(IdRNodeSW.size()>=1)
	{
		MultiBlock_->Send(&f->f[7][IdRNodeSW[0]],1,7,tag_x_l);
	}
	if(IdGNodeNE.size()>=1)
	{
		MultiBlock_->Recv(&f->f[7][IdGNodeNE[0]],1,4,tag_x_l,status);
	}
	if(IdRNodeNE.size()>=1)
	{
		MultiBlock_->Send(&f->f[5][IdRNodeNE[0]],1,4,tag_x_l);
	}
	if(IdGNodeSW.size()>=1)
	{
		MultiBlock_->Recv(&f->f[5][IdGNodeSW[0]],1,7,tag_x_l,status);
	}
	if(IdRNodeNW.size()>=1)
	{
		MultiBlock_->Send(&f->f[6][IdRNodeNW[0]],1,6,tag_d_bl);
	}
	if(IdGNodeSE.size()>=1)
	{
		MultiBlock_->Recv(&f->f[6][IdGNodeSE[0]],1,5,tag_d_bl,status);
	}
}
void D2Q9::SyncFromGhost(){
	GhostNodesSyncFromGhost();
	CornerNodesSyncFromGhost();
}
void D2Q9::SyncToGhost(){
	GhostNodesSyncToGhost();
	CornerNodesSyncToGhost();
}

void D2Q9::SyncMacroVarToGhost(){

	for (unsigned int i=0;i<IdRNodeW.size();i++)
	{
		for(int j=0;j<Nd_MacroVariables_sync;j++)
			buf_MacroSend[j][1][i]=SyncVar[j][IdRNodeW[i]];
	}
	for (unsigned int i=0;i<IdRNodeS.size();i++)
	{
		for(int j=0;j<Nd_MacroVariables_sync;j++)
			buf_MacroSend[j][2][i]=SyncVar[j][IdRNodeS[i]];
	}
	for (unsigned int i=0;i<IdRNodeE.size();i++)
	{
		for(int j=0;j<Nd_MacroVariables_sync;j++)
			buf_MacroSend[j][0][i]=SyncVar[j][IdRNodeE[i]];
	}
	for (unsigned int i=0;i<IdRNodeN.size();i++)
	{
		for(int j=0;j<Nd_MacroVariables_sync;j++)
			buf_MacroSend[j][3][i]=SyncVar[j][IdRNodeN[i]];
	}
	int tag_x_r=1;
	int tag_x_l=2;
	int tag_y_t=3;
	int tag_y_b=4;
	int tag_d_bl=5;
	MPI_Status status;
	if(IdRNodeE.size()>=1)
	{
		for(int i=0;i<Nd_MacroVariables_sync;i++)
			MultiBlock_->Send(&buf_MacroSend[i][0][0],IdRNodeE.size(),1,tag_x_r);
	}
	if(IdGNodeW.size()>=1)
	{
		for(int i=0;i<Nd_MacroVariables_sync;i++)
			MultiBlock_->Recv(&buf_MacroRecv[i][0][0],IdGNodeW.size(),3,tag_x_r,status);
	}

	if(IdRNodeW.size()>=1)
	{
		for(int i=0;i<Nd_MacroVariables_sync;i++)
			MultiBlock_->Send(&buf_MacroSend[i][1][0],IdRNodeW.size(),3,tag_x_l);
	}
	if(IdGNodeE.size()>=1)
	{
		for(int i=0;i<Nd_MacroVariables_sync;i++)
			MultiBlock_->Recv(&buf_MacroRecv[i][1][0],IdGNodeE.size(),1,tag_x_l,status);
	}
		if(IdRNodeN.size()>=1)
	{
		for(int i=0;i<Nd_MacroVariables_sync;i++)
			MultiBlock_->Send(&buf_MacroSend[i][3][0],IdRNodeN.size(),0,tag_y_t);
	}
	if(IdGNodeS.size()>=1)
	{
		for(int i=0;i<Nd_MacroVariables_sync;i++)
			MultiBlock_->Recv(&buf_MacroRecv[i][3][0],IdGNodeS.size(),2,tag_y_t,status);
	}
	if(IdRNodeS.size()>=1)
	{
		for(int i=0;i<Nd_MacroVariables_sync;i++)
			MultiBlock_->Send(&buf_MacroSend[i][2][0],IdRNodeS.size(),2,tag_y_b);
	}
	if(IdGNodeN.size()>=1)
	{
		for(int i=0;i<Nd_MacroVariables_sync;i++)
			MultiBlock_->Recv(&buf_MacroRecv[i][2][0],IdGNodeN.size(),0,tag_y_b,status);
	}
	for (unsigned int i=0;i<IdGNodeE.size();i++)
	{
		for(int j=0;j<Nd_MacroVariables_sync;j++)
			SyncVar[j][IdGNodeE[i]]=buf_MacroRecv[j][1][i];
	}
	for (unsigned int i=0;i<IdGNodeN.size();i++)
	{
		for(int j=0;j<Nd_MacroVariables_sync;j++)
			SyncVar[j][IdGNodeN[i]]=buf_MacroRecv[j][2][i];
	}
	for (unsigned int i=0;i<IdGNodeW.size();i++)
	{
		for(int j=0;j<Nd_MacroVariables_sync;j++)
			SyncVar[j][IdGNodeW[i]]=buf_MacroRecv[j][0][i];
	}

	for (unsigned int i=0;i<IdGNodeS.size();i++)
	{
		for(int j=0;j<Nd_MacroVariables_sync;j++)
			SyncVar[j][IdGNodeS[i]]=buf_MacroRecv[j][3][i];
	}
	//Corner
		if(IdRNodeSE.size()>=1)
		{
			for(int j=0;j<Nd_MacroVariables_sync;j++)
				MultiBlock_->Send(&SyncVar[j][IdRNodeSE[0]],1,5,tag_x_l);
		}
		if(IdGNodeNW.size()>=1)
		{
			for(int j=0;j<Nd_MacroVariables_sync;j++)
				MultiBlock_->Recv(&SyncVar[j][IdGNodeNW[0]],1,6,tag_x_l,status);
		}
		if(IdRNodeSW.size()>=1)
		{
			for(int j=0;j<Nd_MacroVariables_sync;j++)
				MultiBlock_->Send(&SyncVar[j][IdRNodeSW[0]],1,7,tag_x_l);
		}
		if(IdGNodeNE.size()>=1)
		{
			for(int j=0;j<Nd_MacroVariables_sync;j++)
				MultiBlock_->Recv(&SyncVar[j][IdGNodeNE[0]],1,4,tag_x_l,status);
		}
		if(IdRNodeNE.size()>=1)
		{
			for(int j=0;j<Nd_MacroVariables_sync;j++)
				MultiBlock_->Send(&SyncVar[j][IdRNodeNE[0]],1,4,tag_x_l);
		}
		if(IdGNodeSW.size()>=1)
		{
			for(int j=0;j<Nd_MacroVariables_sync;j++)
				MultiBlock_->Recv(&SyncVar[j][IdGNodeSW[0]],1,7,tag_x_l,status);
		}
		if(IdRNodeNW.size()>=1)
		{
			for(int j=0;j<Nd_MacroVariables_sync;j++)
				MultiBlock_->Send(&SyncVar[j][IdRNodeNW[0]],1,6,tag_d_bl);
		}
		if(IdGNodeSE.size()>=1)
		{
			for(int j=0;j<Nd_MacroVariables_sync;j++)
				MultiBlock_->Recv(&SyncVar[j][IdGNodeSE[0]],1,5,tag_d_bl,status);
		}
}
void D2Q9::UpdatePressure(){
	for (int j=0;j<NodeArrays->NodeInterior.size();j++)
	{
		CalculatePressure(NodeArrays->NodeInterior[j].Get_index());
	}
	for (int j=0;j<NodeArrays->NodeCorner.size();j++)
	{
		CalculatePressure(NodeArrays->NodeCorner[j].Get_index());
	}
	for (int j=0;j<NodeArrays->NodeGlobalCorner.size();j++)
	{
		CalculatePressure(NodeArrays->NodeGlobalCorner[j].Get_index());
	}
	for (int j=0;j<NodeArrays->NodeVelocity.size();j++)
	{
		CalculatePressure(NodeArrays->NodeVelocity[j].Get_index());
	}
	for (int j=0;j<NodeArrays->NodePressure.size();j++)
	{
		CalculatePressure(NodeArrays->NodePressure[j].Get_index());
	}
	for (int j=0;j<NodeArrays->NodeWall.size();j++)
	{
		CalculatePressure(NodeArrays->NodeWall[j].Get_index());
	}
	for (int j=0;j<NodeArrays->NodeSpecialWall.size();j++)
	{
		CalculatePressure(NodeArrays->NodeSpecialWall[j].Get_index());
	}
	for (int j=0;j<NodeArrays->NodeSymmetry.size();j++)
	{
		CalculatePressure(NodeArrays->NodeSymmetry[j].Get_index());
	}
	for (int j=0;j<NodeArrays->NodePeriodic.size();j++)
	{
		CalculatePressure(NodeArrays->NodePeriodic[j].Get_index());
	}

}
void D2Q9::CalculatePressure(int const &idx){
	P[idx]=IdealGazIsothermalPressure(Rho[idx]);
}
double D2Q9::IdealGazIsothermalPressure(double const &Rho){
	return Rho/3.0;
}
