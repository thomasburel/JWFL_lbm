/*
 * D2Q9TwoPhases.cpp
 *
 *  Created on: 9 Jun 2015
 *      Author: thomas
 */

#include "D2Q9_TwoPhases.h"
#include <iomanip>
D2Q9TwoPhases::D2Q9TwoPhases() {
	MultiBlock_=0;
	parallel=0;
	Writer=0;
	PtrParameters=0;
	f=0;
	nbvelo=9;
	Nb_VelocityCollide=nbvelo;
	nbnode=0;
	ftmp=0;
	tmpDistribution=0;
	f_tmp=0;
	tmpreturn=0;
	intTmpReturn=0;
	doubleTmpReturn=0;
	buf_recv=0;
	buf_send=0;
	size_buf=0;
	DiagConnect=0;
	Nd_variables_sync=0;
	buf_MacroSend=0;buf_MacroSendSolid=0;
	buf_MacroRecv=0;buf_MacroRecvSolid=0;
	size_MacroBuf=0;size_MacroBufSolid=0;
	Nd_MacroVariables_sync=0;
	Rho1=0; Rho2=0;
	Opposite[0]=0;
	Opposite[1]=3;
	Opposite[2]=4;
	Opposite[3]=1;
	Opposite[4]=2;
	Opposite[5]=7;
	Opposite[6]=8;
	Opposite[7]=5;
	Opposite[8]=6;
}
void D2Q9TwoPhases::InitD2Q9TwoPhases(MultiBlock* MultiBlock__,ParallelManager* parallel__,WriterManager* Writer__, Parameters* Parameters_ ,InitLBM& ini) {
	Opposite[0]=0;
	Opposite[1]=3;
	Opposite[2]=4;
	Opposite[3]=1;
	Opposite[4]=2;
	Opposite[5]=7;
	Opposite[6]=8;
	Opposite[7]=5;
	Opposite[8]=6;


	f_tmp=new double[9];
	f=new DistriFunct*[2];
	f[0]=new DistriFunct(MultiBlock__->Get_nnodes(),Parameters_->Get_NbVelocities());
	f[1]=new DistriFunct(MultiBlock__->Get_nnodes(),Parameters_->Get_NbVelocities());
	Set_Solver(MultiBlock__,parallel__,Writer__,Parameters_);
	ftmp=new double[nbnode];
	PtrFiStream=f[0];
	PtrFiCollide=PtrFiStream;
	InvTau=1.0/PtrParameters->Get_Tau();
	intTmpReturn=0;
	doubleTmpReturn=0;


	EiCollide=Ei;
	omegaCollide=omega;
	Nb_VelocityCollide=nbvelo;


	init(ini);
}
D2Q9TwoPhases::~D2Q9TwoPhases() {
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

}
void D2Q9TwoPhases::init(InitLBM& ini){
	InitAllDomain(ini);
	Set_BcType();
	parallel->barrier();

}
void D2Q9TwoPhases::InitAllDomain(InitLBM& ini){

	InitDomainBc(ini);
	InitWall(ini);
	InitInterior(ini);

	double alpha=0;
	double* pos =new double[2];
	double* U_=new double[2];
	int idx=0;
	for (unsigned int j=0;j<NodeArrays->NodeSolid.size();j++)
	{
		idx=NodeArrays->NodeSolid[j].Get_index();
		pos[0]=NodeArrays->NodeSolid[j].get_x();
		pos[1]=NodeArrays->NodeSolid[j].get_y();
		ini.IniDomainTwoPhases(parallel->getRank(),NodeArrays->NodeSolid[j],0, idx,pos,Rho[idx],U_,alpha);
		U[0][idx]=U_[0];
		U[1][idx]=U_[1];
	}
	delete [] pos;
	delete [] U_;
}
void D2Q9TwoPhases::InitDomainBc(InitLBM& ini){
	double alpha=0;
	double* pos =new double[2];
	double* U_=new double[2];
	int idx=0;
	for (unsigned int j=0;j<NodeArrays->NodeGlobalCorner.size();j++)
	{
		idx=NodeArrays->NodeGlobalCorner[j].Get_index();
		pos[0]=NodeArrays->NodeGlobalCorner[j].get_x();
		pos[1]=NodeArrays->NodeGlobalCorner[j].get_y();
		ini.IniDomainTwoPhases(parallel->getRank(),NodeArrays->NodeGlobalCorner[j],0, idx,pos,Rho[idx],U_,alpha);
		U[0][idx]=U_[0];
		U[1][idx]=U_[1];
		NodeArrays->NodeGlobalCorner[j].Set_UDef(U_[0],U_[1]);
		NodeArrays->NodeGlobalCorner[j].Set_RhoDef(Rho[idx]);
		NodeArrays->NodeGlobalCorner[j].Set_AlphaDef(alpha);
	}
	for (unsigned int j=0;j<NodeArrays->NodeVelocity.size();j++)
	{
		idx=NodeArrays->NodeVelocity[j].Get_index();
		pos[0]=NodeArrays->NodeVelocity[j].get_x();
		pos[1]=NodeArrays->NodeVelocity[j].get_y();
		ini.IniDomainTwoPhases(parallel->getRank(),NodeArrays->NodeVelocity[j],0, idx,pos,Rho[idx],U_,alpha);
		U[0][idx]=U_[0];
		U[1][idx]=U_[1];
		NodeArrays->NodeVelocity[j].Set_UDef(U_[0],U_[1]);
		NodeArrays->NodeVelocity[j].Set_AlphaDef(alpha);;
	}

	for (unsigned int j=0;j<NodeArrays->NodePressure.size();j++)
	{
		idx=NodeArrays->NodePressure[j].Get_index();
		pos[0]=NodeArrays->NodePressure[j].get_x();
		pos[1]=NodeArrays->NodePressure[j].get_y();
		ini.IniDomainTwoPhases(parallel->getRank(),NodeArrays->NodePressure[j],0, idx,pos,Rho[idx],U_,alpha);
		U[0][idx]=U_[0];
		U[1][idx]=U_[1];
		NodeArrays->NodePressure[j].Set_RhoDef(Rho[idx]);
		NodeArrays->NodePressure[j].Set_AlphaDef(alpha);
	}
	for (unsigned int j=0;j<NodeArrays->NodeSymmetry.size();j++)
	{
		idx=NodeArrays->NodeSymmetry[j].Get_index();
		pos[0]=NodeArrays->NodeSymmetry[j].get_x();
		pos[1]=NodeArrays->NodeSymmetry[j].get_y();
		ini.IniDomainTwoPhases(parallel->getRank(),NodeArrays->NodeSymmetry[j],0, idx,pos,Rho[idx],U_,alpha);
		U[0][idx]=U_[0];
		U[1][idx]=U_[1];
	}
	for (unsigned int j=0;j<NodeArrays->NodePeriodic.size();j++)
	{
		idx=NodeArrays->NodePeriodic[j].Get_index();
		pos[0]=NodeArrays->NodePeriodic[j].get_x();
		pos[1]=NodeArrays->NodePeriodic[j].get_y();
		ini.IniDomainTwoPhases(parallel->getRank(),NodeArrays->NodePeriodic[j],0, idx,pos,Rho[idx],U_,alpha);
		U[0][idx]=U_[0];
		U[1][idx]=U_[1];
	}
	delete [] pos;
	delete [] U_;
}
void D2Q9TwoPhases::InitWall(InitLBM& ini){
	double alpha=0;
	double* pos =new double[2];
	double* U_=new double[2];
	int idx=0;
	for (unsigned int j=0;j<NodeArrays->NodeCorner.size();j++)
	{
		idx=NodeArrays->NodeCorner[j].Get_index();
		pos[0]=NodeArrays->NodeCorner[j].get_x();
		pos[1]=NodeArrays->NodeCorner[j].get_y();
		ini.IniDomainTwoPhases(parallel->getRank(),NodeArrays->NodeCorner[j],0, idx,pos,Rho[idx],U_,alpha);
		U[0][idx]=U_[0];
		U[1][idx]=U_[1];
		NodeArrays->NodeCorner[j].Set_UDef(U_[0],U_[1]);
		NodeArrays->NodeCorner[j].Set_RhoDef(Rho[idx]);
		NodeArrays->NodeCorner[j].Set_AlphaDef(alpha);
	}
	for (unsigned int j=0;j<NodeArrays->NodeWall.size();j++)
	{
		idx=NodeArrays->NodeWall[j].Get_index();
		pos[0]=NodeArrays->NodeWall[j].get_x();
		pos[1]=NodeArrays->NodeWall[j].get_y();
		ini.IniDomainTwoPhases(parallel->getRank(),NodeArrays->NodeWall[j],0, idx,pos,Rho[idx],U_,alpha);
		U[0][idx]=U_[0];
		U[1][idx]=U_[1];
	}
	for (unsigned int j=0;j<NodeArrays->NodeSpecialWall.size();j++)
	{
		idx=NodeArrays->NodeSpecialWall[j].Get_index();
		pos[0]=NodeArrays->NodeSpecialWall[j].get_x();
		pos[1]=NodeArrays->NodeSpecialWall[j].get_y();
		ini.IniDomainTwoPhases(parallel->getRank(),NodeArrays->NodeSpecialWall[j],0, idx,pos,Rho[idx],U_,alpha);
		U[0][idx]=U_[0];
		U[1][idx]=U_[1];
		NodeArrays->NodeSpecialWall[j].Set_UDef(U_[0],U_[1]);
		NodeArrays->NodeSpecialWall[j].Set_RhoDef(Rho[idx]);
		NodeArrays->NodeSpecialWall[j].Set_AlphaDef(alpha);
	}

	delete [] pos;
	delete [] U_;
}
void D2Q9TwoPhases::InitInterior(InitLBM& ini){
	double alpha=0;
	double* pos =new double[2];
	double* U_=new double[2];
	int idx=0;

	for (unsigned int j=0;j<NodeArrays->NodeInterior.size();j++)
	{
		idx=NodeArrays->NodeInterior[j].Get_index();
		pos[0]=NodeArrays->NodeInterior[j].get_x();
		pos[1]=NodeArrays->NodeInterior[j].get_y();
		ini.IniDomainTwoPhases(parallel->getRank(),NodeArrays->NodeInterior[j],0, idx,pos,Rho[idx],U_,alpha);
		U[0][idx]=U_[0];
		U[1][idx]=U_[1];
	}
	for (unsigned int j=0;j<NodeArrays->NodeGhost.size();j++)
	{
		idx=NodeArrays->NodeGhost[j].Get_index();
		pos[0]=NodeArrays->NodeGhost[j].get_x();
		pos[1]=NodeArrays->NodeGhost[j].get_y();
		ini.IniDomainTwoPhases(parallel->getRank(),NodeArrays->NodeGhost[j],0, idx,pos,Rho[idx],U_,alpha);
		U[0][idx]=U_[0];
		U[1][idx]=U_[1];
	}
	delete [] pos;
	delete [] U_;
}
void D2Q9TwoPhases::StreamD2Q9() {


	for (int k=0;k<2;k++)
	{
		for (unsigned int i=1;i<(unsigned int)nbvelo;i++) //No need to stream direction 0
		{
			for (unsigned int j=0;j<NodeArrays->NodeInterior.size();j++)
			{
				if (NodeArrays->NodeInterior[j].Get_connect()[i]!=NodeArrays->NodeInterior[j].Get_index())
				{
					ftmp[NodeArrays->NodeInterior[j].Get_connect()[i]]=f[k]->f[i][NodeArrays->NodeInterior[j].Get_index()];
				}
			}

			for (unsigned int j=0;j<NodeArrays->NodeCorner.size();j++)
			{
				if (NodeArrays->NodeCorner[j].stream()[i])
				{
					ftmp[NodeArrays->NodeCorner[j].Get_connect()[i]]=f[k]->f[i][NodeArrays->NodeCorner[j].Get_index()];
				}		}
			for (unsigned int j=0;j<NodeArrays->NodeGlobalCorner.size();j++)
			{
				if (NodeArrays->NodeGlobalCorner[j].stream()[i])
				{
					ftmp[NodeArrays->NodeGlobalCorner[j].Get_connect()[i]]=f[k]->f[i][NodeArrays->NodeGlobalCorner[j].Get_index()];
				}

			}
			for (unsigned int j=0;j<NodeArrays->NodeVelocity.size();j++)
			{
				if (NodeArrays->NodeVelocity[j].stream()[i])
						{
							ftmp[NodeArrays->NodeVelocity[j].Get_connect()[i]]=f[k]->f[i][NodeArrays->NodeVelocity[j].Get_index()];
						}
			}

			for (unsigned int j=0;j<NodeArrays->NodePressure.size();j++)
			{
				if (NodeArrays->NodePressure[j].stream()[i])
						{
							ftmp[NodeArrays->NodePressure[j].Get_connect()[i]]=f[k]->f[i][NodeArrays->NodePressure[j].Get_index()];
						}
			}
			for (unsigned int j=0;j<NodeArrays->NodeWall.size();j++)
			{
				if (NodeArrays->NodeWall[j].stream()[i])
				{
					ftmp[NodeArrays->NodeWall[j].Get_connect()[i]]=f[k]->f[i][NodeArrays->NodeWall[j].Get_index()];
				}
			}
			for (unsigned int j=0;j<NodeArrays->NodeSpecialWall.size();j++)
			{
				if (NodeArrays->NodeSpecialWall[j].stream()[i])
				{
					ftmp[NodeArrays->NodeSpecialWall[j].Get_connect()[i]]=f[k]->f[i][NodeArrays->NodeSpecialWall[j].Get_index()];
				}
			}
			for (unsigned int j=0;j<NodeArrays->NodeSymmetry.size();j++)
			{
				if (NodeArrays->NodeSymmetry[j].stream()[i])
				{
					ftmp[NodeArrays->NodeSymmetry[j].Get_connect()[i]]=f[k]->f[i][NodeArrays->NodeSymmetry[j].Get_index()];
				}
			}
			for (unsigned int j=0;j<NodeArrays->NodePeriodic.size();j++)
			{
				if (NodeArrays->NodePeriodic[j].stream()[i])
				{
					ftmp[NodeArrays->NodePeriodic[j].Get_connect()[i]]=f[k]->f[i][NodeArrays->NodePeriodic[j].Get_index()];
				}
			}
			for (unsigned int j=0;j<NodeArrays->NodeGhost.size();j++)
			{
				if (NodeArrays->NodeGhost[j].stream()[i])
				{
					ftmp[NodeArrays->NodeGhost[j].Get_connect()[i]]=f[k]->f[i][NodeArrays->NodeGhost[j].Get_index()];
				}
			}
			D2Q9TwoPhases::TmptoDistri(i,k);
		}
	}
	if(!PtrParameters->Get_CollisionOnWalls())
	{
//Force no slip condition
		for (int k=0;k<2;k++)
		{
			for (unsigned int j=0;j<NodeArrays->NodeWall.size();j++)
			{
				if(NodeArrays->NodeWall[j].Get_BcNormal()+Opposite[NodeArrays->NodeWall[j].Get_BcNormal()]==4)
				{
					double density_tmp=f[k]->f[2][NodeArrays->NodeWall[j].Get_index()]+f[k]->f[4][NodeArrays->NodeWall[j].Get_index()];
					density_tmp*=0.5;
					f[k]->f[2][NodeArrays->NodeWall[j].Get_index()]=density_tmp;
					f[k]->f[4][NodeArrays->NodeWall[j].Get_index()]=density_tmp;
				}
				else
				{
					double density_tmp=f[k]->f[1][NodeArrays->NodeWall[j].Get_index()]+f[k]->f[3][NodeArrays->NodeWall[j].Get_index()];
					density_tmp*=0.5;
					f[k]->f[1][NodeArrays->NodeWall[j].Get_index()]=density_tmp;
					f[k]->f[3][NodeArrays->NodeWall[j].Get_index()]=density_tmp;
				}
			}
		}
	}
}

void D2Q9TwoPhases::TmptoDistri(unsigned int& direction, int& IdDistri){
	tmpDistribution= f[IdDistri]->f[direction];
	f[IdDistri]->f[direction]=ftmp;
	ftmp=tmpDistribution;
}
void D2Q9TwoPhases::Set_BcType(){
	for (unsigned int j=0;j<NodeArrays->NodeCorner.size();j++)
	{
		Set_CornerType(NodeArrays->NodeCorner[j]);
	}
	for (unsigned int j=0;j<NodeArrays->NodeGlobalCorner.size();j++)
	{
		Set_CornerType(NodeArrays->NodeGlobalCorner[j]);
	}
	for (unsigned int j=0;j<NodeArrays->NodeVelocity.size();j++)
	{
		Set_VelocityType(NodeArrays->NodeVelocity[j]);
	}

	for (unsigned int j=0;j<NodeArrays->NodePressure.size();j++)
	{
		Set_PressureType(NodeArrays->NodePressure[j]);
	}
	for (unsigned int j=0;j<NodeArrays->NodeWall.size();j++)
	{
		Set_WallType(NodeArrays->NodeWall[j]);
	}
	for (unsigned int j=0;j<NodeArrays->NodeSpecialWall.size();j++)
	{
		Set_WallType(NodeArrays->NodeSpecialWall[j]);
	}
	for (unsigned int j=0;j<NodeArrays->NodeSymmetry.size();j++)
	{
		Set_SymmetryType(NodeArrays->NodeSymmetry[j]);
	}
	for (unsigned int j=0;j<NodeArrays->NodePeriodic.size();j++)
	{
		Set_PeriodicType(NodeArrays->NodePeriodic[j]);
	}
	for (unsigned int j=0;j<NodeArrays->NodeGhost.size();j++)
	{
		Set_GhostType(NodeArrays->NodeGhost[j]);
	}
}

void D2Q9TwoPhases::Set_GhostType(NodeGhost2D& NodeIn){
	bool GhostStreaming[9];
	StreamingOrientation(NodeIn,GhostStreaming);
	NodeIn.Set_stream(GhostStreaming,nbvelo);

}
void D2Q9TwoPhases::Set_WallType(NodeWall2D& NodeIn){
	bool WallStreaming[9];
	StreamingOrientation(NodeIn,WallStreaming);
	NodeIn.Set_stream(WallStreaming,nbvelo);
	NodeIn.Set_RhoDef(Rho[NodeIn.Get_index()]);
	NodeIn.Set_UDef(U[0][NodeIn.Get_index()],U[1][NodeIn.Get_index()]);
}
void D2Q9TwoPhases::Set_SymmetryType(NodeSymmetry2D& NodeIn){
	bool SymmetryStreaming[9];
	StreamingOrientation(NodeIn,SymmetryStreaming);
	NodeIn.Set_stream(SymmetryStreaming,nbvelo);

}
void D2Q9TwoPhases::Set_PeriodicType(NodePeriodic2D& NodeIn){
	bool PeriodicStreaming[9];
	StreamingOrientation(NodeIn,PeriodicStreaming);
	NodeIn.Set_stream(PeriodicStreaming,nbvelo);

}
void D2Q9TwoPhases::Set_CornerType(NodeCorner2D& NodeIn){
	bool CornerStreaming[9];
	StreamingOrientation(NodeIn,CornerStreaming);
	NodeIn.Set_stream(CornerStreaming,nbvelo);
	NodeIn.Set_UDef(U[0][NodeIn.Get_index()],U[1][NodeIn.Get_index()]);
	NodeIn.Set_RhoDef(Rho[NodeIn.Get_index()]);
}
void D2Q9TwoPhases::Set_VelocityType(NodeVelocity2D& NodeIn){
	bool VelocityStreaming[9];
	StreamingOrientation(NodeIn,VelocityStreaming);
	NodeIn.Set_stream(VelocityStreaming,nbvelo);
	NodeIn.Set_UDef(U[0][NodeIn.Get_index()],U[1][NodeIn.Get_index()]);
}
void D2Q9TwoPhases::Set_PressureType(NodePressure2D& NodeIn){
	bool PressureStreaming[9];
	StreamingOrientation(NodeIn,PressureStreaming);
	NodeIn.Set_stream(PressureStreaming,nbvelo);
	NodeIn.Set_RhoDef(Rho[NodeIn.Get_index()]);
}
/*
/// Impose velocity by HeZou
void D2Q9TwoPhases::ApplyHeZou_U(NodeVelocity2D& NodeIn, int distID, double &U, double &V){


	for (int i=0;i<9;i++)
	f_tmp[i]=f[distID]->f[i][NodeIn.Get_index()];

	BC_HeZou_U(NodeIn.Get_BcNormal(),f_tmp, U,V);

	for (int i=0;i<9;i++)
	f[distID]->f[i][NodeIn.Get_index()]=f_tmp[i];


}
void D2Q9TwoPhases::ApplyHeZou_U(NodeCorner2D& NodeIn, int normal, int distID, double &U, double &V){

		for (int i=0;i<9;i++)
		f_tmp[i]=f[distID]->f[i][NodeIn.Get_index()];

		BC_HeZou_U(normal,f_tmp, U,V);

		for (int i=0;i<9;i++)
		f[distID]->f[i][NodeIn.Get_index()]=f_tmp[i];

}
/// Impose pressure by HeZou
void D2Q9TwoPhases::ApplyHeZou_P(NodePressure2D& NodeIn, int distID, double Rho, double &U, double &V){

	if(Rho==0)
	{
		for (int i=0;i<9;i++)
					f[distID]->f[i][NodeIn.Get_index()]=0;
	}
	else
	{
		for (int i=0;i<9;i++)
			f_tmp[i]=f[distID]->f[i][NodeIn.Get_index()];
		BC_HeZou_P(NodeIn.Get_BcNormal(),f_tmp,Rho, U,V);
		for (int i=0;i<9;i++)
			f[distID]->f[i][NodeIn.Get_index()]=f_tmp[i];
	}

}
/// Impose pressure by HeZou
void D2Q9TwoPhases::ApplyHeZou_P(NodeCorner2D& NodeIn, int normal, int distID, double Mass, double &U, double &V){

	if(Mass==0)
	{
		for (int i=0;i<9;i++)
					f[distID]->f[i][NodeIn.Get_index()]=0;
	}
	else
	{
		for (int i=0;i<9;i++)
			f_tmp[i]=f[distID]->f[i][NodeIn.Get_index()];
		BC_HeZou_P(normal,f_tmp,Mass, U,V);
		for (int i=0;i<9;i++)
			f[distID]->f[i][NodeIn.Get_index()]=f_tmp[i];
	}

}
//Symmetry treatment on node
void D2Q9TwoPhases::ApplySymmetryOnNode(NodeSymmetry2D& NodeIn){
	switch(NodeIn.Get_BcNormal())
			{
			case 2:
				f[0]->f[2][NodeIn.Get_index()]=f[0]->f[4][NodeIn.Get_index()];
				f[0]->f[5][NodeIn.Get_index()]=f[0]->f[8][NodeIn.Get_index()];
				f[0]->f[6][NodeIn.Get_index()]=f[0]->f[7][NodeIn.Get_index()];

				f[1]->f[2][NodeIn.Get_index()]=f[1]->f[4][NodeIn.Get_index()];
				f[1]->f[5][NodeIn.Get_index()]=f[1]->f[8][NodeIn.Get_index()];
				f[1]->f[6][NodeIn.Get_index()]=f[1]->f[7][NodeIn.Get_index()];
				break;
			case 4:
				f[0]->f[4][NodeIn.Get_index()]=f[0]->f[2][NodeIn.Get_index()];
				f[0]->f[7][NodeIn.Get_index()]=f[0]->f[6][NodeIn.Get_index()];
				f[0]->f[8][NodeIn.Get_index()]=f[0]->f[5][NodeIn.Get_index()];

				f[1]->f[4][NodeIn.Get_index()]=f[1]->f[2][NodeIn.Get_index()];
				f[1]->f[7][NodeIn.Get_index()]=f[1]->f[6][NodeIn.Get_index()];
				f[1]->f[8][NodeIn.Get_index()]=f[1]->f[5][NodeIn.Get_index()];
				break;
			case 1:
				f[0]->f[1][NodeIn.Get_index()]=f[0]->f[3][NodeIn.Get_index()];
				f[0]->f[5][NodeIn.Get_index()]=f[0]->f[6][NodeIn.Get_index()];
				f[0]->f[8][NodeIn.Get_index()]=f[0]->f[7][NodeIn.Get_index()];

				f[1]->f[1][NodeIn.Get_index()]=f[1]->f[3][NodeIn.Get_index()];
				f[1]->f[5][NodeIn.Get_index()]=f[1]->f[6][NodeIn.Get_index()];
				f[1]->f[8][NodeIn.Get_index()]=f[1]->f[7][NodeIn.Get_index()];
				break;
			case 3:
				f[0]->f[3][NodeIn.Get_index()]=f[0]->f[1][NodeIn.Get_index()];
				f[0]->f[6][NodeIn.Get_index()]=f[0]->f[5][NodeIn.Get_index()];
				f[0]->f[7][NodeIn.Get_index()]=f[0]->f[8][NodeIn.Get_index()];

				f[1]->f[3][NodeIn.Get_index()]=f[1]->f[1][NodeIn.Get_index()];
				f[1]->f[6][NodeIn.Get_index()]=f[1]->f[5][NodeIn.Get_index()];
				f[1]->f[7][NodeIn.Get_index()]=f[1]->f[8][NodeIn.Get_index()];
				break;
			default:
				std::cerr<<"Direction symmetry not found"<<std::endl;
				break;
			}
}

void D2Q9TwoPhases::ApplyBounceBackSymmetry(NodeWall2D& NodeIn){

	switch(NodeIn.Get_BcNormal())
			{
			case 2:
				f[0]->f[2][NodeIn.Get_index()]=f[0]->f[Opposite[2]][NodeIn.Get_index()];
				f[0]->f[5][NodeIn.Get_index()]=f[0]->f[Opposite[5]][NodeIn.Get_index()];
				f[0]->f[6][NodeIn.Get_index()]=f[0]->f[Opposite[6]][NodeIn.Get_index()];
				break;
			case 4:
				f[0]->f[4][NodeIn.Get_index()]=f[0]->f[Opposite[4]][NodeIn.Get_index()];
				f[0]->f[7][NodeIn.Get_index()]=f[0]->f[Opposite[7]][NodeIn.Get_index()];
				f[0]->f[8][NodeIn.Get_index()]=f[0]->f[Opposite[8]][NodeIn.Get_index()];
				break;
			case 1:
				//apply Symmetry
				if(NodeArrays->TypeOfNode[NodeIn.Get_connect(0)]==Solid)//Bottom
				{
					f[0]->f[2][NodeIn.Get_index()]=f[0]->f[Opposite[2]][NodeIn.Get_index()];
					f[0]->f[6][NodeIn.Get_index()]=f[0]->f[7][NodeIn.Get_index()];

					f[1]->f[2][NodeIn.Get_index()]=f[1]->f[Opposite[2]][NodeIn.Get_index()];
					f[1]->f[6][NodeIn.Get_index()]=f[1]->f[7][NodeIn.Get_index()];
				}
				else //Top
				{
					f[0]->f[4][NodeIn.Get_index()]=f[0]->f[Opposite[4]][NodeIn.Get_index()];
					f[0]->f[7][NodeIn.Get_index()]=f[0]->f[6][NodeIn.Get_index()];

					f[1]->f[4][NodeIn.Get_index()]=f[1]->f[Opposite[4]][NodeIn.Get_index()];
					f[1]->f[7][NodeIn.Get_index()]=f[1]->f[6][NodeIn.Get_index()];
				}
				f[0]->f[1][NodeIn.Get_index()]=f[0]->f[Opposite[1]][NodeIn.Get_index()];
				f[0]->f[5][NodeIn.Get_index()]=f[0]->f[Opposite[5]][NodeIn.Get_index()];
				f[0]->f[8][NodeIn.Get_index()]=f[0]->f[Opposite[8]][NodeIn.Get_index()];

				f[1]->f[1][NodeIn.Get_index()]=f[1]->f[Opposite[1]][NodeIn.Get_index()];
				f[1]->f[5][NodeIn.Get_index()]=f[1]->f[Opposite[5]][NodeIn.Get_index()];
				f[1]->f[8][NodeIn.Get_index()]=f[1]->f[Opposite[8]][NodeIn.Get_index()];
				break;
			case 3:
				//apply Symmetry
				if(NodeArrays->TypeOfNode[NodeIn.Get_connect(0)]==Solid)//Bottom
				{
					f[0]->f[2][NodeIn.Get_index()]=f[0]->f[Opposite[2]][NodeIn.Get_index()];
					f[0]->f[5][NodeIn.Get_index()]=f[0]->f[8][NodeIn.Get_index()];

					f[1]->f[2][NodeIn.Get_index()]=f[1]->f[Opposite[2]][NodeIn.Get_index()];
					f[1]->f[5][NodeIn.Get_index()]=f[1]->f[8][NodeIn.Get_index()];
				}
				else //Top
				{
					f[0]->f[4][NodeIn.Get_index()]=f[0]->f[Opposite[4]][NodeIn.Get_index()];
					f[0]->f[8][NodeIn.Get_index()]=f[0]->f[5][NodeIn.Get_index()];

					f[1]->f[4][NodeIn.Get_index()]=f[1]->f[Opposite[4]][NodeIn.Get_index()];
					f[1]->f[8][NodeIn.Get_index()]=f[1]->f[5][NodeIn.Get_index()];
				}
				f[0]->f[3][NodeIn.Get_index()]=f[0]->f[Opposite[3]][NodeIn.Get_index()];
				f[0]->f[6][NodeIn.Get_index()]=f[0]->f[Opposite[6]][NodeIn.Get_index()];
				f[0]->f[7][NodeIn.Get_index()]=f[0]->f[Opposite[7]][NodeIn.Get_index()];

				f[1]->f[3][NodeIn.Get_index()]=f[1]->f[Opposite[3]][NodeIn.Get_index()];
				f[1]->f[6][NodeIn.Get_index()]=f[1]->f[Opposite[6]][NodeIn.Get_index()];
				f[1]->f[7][NodeIn.Get_index()]=f[1]->f[Opposite[7]][NodeIn.Get_index()];
				break;
			default:
				std::cerr<<"Direction wall bounce back not found. Index: "<<NodeIn.Get_index()<<std::endl;
				break;
			}
}

///Diffuse Wall treatment
void D2Q9TwoPhases::ApplyDiffuseWall(NodeWall2D& NodeIn){

	switch(NodeIn.Get_BcNormal())
			{
			case 2:
				rhodiff=(f[0]->f[4][NodeIn.Get_index()]+f[0]->f[7][NodeIn.Get_index()]+f[0]->f[8][NodeIn.Get_index()])/SumWeightS;
				f[0]->f[2][NodeIn.Get_index()]=omega[2]*rhodiff;
				f[0]->f[5][NodeIn.Get_index()]=omega[5]*rhodiff;
				f[0]->f[6][NodeIn.Get_index()]=omega[6]*rhodiff;

				rhodiff=(f[1]->f[4][NodeIn.Get_index()]+f[1]->f[7][NodeIn.Get_index()]+f[1]->f[8][NodeIn.Get_index()])/SumWeightS;
				f[1]->f[2][NodeIn.Get_index()]=omega[2]*rhodiff;
				f[1]->f[5][NodeIn.Get_index()]=omega[5]*rhodiff;
				f[1]->f[6][NodeIn.Get_index()]=omega[6]*rhodiff;
				break;
			case 4:
				rhodiff=(f[0]->f[2][NodeIn.Get_index()]+f[0]->f[5][NodeIn.Get_index()]+f[0]->f[6][NodeIn.Get_index()])/SumWeightN;
				f[0]->f[4][NodeIn.Get_index()]=omega[4]*rhodiff;
				f[0]->f[7][NodeIn.Get_index()]=omega[7]*rhodiff;
				f[0]->f[8][NodeIn.Get_index()]=omega[8]*rhodiff;

				rhodiff=(f[1]->f[2][NodeIn.Get_index()]+f[1]->f[5][NodeIn.Get_index()]+f[1]->f[6][NodeIn.Get_index()])/SumWeightN;
				f[1]->f[4][NodeIn.Get_index()]=omega[4]*rhodiff;
				f[1]->f[7][NodeIn.Get_index()]=omega[7]*rhodiff;
				f[1]->f[8][NodeIn.Get_index()]=omega[8]*rhodiff;
				break;
			case 1:
				rhodiff=(f[0]->f[3][NodeIn.Get_index()]+f[0]->f[6][NodeIn.Get_index()]+f[0]->f[7][NodeIn.Get_index()])/SumWeightW;
				f[0]->f[1][NodeIn.Get_index()]=omega[1]*rhodiff;
				f[0]->f[5][NodeIn.Get_index()]=omega[5]*rhodiff;
				f[0]->f[8][NodeIn.Get_index()]=omega[8]*rhodiff;

				rhodiff=(f[1]->f[3][NodeIn.Get_index()]+f[1]->f[6][NodeIn.Get_index()]+f[1]->f[7][NodeIn.Get_index()])/SumWeightW;
				f[1]->f[1][NodeIn.Get_index()]=omega[1]*rhodiff;
				f[1]->f[5][NodeIn.Get_index()]=omega[5]*rhodiff;
				f[1]->f[8][NodeIn.Get_index()]=omega[8]*rhodiff;
				break;
			case 3:
				rhodiff=(f[0]->f[1][NodeIn.Get_index()]+f[0]->f[5][NodeIn.Get_index()]+f[0]->f[8][NodeIn.Get_index()])/SumWeightE;
				f[0]->f[3][NodeIn.Get_index()]=omega[3]*rhodiff;
				f[0]->f[6][NodeIn.Get_index()]=omega[6]*rhodiff;
				f[0]->f[7][NodeIn.Get_index()]=omega[7]*rhodiff;

				rhodiff=(f[1]->f[1][NodeIn.Get_index()]+f[1]->f[5][NodeIn.Get_index()]+f[1]->f[8][NodeIn.Get_index()])/SumWeightE;
				f[1]->f[3][NodeIn.Get_index()]=omega[3]*rhodiff;
				f[1]->f[6][NodeIn.Get_index()]=omega[6]*rhodiff;
				f[1]->f[7][NodeIn.Get_index()]=omega[7]*rhodiff;
				break;
			default:
				std::cerr<<"Direction: "<< NodeIn.Get_BcNormal()<<" (Wall diffuse boundary condition) not found"<<" x: "<<NodeIn.get_x()<<" y: "<<NodeIn.get_y()<<" processor: "<<parallel->getRank()<<std::endl;
				break;
			}
}
void D2Q9TwoPhases::ApplyDiffuseWallSymmetry(NodeWall2D& NodeIn){

	switch(NodeIn.Get_BcNormal())
			{
			case 2:
				rhodiff=(f[0]->f[4][NodeIn.Get_index()]+f[0]->f[7][NodeIn.Get_index()]+f[0]->f[8][NodeIn.Get_index()])/SumWeightS;
				f[0]->f[2][NodeIn.Get_index()]=omega[2]*rhodiff;
				f[0]->f[5][NodeIn.Get_index()]=omega[5]*rhodiff;
				f[0]->f[6][NodeIn.Get_index()]=omega[6]*rhodiff;

				rhodiff=(f[1]->f[4][NodeIn.Get_index()]+f[1]->f[7][NodeIn.Get_index()]+f[1]->f[8][NodeIn.Get_index()])/SumWeightS;
				f[1]->f[2][NodeIn.Get_index()]=omega[2]*rhodiff;
				f[1]->f[5][NodeIn.Get_index()]=omega[5]*rhodiff;
				f[1]->f[6][NodeIn.Get_index()]=omega[6]*rhodiff;
				break;
			case 4:
				rhodiff=(f[0]->f[2][NodeIn.Get_index()]+f[0]->f[5][NodeIn.Get_index()]+f[0]->f[6][NodeIn.Get_index()])/SumWeightN;
				f[0]->f[4][NodeIn.Get_index()]=omega[4]*rhodiff;
				f[0]->f[7][NodeIn.Get_index()]=omega[7]*rhodiff;
				f[0]->f[8][NodeIn.Get_index()]=omega[8]*rhodiff;

				rhodiff=(f[1]->f[2][NodeIn.Get_index()]+f[1]->f[5][NodeIn.Get_index()]+f[1]->f[6][NodeIn.Get_index()])/SumWeightN;
				f[1]->f[4][NodeIn.Get_index()]=omega[4]*rhodiff;
				f[1]->f[7][NodeIn.Get_index()]=omega[7]*rhodiff;
				f[1]->f[8][NodeIn.Get_index()]=omega[8]*rhodiff;
				break;
			case 1:
				//apply Symmetry
				if(NodeArrays->TypeOfNode[NodeIn.Get_connect(0)]==Solid)//Bottom
				{
					f[0]->f[2][NodeIn.Get_index()]=f[0]->f[Opposite[2]][NodeIn.Get_index()];
					f[0]->f[6][NodeIn.Get_index()]=f[0]->f[7][NodeIn.Get_index()];

					f[1]->f[2][NodeIn.Get_index()]=f[1]->f[Opposite[2]][NodeIn.Get_index()];
					f[1]->f[6][NodeIn.Get_index()]=f[1]->f[7][NodeIn.Get_index()];
				}
				else //Top
				{
					f[0]->f[4][NodeIn.Get_index()]=f[0]->f[Opposite[4]][NodeIn.Get_index()];
					f[0]->f[7][NodeIn.Get_index()]=f[0]->f[6][NodeIn.Get_index()];

					f[1]->f[4][NodeIn.Get_index()]=f[1]->f[Opposite[4]][NodeIn.Get_index()];
					f[1]->f[7][NodeIn.Get_index()]=f[1]->f[6][NodeIn.Get_index()];
				}
				rhodiff=(f[0]->f[3][NodeIn.Get_index()]+f[0]->f[6][NodeIn.Get_index()]+f[0]->f[7][NodeIn.Get_index()])/SumWeightW;
				f[0]->f[1][NodeIn.Get_index()]=omega[1]*rhodiff;
				f[0]->f[5][NodeIn.Get_index()]=omega[5]*rhodiff;
				f[0]->f[8][NodeIn.Get_index()]=omega[8]*rhodiff;

				rhodiff=(f[1]->f[3][NodeIn.Get_index()]+f[1]->f[6][NodeIn.Get_index()]+f[1]->f[7][NodeIn.Get_index()])/SumWeightW;
				f[1]->f[1][NodeIn.Get_index()]=omega[1]*rhodiff;
				f[1]->f[5][NodeIn.Get_index()]=omega[5]*rhodiff;
				f[1]->f[8][NodeIn.Get_index()]=omega[8]*rhodiff;
				break;
			case 3:
				//apply Symmetry
				if(NodeArrays->TypeOfNode[NodeIn.Get_connect(0)]==Solid)//Bottom
				{
					f[0]->f[2][NodeIn.Get_index()]=f[0]->f[Opposite[2]][NodeIn.Get_index()];
					f[0]->f[5][NodeIn.Get_index()]=f[0]->f[8][NodeIn.Get_index()];

					f[1]->f[2][NodeIn.Get_index()]=f[1]->f[Opposite[2]][NodeIn.Get_index()];
					f[1]->f[5][NodeIn.Get_index()]=f[1]->f[8][NodeIn.Get_index()];
				}
				else //Top
				{
					f[0]->f[4][NodeIn.Get_index()]=f[0]->f[Opposite[4]][NodeIn.Get_index()];
					f[0]->f[8][NodeIn.Get_index()]=f[0]->f[5][NodeIn.Get_index()];

					f[1]->f[4][NodeIn.Get_index()]=f[1]->f[Opposite[4]][NodeIn.Get_index()];
					f[1]->f[8][NodeIn.Get_index()]=f[1]->f[5][NodeIn.Get_index()];
				}
				rhodiff=(f[0]->f[1][NodeIn.Get_index()]+f[0]->f[5][NodeIn.Get_index()]+f[0]->f[8][NodeIn.Get_index()])/SumWeightE;
				f[0]->f[3][NodeIn.Get_index()]=omega[3]*rhodiff;
				f[0]->f[6][NodeIn.Get_index()]=omega[6]*rhodiff;
				f[0]->f[7][NodeIn.Get_index()]=omega[7]*rhodiff;

				rhodiff=(f[1]->f[1][NodeIn.Get_index()]+f[1]->f[5][NodeIn.Get_index()]+f[1]->f[8][NodeIn.Get_index()])/SumWeightE;
				f[1]->f[3][NodeIn.Get_index()]=omega[3]*rhodiff;
				f[1]->f[6][NodeIn.Get_index()]=omega[6]*rhodiff;
				f[1]->f[7][NodeIn.Get_index()]=omega[7]*rhodiff;
				break;
			default:
				std::cerr<<"Direction: "<< NodeIn.Get_BcNormal()<<" (Wall diffuse boundary condition) not found"<<std::endl;
				break;
			}
}

///Diffuse boundary condition for Concave and Convex corners
///For Concave corners, the incoming distributions go to outcoming so divide by 0.5 and the diagonal distribution doesn't participate for the streaming is approximate by the density at the neibourgh nodes
///For Convex corners, the diffusion goes in all directions.
void D2Q9TwoPhases::ApplyDiffuseWall(NodeCorner2D& NodeIn){
	double tmp=0;
	switch(NodeIn.Get_BcNormal())
			{
			case 5:
					rhodiff=(f[0]->f[3][NodeIn.Get_index()]+f[0]->f[4][NodeIn.Get_index()]+f[0]->f[7][NodeIn.Get_index()])/SumWeightConvexNE;
					f[0]->f[5][NodeIn.Get_index()]=omega[5]*rhodiff;
					f[0]->f[1][NodeIn.Get_index()]=omega[1]*rhodiff;
					f[0]->f[2][NodeIn.Get_index()]=omega[2]*rhodiff;
					f[0]->f[6][NodeIn.Get_index()]=omega[6]*rhodiff;
					f[0]->f[8][NodeIn.Get_index()]=omega[8]*rhodiff;

					rhodiff=(f[1]->f[3][NodeIn.Get_index()]+f[1]->f[4][NodeIn.Get_index()]+f[1]->f[7][NodeIn.Get_index()])/SumWeightConvexNE;
					f[1]->f[5][NodeIn.Get_index()]=omega[5]*rhodiff;
					f[1]->f[1][NodeIn.Get_index()]=omega[1]*rhodiff;
					f[1]->f[2][NodeIn.Get_index()]=omega[2]*rhodiff;
					f[1]->f[6][NodeIn.Get_index()]=omega[6]*rhodiff;
					f[1]->f[8][NodeIn.Get_index()]=omega[8]*rhodiff;

				break;
			case 6:

					rhodiff=(f[0]->f[1][NodeIn.Get_index()]+f[0]->f[4][NodeIn.Get_index()]+f[0]->f[8][NodeIn.Get_index()])/SumWeightConvexNW;
					f[0]->f[6][NodeIn.Get_index()]=omega[6]*rhodiff;
					f[0]->f[2][NodeIn.Get_index()]=omega[2]*rhodiff;
					f[0]->f[3][NodeIn.Get_index()]=omega[3]*rhodiff;
					f[0]->f[5][NodeIn.Get_index()]=omega[5]*rhodiff;
					f[0]->f[7][NodeIn.Get_index()]=omega[7]*rhodiff;

					rhodiff=(f[1]->f[1][NodeIn.Get_index()]+f[1]->f[4][NodeIn.Get_index()]+f[1]->f[8][NodeIn.Get_index()])/SumWeightConvexNW;
					f[1]->f[6][NodeIn.Get_index()]=omega[6]*rhodiff;
					f[1]->f[2][NodeIn.Get_index()]=omega[2]*rhodiff;
					f[1]->f[3][NodeIn.Get_index()]=omega[3]*rhodiff;
					f[1]->f[5][NodeIn.Get_index()]=omega[5]*rhodiff;
					f[1]->f[7][NodeIn.Get_index()]=omega[7]*rhodiff;

				break;
			case 7:

					rhodiff=(f[0]->f[1][NodeIn.Get_index()]+f[0]->f[2][NodeIn.Get_index()]+f[0]->f[5][NodeIn.Get_index()])/SumWeightConvexSW;
					f[0]->f[7][NodeIn.Get_index()]=omega[7]*rhodiff;
					f[0]->f[3][NodeIn.Get_index()]=omega[3]*rhodiff;
					f[0]->f[4][NodeIn.Get_index()]=omega[4]*rhodiff;
					f[0]->f[6][NodeIn.Get_index()]=omega[6]*rhodiff;
					f[0]->f[8][NodeIn.Get_index()]=omega[8]*rhodiff;

					rhodiff=(f[1]->f[1][NodeIn.Get_index()]+f[1]->f[2][NodeIn.Get_index()]+f[1]->f[5][NodeIn.Get_index()])/SumWeightConvexSW;
					f[1]->f[7][NodeIn.Get_index()]=omega[7]*rhodiff;
					f[1]->f[3][NodeIn.Get_index()]=omega[3]*rhodiff;
					f[1]->f[4][NodeIn.Get_index()]=omega[4]*rhodiff;
					f[1]->f[6][NodeIn.Get_index()]=omega[6]*rhodiff;
					f[1]->f[8][NodeIn.Get_index()]=omega[8]*rhodiff;
				//}
				break;
			case 8:

					rhodiff=(f[0]->f[2][NodeIn.Get_index()]+f[0]->f[3][NodeIn.Get_index()]+f[0]->f[6][NodeIn.Get_index()])/SumWeightConvexSE;
					f[0]->f[8][NodeIn.Get_index()]=omega[8]*rhodiff;
					f[0]->f[1][NodeIn.Get_index()]=omega[1]*rhodiff;
					f[0]->f[4][NodeIn.Get_index()]=omega[4]*rhodiff;
					f[0]->f[5][NodeIn.Get_index()]=omega[5]*rhodiff;
					f[0]->f[7][NodeIn.Get_index()]=omega[7]*rhodiff;

					rhodiff=(f[1]->f[2][NodeIn.Get_index()]+f[1]->f[3][NodeIn.Get_index()]+f[1]->f[6][NodeIn.Get_index()])/SumWeightConvexSE;
					f[1]->f[8][NodeIn.Get_index()]=omega[8]*rhodiff;
					f[1]->f[1][NodeIn.Get_index()]=omega[1]*rhodiff;
					f[1]->f[4][NodeIn.Get_index()]=omega[4]*rhodiff;
					f[1]->f[5][NodeIn.Get_index()]=omega[5]*rhodiff;
					f[1]->f[7][NodeIn.Get_index()]=omega[7]*rhodiff;

				break;
			default:
				std::cerr<<"Direction: "<< NodeIn.Get_BcNormal()<<" (Corner diffuse boundary conditions) not found"<<std::endl;
				break;
			}
}

*/

double D2Q9TwoPhases::Cal_RhoCorner(NodeCorner2D& nodeIn){

	unsigned int direction1,direction2;
	switch(nodeIn.Get_BcNormal())
	{
	case 5:
		direction1=1;
		direction2=2;
		doubleTmpReturn=(Rho[nodeIn.Get_connect()[direction1]]+Rho[nodeIn.Get_connect()[direction2]])*0.5;
		break;
	case 6:
		direction1=2;
		direction2=3;
		doubleTmpReturn=(Rho[nodeIn.Get_connect()[direction1]]+Rho[nodeIn.Get_connect()[direction2]])*0.5;
		break;
	case 7:
		direction1=3;
		direction2=4;
		doubleTmpReturn=(Rho[nodeIn.Get_connect()[direction1]]+Rho[nodeIn.Get_connect()[direction2]])*0.5;
		break;
	case 8:
		direction1=1;
		direction2=4;
		doubleTmpReturn=(Rho[nodeIn.Get_connect()[direction1]]+Rho[nodeIn.Get_connect()[direction2]])*0.5;
		break;
	}

	return doubleTmpReturn;
}

void D2Q9TwoPhases::StreamingOrientation(NodeGhost2D& nodeIn, bool Streaming[9]){
	Streaming[0]=false;
	for (unsigned int i=1;i<(unsigned int)nbvelo;i++)
		if((nodeIn.Get_connect()[i]==nodeIn.Get_index())||(nodeIn.get_NodeType()==Solid)||(NodeArrays->TypeOfNode[nodeIn.Get_connect()[i]]==Solid)||(NodeArrays->TypeOfNode[nodeIn.Get_connect()[i]]==Ghost))
			Streaming[i]=false;
		else
			Streaming[i]=true;
}
void D2Q9TwoPhases::StreamingOrientation(NodeCorner2D& nodeIn, bool Streaming[9]){

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
void D2Q9TwoPhases::StreamingOrientation(NodeWall2D& nodeIn, bool Streaming[9]){
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
/*	if(nodeIn.get_x()==206 && (nodeIn.get_y()<182||nodeIn.get_y()<187))
		std::cout<<"processor: "<<parallel->getRank()<<" X: "<<nodeIn.get_x()<<" y: "<<nodeIn.get_y()<<" normal: "<<nodeIn.Get_BcNormal()<<std::endl;
	sleep(0.5);*/

}
void D2Q9TwoPhases::StreamingOrientation(NodeSymmetry2D& nodeIn, bool Streaming[9]){
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
void D2Q9TwoPhases::StreamingOrientation(NodePeriodic2D& nodeIn, bool Streaming[9]){
	Streaming[0]=false;
	for (unsigned int i=1;i<(unsigned int)nbvelo;i++)
		if((nodeIn.Get_connect()[i]==nodeIn.Get_index())||(nodeIn.get_NodeType()==Solid)||(NodeArrays->TypeOfNode[nodeIn.Get_connect()[i]]==Solid))
			Streaming[i]=false;
		else
			Streaming[i]=true;
}
void D2Q9TwoPhases::StreamingOrientation(NodeVelocity2D& nodeIn, bool Streaming[9]){
	Streaming[0]=false;
	for (unsigned int i=1;i<(unsigned int)nbvelo;i++)
		if((nodeIn.Get_connect()[i]==nodeIn.Get_index())||(nodeIn.get_NodeType()==Solid)||(NodeArrays->TypeOfNode[nodeIn.Get_connect()[i]]==Solid))
			Streaming[i]=false;
		else
			Streaming[i]=true;
}
void D2Q9TwoPhases::StreamingOrientation(NodePressure2D& nodeIn, bool Streaming[9]){
	Streaming[0]=false;
	for (unsigned int i=1;i<(unsigned int)nbvelo;i++)
		if((nodeIn.Get_connect()[i]==nodeIn.Get_index())||(nodeIn.get_NodeType()==Solid)||(NodeArrays->TypeOfNode[nodeIn.Get_connect()[i]]==Solid))
			Streaming[i]=false;
		else
			Streaming[i]=true;
}




//Communication Functions
void D2Q9TwoPhases::IniComVariables(){
	MultiBlock_->Get_Connect_Node(IdNodeN,IdNodeE,IdNodeS,IdNodeW,IdNodeSW,IdNodeSE,IdNodeNW,IdNodeNE);
	MultiBlock_->Get_Connect_Node(IdRNodeN,IdRNodeE,IdRNodeS,IdRNodeW,IdGNodeN,IdGNodeE,IdGNodeS,IdGNodeW,
			IdRNodeSW,IdRNodeSE,IdRNodeNW,IdRNodeNE,IdGNodeSW,IdGNodeSE,IdGNodeNW,IdGNodeNE);
	MultiBlock_->Get_Connect_SolidNode(SolidIdRNodeN,SolidIdRNodeE,SolidIdRNodeS,SolidIdRNodeW,SolidIdGNodeN,SolidIdGNodeE,SolidIdGNodeS,SolidIdGNodeW,
			SolidIdRNodeSW,SolidIdRNodeSE,SolidIdRNodeNW,SolidIdRNodeNE,SolidIdGNodeSW,SolidIdGNodeSE,SolidIdGNodeNW,SolidIdGNodeNE);
	Nd_variables_sync=3*2;
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
	size_buf[0]=IdGNodeE.size();
	size_buf[1]=IdGNodeW.size();
	size_buf[2]=IdGNodeS.size();
	size_buf[3]=IdGNodeN.size();
// Macro sync
	Nd_MacroVariables_sync=Dic->Get_NbSyncVar();//6;
	if(parallel->isMainProcessor())
	{
		std::cout<<"Synchromisation ; number of variable: "<<Nd_MacroVariables_sync<<std::endl<<"Synchronise Variables are: ";
		std::vector<std::string> syncname=Dic->Get_SyncVarName();
		for(int i=0;i<Nd_MacroVariables_sync;i++)
			std::cout<<syncname[i]<<" ";
		std::cout<<std::endl;
	}

	SyncVar=Dic->Get_SyncVar();
	buf_MacroSend=new double** [Nd_MacroVariables_sync];
	buf_MacroRecv=new double** [Nd_MacroVariables_sync];
	buf_MacroSendSolid=new double** [Nd_MacroVariables_sync];
	buf_MacroRecvSolid=new double** [Nd_MacroVariables_sync];
	for (int i=0;i<Nd_MacroVariables_sync;i++)
	{
		buf_MacroSend[i]=new double* [4];
		buf_MacroRecv[i]=new double* [4];
		buf_MacroSendSolid[i]=new double* [4];
		buf_MacroRecvSolid[i]=new double* [4];
	}

	size_MacroBuf=new int [8];size_MacroBufSolid=new int [8];

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

		buf_MacroSendSolid[i][0]=new double[SolidIdGNodeE.size()];
		buf_MacroSendSolid[i][1]=new double[SolidIdGNodeW.size()];
		buf_MacroSendSolid[i][2]=new double[SolidIdGNodeS.size()];
		buf_MacroSendSolid[i][3]=new double[SolidIdGNodeN.size()];
		buf_MacroRecvSolid[i][0]=new double[SolidIdRNodeW.size()];
		buf_MacroRecvSolid[i][1]=new double[SolidIdRNodeE.size()];
		buf_MacroRecvSolid[i][2]=new double[SolidIdRNodeN.size()];
		buf_MacroRecvSolid[i][3]=new double[SolidIdRNodeS.size()];
	}
	size_MacroBuf[0]=IdRNodeE.size();
	size_MacroBuf[1]=IdRNodeW.size();
	size_MacroBuf[2]=IdRNodeS.size();
	size_MacroBuf[3]=IdRNodeN.size();
	size_MacroBuf[4]=IdGNodeW.size();
	size_MacroBuf[5]=IdGNodeE.size();
	size_MacroBuf[6]=IdGNodeN.size();
	size_MacroBuf[7]=IdGNodeS.size();

	size_MacroBufSolid[0]=SolidIdGNodeE.size();
	size_MacroBufSolid[1]=SolidIdGNodeW.size();
	size_MacroBufSolid[2]=SolidIdGNodeS.size();
	size_MacroBufSolid[3]=SolidIdGNodeN.size();
	size_MacroBufSolid[4]=SolidIdRNodeW.size();
	size_MacroBufSolid[5]=SolidIdRNodeE.size();
	size_MacroBufSolid[6]=SolidIdRNodeN.size();
	size_MacroBufSolid[7]=SolidIdRNodeS.size();
}
void D2Q9TwoPhases::GhostNodesSyncFromGhost(){

	for (unsigned int i=0;i<IdNodeE.size();i++)
	{
		buf_send[0][0][i]=f[0]->f[1][IdNodeE[i]];
		buf_send[1][0][i]=f[0]->f[5][IdNodeE[i]];
		buf_send[2][0][i]=f[0]->f[8][IdNodeE[i]];

		buf_send[3][0][i]=f[1]->f[1][IdNodeE[i]];
		buf_send[4][0][i]=f[1]->f[5][IdNodeE[i]];
		buf_send[5][0][i]=f[1]->f[8][IdNodeE[i]];
	}

	for (unsigned int i=0;i<IdNodeN.size();i++)
	{
		buf_send[0][3][i]=f[0]->f[2][IdNodeN[i]];
		buf_send[1][3][i]=f[0]->f[5][IdNodeN[i]];
		buf_send[2][3][i]=f[0]->f[6][IdNodeN[i]];

		buf_send[3][3][i]=f[1]->f[2][IdNodeN[i]];
		buf_send[4][3][i]=f[1]->f[5][IdNodeN[i]];
		buf_send[5][3][i]=f[1]->f[6][IdNodeN[i]];
	}

//	int itmp=0;
	for (int i=0;i<Nd_variables_sync;i++)
	{
		MultiBlock_->CommunicationFromGhost(buf_send[i],buf_recv[i],size_buf);
		//MultiBlock_->CommunicationFromGhost(buf_send[i],buf_recv[i],size_buf,&status[itmp],&request[itmp]);
		//itmp+=2;
	}
	for (unsigned int i=0;i<IdNodeW.size();i++)
	{
		f[0]->f[1][IdNodeW[i]]=buf_recv[0][0][i];
		f[0]->f[5][IdNodeW[i]]=buf_recv[1][0][i];
		f[0]->f[8][IdNodeW[i]]=buf_recv[2][0][i];
	}

	for (unsigned int i=0;i<IdNodeS.size();i++)
	{
		f[0]->f[2][IdNodeS[i]]=buf_recv[0][3][i];
		f[0]->f[5][IdNodeS[i]]=buf_recv[1][3][i];
		f[0]->f[6][IdNodeS[i]]=buf_recv[2][3][i];
	}

}
void D2Q9TwoPhases::GhostNodesSyncToGhost(){
/*	for (unsigned int i=0;i<IdNodeW.size();i++)
	{
		buf_send[0][1][i]=f[0]->f[3][IdNodeW[i]];
		buf_send[1][1][i]=f[0]->f[6][IdNodeW[i]];
		buf_send[2][1][i]=f[0]->f[7][IdNodeW[i]];
	}
	for (unsigned int i=0;i<IdNodeS.size();i++)
	{
		buf_send[0][2][i]=f[0]->f[4][IdNodeS[i]];
		buf_send[1][2][i]=f[0]->f[7][IdNodeS[i]];
		buf_send[2][2][i]=f[0]->f[8][IdNodeS[i]];
	}
	//MPI_Barrier(MPI_COMM_WORLD);
	for (int i=0;i<Nd_variables_sync;i++)
	{
		//MultiBlock_->Communication(buf_send[i],buf_recv[i],size_buf);
		MultiBlock_->CommunicationToGhost(buf_send[i],buf_recv[i],size_buf);

	}
	for (unsigned int i=0;i<IdNodeE.size();i++)
	{
		f[0]->f[3][IdNodeE[i]]=buf_recv[0][1][i];
		f[0]->f[6][IdNodeE[i]]=buf_recv[1][1][i];
		f[0]->f[7][IdNodeE[i]]=buf_recv[2][1][i];
	}
	for (unsigned int i=0;i<IdNodeN.size();i++)
	{
		f[0]->f[4][IdNodeN[i]]=buf_recv[0][2][i];
		f[0]->f[7][IdNodeN[i]]=buf_recv[1][2][i];
		f[0]->f[8][IdNodeN[i]]=buf_recv[2][2][i];
	}*/

	for (unsigned int i=0;i<IdRNodeW.size();i++)
	{
		buf_send[0][1][i]=f[0]->f[3][IdRNodeW[i]];
		buf_send[1][1][i]=f[0]->f[6][IdRNodeW[i]];
		buf_send[2][1][i]=f[0]->f[7][IdRNodeW[i]];

		buf_send[3][1][i]=f[1]->f[3][IdRNodeW[i]];
		buf_send[4][1][i]=f[1]->f[6][IdRNodeW[i]];
		buf_send[5][1][i]=f[1]->f[7][IdRNodeW[i]];

	}
	for (unsigned int i=0;i<IdRNodeS.size();i++)
	{
		buf_send[0][2][i]=f[0]->f[4][IdRNodeS[i]];
		buf_send[1][2][i]=f[0]->f[7][IdRNodeS[i]];
		buf_send[2][2][i]=f[0]->f[8][IdRNodeS[i]];

		buf_send[3][2][i]=f[1]->f[4][IdRNodeS[i]];
		buf_send[4][2][i]=f[1]->f[7][IdRNodeS[i]];
		buf_send[5][2][i]=f[1]->f[8][IdRNodeS[i]];
	}
	for (unsigned int i=0;i<IdRNodeE.size();i++)
	{
		buf_send[0][0][i]=f[0]->f[1][IdRNodeE[i]];
		buf_send[1][0][i]=f[0]->f[5][IdRNodeE[i]];
		buf_send[2][0][i]=f[0]->f[8][IdRNodeE[i]];

		buf_send[3][0][i]=f[1]->f[1][IdRNodeE[i]];
		buf_send[4][0][i]=f[1]->f[5][IdRNodeE[i]];
		buf_send[5][0][i]=f[1]->f[8][IdRNodeE[i]];
	}
	for (unsigned int i=0;i<IdRNodeN.size();i++)
	{
		buf_send[0][3][i]=f[0]->f[2][IdRNodeN[i]];
		buf_send[1][3][i]=f[0]->f[5][IdRNodeN[i]];
		buf_send[2][3][i]=f[0]->f[6][IdRNodeN[i]];

		buf_send[3][3][i]=f[1]->f[2][IdRNodeN[i]];
		buf_send[4][3][i]=f[1]->f[5][IdRNodeN[i]];
		buf_send[5][3][i]=f[1]->f[6][IdRNodeN[i]];
	}
/*	for (int i=0;i<Nd_variables_sync;i++)
	{
		MultiBlock_->Communication(buf_MacroSend[i],buf_MacroRecv[i],size_buf);
		//MultiBlock_->CommunicationToGhost(buf_MacroSend[i],buf_MacroRecv[i],size_MacroBuf);

	}*/
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

		MultiBlock_->Send(&buf_send[3][0][0],IdRNodeE.size(),1,tag_x_r);
		MultiBlock_->Send(&buf_send[4][0][0],IdRNodeE.size(),1,tag_x_r);
		MultiBlock_->Send(&buf_send[5][0][0],IdRNodeE.size(),1,tag_x_r);
	}
	if(IdGNodeW.size()>=1)
	{
		MultiBlock_->Recv(&buf_recv[0][0][0],IdGNodeW.size(),3,tag_x_r,status);
		MultiBlock_->Recv(&buf_recv[1][0][0],IdGNodeW.size(),3,tag_x_r,status);
		MultiBlock_->Recv(&buf_recv[2][0][0],IdGNodeW.size(),3,tag_x_r,status);

		MultiBlock_->Recv(&buf_recv[3][0][0],IdGNodeW.size(),3,tag_x_r,status);
		MultiBlock_->Recv(&buf_recv[4][0][0],IdGNodeW.size(),3,tag_x_r,status);
		MultiBlock_->Recv(&buf_recv[5][0][0],IdGNodeW.size(),3,tag_x_r,status);
	}

	if(IdRNodeW.size()>=1)
	{
		MultiBlock_->Send(&buf_send[0][1][0],IdRNodeW.size(),3,tag_x_l);
		MultiBlock_->Send(&buf_send[1][1][0],IdRNodeW.size(),3,tag_x_l);
		MultiBlock_->Send(&buf_send[2][1][0],IdRNodeW.size(),3,tag_x_l);

		MultiBlock_->Send(&buf_send[3][1][0],IdRNodeW.size(),3,tag_x_l);
		MultiBlock_->Send(&buf_send[4][1][0],IdRNodeW.size(),3,tag_x_l);
		MultiBlock_->Send(&buf_send[5][1][0],IdRNodeW.size(),3,tag_x_l);
	}
	if(IdGNodeE.size()>=1)
	{
		MultiBlock_->Recv(&buf_recv[0][1][0],IdGNodeE.size(),1,tag_x_l,status);
		MultiBlock_->Recv(&buf_recv[1][1][0],IdGNodeE.size(),1,tag_x_l,status);
		MultiBlock_->Recv(&buf_recv[2][1][0],IdGNodeE.size(),1,tag_x_l,status);

		MultiBlock_->Recv(&buf_recv[3][1][0],IdGNodeE.size(),1,tag_x_l,status);
		MultiBlock_->Recv(&buf_recv[4][1][0],IdGNodeE.size(),1,tag_x_l,status);
		MultiBlock_->Recv(&buf_recv[5][1][0],IdGNodeE.size(),1,tag_x_l,status);
	}
		if(IdRNodeN.size()>=1)
	{
		MultiBlock_->Send(&buf_send[0][3][0],IdRNodeN.size(),0,tag_y_t);
		MultiBlock_->Send(&buf_send[1][3][0],IdRNodeN.size(),0,tag_y_t);
		MultiBlock_->Send(&buf_send[2][3][0],IdRNodeN.size(),0,tag_y_t);

		MultiBlock_->Send(&buf_send[3][3][0],IdRNodeN.size(),0,tag_y_t);
		MultiBlock_->Send(&buf_send[4][3][0],IdRNodeN.size(),0,tag_y_t);
		MultiBlock_->Send(&buf_send[5][3][0],IdRNodeN.size(),0,tag_y_t);
	}
	if(IdGNodeS.size()>=1)
	{
		MultiBlock_->Recv(&buf_recv[0][3][0],IdGNodeS.size(),2,tag_y_t,status);
		MultiBlock_->Recv(&buf_recv[1][3][0],IdGNodeS.size(),2,tag_y_t,status);
		MultiBlock_->Recv(&buf_recv[2][3][0],IdGNodeS.size(),2,tag_y_t,status);

		MultiBlock_->Recv(&buf_recv[3][3][0],IdGNodeS.size(),2,tag_y_t,status);
		MultiBlock_->Recv(&buf_recv[4][3][0],IdGNodeS.size(),2,tag_y_t,status);
		MultiBlock_->Recv(&buf_recv[5][3][0],IdGNodeS.size(),2,tag_y_t,status);
	}
	if(IdRNodeS.size()>=1)
	{
		MultiBlock_->Send(&buf_send[0][2][0],IdRNodeS.size(),2,tag_y_b);
		MultiBlock_->Send(&buf_send[1][2][0],IdRNodeS.size(),2,tag_y_b);
		MultiBlock_->Send(&buf_send[2][2][0],IdRNodeS.size(),2,tag_y_b);

		MultiBlock_->Send(&buf_send[3][2][0],IdRNodeS.size(),2,tag_y_b);
		MultiBlock_->Send(&buf_send[4][2][0],IdRNodeS.size(),2,tag_y_b);
		MultiBlock_->Send(&buf_send[5][2][0],IdRNodeS.size(),2,tag_y_b);
	}
	if(IdGNodeN.size()>=1)
	{
		MultiBlock_->Recv(&buf_recv[0][2][0],IdGNodeN.size(),0,tag_y_b,status);
		MultiBlock_->Recv(&buf_recv[1][2][0],IdGNodeN.size(),0,tag_y_b,status);
		MultiBlock_->Recv(&buf_recv[2][2][0],IdGNodeN.size(),0,tag_y_b,status);

		MultiBlock_->Recv(&buf_recv[3][2][0],IdGNodeN.size(),0,tag_y_b,status);
		MultiBlock_->Recv(&buf_recv[4][2][0],IdGNodeN.size(),0,tag_y_b,status);
		MultiBlock_->Recv(&buf_recv[5][2][0],IdGNodeN.size(),0,tag_y_b,status);
	}
	//MultiBlock_->Communication(buf_MacroSend[0],buf_MacroRecv[0],size_buf);
	//int* block=MultiBlock_->get_Block_Connect();
	//if(block[1]>=0)
	for (unsigned int i=0;i<IdGNodeE.size();i++)
	{
		f[0]->f[3][IdGNodeE[i]]=buf_recv[0][1][i];
		f[0]->f[6][IdGNodeE[i]]=buf_recv[1][1][i];
		f[0]->f[7][IdGNodeE[i]]=buf_recv[2][1][i];

		f[1]->f[3][IdGNodeE[i]]=buf_recv[3][1][i];
		f[1]->f[6][IdGNodeE[i]]=buf_recv[4][1][i];
		f[1]->f[7][IdGNodeE[i]]=buf_recv[5][1][i];
	}
	//if(block[0]>=0)
	for (unsigned int i=0;i<IdGNodeN.size();i++)
	{
		f[0]->f[4][IdGNodeN[i]]=buf_recv[0][2][i];
		f[0]->f[7][IdGNodeN[i]]=buf_recv[1][2][i];
		f[0]->f[8][IdGNodeN[i]]=buf_recv[2][2][i];

		f[1]->f[4][IdGNodeN[i]]=buf_recv[3][2][i];
		f[1]->f[7][IdGNodeN[i]]=buf_recv[4][2][i];
		f[1]->f[8][IdGNodeN[i]]=buf_recv[5][2][i];
	}
	//if(block[3]>=0)
	for (unsigned int i=0;i<IdGNodeW.size();i++)
	{
		f[0]->f[1][IdGNodeW[i]]=buf_recv[0][0][i];
		f[0]->f[5][IdGNodeW[i]]=buf_recv[1][0][i];
		f[0]->f[8][IdGNodeW[i]]=buf_recv[2][0][i];

		f[1]->f[1][IdGNodeW[i]]=buf_recv[3][0][i];
		f[1]->f[5][IdGNodeW[i]]=buf_recv[4][0][i];
		f[1]->f[8][IdGNodeW[i]]=buf_recv[5][0][i];
	}
	//if(block[2]>=0)
	for (unsigned int i=0;i<IdGNodeS.size();i++)
	{
		f[0]->f[2][IdGNodeS[i]]=buf_recv[0][3][i];
		f[0]->f[5][IdGNodeS[i]]=buf_recv[1][3][i];
		f[0]->f[6][IdGNodeS[i]]=buf_recv[2][3][i];

		f[1]->f[2][IdGNodeS[i]]=buf_recv[3][3][i];
		f[1]->f[5][IdGNodeS[i]]=buf_recv[4][3][i];
		f[1]->f[6][IdGNodeS[i]]=buf_recv[5][3][i];
	}
}
void D2Q9TwoPhases::CornerNodesSyncFromGhost(){
	 int tag_x_r=1;
	 int tag_y_t=2;
	 int tag_d_tr=3;

	//N=0,E=1,S=2,W=3,NE=4, SE=5, NW=6, SW=7;
	MPI_Status status;
	if(IdNodeNE.size()>=1)
	{
		MultiBlock_->Send(&f[0]->f[5][IdNodeNE[0]],1,4,tag_d_tr);
		MultiBlock_->Send(&f[1]->f[5][IdNodeNE[0]],1,4,tag_d_tr);
	}
	if(IdNodeSW.size()>=1)
	{
		MultiBlock_->Recv(&f[0]->f[5][IdNodeSW[0]],1,7,tag_d_tr,status);
		MultiBlock_->Recv(&f[1]->f[5][IdNodeSW[0]],1,7,tag_d_tr,status);
	}
	if(IdNodeNW.size()>=1)
	{
		MultiBlock_->Send(&f[0]->f[2][IdNodeNW[0]],1,0,tag_y_t);
		MultiBlock_->Send(&f[1]->f[2][IdNodeNW[0]],1,0,tag_y_t);
	}
	if(IdNodeSW.size()>=1)
	{
		MultiBlock_->Recv(&f[0]->f[2][IdNodeSW[0]],1,2,tag_y_t,status);
		MultiBlock_->Recv(&f[1]->f[2][IdNodeSW[0]],1,2,tag_y_t,status);
	}
	if(IdNodeNW.size()>=1)
	{
		MultiBlock_->Send(&f[0]->f[6][IdNodeNW[0]],1,0,tag_y_t);
		MultiBlock_->Send(&f[1]->f[6][IdNodeNW[0]],1,0,tag_y_t);
	}
	if(IdNodeSW.size()>=1)
	{
		MultiBlock_->Recv(&f[0]->f[6][IdNodeSW[0]],1,2,tag_y_t,status);
		MultiBlock_->Recv(&f[1]->f[6][IdNodeSW[0]],1,2,tag_y_t,status);
	}
	if(IdNodeSE.size()>=1)
	{
		MultiBlock_->Send(&f[0]->f[1][IdNodeSE[0]],1,1,tag_x_r);
		MultiBlock_->Send(&f[1]->f[1][IdNodeSE[0]],1,1,tag_x_r);
	}
	if(IdNodeSW.size()>=1)
	{
		MultiBlock_->Recv(&f[0]->f[1][IdNodeSW[0]],1,3,tag_x_r,status);
		MultiBlock_->Recv(&f[1]->f[1][IdNodeSW[0]],1,3,tag_x_r,status);
	}
	if(IdNodeSE.size()>=1)
	{
		MultiBlock_->Send(&f[0]->f[8][IdNodeSE[0]],1,1,tag_x_r);
		MultiBlock_->Send(&f[1]->f[8][IdNodeSE[0]],1,1,tag_x_r);
	}
	if(IdNodeSW.size()>=1)
	{
		MultiBlock_->Recv(&f[0]->f[8][IdNodeSW[0]],1,3,tag_x_r,status);
		MultiBlock_->Recv(&f[1]->f[8][IdNodeSW[0]],1,3,tag_x_r,status);
	}
/*	if(IdGNodeNE.size()>=1)
		MultiBlock_->Send(&f[0]->f[5][IdGNodeNE[0]],1,4,tag_d_tr);
	if(IdRNodeSW.size()>=1)
		MultiBlock_->Recv(&f[0]->f[5][IdRNodeSW[0]],1,7,tag_d_tr,status);
	if(IdGNodeNW.size()>=1)
		MultiBlock_->Send(&f[0]->f[2][IdGNodeNW[0]],1,0,tag_y_t);
	if(IdRNodeSW.size()>=1)
		MultiBlock_->Recv(&f[0]->f[2][IdRNodeSW[0]],1,2,tag_y_t,status);
	if(IdGNodeNW.size()>=1)
		MultiBlock_->Send(&f[0]->f[6][IdGNodeNW[0]],1,0,tag_y_t);
	if(IdRNodeSW.size()>=1)
		MultiBlock_->Recv(&f[0]->f[6][IdRNodeSW[0]],1,2,tag_y_t,status);
	if(IdGNodeSE.size()>=1)
		MultiBlock_->Send(&f[0]->f[1][IdGNodeSE[0]],1,1,tag_x_r);
	if(IdRNodeSW.size()>=1)
		MultiBlock_->Recv(&f[0]->f[1][IdRNodeSW[0]],1,3,tag_x_r,status);
	if(IdGNodeSE.size()>=1)
		MultiBlock_->Send(&f[0]->f[8][IdGNodeSE[0]],1,1,tag_x_r);
	if(IdRNodeSW.size()>=1)
		MultiBlock_->Recv(&f[0]->f[8][IdRNodeSW[0]],1,3,tag_x_r,status);*/
}
void D2Q9TwoPhases::CornerNodesSyncToGhost(){
	 int tag_x_l=1;
//	 int tag_y_b=2;
	 int tag_d_bl=3;
	//N=0,E=1,S=2,W=3,NE=4, SE=5, NW=6, SW=7;
	MPI_Status status;

	if(IdRNodeSE.size()>=1)
	{
		MultiBlock_->Send(&f[0]->f[8][IdRNodeSE[0]],1,5,tag_x_l);
		MultiBlock_->Send(&f[1]->f[8][IdRNodeSE[0]],1,5,tag_x_l);
	}
	if(IdGNodeNW.size()>=1)
	{
		MultiBlock_->Recv(&f[0]->f[8][IdGNodeNW[0]],1,6,tag_x_l,status);
		MultiBlock_->Recv(&f[1]->f[8][IdGNodeNW[0]],1,6,tag_x_l,status);
	}
	if(IdRNodeSW.size()>=1)
	{
		MultiBlock_->Send(&f[0]->f[7][IdRNodeSW[0]],1,7,tag_x_l);
		MultiBlock_->Send(&f[1]->f[7][IdRNodeSW[0]],1,7,tag_x_l);
	}
	if(IdGNodeNE.size()>=1)
	{
		MultiBlock_->Recv(&f[0]->f[7][IdGNodeNE[0]],1,4,tag_x_l,status);
		MultiBlock_->Recv(&f[1]->f[7][IdGNodeNE[0]],1,4,tag_x_l,status);
	}
	if(IdRNodeNE.size()>=1)
	{
		MultiBlock_->Send(&f[0]->f[5][IdRNodeNE[0]],1,4,tag_x_l);
		MultiBlock_->Send(&f[1]->f[5][IdRNodeNE[0]],1,4,tag_x_l);
	}
	if(IdGNodeSW.size()>=1)
	{
		MultiBlock_->Recv(&f[0]->f[5][IdGNodeSW[0]],1,7,tag_x_l,status);
		MultiBlock_->Recv(&f[1]->f[5][IdGNodeSW[0]],1,7,tag_x_l,status);
	}
	if(IdRNodeNW.size()>=1)
	{
		MultiBlock_->Send(&f[0]->f[6][IdRNodeNW[0]],1,6,tag_d_bl);
		MultiBlock_->Send(&f[1]->f[6][IdRNodeNW[0]],1,6,tag_d_bl);
	}
	if(IdGNodeSE.size()>=1)
	{
		MultiBlock_->Recv(&f[0]->f[6][IdGNodeSE[0]],1,5,tag_d_bl,status);
		MultiBlock_->Recv(&f[1]->f[6][IdGNodeSE[0]],1,5,tag_d_bl,status);
	}
}
void D2Q9TwoPhases::SyncFromGhost(){
	GhostNodesSyncFromGhost();
	CornerNodesSyncFromGhost();
}
void D2Q9TwoPhases::SyncToGhost(){
	GhostNodesSyncToGhost();
	CornerNodesSyncToGhost();
}

void D2Q9TwoPhases::SyncMacroVarToGhost(){
//Put variables in buffer arrays
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
//		MultiBlock_->Communication(buf_MacroSend[0],buf_MacroRecv[0],size_MacroBuf);
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
//Set variables from buffer to real variables
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

void D2Q9TwoPhases::SyncVarSolidGhost(double *&VarIn){
	//Put variables in buffer arrays
		for (unsigned int i=0;i<SolidIdGNodeW.size();i++)
		{
				buf_MacroSendSolid[0][1][i]=VarIn[SolidIdGNodeW[i]];
		}
		for (unsigned int i=0;i<SolidIdRNodeS.size();i++)
		{
				buf_MacroRecvSolid[0][3][i]=VarIn[SolidIdRNodeS[i]];
		}
		for (unsigned int i=0;i<SolidIdRNodeE.size();i++)
		{
				buf_MacroRecvSolid[0][1][i]=VarIn[SolidIdRNodeE[i]];
		}
		for (unsigned int i=0;i<SolidIdGNodeN.size();i++)
		{
				buf_MacroSendSolid[0][3][i]=VarIn[SolidIdGNodeN[i]];
		}
	//		MultiBlock_->Communication(buf_MacroSend[0],buf_MacroRecv[0],size_MacroBuf);
		int tag_x_r=1;
		int tag_x_l=2;
		int tag_y_t=3;
		int tag_y_b=4;
		int tag_d_bl=5;
		MPI_Status status;
		if(SolidIdRNodeE.size()>=1)
		{
				MultiBlock_->Send(&buf_MacroRecvSolid[0][1][0],SolidIdRNodeE.size(),1,tag_x_r);
		}
		if(SolidIdGNodeW.size()>=1)
		{
				MultiBlock_->Recv(&buf_MacroSendSolid[0][1][0],SolidIdGNodeW.size(),3,tag_x_r,status);
		}

		if(SolidIdGNodeW.size()>=1)
		{
				MultiBlock_->Send(&buf_MacroSendSolid[0][1][0],SolidIdGNodeW.size(),3,tag_x_l);
		}
		if(SolidIdRNodeE.size()>=1)
		{
				MultiBlock_->Recv(&buf_MacroRecvSolid[0][1][0],SolidIdRNodeE.size(),1,tag_x_l,status);
		}
		if(SolidIdRNodeN.size()>=1)
		{
				MultiBlock_->Send(&buf_MacroRecvSolid[0][2][0],SolidIdRNodeN.size(),0,tag_y_t);
		}
		if(SolidIdGNodeS.size()>=1)
		{
				MultiBlock_->Recv(&buf_MacroSendSolid[0][2][0],SolidIdGNodeS.size(),2,tag_y_t,status);
		}
		if(SolidIdRNodeS.size()>=1)
		{
				MultiBlock_->Send(&buf_MacroRecvSolid[0][3][0],SolidIdRNodeS.size(),2,tag_y_b);
		}
		if(SolidIdGNodeN.size()>=1)
		{
				MultiBlock_->Recv(&buf_MacroSendSolid[0][3][0],SolidIdGNodeN.size(),0,tag_y_b,status);
		}
	//Set variables from buffer to real variables
		for (unsigned int i=0;i<SolidIdGNodeE.size();i++)
		{
			VarIn[SolidIdGNodeE[i]]=buf_MacroSendSolid[0][0][i];
		}
		for (unsigned int i=0;i<SolidIdGNodeN.size();i++)
		{
			VarIn[SolidIdGNodeN[i]]=buf_MacroSendSolid[0][3][i];
		}
		for (unsigned int i=0;i<SolidIdGNodeW.size();i++)
		{
			VarIn[SolidIdGNodeW[i]]=buf_MacroSendSolid[0][1][i];
		}
		for (unsigned int i=0;i<SolidIdGNodeS.size();i++)
		{
			VarIn[SolidIdGNodeS[i]]=buf_MacroSendSolid[0][2][i];
		}

		if(SolidIdRNodeSE.size()>=1)
		{
				MultiBlock_->Send(&VarIn[SolidIdRNodeSE[0]],1,5,tag_x_l);
		}
		if(SolidIdGNodeNW.size()>=1)
		{
				MultiBlock_->Recv(&VarIn[SolidIdGNodeNW[0]],1,6,tag_x_l,status);
		}
		if(SolidIdRNodeSW.size()>=1)
		{
				MultiBlock_->Send(&VarIn[SolidIdRNodeSW[0]],1,7,tag_x_l);
		}
		if(SolidIdGNodeNE.size()>=1)
		{
				MultiBlock_->Recv(&VarIn[SolidIdGNodeNE[0]],1,4,tag_x_l,status);
		}
		if(SolidIdRNodeNE.size()>=1)
		{
				MultiBlock_->Send(&VarIn[SolidIdRNodeNE[0]],1,4,tag_x_l);
		}
		if(SolidIdGNodeSW.size()>=1)
		{
				MultiBlock_->Recv(&VarIn[SolidIdGNodeSW[0]],1,7,tag_x_l,status);
		}
		if(SolidIdRNodeNW.size()>=1)
		{
				MultiBlock_->Send(&VarIn[SolidIdRNodeNW[0]],1,6,tag_d_bl);
		}
		if(SolidIdGNodeSE.size()>=1)
		{
				MultiBlock_->Recv(&VarIn[SolidIdGNodeSE[0]],1,5,tag_d_bl,status);
		}
}
void D2Q9TwoPhases::SyncVarToSolidGhost(double *&VarIn){
	//Put variables in buffer arrays
		for (unsigned int i=0;i<SolidIdRNodeW.size();i++)
		{
				buf_MacroRecvSolid[0][0][i]=VarIn[SolidIdRNodeW[i]];
		}
		for (unsigned int i=0;i<SolidIdRNodeS.size();i++)
		{
				buf_MacroRecvSolid[0][3][i]=VarIn[SolidIdRNodeS[i]];
		}
		for (unsigned int i=0;i<SolidIdRNodeE.size();i++)
		{
				buf_MacroRecvSolid[0][1][i]=VarIn[SolidIdRNodeE[i]];
		}
		for (unsigned int i=0;i<SolidIdRNodeN.size();i++)
		{
				buf_MacroRecvSolid[0][2][i]=VarIn[SolidIdRNodeN[i]];
		}
	//		MultiBlock_->Communication(buf_MacroSend[0],buf_MacroRecv[0],size_MacroBuf);
		int tag_x_r=1;
		int tag_x_l=2;
		int tag_y_t=3;
		int tag_y_b=4;
		int tag_d_bl=5;
		MPI_Status status;
		if(SolidIdRNodeE.size()>=1)
		{
				MultiBlock_->Send(&buf_MacroRecvSolid[0][1][0],SolidIdRNodeE.size(),1,tag_x_r);
		}
		if(SolidIdGNodeW.size()>=1)
		{
				MultiBlock_->Recv(&buf_MacroSendSolid[0][1][0],SolidIdGNodeW.size(),3,tag_x_r,status);
		}

		if(SolidIdRNodeW.size()>=1)
		{
				MultiBlock_->Send(&buf_MacroRecvSolid[0][0][0],SolidIdRNodeW.size(),3,tag_x_l);
		}
		if(SolidIdGNodeE.size()>=1)
		{
				MultiBlock_->Recv(&buf_MacroSendSolid[0][0][0],SolidIdGNodeE.size(),1,tag_x_l,status);
		}
		if(SolidIdRNodeN.size()>=1)
		{
				MultiBlock_->Send(&buf_MacroRecvSolid[0][2][0],SolidIdRNodeN.size(),0,tag_y_t);
		}
		if(SolidIdGNodeS.size()>=1)
		{
				MultiBlock_->Recv(&buf_MacroSendSolid[0][2][0],SolidIdGNodeS.size(),2,tag_y_t,status);
		}
		if(SolidIdRNodeS.size()>=1)
		{
				MultiBlock_->Send(&buf_MacroRecvSolid[0][3][0],SolidIdRNodeS.size(),2,tag_y_b);
		}
		if(SolidIdGNodeN.size()>=1)
		{
				MultiBlock_->Recv(&buf_MacroSendSolid[0][3][0],SolidIdGNodeN.size(),0,tag_y_b,status);
		}
	//Set variables from buffer to real variables
		for (unsigned int i=0;i<SolidIdGNodeE.size();i++)
		{
			VarIn[SolidIdGNodeE[i]]=buf_MacroSendSolid[0][0][i];
		}
		for (unsigned int i=0;i<SolidIdGNodeN.size();i++)
		{
			VarIn[SolidIdGNodeN[i]]=buf_MacroSendSolid[0][3][i];
		}
		for (unsigned int i=0;i<SolidIdGNodeW.size();i++)
		{
			VarIn[SolidIdGNodeW[i]]=buf_MacroSendSolid[0][1][i];
		}
		for (unsigned int i=0;i<SolidIdGNodeS.size();i++)
		{
			VarIn[SolidIdGNodeS[i]]=buf_MacroSendSolid[0][2][i];
		}

		if(SolidIdRNodeSE.size()>=1)
		{
				MultiBlock_->Send(&VarIn[SolidIdRNodeSE[0]],1,5,tag_x_l);
		}
		if(SolidIdGNodeNW.size()>=1)
		{
				MultiBlock_->Recv(&VarIn[SolidIdGNodeNW[0]],1,6,tag_x_l,status);
		}
		if(SolidIdRNodeSW.size()>=1)
		{
				MultiBlock_->Send(&VarIn[SolidIdRNodeSW[0]],1,7,tag_x_l);
		}
		if(SolidIdGNodeNE.size()>=1)
		{
				MultiBlock_->Recv(&VarIn[SolidIdGNodeNE[0]],1,4,tag_x_l,status);
		}
		if(SolidIdRNodeNE.size()>=1)
		{
				MultiBlock_->Send(&VarIn[SolidIdRNodeNE[0]],1,4,tag_x_l);
		}
		if(SolidIdGNodeSW.size()>=1)
		{
				MultiBlock_->Recv(&VarIn[SolidIdGNodeSW[0]],1,7,tag_x_l,status);
		}
		if(SolidIdRNodeNW.size()>=1)
		{
				MultiBlock_->Send(&VarIn[SolidIdRNodeNW[0]],1,6,tag_d_bl);
		}
		if(SolidIdGNodeSE.size()>=1)
		{
				MultiBlock_->Recv(&VarIn[SolidIdGNodeSE[0]],1,5,tag_d_bl,status);
		}
}
void D2Q9TwoPhases::SyncVarFromSolidGhost(double *&VarIn){
	//Put variables in buffer arrays
		for (unsigned int i=0;i<SolidIdGNodeW.size();i++)
		{
				buf_MacroSendSolid[0][1][i]=VarIn[SolidIdGNodeW[i]];
		}
		for (unsigned int i=0;i<SolidIdGNodeS.size();i++)
		{
				buf_MacroSendSolid[0][2][i]=VarIn[SolidIdGNodeS[i]];
		}
		for (unsigned int i=0;i<SolidIdGNodeE.size();i++)
		{
				buf_MacroSendSolid[0][0][i]=VarIn[SolidIdGNodeE[i]];
		}
		for (unsigned int i=0;i<SolidIdGNodeN.size();i++)
		{
				buf_MacroSendSolid[0][3][i]=VarIn[SolidIdGNodeN[i]];
		}
	//		MultiBlock_->Communication(buf_MacroSend[0],buf_MacroRecv[0],size_MacroBuf);
		int tag_x_r=1;
		int tag_x_l=2;
		int tag_y_t=3;
		int tag_y_b=4;
		int tag_d_bl=5;
		MPI_Status status;
		if(SolidIdGNodeE.size()>=1)
		{
				MultiBlock_->Send(&buf_MacroSendSolid[0][0][0],SolidIdGNodeE.size(),1,tag_x_r);
		}
		if(SolidIdRNodeW.size()>=1)
		{
				MultiBlock_->Recv(&buf_MacroRecvSolid[0][0][0],SolidIdRNodeW.size(),3,tag_x_r,status);
		}

		if(SolidIdGNodeW.size()>=1)
		{
				MultiBlock_->Send(&buf_MacroSendSolid[0][1][0],SolidIdGNodeW.size(),3,tag_x_l);
		}
		if(SolidIdRNodeE.size()>=1)
		{
				MultiBlock_->Recv(&buf_MacroRecvSolid[0][1][0],SolidIdRNodeE.size(),1,tag_x_l,status);
		}
		if(SolidIdGNodeN.size()>=1)
		{
				MultiBlock_->Send(&buf_MacroSendSolid[0][3][0],SolidIdGNodeN.size(),0,tag_y_t);
		}
		if(SolidIdRNodeS.size()>=1)
		{
				MultiBlock_->Recv(&buf_MacroRecvSolid[0][3][0],SolidIdRNodeS.size(),2,tag_y_t,status);
		}
		if(SolidIdGNodeS.size()>=1)
		{
				MultiBlock_->Send(&buf_MacroSendSolid[0][2][0],SolidIdGNodeS.size(),2,tag_y_b);
		}
		if(SolidIdRNodeN.size()>=1)
		{
				MultiBlock_->Recv(&buf_MacroRecvSolid[0][2][0],SolidIdRNodeN.size(),0,tag_y_b,status);
		}
	//Set variables from buffer to real variables
		for (unsigned int i=0;i<SolidIdRNodeE.size();i++)
		{
			VarIn[SolidIdRNodeE[i]]=buf_MacroRecvSolid[0][1][i];
		}
		for (unsigned int i=0;i<SolidIdRNodeN.size();i++)
		{
			VarIn[SolidIdRNodeN[i]]=buf_MacroRecvSolid[0][2][i];
		}
		for (unsigned int i=0;i<SolidIdRNodeW.size();i++)
		{
			VarIn[SolidIdRNodeW[i]]=buf_MacroRecvSolid[0][0][i];
		}

		for (unsigned int i=0;i<SolidIdRNodeS.size();i++)
		{
			VarIn[SolidIdRNodeS[i]]=buf_MacroRecvSolid[0][3][i];
		}
		if(SolidIdGNodeSE.size()>=1)
		{
				MultiBlock_->Send(&VarIn[SolidIdGNodeSE[0]],1,5,tag_x_l);
		}
		if(SolidIdRNodeNW.size()>=1)
		{
				MultiBlock_->Recv(&VarIn[SolidIdRNodeNW[0]],1,6,tag_x_l,status);
		}
		if(SolidIdGNodeSW.size()>=1)
		{
				MultiBlock_->Send(&VarIn[SolidIdGNodeSW[0]],1,7,tag_x_l);
		}
		if(SolidIdRNodeNE.size()>=1)
		{
				MultiBlock_->Recv(&VarIn[SolidIdRNodeNE[0]],1,4,tag_x_l,status);
		}
		if(SolidIdGNodeNE.size()>=1)
		{
				MultiBlock_->Send(&VarIn[SolidIdGNodeNE[0]],1,4,tag_x_l);
		}
		if(SolidIdRNodeSW.size()>=1)
		{
				MultiBlock_->Recv(&VarIn[SolidIdRNodeSW[0]],1,7,tag_x_l,status);
		}
		if(SolidIdGNodeNW.size()>=1)
		{
				MultiBlock_->Send(&VarIn[SolidIdGNodeNW[0]],1,6,tag_d_bl);
		}
		if(SolidIdRNodeSE.size()>=1)
		{
				MultiBlock_->Recv(&VarIn[SolidIdRNodeSE[0]],1,5,tag_d_bl,status);
		}
}
double D2Q9TwoPhases::Convert_Alpha_To_Rho(double alpha)
{

	return PtrParameters->Get_Rho_1()*alpha + PtrParameters->Get_Rho_2()*(1-alpha);
}
double D2Q9TwoPhases::Convert_Rho_To_Alpha(double Rho)
{

	return (Rho-PtrParameters->Get_Rho_2()) /(PtrParameters->Get_Rho_1() - PtrParameters->Get_Rho_2());
}
void D2Q9TwoPhases::InitialiseFromFile(){
	//Read From file
	if(PtrParameters->IsInitFromFile())
	{
		for(int i=0;i<PtrParameters->Get_NumberVariableToInit();i++)
		{
			Read_Variable(PtrParameters->Get_VariableNameToInit(i),PtrParameters->Get_FileNameToInit(i));
		}
	}
}
