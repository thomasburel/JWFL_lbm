/*
 * D2Q9ColourFluid.cpp
 *
 *  Created on: 31 Jul 2016
 *      Author: Thomas Burel
 */

#include "D2Q9ColourFluid.h"

D2Q9ColourFluid::D2Q9ColourFluid() {
	Rhor=0;
	Rhob=0;
}

D2Q9ColourFluid::~D2Q9ColourFluid() {
	delete [] Rhor;
	delete [] Rhob;
}

D2Q9ColourFluid::D2Q9ColourFluid(MultiBlock* MultiBlock__,ParallelManager* parallel__,WriterManager* Writer__, Parameters* Parameters__,InitLBM& ini){
	InitD2Q9TwoPhases(MultiBlock__,parallel__,Writer__, Parameters__, ini);
	InitMultiphase(ini);
	Set_Colour_gradient();
}

void D2Q9ColourFluid::InitMultiphase(InitLBM& ini){
// Initialise only the Colour Fluid part.
	Rhor=new double [nbnodes_total];
	Rhob=new double [nbnodes_total];
	double alpha=0;
	double* pos =new double[2];
	double* U_=new double[2];
// Loops for all kind of nodes
	for (int j=0;j<NodeArrays->NodeInterior.size();j++)
	{
// Get position
		pos[0]=NodeArrays->NodeInterior[j].get_x();
		pos[1]=NodeArrays->NodeInterior[j].get_y();
// Get initialise value from the user
		ini.IniDomain(parallel->getRank(),NodeArrays->NodeInterior[j],0, NodeArrays->NodeInterior[j].Get_index(),pos,Rho[NodeArrays->NodeInterior[j].Get_index()],U_,alpha);
// Initialise the blue and red densities
		Rhor[NodeArrays->NodeInterior[j].Get_index()]=alpha*Rho[NodeArrays->NodeInterior[j].Get_index()];//PtrParameters->Get_Rho_1();
		Rhob[NodeArrays->NodeInterior[j].Get_index()]=(1- alpha) *Rho[NodeArrays->NodeInterior[j].Get_index()];//PtrParameters->Get_Rho_2();
	}
	for (int j=0;j<NodeArrays->NodeCorner.size();j++)
	{
// Get position
		pos[0]=NodeArrays->NodeCorner[j].get_x();
		pos[1]=NodeArrays->NodeCorner[j].get_y();
// Get initialise value from the user
		ini.IniDomain(parallel->getRank(),NodeArrays->NodeCorner[j],0, NodeArrays->NodeCorner[j].Get_index(),pos,Rho[NodeArrays->NodeCorner[j].Get_index()],U_,alpha);
// Initialise the blue and red densities
		Rhor[NodeArrays->NodeCorner[j].Get_index()]=alpha*Rho[NodeArrays->NodeCorner[j].Get_index()];//PtrParameters->Get_Rho_1();
		Rhob[NodeArrays->NodeCorner[j].Get_index()]=(1- alpha) *Rho[NodeArrays->NodeCorner[j].Get_index()];//PtrParameters->Get_Rho_2();
	}
	for (int j=0;j<NodeArrays->NodeGlobalCorner.size();j++)
	{
// Get position
		pos[0]=NodeArrays->NodeGlobalCorner[j].get_x();
		pos[1]=NodeArrays->NodeGlobalCorner[j].get_y();
// Get initialise values from the user
		ini.IniDomain(parallel->getRank(),NodeArrays->NodeGlobalCorner[j],0, NodeArrays->NodeGlobalCorner[j].Get_index(),pos,Rho[NodeArrays->NodeGlobalCorner[j].Get_index()],U_,alpha);
// Initialise the blue and red densities
		Rhor[NodeArrays->NodeGlobalCorner[j].Get_index()]=alpha*Rho[NodeArrays->NodeGlobalCorner[j].Get_index()];//PtrParameters->Get_Rho_1();
		Rhob[NodeArrays->NodeGlobalCorner[j].Get_index()]=(1- alpha) *Rho[NodeArrays->NodeGlobalCorner[j].Get_index()];//PtrParameters->Get_Rho_2();
	}
	for (int j=0;j<NodeArrays->NodeVelocity.size();j++)
	{
// Get position
		pos[0]=NodeArrays->NodeVelocity[j].get_x();
		pos[1]=NodeArrays->NodeVelocity[j].get_y();
// Get initialise values from the user
		ini.IniDomain(parallel->getRank(),NodeArrays->NodeVelocity[j],0, NodeArrays->NodeVelocity[j].Get_index(),pos,Rho[NodeArrays->NodeVelocity[j].Get_index()],U_,alpha);
// Initialise the blue and red densities
		Rhor[NodeArrays->NodeVelocity[j].Get_index()]=alpha*Rho[NodeArrays->NodeVelocity[j].Get_index()];//PtrParameters->Get_Rho_1();
		Rhob[NodeArrays->NodeVelocity[j].Get_index()]=(1- alpha) *Rho[NodeArrays->NodeVelocity[j].Get_index()];//PtrParameters->Get_Rho_2();
	}

	for (int j=0;j<NodeArrays->NodePressure.size();j++)
	{
// Get position
		pos[0]=NodeArrays->NodePressure[j].get_x();
		pos[1]=NodeArrays->NodePressure[j].get_y();
// Get initialise values from the user
		ini.IniDomain(parallel->getRank(),NodeArrays->NodePressure[j],0, NodeArrays->NodePressure[j].Get_index(),pos,Rho[NodeArrays->NodePressure[j].Get_index()],U_,alpha);
// Initialise the blue and red densities
		Rhor[NodeArrays->NodePressure[j].Get_index()]=alpha*Rho[NodeArrays->NodePressure[j].Get_index()];//PtrParameters->Get_Rho_1();
		Rhob[NodeArrays->NodePressure[j].Get_index()]=(1- alpha) *Rho[NodeArrays->NodePressure[j].Get_index()];//PtrParameters->Get_Rho_2();
	}
	for (int j=0;j<NodeArrays->NodeWall.size();j++)
	{
// Get position
		pos[0]=NodeArrays->NodeWall[j].get_x();
		pos[1]=NodeArrays->NodeWall[j].get_y();
// Get initialise values from the user
		ini.IniDomain(parallel->getRank(),NodeArrays->NodeWall[j],0, NodeArrays->NodeWall[j].Get_index(),pos,Rho[NodeArrays->NodeWall[j].Get_index()],U_,alpha);
// Initialise the blue and red densities
		Rhor[NodeArrays->NodeWall[j].Get_index()]=alpha*Rho[NodeArrays->NodeWall[j].Get_index()];//PtrParameters->Get_Rho_1();
		Rhob[NodeArrays->NodeWall[j].Get_index()]=(1- alpha) *Rho[NodeArrays->NodeWall[j].Get_index()];//PtrParameters->Get_Rho_2();
	}
	for (int j=0;j<NodeArrays->NodeSpecialWall.size();j++)
	{
// Get position
		pos[0]=NodeArrays->NodeSpecialWall[j].get_x();
		pos[1]=NodeArrays->NodeSpecialWall[j].get_y();
// Get initialise values from the user
		ini.IniDomain(parallel->getRank(),NodeArrays->NodeSpecialWall[j],0, NodeArrays->NodeSpecialWall[j].Get_index(),pos,Rho[NodeArrays->NodeSpecialWall[j].Get_index()],U_,alpha);
// Initialise the blue and red densities
		Rhor[NodeArrays->NodeSpecialWall[j].Get_index()]=alpha*Rho[NodeArrays->NodeSpecialWall[j].Get_index()];//PtrParameters->Get_Rho_1();
		Rhob[NodeArrays->NodeSpecialWall[j].Get_index()]=(1- alpha) *Rho[NodeArrays->NodeSpecialWall[j].Get_index()];//PtrParameters->Get_Rho_2();
	}
	for (int j=0;j<NodeArrays->NodeSymmetry.size();j++)
	{
// Get position
		pos[0]=NodeArrays->NodeSymmetry[j].get_x();
		pos[1]=NodeArrays->NodeSymmetry[j].get_y();
// Get initialise values from the user
		ini.IniDomain(parallel->getRank(),NodeArrays->NodeSymmetry[j],0, NodeArrays->NodeSymmetry[j].Get_index(),pos,Rho[NodeArrays->NodeSymmetry[j].Get_index()],U_,alpha);
// Initialise the blue and red densities
		Rhor[NodeArrays->NodeSymmetry[j].Get_index()]=alpha*Rho[NodeArrays->NodeSymmetry[j].Get_index()];//PtrParameters->Get_Rho_1();
		Rhob[NodeArrays->NodeSymmetry[j].Get_index()]=(1- alpha) *Rho[NodeArrays->NodeSymmetry[j].Get_index()];//PtrParameters->Get_Rho_2();
	}
	for (int j=0;j<NodeArrays->NodeGhost.size();j++)
	{
// Get position
		pos[0]=NodeArrays->NodeGhost[j].get_x();
		pos[1]=NodeArrays->NodeGhost[j].get_y();
// Get initialise values from the user
		ini.IniDomain(parallel->getRank(),NodeArrays->NodeGhost[j],0, NodeArrays->NodeGhost[j].Get_index(),pos,Rho[NodeArrays->NodeGhost[j].Get_index()],U_,alpha);
// Initialise the blue and red densities
		Rhor[NodeArrays->NodeGhost[j].Get_index()]=alpha*Rho[NodeArrays->NodeGhost[j].Get_index()];//PtrParameters->Get_Rho_1();
		Rhob[NodeArrays->NodeGhost[j].Get_index()]=(1- alpha) *Rho[NodeArrays->NodeGhost[j].Get_index()];//PtrParameters->Get_Rho_2();
	}
	for (int j=0;j<NodeArrays->NodeSolid.size();j++)
	{
// Get position
		pos[0]=NodeArrays->NodeSolid[j].get_x();
		pos[1]=NodeArrays->NodeSolid[j].get_y();
// Get initialise values from the user
		ini.IniDomain(parallel->getRank(),NodeArrays->NodeSolid[j],0, NodeArrays->NodeSolid[j].Get_index(),pos,Rho[NodeArrays->NodeSolid[j].Get_index()],U_,alpha);
// Initialise the blue and red densities
		Rhor[NodeArrays->NodeSolid[j].Get_index()]=alpha*Rho[NodeArrays->NodeSolid[j].Get_index()];//PtrParameters->Get_Rho_1();
		Rhob[NodeArrays->NodeSolid[j].Get_index()]=(1- alpha) *Rho[NodeArrays->NodeSolid[j].Get_index()];//PtrParameters->Get_Rho_2();
	}
	for (int i=0;i<nbvelo;i++)
	{
		for (int j=0;j<nbnode;j++)
		{
			f[0]->f[i][j]=CollideLowOrder::EquiDistriFunct2D(Rhor[j], U[0][j], U[1][j], &Ei[i][0], omega[i]);
			f[1]->f[i][j]=CollideLowOrder::EquiDistriFunct2D(Rhob[j], U[0][j], U[1][j], &Ei[i][0], omega[i]);
		}
	}

	delete [] pos;
	delete [] U_;
}
void D2Q9ColourFluid::run(){
	int NbStep=PtrParameters->Get_NbStep();
	int OutPutNStep=PtrParameters->Get_OutPutNSteps();
	int listing=PtrParameters->Get_listing();
	double time_inirun=parallel->getTime();
	double time_run=0;
	int it=0;

	if(parallel->getSize()>1)
		SyncToGhost();
	ApplyBc();
	UpdateMacroVariables();
	Writer->Write_Output(it);

	if(parallel->getSize()>1)
	{
		for (int i=1;i<NbStep+1;i++)
		{

			ColourFluid_Collision();
			SyncToGhost();
			StreamD2Q9();;
			ApplyBc();
			UpdateMacroVariables();
			SyncMacroVarToGhost();
			if(i%OutPutNStep==0)
			{
				Writer->Write_Output(i);
			}
			if(i%listing==0 && parallel->isMainProcessor() )
			{
				time_run=parallel->getTime()-time_inirun;

				std::cout<<"Iteration number: "<<i<< " Running time: ";
				if(time_run>60.0)
					std::cout<<trunc(time_run/60)<<"Minutes ";
				else //less than 1 minute
				std::cout<<trunc(time_run)<<"Seconds ";
				std::cout<< " Time per iteration: "<<time_run/i<<"s"<<std::endl;
			}
		}
	}
	else
	{
		for (int i=1;i<NbStep+1;i++)
		{

			ColourFluid_Collision();
			StreamD2Q9();
			ApplyBc();
			UpdateMacroVariables();
			if(i%OutPutNStep==0)
			{
				Writer->Write_Output(i);
			}
			if(i%listing==0)
			{
				time_run=parallel->getTime()-time_inirun;
				std::cout<<"Iteration number: "<<i<< " Running time: ";
				if(time_run>60.0)
					std::cout<<trunc(time_run/60)<<"Minutes ";
				else //less than 1 minute
				std::cout<<trunc(time_run)<<"Seconds ";
				std::cout<< " Time per iteration: "<<time_run/i<<"s"<<std::endl;			}
		}
	}
	Writer->Write_breakpoint(*PtrParameters);
}
void D2Q9ColourFluid::ColourFluid_Collision()
{
	double Ak=0.65;
	double F[2];
	double F_Norm=0;
	int idx_tmp;
	double wtmp=0;
//	double *U_tmp, *V_tmp, *F_tmp;
	double  InvTau_,fi_tmp;
	InvTau_=InvTau; //Tmp

	for (int j=0;j<NodeArrays->NodeInterior.size();j++)
	{
	// Common variables
		idx_tmp=NodeArrays->NodeInterior[j].Get_index();
	// Calculate gradients
		(this->*PtrColourGrad)(idx_tmp,&F[0]);
	//Calculate Norms
		F_Norm=std::sqrt(F[0]*F[0]+F[1]*F[1]);
	//Model
		for (int i=0;i<9;i++)
		{
			fi_tmp=f[0]->f[i][idx_tmp]+f[1]->f[i][idx_tmp];
			Collide_ColorFluid(i,fi_tmp,Rho[idx_tmp],&F[0],F_Norm, InvTau_, U[0][idx_tmp], U[1][idx_tmp]);
			Recoloring(fi_tmp, f[0]->f[i][idx_tmp], f[1]->f[i][idx_tmp], Rho[idx_tmp], Rhor[idx_tmp], Rhob[idx_tmp]);
		}
	}
	for (int j=0;j<NodeArrays->NodeCorner.size();j++)
	{
	// Common variables
		idx_tmp=NodeArrays->NodeCorner[j].Get_index();
	// Calculate gradients
		Colour_gradient(idx_tmp,&F[0]);
	//Calculate Norms
		F_Norm=std::sqrt(F[0]*F[0]+F[1]*F[1]);
	//Model
		for (int i=0;i<9;i++)
		{
			fi_tmp=f[0]->f[i][idx_tmp]+f[1]->f[i][idx_tmp];
			Collide_ColorFluid(i,fi_tmp,Rho[idx_tmp],&F[0],F_Norm, InvTau_, U[0][idx_tmp], U[1][idx_tmp]);
			Recoloring(fi_tmp, f[0]->f[i][idx_tmp], f[1]->f[i][idx_tmp], Rho[idx_tmp], Rhor[idx_tmp], Rhob[idx_tmp]);
		}
	}
	for (int j=0;j<NodeArrays->NodeGlobalCorner.size();j++)
	{
	// Common variables
		idx_tmp=NodeArrays->NodeGlobalCorner[j].Get_index();
	// Calculate gradients
		Colour_gradient(idx_tmp,&F[0]);
	//Calculate Norms
		F_Norm=std::sqrt(F[0]*F[0]+F[1]*F[1]);
	//Model
		for (int i=0;i<9;i++)
		{
			fi_tmp=f[0]->f[i][idx_tmp]+f[1]->f[i][idx_tmp];
			Collide_ColorFluid(i,fi_tmp,Rho[idx_tmp],&F[0],F_Norm, InvTau_, U[0][idx_tmp], U[1][idx_tmp]);
			Recoloring(fi_tmp, f[0]->f[i][idx_tmp], f[1]->f[i][idx_tmp], Rho[idx_tmp], Rhor[idx_tmp], Rhob[idx_tmp]);
		}
	}
	for (int j=0;j<NodeArrays->NodeVelocity.size();j++)
	{
	// Common variables
		idx_tmp=NodeArrays->NodeVelocity[j].Get_index();
	// Calculate gradients
		Colour_gradient(idx_tmp,&F[0]);
	//Calculate Norms
		F_Norm=std::sqrt(F[0]*F[0]+F[1]*F[1]);
	//Model
		for (int i=0;i<9;i++)
		{
			fi_tmp=f[0]->f[i][idx_tmp]+f[1]->f[i][idx_tmp];
			Collide_ColorFluid(i,fi_tmp,Rho[idx_tmp],&F[0],F_Norm, InvTau_, U[0][idx_tmp], U[1][idx_tmp]);
			Recoloring(fi_tmp, f[0]->f[i][idx_tmp], f[1]->f[i][idx_tmp], Rho[idx_tmp], Rhor[idx_tmp], Rhob[idx_tmp]);
		}
	}
	for (int j=0;j<NodeArrays->NodePressure.size();j++)
	{
	// Common variables
		idx_tmp=NodeArrays->NodePressure[j].Get_index();
	// Calculate gradients
		Colour_gradient(idx_tmp,&F[0]);
	//Calculate Norms
		F_Norm=std::sqrt(F[0]*F[0]+F[1]*F[1]);
	//Model
		for (int i=0;i<9;i++)
		{
			fi_tmp=f[0]->f[i][idx_tmp]+f[1]->f[i][idx_tmp];
			Collide_ColorFluid(i,fi_tmp,Rho[idx_tmp],&F[0],F_Norm, InvTau_, U[0][idx_tmp], U[1][idx_tmp]);
			Recoloring(fi_tmp, f[0]->f[i][idx_tmp], f[1]->f[i][idx_tmp], Rho[idx_tmp], Rhor[idx_tmp], Rhob[idx_tmp]);
		}
	}
	for (int j=0;j<NodeArrays->NodeWall.size();j++)
	{
	// Common variables
		idx_tmp=NodeArrays->NodeWall[j].Get_index();
	// Calculate gradients
		Colour_gradient(idx_tmp,&F[0]);
	//Calculate Norms
		F_Norm=std::sqrt(F[0]*F[0]+F[1]*F[1]);
	//Model
		for (int i=0;i<9;i++)
		{
			fi_tmp=f[0]->f[i][idx_tmp]+f[1]->f[i][idx_tmp];
			Collide_ColorFluid(i,fi_tmp,Rho[idx_tmp],&F[0],F_Norm, InvTau_, U[0][idx_tmp], U[1][idx_tmp]);
			Recoloring(fi_tmp, f[0]->f[i][idx_tmp], f[1]->f[i][idx_tmp], Rho[idx_tmp], Rhor[idx_tmp], Rhob[idx_tmp]);
		}
	}
	for (int j=0;j<NodeArrays->NodeSymmetry.size();j++)
	{
	// Common variables
		idx_tmp=NodeArrays->NodeSymmetry[j].Get_index();
	// Calculate gradients
		Colour_gradient(idx_tmp,&F[0]);
	//Calculate Norms
		F_Norm=std::sqrt(F[0]*F[0]+F[1]*F[1]);
	//Model
		for (int i=0;i<9;i++)
		{
			fi_tmp=f[0]->f[i][idx_tmp]+f[1]->f[i][idx_tmp];
			Collide_ColorFluid(i,fi_tmp,Rho[idx_tmp],&F[0],F_Norm, InvTau_, U[0][idx_tmp], U[1][idx_tmp]);
			Recoloring(fi_tmp, f[0]->f[i][idx_tmp], f[1]->f[i][idx_tmp], Rho[idx_tmp], Rhor[idx_tmp], Rhob[idx_tmp]);
		}
	}
}

void D2Q9ColourFluid::Set_Colour_gradient(){
	switch(PtrParameters->Get_ColourGradType())
	{
	case Gunstensen:
		PtrColourGrad =&D2Q9ColourFluid::Colour_gradient_Gunstensen;
		break;
	case DensityGrad:
		PtrColourGrad =&D2Q9ColourFluid::Colour_gradient_DensityGrad;
		break;
	case DensityNormalGrad:
		PtrColourGrad =&D2Q9ColourFluid::Colour_gradient_DensityNormalGrad;
		break;
	}
}
void D2Q9ColourFluid::Colour_gradient(int & nodenumber, double* F){
	for (int k=0; k<nbvelo;k++)
	{
		for (int j=0;j<2;j++)
		{
			F[j]=(Rhor[nodenumber]-Rhob[nodenumber])*Ei[k][j];
		}
	}
}
void D2Q9ColourFluid::Colour_gradient_Gunstensen(int & nodenumber, double* F){
	for (int k=0; k<nbvelo;k++)
	{
		for (int j=0;j<2;j++)
		{
			F[j]=(Rhor[nodenumber]-Rhob[nodenumber])*Ei[k][j];
		}
	}
}
void D2Q9ColourFluid::Colour_gradient_DensityGrad(int & nodenumber, double* F){

}
void D2Q9ColourFluid::Colour_gradient_DensityNormalGrad(int & nodenumber, double* F){

}
double D2Q9ColourFluid::TwoPhase_Collision_operator(int & nodenumber, int & i, double & Ak, double* F, double & F_Norm){
 return Ak*0.5*F_Norm*(((F[0]*Ei[i][0]+F[1]*Ei[i][1])/F_Norm)*((F[0]*Ei[i][0]+F[1]*Ei[i][1])/F_Norm)-3/4);
}
///Select and apply boundary conditions
void D2Q9ColourFluid::ApplyBc(){

	for (int j=0;j<NodeArrays->NodeVelocity.size();j++)
	{
		ApplyHeZou_U(NodeArrays->NodeVelocity[j],0,NodeArrays->NodeVelocity[j].Get_UDef()[0],NodeArrays->NodeVelocity[j].Get_UDef()[1]);
		ApplyHeZou_U(NodeArrays->NodeVelocity[j],1,NodeArrays->NodeVelocity[j].Get_UDef()[0],NodeArrays->NodeVelocity[j].Get_UDef()[1]);

	}

	for (int j=0;j<NodeArrays->NodePressure.size();j++)
	{
		ApplyHeZou_P(NodeArrays->NodePressure[j],0,NodeArrays->NodePressure[j].Get_RhoDef(),U[0][NodeArrays->NodePressure[j].Get_index()],U[1][NodeArrays->NodePressure[j].Get_index()]);
		ApplyHeZou_P(NodeArrays->NodePressure[j],1,NodeArrays->NodePressure[j].Get_RhoDef(),U[0][NodeArrays->NodePressure[j].Get_index()],U[1][NodeArrays->NodePressure[j].Get_index()]);
	}
	switch(PtrParameters->Get_WallType())
	{
		case BounceBack:
		for (int j=0;j<NodeArrays->NodeWall.size();j++)
		{
			ApplyBounceBack(NodeArrays->NodeWall[j]);
		}
		for (int j=0;j<NodeArrays->NodeSpecialWall.size();j++)
		{
			ApplyBounceBackSymmetry(NodeArrays->NodeSpecialWall[j]);
		}
		for (int j=0;j<NodeArrays->NodeCorner.size();j++)
		{
	//		ApplyCorner(NodeArrays->NodeCorner[j]);
			ApplyBounceBack(NodeArrays->NodeCorner[j]);
		}
		break;
		case Diffuse:
		for (int j=0;j<NodeArrays->NodeWall.size();j++)
		{
			ApplyDiffuseWall(NodeArrays->NodeWall[j]);
		}
		for (int j=0;j<NodeArrays->NodeSpecialWall.size();j++)
		{
			ApplyDiffuseWallSymmetry(NodeArrays->NodeSpecialWall[j]);
		}
		for (int j=0;j<NodeArrays->NodeCorner.size();j++)
		{
			ApplyDiffuseWall(NodeArrays->NodeCorner[j]);
		}
		break;
	}
	for (int j=0;j<NodeArrays->NodeGlobalCorner.size();j++)
	{
		ApplyGlobalCorner(NodeArrays->NodeGlobalCorner[j]);
	}
	for (int j=0;j<NodeArrays->NodeSymmetry.size();j++)
	{
		if(NodeArrays->NodeSymmetry[j].Get_SymmetryType()==OnNode)
			ApplySymmetryOnNode(NodeArrays->NodeSymmetry[j]);
		if(NodeArrays->NodeSymmetry[j].Get_SymmetryType()==SymmetryPressureOnNode)
			ApplySymmetryPressureOnNode(NodeArrays->NodeSymmetry[j]);
	}
   parallel->barrier();

}


void D2Q9ColourFluid::ApplyGlobalCorner(NodeCorner2D& NodeIn){
	NodeType NodeTypeTmp1;
	NodeType NodeTypeTmp2;
	int normal=0;
	for(int j=0;j<2;j++)//Not an efficient way to treat the two distribution function but that affects only 4 nodes maximum
	{
		switch(NodeIn.Get_BcNormal())
		{
		case 5:
			NodeTypeTmp1=NodeArrays->TypeOfNode[NodeIn.Get_connect()[2]];
			NodeTypeTmp2=NodeArrays->TypeOfNode[NodeIn.Get_connect()[1]];
			if(NodeTypeTmp1==Symmetry||NodeTypeTmp2==Symmetry)
			{
				if(NodeTypeTmp1==NodeTypeTmp2)
				{
					f[j]->f[1][NodeIn.Get_index()]=f[j]->f[3][NodeIn.Get_index()];//+0.00007/3;
					f[j]->f[2][NodeIn.Get_index()]=f[j]->f[4][NodeIn.Get_index()];
					f[j]->f[6][NodeIn.Get_index()]=f[j]->f[7][NodeIn.Get_index()];
					f[j]->f[8][NodeIn.Get_index()]=f[j]->f[7][NodeIn.Get_index()];//+0.00007/3;
					f[j]->f[5][NodeIn.Get_index()]=f[j]->f[7][NodeIn.Get_index()];//+0.00007/3;
//					for (int i=0;i<9;i++)
//						f[j]->f[i][NodeIn.Get_index()]=f[j]->f[i][NodeIn.Get_index()]+omega[i]*Ei[i][0]*0.00007;
					std::cerr<<"Symmetry-Symmetry BC not yet fully implemented"<<std::endl;
				}

				else
				{
				if(NodeTypeTmp1!=Symmetry)
					{
						normal=1;
						f[j]->f[2][NodeIn.Get_index()]=f[j]->f[4][NodeIn.Get_index()];
						f[j]->f[6][NodeIn.Get_index()]=f[j]->f[7][NodeIn.Get_index()];
						if(NodeTypeTmp1==Pressure)
						{
							ApplyHeZou_P(NodeIn,normal,0,NodeIn.Get_RhoDef(),U[0][NodeIn.Get_index()],U[1][NodeIn.Get_index()]);
							ApplyHeZou_P(NodeIn,normal,1,NodeIn.Get_RhoDef(),U[0][NodeIn.Get_index()],U[1][NodeIn.Get_index()]);
						}
						else
						{
							ApplyHeZou_U(NodeIn,normal,0,U[0][NodeIn.Get_index()],U[1][NodeIn.Get_index()]);
							ApplyHeZou_U(NodeIn,normal,1,U[0][NodeIn.Get_index()],U[1][NodeIn.Get_index()]);
						}

					}
					else
					{
						normal=2;
						f[j]->f[8][NodeIn.Get_index()]=f[j]->f[7][NodeIn.Get_index()];
						f[j]->f[1][NodeIn.Get_index()]=f[j]->f[3][NodeIn.Get_index()];
						if(NodeTypeTmp2==Pressure)
						{
							ApplyHeZou_P(NodeIn,normal,0,NodeIn.Get_RhoDef(),U[0][NodeIn.Get_index()],U[1][NodeIn.Get_index()]);
							ApplyHeZou_P(NodeIn,normal,1,NodeIn.Get_RhoDef(),U[0][NodeIn.Get_index()],U[1][NodeIn.Get_index()]);
						}
						else
						{
							ApplyHeZou_U(NodeIn,normal,0,U[0][NodeIn.Get_index()],U[1][NodeIn.Get_index()]);
							ApplyHeZou_U(NodeIn,normal,1,U[0][NodeIn.Get_index()],U[1][NodeIn.Get_index()]);
						}
					}
				}
			}
			else if(NodeTypeTmp1==Wall||NodeTypeTmp2==Wall)
			{
				switch(PtrParameters->Get_WallType())
				{
					case BounceBack:
						ApplyBounceBack(NodeIn);
						break;
					case Diffuse:
						ApplyDiffuseWall(NodeIn);
						break;
				}
				if(NodeTypeTmp1==Pressure)
				{
					normal=1;
					ApplyHeZou_P(NodeIn,normal,0,NodeIn.Get_RhoDef(),U[0][NodeIn.Get_index()],U[1][NodeIn.Get_index()]);
					ApplyHeZou_P(NodeIn,normal,1,NodeIn.Get_RhoDef(),U[0][NodeIn.Get_index()],U[1][NodeIn.Get_index()]);
				}
				if(NodeTypeTmp2==Pressure)
				{
					normal=2;
					ApplyHeZou_P(NodeIn,normal,0,NodeIn.Get_RhoDef(),U[0][NodeIn.Get_index()],U[1][NodeIn.Get_index()]);
					ApplyHeZou_P(NodeIn,normal,1,NodeIn.Get_RhoDef(),U[0][NodeIn.Get_index()],U[1][NodeIn.Get_index()]);
				}
				if(NodeTypeTmp1==Velocity)
				{
					normal=1;
					ApplyHeZou_U(NodeIn,normal,0,U[0][NodeIn.Get_index()],U[1][NodeIn.Get_index()]);
					ApplyHeZou_U(NodeIn,normal,1,U[0][NodeIn.Get_index()],U[1][NodeIn.Get_index()]);
				}
				if(NodeTypeTmp2==Velocity)
				{
					normal=2;
					ApplyHeZou_U(NodeIn,normal,0,U[0][NodeIn.Get_index()],U[1][NodeIn.Get_index()]);
					ApplyHeZou_U(NodeIn,normal,1,U[0][NodeIn.Get_index()],U[1][NodeIn.Get_index()]);
				}
			}
			break;
		case 6:
			NodeTypeTmp1=NodeArrays->TypeOfNode[NodeIn.Get_connect()[2]];
			NodeTypeTmp2=NodeArrays->TypeOfNode[NodeIn.Get_connect()[3]];
			if(NodeTypeTmp1==Symmetry||NodeTypeTmp2==Symmetry)
			{
				if(NodeTypeTmp1==NodeTypeTmp2)
				{
					f[j]->f[2][NodeIn.Get_index()]=f[j]->f[4][NodeIn.Get_index()];
					f[j]->f[3][NodeIn.Get_index()]=f[j]->f[1][NodeIn.Get_index()];//-0.00007/3;
					f[j]->f[6][NodeIn.Get_index()]=f[j]->f[8][NodeIn.Get_index()];//-0.00007/3;
					f[j]->f[5][NodeIn.Get_index()]=f[j]->f[8][NodeIn.Get_index()];
					f[j]->f[7][NodeIn.Get_index()]=f[j]->f[8][NodeIn.Get_index()];//-0.00007/3;
//					for (int i=0;i<9;i++)
//						f[j]->f[i][NodeIn.Get_index()]=f[j]->f[i][NodeIn.Get_index()]-omega[i]*Ei[i][0]*0.00007;
					std::cerr<<"Symmetry-Symmetry BC not yet fully implemented"<<std::endl;
				}
				else
				{
				if(NodeTypeTmp1!=Symmetry)
					{
						normal=3;
						f[j]->f[2][NodeIn.Get_index()]=f[j]->f[4][NodeIn.Get_index()];
						f[j]->f[5][NodeIn.Get_index()]=f[j]->f[8][NodeIn.Get_index()];
						if(NodeTypeTmp1==Pressure)
						{
							ApplyHeZou_P(NodeIn,normal,0,NodeIn.Get_RhoDef(),U[0][NodeIn.Get_index()],U[1][NodeIn.Get_index()]);
							ApplyHeZou_P(NodeIn,normal,1,NodeIn.Get_RhoDef(),U[0][NodeIn.Get_index()],U[1][NodeIn.Get_index()]);
						}
						else
						{
							ApplyHeZou_U(NodeIn,normal,0,U[0][NodeIn.Get_index()],U[1][NodeIn.Get_index()]);
							ApplyHeZou_U(NodeIn,normal,1,U[0][NodeIn.Get_index()],U[1][NodeIn.Get_index()]);
						}

					}
					else
					{
						normal=2;
						f[j]->f[7][NodeIn.Get_index()]=f[j]->f[8][NodeIn.Get_index()];
						f[j]->f[3][NodeIn.Get_index()]=f[j]->f[1][NodeIn.Get_index()];
						if(NodeTypeTmp2==Pressure)
						{
							ApplyHeZou_P(NodeIn,normal,0,NodeIn.Get_RhoDef(),U[0][NodeIn.Get_index()],U[1][NodeIn.Get_index()]);
							ApplyHeZou_P(NodeIn,normal,1,NodeIn.Get_RhoDef(),U[0][NodeIn.Get_index()],U[1][NodeIn.Get_index()]);
						}
						else
						{
							ApplyHeZou_U(NodeIn,normal,0,U[0][NodeIn.Get_index()],U[1][NodeIn.Get_index()]);
							ApplyHeZou_U(NodeIn,normal,1,U[0][NodeIn.Get_index()],U[1][NodeIn.Get_index()]);
						}
					}
				}
			}
			else if(NodeTypeTmp1==Wall||NodeTypeTmp2==Wall)
			{
				switch(PtrParameters->Get_WallType())
				{
					case BounceBack:
						ApplyBounceBack(NodeIn);
						break;
					case Diffuse:
						ApplyDiffuseWall(NodeIn);
						break;
				}
				if(NodeTypeTmp1==Pressure)
				{
					normal=3;
					ApplyHeZou_P(NodeIn,normal,0,NodeIn.Get_RhoDef(),U[0][NodeIn.Get_index()],U[1][NodeIn.Get_index()]);
					ApplyHeZou_P(NodeIn,normal,1,NodeIn.Get_RhoDef(),U[0][NodeIn.Get_index()],U[1][NodeIn.Get_index()]);
				}
				if(NodeTypeTmp2==Pressure)
				{
					normal=2;
					ApplyHeZou_P(NodeIn,normal,0,NodeIn.Get_RhoDef(),U[0][NodeIn.Get_index()],U[1][NodeIn.Get_index()]);
					ApplyHeZou_P(NodeIn,normal,1,NodeIn.Get_RhoDef(),U[0][NodeIn.Get_index()],U[1][NodeIn.Get_index()]);
				}
				if(NodeTypeTmp1==Velocity)
				{
					normal=3;
					ApplyHeZou_U(NodeIn,normal,0,U[0][NodeIn.Get_index()],U[1][NodeIn.Get_index()]);
					ApplyHeZou_U(NodeIn,normal,1,U[0][NodeIn.Get_index()],U[1][NodeIn.Get_index()]);
				}
				if(NodeTypeTmp2==Velocity)
				{
					normal=2;
					ApplyHeZou_U(NodeIn,normal,0,U[0][NodeIn.Get_index()],U[1][NodeIn.Get_index()]);
					ApplyHeZou_U(NodeIn,normal,1,U[0][NodeIn.Get_index()],U[1][NodeIn.Get_index()]);
				}
			}
			break;
		case 7:
			NodeTypeTmp1=NodeArrays->TypeOfNode[NodeIn.Get_connect()[4]];
			NodeTypeTmp2=NodeArrays->TypeOfNode[NodeIn.Get_connect()[3]];
			if(NodeTypeTmp1==Symmetry||NodeTypeTmp2==Symmetry)
			{
				if(NodeTypeTmp1==NodeTypeTmp2)
				{
					f[j]->f[3][NodeIn.Get_index()]=f[j]->f[1][NodeIn.Get_index()];//-0.00007/3;
					f[j]->f[4][NodeIn.Get_index()]=f[j]->f[2][NodeIn.Get_index()];
					f[j]->f[7][NodeIn.Get_index()]=f[j]->f[5][NodeIn.Get_index()];//-0.00007/3;
					f[j]->f[6][NodeIn.Get_index()]=f[j]->f[5][NodeIn.Get_index()];//-0.00007/3;
					f[j]->f[8][NodeIn.Get_index()]=f[j]->f[5][NodeIn.Get_index()];
//					for (int i=0;i<9;i++)
//						f[j]->f[i][NodeIn.Get_index()]=f[j]->f[i][NodeIn.Get_index()]-omega[i]*Ei[i][0]*0.00007;
					std::cerr<<"Symmetry-Symmetry BC not yet fully implemented"<<std::endl;
				}
				else
				{
				if(NodeTypeTmp1!=Symmetry)
					{
						normal=3;
						f[j]->f[4][NodeIn.Get_index()]=f[j]->f[2][NodeIn.Get_index()];
						f[j]->f[8][NodeIn.Get_index()]=f[j]->f[5][NodeIn.Get_index()];
						if(NodeTypeTmp1==Pressure)
						{
							ApplyHeZou_P(NodeIn,normal,0,NodeIn.Get_RhoDef(),U[0][NodeIn.Get_index()],U[1][NodeIn.Get_index()]);
							ApplyHeZou_P(NodeIn,normal,1,NodeIn.Get_RhoDef(),U[0][NodeIn.Get_index()],U[1][NodeIn.Get_index()]);
						}
						else
						{
							ApplyHeZou_U(NodeIn,normal,0,U[0][NodeIn.Get_index()],U[1][NodeIn.Get_index()]);
							ApplyHeZou_U(NodeIn,normal,1,U[0][NodeIn.Get_index()],U[1][NodeIn.Get_index()]);
						}
					}
					else
					{
						normal=4;
						f[j]->f[3][NodeIn.Get_index()]=f[j]->f[1][NodeIn.Get_index()];
						f[j]->f[6][NodeIn.Get_index()]=f[j]->f[5][NodeIn.Get_index()];
						if(NodeTypeTmp2==Pressure)
						{
							ApplyHeZou_P(NodeIn,normal,0,NodeIn.Get_RhoDef(),U[0][NodeIn.Get_index()],U[1][NodeIn.Get_index()]);
							ApplyHeZou_P(NodeIn,normal,1,NodeIn.Get_RhoDef(),U[0][NodeIn.Get_index()],U[1][NodeIn.Get_index()]);
						}
						else
						{
							ApplyHeZou_U(NodeIn,normal,0,U[0][NodeIn.Get_index()],U[1][NodeIn.Get_index()]);
							ApplyHeZou_U(NodeIn,normal,1,U[0][NodeIn.Get_index()],U[1][NodeIn.Get_index()]);
						}
					}
				}
			}
			else if(NodeTypeTmp1==Wall||NodeTypeTmp2==Wall)
			{
				switch(PtrParameters->Get_WallType())
				{
					case BounceBack:
						ApplyBounceBack(NodeIn);
						break;
					case Diffuse:
						ApplyDiffuseWall(NodeIn);
						break;
				}
				if(NodeTypeTmp1==Pressure)
				{
					normal=3;
					ApplyHeZou_P(NodeIn,normal,0,NodeIn.Get_RhoDef(),U[0][NodeIn.Get_index()],U[1][NodeIn.Get_index()]);
					ApplyHeZou_P(NodeIn,normal,1,NodeIn.Get_RhoDef(),U[0][NodeIn.Get_index()],U[1][NodeIn.Get_index()]);
				}
				if(NodeTypeTmp2==Pressure)
				{
					normal=4;
					ApplyHeZou_P(NodeIn,normal,0,NodeIn.Get_RhoDef(),U[0][NodeIn.Get_index()],U[1][NodeIn.Get_index()]);
					ApplyHeZou_P(NodeIn,normal,1,NodeIn.Get_RhoDef(),U[0][NodeIn.Get_index()],U[1][NodeIn.Get_index()]);
				}
				if(NodeTypeTmp1==Velocity)
				{
					normal=3;
					ApplyHeZou_U(NodeIn,normal,0,U[0][NodeIn.Get_index()],U[1][NodeIn.Get_index()]);
					ApplyHeZou_U(NodeIn,normal,1,U[0][NodeIn.Get_index()],U[1][NodeIn.Get_index()]);
				}
				if(NodeTypeTmp2==Velocity)
				{
					normal=4;
					ApplyHeZou_U(NodeIn,normal,0,U[0][NodeIn.Get_index()],U[1][NodeIn.Get_index()]);
					ApplyHeZou_U(NodeIn,normal,1,U[0][NodeIn.Get_index()],U[1][NodeIn.Get_index()]);
				}
			}
			break;
		case 8:
			NodeTypeTmp1=NodeArrays->TypeOfNode[NodeIn.Get_connect()[4]];
			NodeTypeTmp2=NodeArrays->TypeOfNode[NodeIn.Get_connect()[1]];
			if(NodeTypeTmp1==Symmetry||NodeTypeTmp2==Symmetry)
			{
				if(NodeTypeTmp1==NodeTypeTmp2)
				{
					f[j]->f[1][NodeIn.Get_index()]=f[j]->f[3][NodeIn.Get_index()];//+0.00007/3;
					f[j]->f[4][NodeIn.Get_index()]=f[j]->f[2][NodeIn.Get_index()];
					f[j]->f[8][NodeIn.Get_index()]=f[j]->f[6][NodeIn.Get_index()];//+0.00007/3;
					f[j]->f[5][NodeIn.Get_index()]=f[j]->f[6][NodeIn.Get_index()];//+0.00007/3;
					f[j]->f[7][NodeIn.Get_index()]=f[j]->f[6][NodeIn.Get_index()];
//					for (int i=0;i<9;i++)
//						f[j]->f[i][NodeIn.Get_index()]=f[j]->f[i][NodeIn.Get_index()]+omega[i]*Ei[i][0]*0.00007;
					std::cerr<<"Symmetry-Symmetry BC not yet fully implemented"<<std::endl;
				}
				else
				{
				if(NodeTypeTmp1!=Symmetry)
					{
						normal=1;
						f[j]->f[4][NodeIn.Get_index()]=f[j]->f[2][NodeIn.Get_index()];
						f[j]->f[7][NodeIn.Get_index()]=f[j]->f[6][NodeIn.Get_index()];
						if(NodeTypeTmp1==Pressure)
						{
							ApplyHeZou_P(NodeIn,normal,0,NodeIn.Get_RhoDef(),U[0][NodeIn.Get_index()],U[1][NodeIn.Get_index()]);
							ApplyHeZou_P(NodeIn,normal,1,NodeIn.Get_RhoDef(),U[0][NodeIn.Get_index()],U[1][NodeIn.Get_index()]);
						}
						else
						{
							ApplyHeZou_U(NodeIn,normal,0,U[0][NodeIn.Get_index()],U[1][NodeIn.Get_index()]);
							ApplyHeZou_U(NodeIn,normal,1,U[0][NodeIn.Get_index()],U[1][NodeIn.Get_index()]);
						}

					}
					else
					{
						normal=4;
						f[j]->f[5][NodeIn.Get_index()]=f[j]->f[6][NodeIn.Get_index()];
						f[j]->f[1][NodeIn.Get_index()]=f[j]->f[3][NodeIn.Get_index()];
						if(NodeTypeTmp2==Pressure)
						{
							ApplyHeZou_P(NodeIn,normal,0,NodeIn.Get_RhoDef(),U[0][NodeIn.Get_index()],U[1][NodeIn.Get_index()]);
							ApplyHeZou_P(NodeIn,normal,1,NodeIn.Get_RhoDef(),U[0][NodeIn.Get_index()],U[1][NodeIn.Get_index()]);
						}
						else
						{
							ApplyHeZou_U(NodeIn,normal,0,U[0][NodeIn.Get_index()],U[1][NodeIn.Get_index()]);
							ApplyHeZou_U(NodeIn,normal,1,U[0][NodeIn.Get_index()],U[1][NodeIn.Get_index()]);
						}
					}
				}
			}
			else if(NodeTypeTmp1==Wall||NodeTypeTmp2==Wall)
			{
				switch(PtrParameters->Get_WallType())
				{
					case BounceBack:
						ApplyBounceBack(NodeIn);
						break;
					case Diffuse:
						ApplyDiffuseWall(NodeIn);
						break;
				}
				if(NodeTypeTmp1==Pressure)
				{
					normal=1;
					ApplyHeZou_P(NodeIn,normal,0,NodeIn.Get_RhoDef(),U[0][NodeIn.Get_index()],U[1][NodeIn.Get_index()]);
					ApplyHeZou_P(NodeIn,normal,1,NodeIn.Get_RhoDef(),U[0][NodeIn.Get_index()],U[1][NodeIn.Get_index()]);
				}
				if(NodeTypeTmp2==Pressure)
				{
					normal=4;
					ApplyHeZou_P(NodeIn,normal,0,NodeIn.Get_RhoDef(),U[0][NodeIn.Get_index()],U[1][NodeIn.Get_index()]);
					ApplyHeZou_P(NodeIn,normal,1,NodeIn.Get_RhoDef(),U[0][NodeIn.Get_index()],U[1][NodeIn.Get_index()]);
				}
				if(NodeTypeTmp1==Velocity)
				{
					normal=1;
					ApplyHeZou_U(NodeIn,normal,0,U[0][NodeIn.Get_index()],U[1][NodeIn.Get_index()]);
					ApplyHeZou_U(NodeIn,normal,1,U[0][NodeIn.Get_index()],U[1][NodeIn.Get_index()]);
				}
				if(NodeTypeTmp2==Velocity)
				{
					normal=4;
					ApplyHeZou_U(NodeIn,normal,0,U[0][NodeIn.Get_index()],U[1][NodeIn.Get_index()]);
					ApplyHeZou_U(NodeIn,normal,1,U[0][NodeIn.Get_index()],U[1][NodeIn.Get_index()]);
				}
			}
			break;
		default:
			std::cerr<<"Direction: "<< NodeIn.Get_BcNormal()<<" not found for Global Corner"<<std::endl;
		}
	}
}

void D2Q9ColourFluid::ApplySymmetryPressureOnNode(NodeSymmetry2D& NodeIn){

	double rhosym=0;
	double rhoadd=0;
	switch(NodeIn.Get_BcNormal())
			{
			case 2:
				for(int j=0;j<2;j++)
				{
					//Applying symmetry
					f_tmp[0]=f[j]->f[0][NodeIn.Get_index()];
					f_tmp[1]=f[j]->f[1][NodeIn.Get_index()];
					f_tmp[2]=f[j]->f[4][NodeIn.Get_index()];
					f_tmp[3]=f[j]->f[3][NodeIn.Get_index()];
					f_tmp[4]=f[j]->f[4][NodeIn.Get_index()];
					f_tmp[5]=f[j]->f[8][NodeIn.Get_index()];
					f_tmp[6]=f[j]->f[7][NodeIn.Get_index()];
					f_tmp[7]=f[j]->f[7][NodeIn.Get_index()];
					f_tmp[8]=f[j]->f[8][NodeIn.Get_index()];
					//Applying He Zou Pressure
					BC_HeZou_P(NodeIn.Get_BcNormal(),f_tmp,NodeIn.Get_RhoDef(), U[0][NodeIn.Get_index()],U[1][NodeIn.Get_index()]);
					for (int i=0;i<9;i++)
						f[j]->f[i][NodeIn.Get_index()]=f_tmp[i];
				}
				break;
			case 4:
				for(int j=0;j<2;j++)
				{
					//Applying symmetry
					f_tmp[0]=f[j]->f[0][NodeIn.Get_index()];
					f_tmp[1]=f[j]->f[1][NodeIn.Get_index()];
					f_tmp[2]=f[j]->f[2][NodeIn.Get_index()];
					f_tmp[3]=f[j]->f[3][NodeIn.Get_index()];
					f_tmp[4]=f[j]->f[2][NodeIn.Get_index()];
					f_tmp[5]=f[j]->f[5][NodeIn.Get_index()];
					f_tmp[6]=f[j]->f[6][NodeIn.Get_index()];
					f_tmp[7]=f[j]->f[6][NodeIn.Get_index()];
					f_tmp[8]=f[j]->f[5][NodeIn.Get_index()];
					//Applying He Zou Pressure
					BC_HeZou_P(NodeIn.Get_BcNormal(),f_tmp,NodeIn.Get_RhoDef(), U[0][NodeIn.Get_index()],U[1][NodeIn.Get_index()]);
					for (int i=0;i<9;i++)
						f[j]->f[i][NodeIn.Get_index()]=f_tmp[i];
				}
				break;
			case 1:
				for(int j=0;j<2;j++)
				{
					//Applying symmetry
					f_tmp[0]=f[j]->f[0][NodeIn.Get_index()];
					f_tmp[1]=f[j]->f[3][NodeIn.Get_index()];
					f_tmp[2]=f[j]->f[2][NodeIn.Get_index()];
					f_tmp[3]=f[j]->f[3][NodeIn.Get_index()];
					f_tmp[4]=f[j]->f[4][NodeIn.Get_index()];
					f_tmp[5]=f[j]->f[6][NodeIn.Get_index()];
					f_tmp[6]=f[j]->f[6][NodeIn.Get_index()];
					f_tmp[7]=f[j]->f[7][NodeIn.Get_index()];
					f_tmp[8]=f[j]->f[7][NodeIn.Get_index()];
					//Applying He Zou Pressure
					BC_HeZou_P(NodeIn.Get_BcNormal(),f_tmp,NodeIn.Get_RhoDef(), U[0][NodeIn.Get_index()],U[1][NodeIn.Get_index()]);
					for (int i=0;i<9;i++)
						f[j]->f[i][NodeIn.Get_index()]=f_tmp[i];
				}
				break;
			case 3:
				for(int j=0;j<2;j++)
				{
					//Applying symmetry
					f_tmp[0]=f[j]->f[0][NodeIn.Get_index()];
					f_tmp[1]=f[j]->f[1][NodeIn.Get_index()];
					f_tmp[2]=f[j]->f[2][NodeIn.Get_index()];
					f_tmp[3]=f[j]->f[1][NodeIn.Get_index()];
					f_tmp[4]=f[j]->f[4][NodeIn.Get_index()];
					f_tmp[5]=f[j]->f[5][NodeIn.Get_index()];
					f_tmp[6]=f[j]->f[5][NodeIn.Get_index()];
					f_tmp[7]=f[j]->f[8][NodeIn.Get_index()];
					f_tmp[8]=f[j]->f[8][NodeIn.Get_index()];
					//Applying He Zou Pressure
					BC_HeZou_P(NodeIn.Get_BcNormal(),f_tmp,NodeIn.Get_RhoDef(), U[0][NodeIn.Get_index()],U[1][NodeIn.Get_index()]);
					for (int i=0;i<9;i++)
						f[j]->f[i][NodeIn.Get_index()]=f_tmp[i];
				}
				break;
			default:
				std::cerr<<"Direction symmetry not found"<<std::endl;
				break;
			}

}
void D2Q9ColourFluid::UpdateMacroVariables(){

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
		for (int j=0;j<NodeArrays->NodeSymmetry.size();j++)
		{
			MacroVariables(NodeArrays->NodeSymmetry[j].Get_index());
		}

}
void D2Q9ColourFluid::MacroVariables(int& idx){

		U[0][idx]=0;
		U[1][idx]=0;
		Rhor[idx]=0,Rhob[idx]=0;
		//Rho[idx]=0;
		for (int k=0; k<nbvelo;k++)
		{
			Rhor[idx]+=f[0]->f[k][idx];
			Rhob[idx]+=f[1]->f[k][idx];
			for (int j=0;j<2;j++)
			{
				U[j][idx]+=(f[0]->f[k][idx]+f[1]->f[k][idx])*Ei[k][j];
			}
		}
		Rho[idx]=Rhor[idx]+Rhob[idx];
		U[0][idx]=U[0][idx]/Rho[idx];
		U[1][idx]=U[1][idx]/Rho[idx];

}

double D2Q9ColourFluid::Cal_RhoR_Corner(NodeCorner2D& nodeIn){

	unsigned int direction1,direction2;
	switch(nodeIn.Get_BcNormal())
	{
	case 5:
		direction1=1;
		direction2=2;
		doubleTmpReturn=(Rhor[nodeIn.Get_connect()[direction1]]+Rhor[nodeIn.Get_connect()[direction2]])*0.5;//Connect(nodenumber,direction1)]+Rho[Connect(nodenumber,direction2)])*0.5;
		break;
	case 6:
		direction1=2;
		direction2=3;
		doubleTmpReturn=(Rhor[nodeIn.Get_connect()[direction1]]+Rhor[nodeIn.Get_connect()[direction2]])*0.5;//(Rho[Connect(nodenumber,direction1)]+Rho[Connect(nodenumber,direction2)])*0.5;
		break;
	case 7:
		direction1=3;
		direction2=4;
		doubleTmpReturn=(Rhor[nodeIn.Get_connect()[direction1]]+Rhor[nodeIn.Get_connect()[direction2]])*0.5;//(Rho[Connect(nodenumber,direction1)]+Rho[Connect(nodenumber,direction2)])*0.5;
		break;
	case 8:
		direction1=1;
		direction2=4;
		doubleTmpReturn=(Rhor[nodeIn.Get_connect()[direction1]]+Rhor[nodeIn.Get_connect()[direction2]])*0.5;//(Rho[Connect(nodenumber,direction1)]+Rho[Connect(nodenumber,direction2)])*0.5;
		break;
	}

	return doubleTmpReturn;
}
double D2Q9ColourFluid::Cal_RhoB_Corner(NodeCorner2D& nodeIn){

	unsigned int direction1,direction2;
	switch(nodeIn.Get_BcNormal())
	{
	case 5:
		direction1=1;
		direction2=2;
		doubleTmpReturn=(Rhob[nodeIn.Get_connect()[direction1]]+Rhob[nodeIn.Get_connect()[direction2]])*0.5;//Connect(nodenumber,direction1)]+Rho[Connect(nodenumber,direction2)])*0.5;
		break;
	case 6:
		direction1=2;
		direction2=3;
		doubleTmpReturn=(Rhob[nodeIn.Get_connect()[direction1]]+Rhob[nodeIn.Get_connect()[direction2]])*0.5;//(Rho[Connect(nodenumber,direction1)]+Rho[Connect(nodenumber,direction2)])*0.5;
		break;
	case 7:
		direction1=3;
		direction2=4;
		doubleTmpReturn=(Rhob[nodeIn.Get_connect()[direction1]]+Rhob[nodeIn.Get_connect()[direction2]])*0.5;//(Rho[Connect(nodenumber,direction1)]+Rho[Connect(nodenumber,direction2)])*0.5;
		break;
	case 8:
		direction1=1;
		direction2=4;
		doubleTmpReturn=(Rhob[nodeIn.Get_connect()[direction1]]+Rhob[nodeIn.Get_connect()[direction2]])*0.5;//(Rho[Connect(nodenumber,direction1)]+Rho[Connect(nodenumber,direction2)])*0.5;
		break;
	}

	return doubleTmpReturn;
}

/// Corner treat by Chih-Fung Ho, Cheng Chang, Kuen-Hau Lin and Chao-An Lin
/// Consistent Boundary Conditions for 2D and 3D Lattice Boltzmann Simulations
void D2Q9ColourFluid::ApplyCorner(NodeCorner2D& NodeIn){

	for (int i=0;i<9;i++)
		f_tmp[i]=f[0]->f[i][NodeIn.Get_index()];
	BC_corner(NodeIn.Get_BcNormal(),f_tmp,Cal_RhoR_Corner(NodeIn) , NodeIn.Get_UDef()[0],NodeIn.Get_UDef()[1]);
	for (int i=0;i<9;i++)
	f[0]->f[i][NodeIn.Get_index()]=f_tmp[i];

	for (int i=0;i<9;i++)
		f_tmp[i]=f[1]->f[i][NodeIn.Get_index()];
	BC_corner(NodeIn.Get_BcNormal(),f_tmp,Cal_RhoB_Corner(NodeIn) , NodeIn.Get_UDef()[0],NodeIn.Get_UDef()[1]);
	for (int i=0;i<9;i++)
		f[1]->f[i][NodeIn.Get_index()]=f_tmp[i];

}

void D2Q9ColourFluid::ApplyBounceBack(NodeCorner2D& NodeIn){
	double tmp=0;
	switch(NodeIn.Get_BcNormal())
			{
			case 5:
				f[0]->f[5][NodeIn.Get_index()]=f[0]->f[Opposite[5]][NodeIn.Get_index()];
				f[0]->f[1][NodeIn.Get_index()]=f[0]->f[Opposite[1]][NodeIn.Get_index()];
				f[0]->f[2][NodeIn.Get_index()]=f[0]->f[Opposite[2]][NodeIn.Get_index()];

				f[1]->f[5][NodeIn.Get_index()]=f[1]->f[Opposite[5]][NodeIn.Get_index()];
				f[1]->f[1][NodeIn.Get_index()]=f[1]->f[Opposite[1]][NodeIn.Get_index()];
				f[1]->f[2][NodeIn.Get_index()]=f[1]->f[Opposite[2]][NodeIn.Get_index()];
				if (NodeIn.stream()[3]==false)
				{

					tmp=f[0]->f[0][NodeIn.Get_index()]+f[0]->f[1][NodeIn.Get_index()]+f[0]->f[2][NodeIn.Get_index()]+f[0]->f[3][NodeIn.Get_index()]+f[0]->f[4][NodeIn.Get_index()]+f[0]->f[5][NodeIn.Get_index()]+f[0]->f[7][NodeIn.Get_index()];
					f[0]->f[6][NodeIn.Get_index()]=(Cal_RhoR_Corner(NodeIn)-tmp)/2;
					f[0]->f[8][NodeIn.Get_index()]=(Cal_RhoR_Corner(NodeIn)-tmp)/2;

					tmp=f[1]->f[0][NodeIn.Get_index()]+f[1]->f[1][NodeIn.Get_index()]+f[1]->f[2][NodeIn.Get_index()]+f[1]->f[3][NodeIn.Get_index()]+f[1]->f[4][NodeIn.Get_index()]+f[1]->f[5][NodeIn.Get_index()]+f[1]->f[7][NodeIn.Get_index()];
					f[1]->f[6][NodeIn.Get_index()]=(Cal_RhoB_Corner(NodeIn)-tmp)/2;
					f[1]->f[8][NodeIn.Get_index()]=(Cal_RhoB_Corner(NodeIn)-tmp)/2;
				}
				else
				{
					tmp=f[0]->f[6][NodeIn.Get_index()]+f[0]->f[8][NodeIn.Get_index()];
					f[0]->f[6][NodeIn.Get_index()]=tmp*0.5;
					f[0]->f[8][NodeIn.Get_index()]=tmp*0.5;

					tmp=f[1]->f[6][NodeIn.Get_index()]+f[1]->f[8][NodeIn.Get_index()];
					f[1]->f[6][NodeIn.Get_index()]=tmp*0.5;
					f[1]->f[8][NodeIn.Get_index()]=tmp*0.5;
				}
				break;
			case 6:
				f[0]->f[6][NodeIn.Get_index()]=f[0]->f[Opposite[6]][NodeIn.Get_index()];
				f[0]->f[2][NodeIn.Get_index()]=f[0]->f[Opposite[2]][NodeIn.Get_index()];
				f[0]->f[3][NodeIn.Get_index()]=f[0]->f[Opposite[3]][NodeIn.Get_index()];

				f[1]->f[6][NodeIn.Get_index()]=f[1]->f[Opposite[6]][NodeIn.Get_index()];
				f[1]->f[2][NodeIn.Get_index()]=f[1]->f[Opposite[2]][NodeIn.Get_index()];
				f[1]->f[3][NodeIn.Get_index()]=f[1]->f[Opposite[3]][NodeIn.Get_index()];
				if (NodeIn.stream()[1]==false)
				{
					tmp=f[0]->f[0][NodeIn.Get_index()]+f[0]->f[1][NodeIn.Get_index()]+f[0]->f[2][NodeIn.Get_index()]+f[0]->f[3][NodeIn.Get_index()]+f[0]->f[4][NodeIn.Get_index()]+f[0]->f[6][NodeIn.Get_index()]+f[0]->f[8][NodeIn.Get_index()];
					f[0]->f[5][NodeIn.Get_index()]=(Cal_RhoR_Corner(NodeIn)-tmp)/2;
					f[0]->f[7][NodeIn.Get_index()]=(Cal_RhoR_Corner(NodeIn)-tmp)/2;

					tmp=f[1]->f[0][NodeIn.Get_index()]+f[1]->f[1][NodeIn.Get_index()]+f[1]->f[2][NodeIn.Get_index()]+f[1]->f[3][NodeIn.Get_index()]+f[1]->f[4][NodeIn.Get_index()]+f[1]->f[6][NodeIn.Get_index()]+f[1]->f[8][NodeIn.Get_index()];
					f[1]->f[5][NodeIn.Get_index()]=(Cal_RhoB_Corner(NodeIn)-tmp)/2;
					f[1]->f[7][NodeIn.Get_index()]=(Cal_RhoB_Corner(NodeIn)-tmp)/2;
				}
				else
				{
					tmp=f[0]->f[7][NodeIn.Get_index()]+f[0]->f[5][NodeIn.Get_index()];
					f[0]->f[5][NodeIn.Get_index()]=tmp*0.5;
					f[0]->f[7][NodeIn.Get_index()]=tmp*0.5;

					tmp=f[1]->f[7][NodeIn.Get_index()]+f[1]->f[5][NodeIn.Get_index()];
					f[1]->f[5][NodeIn.Get_index()]=tmp*0.5;
					f[1]->f[7][NodeIn.Get_index()]=tmp*0.5;
				}
				break;
			case 7:
				f[0]->f[7][NodeIn.Get_index()]=f[0]->f[Opposite[7]][NodeIn.Get_index()];
				f[0]->f[3][NodeIn.Get_index()]=f[0]->f[Opposite[3]][NodeIn.Get_index()];
				f[0]->f[4][NodeIn.Get_index()]=f[0]->f[Opposite[4]][NodeIn.Get_index()];

				f[1]->f[7][NodeIn.Get_index()]=f[1]->f[Opposite[7]][NodeIn.Get_index()];
				f[1]->f[3][NodeIn.Get_index()]=f[1]->f[Opposite[3]][NodeIn.Get_index()];
				f[1]->f[4][NodeIn.Get_index()]=f[1]->f[Opposite[4]][NodeIn.Get_index()];
				if (NodeIn.stream()[1]==false)
				{
					tmp=f[0]->f[0][NodeIn.Get_index()]+f[0]->f[1][NodeIn.Get_index()]+f[0]->f[2][NodeIn.Get_index()]+f[0]->f[3][NodeIn.Get_index()]+f[0]->f[4][NodeIn.Get_index()]+f[0]->f[5][NodeIn.Get_index()]+f[0]->f[7][NodeIn.Get_index()];
					f[0]->f[6][NodeIn.Get_index()]=(Cal_RhoR_Corner(NodeIn)-tmp)/2;
					f[0]->f[8][NodeIn.Get_index()]=(Cal_RhoR_Corner(NodeIn)-tmp)/2;

					tmp=f[1]->f[0][NodeIn.Get_index()]+f[1]->f[1][NodeIn.Get_index()]+f[1]->f[2][NodeIn.Get_index()]+f[1]->f[3][NodeIn.Get_index()]+f[1]->f[4][NodeIn.Get_index()]+f[1]->f[5][NodeIn.Get_index()]+f[1]->f[7][NodeIn.Get_index()];
					f[1]->f[6][NodeIn.Get_index()]=(Cal_RhoB_Corner(NodeIn)-tmp)/2;
					f[1]->f[8][NodeIn.Get_index()]=(Cal_RhoB_Corner(NodeIn)-tmp)/2;
				}
				else
				{
					tmp=f[0]->f[6][NodeIn.Get_index()]+f[0]->f[8][NodeIn.Get_index()];
					f[0]->f[6][NodeIn.Get_index()]=tmp*0.5;
					f[0]->f[8][NodeIn.Get_index()]=tmp*0.5;

					tmp=f[1]->f[6][NodeIn.Get_index()]+f[1]->f[8][NodeIn.Get_index()];
					f[1]->f[6][NodeIn.Get_index()]=tmp*0.5;
					f[1]->f[8][NodeIn.Get_index()]=tmp*0.5;
				}
				break;
			case 8:
				f[0]->f[8][NodeIn.Get_index()]=f[0]->f[Opposite[8]][NodeIn.Get_index()];
				f[0]->f[1][NodeIn.Get_index()]=f[0]->f[Opposite[1]][NodeIn.Get_index()];
				f[0]->f[4][NodeIn.Get_index()]=f[0]->f[Opposite[4]][NodeIn.Get_index()];

				f[1]->f[8][NodeIn.Get_index()]=f[1]->f[Opposite[8]][NodeIn.Get_index()];
				f[1]->f[1][NodeIn.Get_index()]=f[1]->f[Opposite[1]][NodeIn.Get_index()];
				f[1]->f[4][NodeIn.Get_index()]=f[1]->f[Opposite[4]][NodeIn.Get_index()];
				if (NodeIn.stream()[2]==false)
				{
					tmp=f[0]->f[0][NodeIn.Get_index()]+f[0]->f[1][NodeIn.Get_index()]+f[0]->f[2][NodeIn.Get_index()]+f[0]->f[3][NodeIn.Get_index()]+f[0]->f[4][NodeIn.Get_index()]+f[0]->f[6][NodeIn.Get_index()]+f[0]->f[8][NodeIn.Get_index()];
					f[0]->f[5][NodeIn.Get_index()]=(Cal_RhoR_Corner(NodeIn)-tmp)/2;
					f[0]->f[7][NodeIn.Get_index()]=(Cal_RhoR_Corner(NodeIn)-tmp)/2;

					tmp=f[1]->f[0][NodeIn.Get_index()]+f[1]->f[1][NodeIn.Get_index()]+f[1]->f[2][NodeIn.Get_index()]+f[1]->f[3][NodeIn.Get_index()]+f[1]->f[4][NodeIn.Get_index()]+f[1]->f[6][NodeIn.Get_index()]+f[1]->f[8][NodeIn.Get_index()];
					f[1]->f[5][NodeIn.Get_index()]=(Cal_RhoB_Corner(NodeIn)-tmp)/2;
					f[1]->f[7][NodeIn.Get_index()]=(Cal_RhoB_Corner(NodeIn)-tmp)/2;
				}
				else
				{
					tmp=f[0]->f[5][NodeIn.Get_index()]+f[0]->f[7][NodeIn.Get_index()];
					f[0]->f[5][NodeIn.Get_index()]=tmp*0.5;
					f[0]->f[7][NodeIn.Get_index()]=tmp*0.5;

					tmp=f[1]->f[5][NodeIn.Get_index()]+f[1]->f[7][NodeIn.Get_index()];
					f[1]->f[5][NodeIn.Get_index()]=tmp*0.5;
					f[1]->f[7][NodeIn.Get_index()]=tmp*0.5;
				}
				break;
			default:
				std::cerr<<"Direction corner bounce back not found. x:"<<NodeIn.get_x()<<" y: "<<NodeIn.get_y()<<std::endl;
				break;
			}
}
///Bounceback Wall treatment
void D2Q9ColourFluid::ApplyBounceBack(NodeWall2D& NodeIn){

	switch(NodeIn.Get_BcNormal())
			{
			case 2:
				f[0]->f[2][NodeIn.Get_index()]=f[0]->f[Opposite[2]][NodeIn.Get_index()];
				f[0]->f[5][NodeIn.Get_index()]=f[0]->f[Opposite[5]][NodeIn.Get_index()];
				f[0]->f[6][NodeIn.Get_index()]=f[0]->f[Opposite[6]][NodeIn.Get_index()];

				f[1]->f[2][NodeIn.Get_index()]=f[1]->f[Opposite[2]][NodeIn.Get_index()];
				f[1]->f[5][NodeIn.Get_index()]=f[1]->f[Opposite[5]][NodeIn.Get_index()];
				f[1]->f[6][NodeIn.Get_index()]=f[1]->f[Opposite[6]][NodeIn.Get_index()];
				break;
			case 4:
				f[0]->f[4][NodeIn.Get_index()]=f[0]->f[Opposite[4]][NodeIn.Get_index()];
				f[0]->f[7][NodeIn.Get_index()]=f[0]->f[Opposite[7]][NodeIn.Get_index()];
				f[0]->f[8][NodeIn.Get_index()]=f[0]->f[Opposite[8]][NodeIn.Get_index()];

				f[1]->f[4][NodeIn.Get_index()]=f[1]->f[Opposite[4]][NodeIn.Get_index()];
				f[1]->f[7][NodeIn.Get_index()]=f[1]->f[Opposite[7]][NodeIn.Get_index()];
				f[1]->f[8][NodeIn.Get_index()]=f[1]->f[Opposite[8]][NodeIn.Get_index()];
				break;
			case 1:
				f[0]->f[1][NodeIn.Get_index()]=f[0]->f[Opposite[1]][NodeIn.Get_index()];
				f[0]->f[5][NodeIn.Get_index()]=f[0]->f[Opposite[5]][NodeIn.Get_index()];
				f[0]->f[8][NodeIn.Get_index()]=f[0]->f[Opposite[8]][NodeIn.Get_index()];

				f[1]->f[1][NodeIn.Get_index()]=f[1]->f[Opposite[1]][NodeIn.Get_index()];
				f[1]->f[5][NodeIn.Get_index()]=f[1]->f[Opposite[5]][NodeIn.Get_index()];
				f[1]->f[8][NodeIn.Get_index()]=f[1]->f[Opposite[8]][NodeIn.Get_index()];
				break;
			case 3:
				f[0]->f[3][NodeIn.Get_index()]=f[0]->f[Opposite[3]][NodeIn.Get_index()];
				f[0]->f[6][NodeIn.Get_index()]=f[0]->f[Opposite[6]][NodeIn.Get_index()];
				f[0]->f[7][NodeIn.Get_index()]=f[0]->f[Opposite[7]][NodeIn.Get_index()];


				f[1]->f[3][NodeIn.Get_index()]=f[1]->f[Opposite[3]][NodeIn.Get_index()];
				f[1]->f[6][NodeIn.Get_index()]=f[1]->f[Opposite[6]][NodeIn.Get_index()];
				f[1]->f[7][NodeIn.Get_index()]=f[1]->f[Opposite[7]][NodeIn.Get_index()];				break;
			default:
				std::cerr<<"Direction wall bounce back not found. Index: "<<NodeIn.Get_index()<<" x: "<<NodeIn.get_x()<<" y: "<<NodeIn.get_y()<<" direction not found: "<<NodeIn.Get_BcNormal()<<std::endl;
				break;
			}
}
