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
	PtrColourGrad=0;
	PtrColourGradBc=0;
	PtrColourGradCorner=0;
	PtrCollision=0;
	D_tmp=0;
	PtrD_tmp=0;
	PtrRecolour=0;
	PtrMacro=0;
	beta=0.7;
	Rho_limiter=1.E-8;
	A1=0;
	A2=0;
	F=0;
	G=0;
	G_Norm=0;
	tension=0;


}

D2Q9ColourFluid::~D2Q9ColourFluid() {
	delete [] Rhor;
	delete [] Rhob;
}

D2Q9ColourFluid::D2Q9ColourFluid(MultiBlock* MultiBlock__,ParallelManager* parallel__,WriterManager* Writer__, Parameters* Parameters__,InitLBM& ini){
	InitD2Q9TwoPhases(MultiBlock__,parallel__,Writer__, Parameters__, ini);
	InitMultiphase(ini);
	Set_PointersOnFunctions();
}
void D2Q9ColourFluid::Set_PointersOnFunctions(){
	Set_Collide();
	Set_Colour_gradient();
	Set_Recolouring();
	Set_Macro();
}
//Std2D,Std2DLocal,Std2DBody,Std2DNonCstTau,Std2DNonCstTauLocal,Std2DNonCstTauBody
void D2Q9ColourFluid::Select_Colour_Operator(ColourOperatorType OperatorType_){
	switch(OperatorType_)
	{
	case ::SurfaceForce:
		PtrCollision=&D2Q9ColourFluid::Collision_SurfaceForce;
		break;
	case ::Gunstensen:
		PtrCollision=&D2Q9ColourFluid::Collision_Gunstensen;
		break;
	default:
		std::cerr<<" Colour operator not found."<<std::endl;
		break;
	}
}
void D2Q9ColourFluid::Set_Collide(){
if(PtrParameters->Get_ColourOperatorType()== ::SurfaceForce && PtrParameters->Get_UserForceType()== ::LocalForce)
{
	std::cout<<" Warming: The local user force will be ignored. The colour model with the surface force is incompatible with a local user force"<<std::endl;
}
	if(PtrParameters->Get_FluidType()==Newtonian)
	{
		if(PtrParameters->Get_ColourOperatorType()== ::SurfaceForce)
			{Select_Collide_2D(Std2DBody);Select_Colour_Operator(PtrParameters->Get_ColourOperatorType());}
		else if(PtrParameters->Get_UserForceType()== ::LocalForce)
			{Select_Collide_2D(Std2DLocal);Select_Colour_Operator(PtrParameters->Get_ColourOperatorType());}
		else if(PtrParameters->Get_UserForceType()== ::BodyForce)
			{Select_Collide_2D(Std2DBody);Select_Colour_Operator(PtrParameters->Get_ColourOperatorType());}
		else
			{Select_Collide_2D(Std2D);Select_Colour_Operator(PtrParameters->Get_ColourOperatorType());}
	}
	else
	{
		if(PtrParameters->Get_ColourOperatorType()== ::SurfaceForce)
			{Select_Collide_2D(Std2DNonCstTauBody);Select_Colour_Operator(PtrParameters->Get_ColourOperatorType());}
		else if(PtrParameters->Get_UserForceType()== ::LocalForce)
			{Select_Collide_2D(Std2DNonCstTauLocal);Select_Colour_Operator(PtrParameters->Get_ColourOperatorType());}
		else if(PtrParameters->Get_UserForceType()== ::BodyForce)
			{Select_Collide_2D(Std2DNonCstTauBody);Select_Colour_Operator(PtrParameters->Get_ColourOperatorType());}
		else
			{Select_Collide_2D(Std2DNonCstTau);Select_Colour_Operator(PtrParameters->Get_ColourOperatorType());}
	}

}
void D2Q9ColourFluid::Set_Colour_gradient(){

		switch(PtrParameters->Get_ColourGradType())
		{
		case Gunstensen:
			if(PtrParameters->Get_ColourOperatorType()== ::SurfaceForce)
			{
				std::cout<<" Force colour gradient to Density Normal gradient due to Surface force model."<<std::endl;
				PtrColourGrad =&D2Q9ColourFluid::Colour_gradient_DensityNormalGrad;
				PtrColourGradBc =&D2Q9ColourFluid::Colour_gradient_DensityNormalGradBc;
				PtrColourGradCorner =&D2Q9ColourFluid::Colour_gradient_DensityNormalGradCorner;
			}
			else
			{
			PtrColourGrad =&D2Q9ColourFluid::Colour_gradient_Gunstensen;
			PtrColourGradBc =&D2Q9ColourFluid::Colour_gradient_Gunstensen;
			PtrColourGradCorner =&D2Q9ColourFluid::Colour_gradient_Gunstensen;
			}
			break;
		case DensityGrad:
			PtrColourGrad =&D2Q9ColourFluid::Colour_gradient_DensityGrad;
			PtrColourGradBc =&D2Q9ColourFluid::Colour_gradient_DensityGradBc;
			PtrColourGradCorner =&D2Q9ColourFluid::Colour_gradient_DensityGradCorner;
			break;
		case DensityNormalGrad:
			PtrColourGrad =&D2Q9ColourFluid::Colour_gradient_DensityNormalGrad;
			PtrColourGradBc =&D2Q9ColourFluid::Colour_gradient_DensityNormalGradBc;
			PtrColourGradCorner =&D2Q9ColourFluid::Colour_gradient_DensityNormalGradCorner;
			break;
		}
}
void D2Q9ColourFluid::Set_Macro(){
	switch(PtrParameters->Get_ColourGradType())
	{
	case Gunstensen:
		if(PtrParameters->Get_ColourOperatorType()== ::SurfaceForce )
			PtrMacro=&D2Q9ColourFluid::MacroVariablesWithNormalDensityAndForce;
		else
			if(PtrParameters->Get_NormalDensityOutput())
				PtrMacro=&D2Q9ColourFluid::MacroVariablesWithNormalDensity;
			else
				PtrMacro=&D2Q9ColourFluid::MacroVariables;
		break;
	case DensityGrad:
		if(PtrParameters->Get_ColourOperatorType()== ::SurfaceForce)
			PtrMacro=&D2Q9ColourFluid::MacroVariablesWithNormalDensityAndForce;
		else
			if(PtrParameters->Get_NormalDensityOutput())
				PtrMacro=&D2Q9ColourFluid::MacroVariablesWithNormalDensity;
			else
				PtrMacro=&D2Q9ColourFluid::MacroVariables;
		break;
	case DensityNormalGrad:
		if(PtrParameters->Get_ColourOperatorType()== ::SurfaceForce)
			PtrMacro=&D2Q9ColourFluid::MacroVariablesWithNormalDensityAndForce;
		else
			PtrMacro=&D2Q9ColourFluid::MacroVariablesWithNormalDensity;
		break;
	}
//	PtrMacro=&D2Q9ColourFluid::MacroVariables;
//	PtrMacro=&D2Q9ColourFluid::MacroVariablesWithNormalDensity;
}
void D2Q9ColourFluid::Set_Recolouring(){
	PtrRecolour=&D2Q9ColourFluid::Recolouring_Latva;
}
void D2Q9ColourFluid::InitMultiphase(InitLBM& ini){
//Initialise parameters
	beta=PtrParameters->Get_Beta();
	A1=PtrParameters->Get_A1();
	A2=PtrParameters->Get_A2();
	tension=PtrParameters->Get_SurfaceTension();

// Initialise only the Colour Fluid part.
/*	G[0]=V1[0];
	G[1]=V1[1];
	F[0]=V2[0];
	F[1]=V2[1];*/
	G=V1;
	F=V2;
	G_Norm=new double [nbnodes_total];
/*	G=new double* [2];
	G[0]=new double [nbnodes_total];
	G[1]=new double [nbnodes_total];*/
	for(int i=0;i<nbnodes_total;i++)
	{
		G[0][i]=0;
		G[1][i]=0;
		G_Norm[i]=0;
	}

//	F=new double* [nbnodes_total];
	for(int i=0;i<nbnodes_total;i++)
	{
//		F[i]=new double [2];
		F[0][i]=0;
		F[1][i]=0;
	}
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
		ini.IniDomainTwoPhases(parallel->getRank(),NodeArrays->NodeInterior[j],0, NodeArrays->NodeInterior[j].Get_index(),pos,Rho[NodeArrays->NodeInterior[j].Get_index()],U_,alpha);
// Initialise the blue and red densities
		Rhor[NodeArrays->NodeInterior[j].Get_index()]=alpha*Rho[NodeArrays->NodeInterior[j].Get_index()];//PtrParameters->Get_Rho_1();
		Rhob[NodeArrays->NodeInterior[j].Get_index()]=(1- alpha) *Rho[NodeArrays->NodeInterior[j].Get_index()];//PtrParameters->Get_Rho_2();
		RhoN[NodeArrays->NodeInterior[j].Get_index()]=(Rhor[NodeArrays->NodeInterior[j].Get_index()]-Rhob[NodeArrays->NodeInterior[j].Get_index()])/Rho[NodeArrays->NodeInterior[j].Get_index()];
	}

	for (int j=0;j<NodeArrays->NodeCorner.size();j++)
	{
// Get position
		pos[0]=NodeArrays->NodeCorner[j].get_x();
		pos[1]=NodeArrays->NodeCorner[j].get_y();
// Get initialise value from the user
		ini.IniDomainTwoPhases(parallel->getRank(),NodeArrays->NodeCorner[j],0, NodeArrays->NodeCorner[j].Get_index(),pos,Rho[NodeArrays->NodeCorner[j].Get_index()],U_,alpha);
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
		ini.IniDomainTwoPhases(parallel->getRank(),NodeArrays->NodeGlobalCorner[j],0, NodeArrays->NodeGlobalCorner[j].Get_index(),pos,Rho[NodeArrays->NodeGlobalCorner[j].Get_index()],U_,alpha);
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
		ini.IniDomainTwoPhases(parallel->getRank(),NodeArrays->NodeVelocity[j],0, NodeArrays->NodeVelocity[j].Get_index(),pos,Rho[NodeArrays->NodeVelocity[j].Get_index()],U_,alpha);
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
		ini.IniDomainTwoPhases(parallel->getRank(),NodeArrays->NodePressure[j],0, NodeArrays->NodePressure[j].Get_index(),pos,Rho[NodeArrays->NodePressure[j].Get_index()],U_,alpha);
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
		ini.IniDomainTwoPhases(parallel->getRank(),NodeArrays->NodeWall[j],0, NodeArrays->NodeWall[j].Get_index(),pos,Rho[NodeArrays->NodeWall[j].Get_index()],U_,alpha);
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
		ini.IniDomainTwoPhases(parallel->getRank(),NodeArrays->NodeSpecialWall[j],0, NodeArrays->NodeSpecialWall[j].Get_index(),pos,Rho[NodeArrays->NodeSpecialWall[j].Get_index()],U_,alpha);
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
		ini.IniDomainTwoPhases(parallel->getRank(),NodeArrays->NodeSymmetry[j],0, NodeArrays->NodeSymmetry[j].Get_index(),pos,Rho[NodeArrays->NodeSymmetry[j].Get_index()],U_,alpha);
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
		ini.IniDomainTwoPhases(parallel->getRank(),NodeArrays->NodeGhost[j],0, NodeArrays->NodeGhost[j].Get_index(),pos,Rho[NodeArrays->NodeGhost[j].Get_index()],U_,alpha);
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
		ini.IniDomainTwoPhases(parallel->getRank(),NodeArrays->NodeSolid[j],0, NodeArrays->NodeSolid[j].Get_index(),pos,Rho[NodeArrays->NodeSolid[j].Get_index()],U_,alpha);
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
	if(parallel->getSize()>1)
		SyncMacroVarToGhost();
	Colour_gradient();
	ColourFluid_Collision();
	Writer->Write_Output(it);

	if(parallel->getSize()>1)
	{
		for (int i=1;i<NbStep+1;i++)
		{
			Colour_gradient();
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
			Colour_gradient();
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
///Calculate the colour gradient by Gunstensen formulation, density gradient or normal density gradient
void D2Q9ColourFluid::Colour_gradient(){
	int idx_tmp;int normal_interior=0;
	for (int j=0;j<NodeArrays->NodeInterior.size();j++)
	{
	// Common variables
		idx_tmp=NodeArrays->NodeInterior[j].Get_index();
	// Calculate gradients
		(this->*PtrColourGrad)(idx_tmp,NodeArrays->NodeInterior[j].Get_connect(),normal_interior);
	}

	for (int j=0;j<NodeArrays->NodeCorner.size();j++)
	{
	// Common variables
		idx_tmp=NodeArrays->NodeCorner[j].Get_index();
	// Calculate gradients
		(this->*PtrColourGradCorner)(idx_tmp,NodeArrays->NodeCorner[j].Get_connect(),NodeArrays->NodeCorner[j].Get_BcNormal());
	}

	for (int j=0;j<NodeArrays->NodeGlobalCorner.size();j++)
	{
	// Common variables
		idx_tmp=NodeArrays->NodeGlobalCorner[j].Get_index();
	// Calculate gradients
		(this->*PtrColourGradBc)(idx_tmp,NodeArrays->NodeGlobalCorner[j].Get_connect(),NodeArrays->NodeGlobalCorner[j].Get_BcNormal());
	}

	for (int j=0;j<NodeArrays->NodeVelocity.size();j++)
	{
	// Common variables
		idx_tmp=NodeArrays->NodeVelocity[j].Get_index();
	// Calculate gradients
		(this->*PtrColourGradBc)(idx_tmp,NodeArrays->NodeVelocity[j].Get_connect(),NodeArrays->NodeVelocity[j].Get_BcNormal());
	}

	for (int j=0;j<NodeArrays->NodePressure.size();j++)
	{
	// Common variables
		idx_tmp=NodeArrays->NodePressure[j].Get_index();
	// Calculate gradients
		(this->*PtrColourGradBc)(idx_tmp,NodeArrays->NodePressure[j].Get_connect(),NodeArrays->NodePressure[j].Get_BcNormal());
	}

	for (int j=0;j<NodeArrays->NodeWall.size();j++)
	{
	// Common variables
		idx_tmp=NodeArrays->NodeWall[j].Get_index();
	// Calculate gradients
		(this->*PtrColourGradBc)(idx_tmp,NodeArrays->NodeWall[j].Get_connect(),NodeArrays->NodeWall[j].Get_BcNormal());
	}

	for (int j=0;j<NodeArrays->NodeSymmetry.size();j++)
	{
	// Common variables
		idx_tmp=NodeArrays->NodeSymmetry[j].Get_index();
	// Calculate gradients
		(this->*PtrColourGradBc)(idx_tmp,NodeArrays->NodeSymmetry[j].Get_connect(),NodeArrays->NodeSymmetry[j].Get_BcNormal());
	}
}
void D2Q9ColourFluid::ColourFluid_Collision()
{
	double Ak=0.65;
//	double F_tmp[2];
//	double G_Norm=0;
	int idx_tmp;
	double wtmp=0;
	int normal_interior=0;
//	double *U_tmp, *V_tmp, *F_tmp;
	double  InvTau_;
	double fi_tmp[9];
	InvTau_=InvTau; //Tmp
//	std::cout<<" Nd Interior: "<<NodeArrays->NodeInterior.size()<<std::endl;
	for (int j=0;j<NodeArrays->NodeInterior.size();j++)
	{
	// Common variables
		idx_tmp=NodeArrays->NodeInterior[j].Get_index();
	//Model
		(this->*PtrCollision)(idx_tmp, NodeArrays->NodeInterior[j].Get_connect(),normal_interior,&fi_tmp[0]);
		(this->*PtrRecolour)(idx_tmp,&fi_tmp[0]);
	}

	for (int j=0;j<NodeArrays->NodeCorner.size();j++)
	{
	// Common variables
		idx_tmp=NodeArrays->NodeCorner[j].Get_index();
	//Model
		(this->*PtrCollision)(idx_tmp, NodeArrays->NodeCorner[j].Get_connect(),NodeArrays->NodeCorner[j].Get_BcNormal(),&fi_tmp[0]);
		(this->*PtrRecolour)(idx_tmp,&fi_tmp[0]);
	}

	for (int j=0;j<NodeArrays->NodeGlobalCorner.size();j++)
	{
	// Common variables
		idx_tmp=NodeArrays->NodeGlobalCorner[j].Get_index();
	//Calculate Norms
		G_Norm[idx_tmp]=0;
//		G_Norm=std::sqrt(G[idx_tmp][0]*G[idx_tmp][0]+G[idx_tmp][1]*G[idx_tmp][1]);
	//Model
		(this->*PtrCollision)(idx_tmp, NodeArrays->NodeGlobalCorner[j].Get_connect(),NodeArrays->NodeGlobalCorner[j].Get_BcNormal(),&fi_tmp[0]);
		(this->*PtrRecolour)(idx_tmp,&fi_tmp[0]);
	}

	for (int j=0;j<NodeArrays->NodeVelocity.size();j++)
	{
	// Common variables
		idx_tmp=NodeArrays->NodeVelocity[j].Get_index();
	//Model
		(this->*PtrCollision)(idx_tmp, NodeArrays->NodeVelocity[j].Get_connect(),NodeArrays->NodeVelocity[j].Get_BcNormal(),&fi_tmp[0]);
		(this->*PtrRecolour)(idx_tmp,&fi_tmp[0]);
	}

	for (int j=0;j<NodeArrays->NodePressure.size();j++)
	{
	// Common variables
		idx_tmp=NodeArrays->NodePressure[j].Get_index();

		(this->*PtrCollision)(idx_tmp, NodeArrays->NodePressure[j].Get_connect(),NodeArrays->NodePressure[j].Get_BcNormal(),&fi_tmp[0]);
	//Model
		(this->*PtrRecolour)(idx_tmp,&fi_tmp[0]);
	}

	for (int j=0;j<NodeArrays->NodeWall.size();j++)
	{
	// Common variables
		idx_tmp=NodeArrays->NodeWall[j].Get_index();
	//Calculate Norms
		G_Norm[idx_tmp]=0;
//		G_Norm=std::sqrt(G[idx_tmp][0]*G[idx_tmp][0]+G[idx_tmp][1]*G[idx_tmp][1]);
	//Model
		(this->*PtrCollision)(idx_tmp, NodeArrays->NodeWall[j].Get_connect(),NodeArrays->NodeWall[j].Get_BcNormal(),&fi_tmp[0]);
		(this->*PtrRecolour)(idx_tmp,&fi_tmp[0]);
	}

	for (int j=0;j<NodeArrays->NodeSymmetry.size();j++)
	{
	// Common variables
		idx_tmp=NodeArrays->NodeSymmetry[j].Get_index();

	//Model
		(this->*PtrCollision)(idx_tmp, NodeArrays->NodeSymmetry[j].Get_connect(),NodeArrays->NodeSymmetry[j].Get_BcNormal(),&fi_tmp[0]);
		(this->*PtrRecolour)(idx_tmp,&fi_tmp[0]);

	}
}
/*void D2Q9ColourFluid::CollideD2Q9ColourFluid(int & direction, double & fi,double &rho,double*  F_tmp,double & F_Norm, double & InvTau_, double &u, double &v){
	CollideLowOrder::Collide_2D(direction,fi,rho, u, v,F_tmp,F_Norm);
	if(F_Norm>0)
		fi+=TwoPhase_Collision_operator(direction, F_tmp, F_Norm);
}*/
/*double D2Q9ColourFluid::TwoPhase_Collision_operator(int & i, double* F_tmp, double & F_Norm){
	double EiGperGNorm=(F_tmp[0]* Ei[i][0]+F_tmp[1]* Ei[i][1])/F_Norm;
	 return A1*0.5*F_Norm*(EiGperGNorm*EiGperGNorm-3/4);
}*/



void D2Q9ColourFluid::Colour_gradient_Gunstensen(int & nodenumber, int* connect,int & normal){

	for (int j=0;j<2;j++)
	{
		//F[j]=0;
		G[j][nodenumber]=(Rhor[nodenumber]-Rhob[nodenumber])*Ei[0][j];

		for (int k=1; k<nbvelo;k++)
		{
			G[j][nodenumber]+=(Rhor[connect[k]]-Rhob[connect[k]])*Ei[k][j];
		}
	}
	G_Norm[nodenumber]=std::sqrt(G[0][nodenumber]*G[0][nodenumber]+G[1][nodenumber]*G[1][nodenumber]);
}
void D2Q9ColourFluid::Colour_gradient_DensityGrad(int & nodenumber, int* connect,int & normal){
	double tmp[2];
	Grad(&tmp[0],&Rho[0],connect,normal);
	G[0][nodenumber]=tmp[0];G[1][nodenumber]=tmp[1];
	G_Norm[nodenumber]=std::sqrt(G[0][nodenumber]*G[0][nodenumber]+G[1][nodenumber]*G[1][nodenumber]);
}
void D2Q9ColourFluid::Colour_gradient_DensityGradBc(int & nodenumber, int* connect,int & normal){
	double tmp[2];
	GradBc(&tmp[0],&Rho[0],connect,normal);
	G[0][nodenumber]=tmp[0];G[1][nodenumber]=tmp[1];
	G_Norm[nodenumber]=std::sqrt(G[0][nodenumber]*G[0][nodenumber]+G[1][nodenumber]*G[1][nodenumber]);
}
void D2Q9ColourFluid::Colour_gradient_DensityGradCorner(int & nodenumber, int* connect,int & normal){
	double tmp[2];
	GradCorner(&tmp[0],&Rho[0],connect,normal);
	G[0][nodenumber]=tmp[0];G[1][nodenumber]=tmp[1];
	G_Norm[nodenumber]=std::sqrt(G[0][nodenumber]*G[0][nodenumber]+G[1][nodenumber]*G[1][nodenumber]);
}
void D2Q9ColourFluid::Colour_gradient_DensityNormalGrad(int & nodenumber, int* connect,int & normal){
	double tmp[2];
//	if(nodenumber==4600)
//		std::cout<<"G Norm is: "<<G_Norm[nodenumber]<<std::endl;
	Grad(&tmp[0],&RhoN[0],connect,normal);
	G[0][nodenumber]=tmp[0];G[1][nodenumber]=tmp[1];
	G_Norm[nodenumber]=std::sqrt(G[0][nodenumber]*G[0][nodenumber]+G[1][nodenumber]*G[1][nodenumber]);

	//Normalise grad(RhoN)
	if(G_Norm[nodenumber]>0)
	{
		G[0][nodenumber]=G[0][nodenumber]/G_Norm[nodenumber];
		G[1][nodenumber]=G[1][nodenumber]/G_Norm[nodenumber];
	}
}
void D2Q9ColourFluid::Colour_gradient_DensityNormalGradBc(int & nodenumber, int* connect,int & normal){
	double tmp[2];
	GradBc(&tmp[0],&RhoN[0],connect,normal);
	G[0][nodenumber]=tmp[0];G[1][nodenumber]=tmp[1];
	G_Norm[nodenumber]=std::sqrt(G[0][nodenumber]*G[0][nodenumber]+G[1][nodenumber]*G[1][nodenumber]);
	//Normalise grad(RhoN)
	if(G_Norm[nodenumber]>0)
	{
		G[0][nodenumber]=G[0][nodenumber]/G_Norm[nodenumber];
		G[1][nodenumber]=G[1][nodenumber]/G_Norm[nodenumber];
	}
}
void D2Q9ColourFluid::Colour_gradient_DensityNormalGradCorner(int & nodenumber, int* connect,int & normal){
	double tmp[2];
	GradCorner(&tmp[0],&RhoN[0],connect,normal);
	G[0][nodenumber]=tmp[0];G[1][nodenumber]=tmp[1];
	G_Norm[nodenumber]=std::sqrt(G[0][nodenumber]*G[0][nodenumber]+G[1][nodenumber]*G[1][nodenumber]);
	//Normalise grad(RhoN)
	if(G_Norm[nodenumber]>0)
	{
		G[0][nodenumber]=G[0][nodenumber]/G_Norm[nodenumber];
		G[1][nodenumber]=G[1][nodenumber]/G_Norm[nodenumber];
	}
}
//(idx_tmp,i,fi_tmp[i],&G[idx_tmp][0]);
void D2Q9ColourFluid::Recolouring_Latva(int & nodenumber, double * fi_tmp){
//(idx_tmp,i,fi_tmp,&F[0],F_Norm)
/*	if(CosPhi(i,F,F_Norm)<-1 ||CosPhi(i,F,F_Norm)>1)
		std::cout<<"CosPhi: "<<CosPhi(i,F,F_Norm)<<std::endl;*/

		f[0]->f[0][nodenumber]=Rhor[nodenumber]*fi_tmp[0]/Rho[nodenumber];
		for(int i=1;i<nbvelo;i++)
		{
		f[0]->f[i][nodenumber]=Rhor[nodenumber]*fi_tmp[i]/Rho[nodenumber]
								+beta*omega[i]*Rhor[nodenumber]*Rhob[nodenumber]*CosPhi(nodenumber,i,G_Norm[nodenumber])/Rho[nodenumber];
		}
		for(int i=0;i<nbvelo;i++)
		{
		if(Rhor[nodenumber]<=Rho_limiter)
			f[0]->f[i][nodenumber]=0;
		if(Rhob[nodenumber]<=Rho_limiter)
		{
			f[0]->f[i][nodenumber]=fi_tmp[i];
			f[1]->f[i][nodenumber]=0;
		}
		else
			f[1]->f[i][nodenumber]=fi_tmp[i]-f[0]->f[i][nodenumber];
		}

}
//CosPhi(nodenumber,i,G_Norm[nodenumber])
double D2Q9ColourFluid::CosPhi(int nodenumber, int & direction,double & F_Norm){
		return (G[0][nodenumber]* Ei[direction][0]+G[1][nodenumber]* Ei[direction][1]);///(F_Norm*std::sqrt(Ei[direction][0]*Ei[direction][0]+Ei[direction][1]*Ei[direction][1]));
}
/*double D2Q9ColourFluid::TwoPhase_Collision_operator(int & nodenumber, int & i, double & Ak, double* F_tmp, double & F_Norm){
 return Ak*0.5*F_Norm*(((F_tmp[0]*Ei[i][0]+F_tmp[1]*Ei[i][1])/F_Norm)*((F_tmp[0]*Ei[i][0]+F_tmp[1]*Ei[i][1])/F_Norm)-3/4);
}*/
double& D2Q9ColourFluid::Collision_operator_Gunstensen(int & i, int & nodenumber, double Ak){
	if(G_Norm[nodenumber]>0)
		D_tmp=Ak*0.5*G_Norm[nodenumber]*(((G[0][nodenumber]*Ei[i][0]+G[1][nodenumber]*Ei[i][1])/G_Norm[nodenumber])*((G[0][nodenumber]*Ei[i][0]+G[1][nodenumber]*Ei[i][1])/G_Norm[nodenumber])-3/4);
	return D_tmp;
}
void D2Q9ColourFluid::Collision_Gunstensen(int & nodenumber, int* connect,int & normal,double* fi){
	double InvTau_tmp=InvTau; int i0=0;
	fi[0]=f[0]->f[0][nodenumber]+f[1]->f[0][nodenumber];
	DVec_2D_tmp[0]=0;DVec_2D_tmp[1]=0;
	Collide_2D(i0, fi[0],Rho[nodenumber], U[0][nodenumber], U[1][nodenumber], DVec_2D_tmp[0],DVec_2D_tmp[1], InvTau_tmp);
	for (int i=1;i<9;i++)
	{
		//Save the mixture distribution for recolouring
		fi[i]=f[0]->f[i][nodenumber]+f[1]->f[i][nodenumber];
		Collide_2D(i, fi[i],Rho[nodenumber], U[0][nodenumber], U[1][nodenumber], Collision_operator_Gunstensen(i, nodenumber, A1),DVec_2D_tmp[1], InvTau_tmp);
	}
}
void D2Q9ColourFluid::Collision_SurfaceForce(int & nodenumber, int* connect,int & normal,double* fi){
	double InvTau_tmp=InvTau;

	if(G_Norm[nodenumber]>0)
	{
		SurfaceForce(nodenumber,connect,normal,F[0][nodenumber],F[1][nodenumber]);
	}
	else
		{F[0][nodenumber]=0;F[1][nodenumber]=0;}
	for (int i=0;i<9;i++)
	{
		//Save the mixture distribution for recolouring
		fi[i]=f[0]->f[i][nodenumber]+f[1]->f[i][nodenumber];
		Collide_2D(i, fi[i],Rho[nodenumber], U[0][nodenumber], U[1][nodenumber], F[0][nodenumber],F[1][nodenumber], InvTau_tmp);
	}
}

void D2Q9ColourFluid::SurfaceForce(int & nodenumber, int* connect,int & normal,double & Fx,double & Fy){
	D_tmp=0.5*tension*Curvature(nodenumber,connect,normal)*G_Norm[nodenumber]; //G_Norm is to get back the density gradient and not the normalise one
	Fx=D_tmp*G[0][nodenumber];
	Fy=D_tmp*G[1][nodenumber];

}
double D2Q9ColourFluid::Curvature(int & nodenumber, int* connect,int & normal){
	double tmp_0[2],tmp_1[2];
	Grad(&tmp_0[0],&G[0][0],connect,normal);//d(Gx)/d(x)
	Grad(&tmp_1[0],&G[1][0],connect,normal);//d(Gy)/d(y)
	return -(tmp_0[0]+tmp_1[1]);//Dot product

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
			(this->*PtrMacro)(NodeArrays->NodeInterior[j].Get_index());
		}
		for (int j=0;j<NodeArrays->NodeCorner.size();j++)
		{
			(this->*PtrMacro)(NodeArrays->NodeCorner[j].Get_index());
		}
		for (int j=0;j<NodeArrays->NodeGlobalCorner.size();j++)
		{
			(this->*PtrMacro)(NodeArrays->NodeGlobalCorner[j].Get_index());
		}
		for (int j=0;j<NodeArrays->NodeVelocity.size();j++)
		{
			(this->*PtrMacro)(NodeArrays->NodeVelocity[j].Get_index());
		}
		for (int j=0;j<NodeArrays->NodePressure.size();j++)
		{
			(this->*PtrMacro)(NodeArrays->NodePressure[j].Get_index());
		}
		for (int j=0;j<NodeArrays->NodeWall.size();j++)
		{
			(this->*PtrMacro)(NodeArrays->NodeWall[j].Get_index());
		}
		for (int j=0;j<NodeArrays->NodeSymmetry.size();j++)
		{
			(this->*PtrMacro)(NodeArrays->NodeSymmetry[j].Get_index());
		}

}

void D2Q9ColourFluid::MacroVariables(int& idx){

		U[0][idx]=0;
		U[1][idx]=0;
		Rhor[idx]=0,Rhob[idx]=0;
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
void D2Q9ColourFluid::MacroVariablesWithNormalDensity(int& idx){
	U[0][idx]=0;
	U[1][idx]=0;
	Rhor[idx]=0,Rhob[idx]=0;
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
	RhoN[idx]=(Rhor[idx]-Rhob[idx])/Rho[idx];
	U[0][idx]=U[0][idx]/Rho[idx];
	U[1][idx]=U[1][idx]/Rho[idx];
}
void D2Q9ColourFluid::MacroVariablesWithNormalDensityAndForce(int& idx){
	U[0][idx]=0;
	U[1][idx]=0;
	Rhor[idx]=0,Rhob[idx]=0;
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
	RhoN[idx]=(Rhor[idx]-Rhob[idx])/Rho[idx];
	U[0][idx]=(U[0][idx]+0.5*F[0][idx])/Rho[idx];
	U[1][idx]=(U[1][idx]+0.5*F[1][idx])/Rho[idx];
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
				f[1]->f[7][NodeIn.Get_index()]=f[1]->f[Opposite[7]][NodeIn.Get_index()];
				break;
			default:
				std::cerr<<"Direction wall bounce back not found. Index: "<<NodeIn.Get_index()<<" x: "<<NodeIn.get_x()<<" y: "<<NodeIn.get_y()<<" direction not found: "<<NodeIn.Get_BcNormal()<<std::endl;
				break;
			}
}
