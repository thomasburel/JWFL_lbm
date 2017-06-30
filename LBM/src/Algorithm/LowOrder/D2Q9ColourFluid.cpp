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
	RhoN=0;
	PtrColourGrad=0;
	PtrColourGradWall=0;
	PtrColourGradBc=0;
	PtrColourGradCorner=0;
	PtrCollision=0;
//	PtrCalNormal=0;
//	PtrCalNormalWall=0;

	Curv=0;
	D_tmp=0;
	PtrD_tmp=0;
	PtrRecolour=0;
	PtrMacro=0;
	beta=0.7;
	Rho_limiter=1.e-5;
	A1=0;
	A2=0;
	F=0;
	G=0;
	Fi=0;
	Normal=0;
	G_Norm=0;
	tension=0;
	teta=0;
	LimitGNorm=0.02;


}

D2Q9ColourFluid::~D2Q9ColourFluid() {
	delete [] Rhor;
	delete [] Rhob;
	delete[]Fi;
}

D2Q9ColourFluid::D2Q9ColourFluid(MultiBlock* MultiBlock__,ParallelManager* parallel__,WriterManager* Writer__, Parameters* Parameters__,InitLBM& ini){
//Initialise the common variables for the two phase models.
	InitD2Q9TwoPhases(MultiBlock__,parallel__,Writer__, Parameters__, ini);
//Set Pointers On Functions for selecting the right model dynamically
	Set_PointersOnFunctions();
//Initialise variables of the colour fluid model.
	InitColourFluid(ini);
//Initialise boundary conditions.
	InitD2Q9Bc(Dic, Parameters__,Ei);
//Set_Convergence
	Set_Convergence();
//Initialise communication between processors
	IniComVariables();
//Set the variables names and the variable pointers for output in solution
	Solution2D::Set_output();
//Set the variables names and the variable pointers for breakpoints in solution
	Solution2D::Set_breakpoint();

}
void D2Q9ColourFluid::Set_PointersOnFunctions(){
// Select the model for two-phase operator in the collision step
	Set_Collide();
// Select the model for the colour gradient
	Set_Colour_gradient();
// Select the method for recolouring
	Set_Recolouring();
// Select the macrocospic function for the colour fluid model depending of the output variables and models
	Set_Macro();
// Select the extrapolations needed
	Set_Extrapolations();
}

void D2Q9ColourFluid::Select_Colour_Operator(ColourFluidEnum::ColourOperatorType OperatorType_){
	switch(OperatorType_)
	{
	case ColourFluidEnum::SurfaceForce:
		PtrCollision=&D2Q9ColourFluid::Collision_SurfaceForce;
		PtrCollisionBc=&D2Q9ColourFluid::Collision_SurfaceForceBc;
		PtrCollisionCorner=&D2Q9ColourFluid::Collision_SurfaceForceCorner;
		break;
	case ColourFluidEnum::Grunau:
		PtrCollision=&D2Q9ColourFluid::Collision_Grunau;
		PtrCollisionBc=PtrCollision;
		PtrCollisionCorner=PtrCollision;
		break;
	case ColourFluidEnum::Reis:
		PtrCollision=&D2Q9ColourFluid::Collision_Reis;
		PtrCollisionBc=PtrCollision;
		PtrCollisionCorner=PtrCollision;
		break;
	default:
		std::cerr<<" Colour operator not found."<<std::endl;
		break;
	}
}


void D2Q9ColourFluid::Set_Extrapolations(){
	if(PtrParameters->Get_ColourExtrapolNoramlDensity())
	{
		ExtrapolDensity.initExtrapolation(2,9,ModelEnum::WeightDistanceExtrapol);
	}
	else
		ExtrapolDensity.initExtrapolation(2,9,ModelEnum::NoExtrapol);

}
void D2Q9ColourFluid::Set_Collide(){
	PtrDicCollide=Dic;


if(PtrParameters->Get_ColourOperatorType()== ColourFluidEnum::SurfaceForce && PtrParameters->Get_UserForceType()== ModelEnum::LocalForce)
{
	std::cout<<" Warming: The local user force will be ignored. The colour model with the surface force is incompatible with a local user force"<<std::endl;
}
	if(PtrParameters->Get_ViscosityType()==ConstViscosity)
	{

		if(PtrParameters->Get_ColourOperatorType()== ColourFluidEnum::SurfaceForce)
			{
			Select_Collide_2D(Std2DBody,PtrParameters->Get_cs2(),PtrParameters->Get_ReferenceDensity(),PtrParameters->Get_ModelOfFluid());Select_Colour_Operator(PtrParameters->Get_ColourOperatorType());
			Select_Collide_2D_V2(*PtrParameters,Std2DBody,PtrParameters->Get_cs2(),PtrParameters->Get_ReferenceDensity(),PtrParameters->Get_ModelOfFluid());Select_Colour_Operator(PtrParameters->Get_ColourOperatorType());
			}
		else if(PtrParameters->Get_UserForceType()== ModelEnum::LocalForce ||PtrParameters->Get_ColourOperatorType()== ColourFluidEnum::Reis||PtrParameters->Get_ColourOperatorType()== ColourFluidEnum::Grunau)
			{
			Select_Collide_2D(Std2DLocal,PtrParameters->Get_cs2(),PtrParameters->Get_ReferenceDensity(),PtrParameters->Get_ModelOfFluid());Select_Colour_Operator(PtrParameters->Get_ColourOperatorType());
			Select_Collide_2D_V2(*PtrParameters,Std2DLocal,PtrParameters->Get_cs2(),PtrParameters->Get_ReferenceDensity(),PtrParameters->Get_ModelOfFluid());Select_Colour_Operator(PtrParameters->Get_ColourOperatorType());
			}
		else if(PtrParameters->Get_UserForceType()== ModelEnum::BodyForce)
			{
			Select_Collide_2D(Std2DBody,PtrParameters->Get_cs2(),PtrParameters->Get_ReferenceDensity(),PtrParameters->Get_ModelOfFluid());Select_Colour_Operator(PtrParameters->Get_ColourOperatorType());
			Select_Collide_2D_V2(*PtrParameters,Std2DBody,PtrParameters->Get_cs2(),PtrParameters->Get_ReferenceDensity(),PtrParameters->Get_ModelOfFluid());Select_Colour_Operator(PtrParameters->Get_ColourOperatorType());
			}
		else
			{
			Select_Collide_2D(Std2D,PtrParameters->Get_cs2(),PtrParameters->Get_ReferenceDensity(),PtrParameters->Get_ModelOfFluid());Select_Colour_Operator(PtrParameters->Get_ColourOperatorType());
			Select_Collide_2D_V2(*PtrParameters,Std2D,PtrParameters->Get_cs2(),PtrParameters->Get_ReferenceDensity(),PtrParameters->Get_ModelOfFluid());Select_Colour_Operator(PtrParameters->Get_ColourOperatorType());
			}
	}
	else
	{
		if(PtrParameters->Get_ColourOperatorType()== ColourFluidEnum::SurfaceForce)
			{
			Select_Collide_2D(Std2DNonCstTauBody,PtrParameters->Get_cs2(),PtrParameters->Get_ReferenceDensity(),PtrParameters->Get_ModelOfFluid());Select_Colour_Operator(PtrParameters->Get_ColourOperatorType());
			Select_Collide_2D_V2(*PtrParameters,Std2DNonCstTauBody,PtrParameters->Get_cs2(),PtrParameters->Get_ReferenceDensity(),PtrParameters->Get_ModelOfFluid());Select_Colour_Operator(PtrParameters->Get_ColourOperatorType());
			}
		else if(PtrParameters->Get_UserForceType()== ModelEnum::LocalForce||PtrParameters->Get_ColourOperatorType()== ColourFluidEnum::Reis||PtrParameters->Get_ColourOperatorType()== ColourFluidEnum::Grunau)
			{
			Select_Collide_2D(Std2DNonCstTauLocal,PtrParameters->Get_cs2(),PtrParameters->Get_ReferenceDensity(),PtrParameters->Get_ModelOfFluid());Select_Colour_Operator(PtrParameters->Get_ColourOperatorType());
			Select_Collide_2D_V2(*PtrParameters,Std2DNonCstTauLocal,PtrParameters->Get_cs2(),PtrParameters->Get_ReferenceDensity(),PtrParameters->Get_ModelOfFluid());Select_Colour_Operator(PtrParameters->Get_ColourOperatorType());
			}
		else if(PtrParameters->Get_UserForceType()== ModelEnum::BodyForce)
			{
			Select_Collide_2D(Std2DNonCstTauBody,PtrParameters->Get_cs2(),PtrParameters->Get_ReferenceDensity(),PtrParameters->Get_ModelOfFluid());Select_Colour_Operator(PtrParameters->Get_ColourOperatorType());
			Select_Collide_2D_V2(*PtrParameters,Std2DNonCstTauBody,PtrParameters->Get_cs2(),PtrParameters->Get_ReferenceDensity(),PtrParameters->Get_ModelOfFluid());Select_Colour_Operator(PtrParameters->Get_ColourOperatorType());
			}
		else
			{
			Select_Collide_2D(Std2DNonCstTau,PtrParameters->Get_cs2(),PtrParameters->Get_ReferenceDensity(),PtrParameters->Get_ModelOfFluid());Select_Colour_Operator(PtrParameters->Get_ColourOperatorType());
			Select_Collide_2D_V2(*PtrParameters,Std2DNonCstTau,PtrParameters->Get_cs2(),PtrParameters->Get_ReferenceDensity(),PtrParameters->Get_ModelOfFluid());Select_Colour_Operator(PtrParameters->Get_ColourOperatorType());
			}
	}

}
void D2Q9ColourFluid::Set_Colour_gradient(){
		switch(PtrParameters->Get_ColourGradType())
		{
		case ColourFluidEnum::Gunstensen:
			if(PtrParameters->Get_ColourOperatorType()== ColourFluidEnum::SurfaceForce)
			{
				std::cout<<" Force colour gradient to Density Normal gradient due to Surface force model."<<std::endl;
				PtrColourGrad =&D2Q9ColourFluid::Colour_gradient_DensityNormalGrad;
				PtrColourGradWall =&D2Q9ColourFluid::Colour_gradient_DensityNormalGrad;
				PtrColourGradBc =&D2Q9ColourFluid::Colour_gradient_DensityNormalGradBc;
				PtrColourGradCorner =&D2Q9ColourFluid::Colour_gradient_DensityNormalGradCorner;
			}
			else
			{
			PtrColourGrad =&D2Q9ColourFluid::Colour_gradient_Gunstensen;
			PtrColourGradWall =&D2Q9ColourFluid::Colour_gradient_Gunstensen;
			PtrColourGradBc =&D2Q9ColourFluid::Colour_gradient_Gunstensen;
			PtrColourGradCorner =&D2Q9ColourFluid::Colour_gradient_Gunstensen;
			}
			break;
		case ColourFluidEnum::DensityGrad:
			PtrColourGrad =&D2Q9ColourFluid::Colour_gradient_DensityGrad;
			PtrColourGradWall =&D2Q9ColourFluid::Colour_gradient_DensityGrad;
			PtrColourGradBc =&D2Q9ColourFluid::Colour_gradient_DensityGradBc;
			PtrColourGradCorner =&D2Q9ColourFluid::Colour_gradient_DensityGradCorner;
			break;
		case ColourFluidEnum::DensityNormalGrad:
			PtrColourGrad =&D2Q9ColourFluid::Colour_gradient_DensityNormalGrad;
			PtrColourGradWall =&D2Q9ColourFluid::Colour_gradient_DensityNormalGrad;
			PtrColourGradBc =&D2Q9ColourFluid::Colour_gradient_DensityNormalGradBc;
			PtrColourGradCorner =&D2Q9ColourFluid::Colour_gradient_DensityNormalGradCorner;
			break;
		}

}
void D2Q9ColourFluid::Set_Macro(){
	switch(PtrParameters->Get_ColourGradType())
	{
	case ColourFluidEnum::Gunstensen:
		if(PtrParameters->Get_ColourOperatorType()== ColourFluidEnum::SurfaceForce )
			PtrMacro=&D2Q9ColourFluid::MacroVariablesWithNormalDensityAndForce;
		else
			if(PtrParameters->Get_NormalDensityOutput())
				PtrMacro=&D2Q9ColourFluid::MacroVariablesWithNormalDensity;
			else
				PtrMacro=&D2Q9ColourFluid::MacroVariables;
		break;
	case ColourFluidEnum::DensityGrad:
		if(PtrParameters->Get_ColourOperatorType()== ColourFluidEnum::SurfaceForce)
			PtrMacro=&D2Q9ColourFluid::MacroVariablesWithNormalDensityAndForce;
		else
			if(PtrParameters->Get_NormalDensityOutput())
				PtrMacro=&D2Q9ColourFluid::MacroVariablesWithNormalDensity;
			else
				PtrMacro=&D2Q9ColourFluid::MacroVariables;
		break;
	case ColourFluidEnum::DensityNormalGrad:
		if(PtrParameters->Get_ColourOperatorType()== ColourFluidEnum::SurfaceForce)
			PtrMacro=&D2Q9ColourFluid::MacroVariablesWithNormalDensityAndForce;
		else
			PtrMacro=&D2Q9ColourFluid::MacroVariablesWithNormalDensity;
		break;
	}

}
void D2Q9ColourFluid::Set_Recolouring(){
	switch(PtrParameters->Get_RecolouringType())
	{
	case ColourFluidEnum::LatvaKokkoRothman:
		PtrRecolour=&D2Q9ColourFluid::Recolouring_Latva;
		break;
	case ColourFluidEnum::LatvaKokkoRothmanEstimator:
		PtrRecolour=&D2Q9ColourFluid::Recolouring_LatvaWithEstimator;
		break;
	default:
		std::cerr<<"Recolouring operator not found. Set default: Latva-Kokko-Rothman operator"<<std::endl;
		PtrRecolour=&D2Q9ColourFluid::Recolouring_Latva;
		break;
	}



}
void D2Q9ColourFluid::InitColourFluid(InitLBM& ini){
//Initialise parameters
	Rho_limiter=PtrParameters->Get_RhoLimiter();
	LimitGNorm=PtrParameters->Get_ColourGradLimiter();
	beta=PtrParameters->Get_Beta();
	A1=PtrParameters->Get_ATau();
	//A1=PtrParameters->Get_A1();
	A2=PtrParameters->Get_A2();
	tension=PtrParameters->Get_SurfaceTension();

// Initialise only the Colour Fluid part.
	Bi[0]=-4.0/27.0;
	for(int i=1;i<5;i++) Bi[i]=2.0/27.0;
	for(int i=5;i<9;i++) Bi[i]=5.0/108.0;
	G=new double*[2];
	F=new double*[2];
	Fi=new double[9];for(int i=0;i<9;i++) Fi[i]=0;
	Normal=new double*[2];
	CAngle.AllocateTeta(NodeArrays,PtrParameters,teta);
	Dic->AddSync("Density",Rho);
	Dic->AddVar(Vector,"ColourGrad",PtrParameters->Get_ColourGradientOutput(),true,false,G[0],G[1]);
	Dic->AddVar(Scalar,"ColourGrad_Norm",PtrParameters->Get_NormColourGradientOutput(),true,false,G_Norm);
	Dic->AddVar(Vector,"Normal",PtrParameters->Get_NormalOutput(),true,false,Normal[0],Normal[1]);
	Dic->AddVar(Scalar,"RhoRed",PtrParameters->Get_RedDensityOutput(),true,true,Rhor);
	Dic->AddVar(Scalar,"RhoBlue",PtrParameters->Get_BlueDensityOutput(),true,true,Rhob);
	Dic->AddVar(Scalar,"RhoN",PtrParameters->Get_NormalDensityOutput(),true,true,RhoN);
	if (PtrParameters->Get_ColourOperatorType()==ColourFluidEnum::SurfaceForce)
	{
		Dic->AddVar(Vector,"SurfaceForce",PtrParameters->Get_SurfaceForceOutput(),true,false,F[0],F[1]);
		Dic->AddVar(Scalar,"Curvature",PtrParameters->Get_CurvatureOutput(),false,false,Curv);
		for(int i=0;i<nbnodes_total;i++)
		{
			F[0][i]=0;
			F[1][i]=0;
			Curv[i]=0;
		}
	}
	for(int i=0;i<nbnodes_total;i++)
	{
		G[0][i]=0;
		G[1][i]=0;
		G_Norm[i]=0;
		Normal[0][i]=0;
		Normal[1][i]=0;
	}
	InitColourFluidAllDomain(ini);
	//Initialise from file or Restart
	InitialiseFromFile();
	//Initialise distribution
	InitDistColourFluidAllDomain();
}
void D2Q9ColourFluid::InitColourFluidAllDomain(InitLBM& ini){
	InitColourFluidDomainBc(ini);
	InitColourFluidWall(ini);
	InitColourFluidInterior(ini);

// Init Solid
	double alpha=0;
	double* pos =new double[2];
	double* U_=new double[2];
	double Rho_tmp;
	int idx=0;

	for (int j=0;j<NodeArrays->NodeSolid.size();j++)
	{
// Set Index
		idx=NodeArrays->NodeSolid[j].Get_index();
// Get position
		pos[0]=NodeArrays->NodeSolid[j].get_x();
		pos[1]=NodeArrays->NodeSolid[j].get_y();
// Get initialise values from the user
		ini.IniDomainTwoPhases(parallel->getRank(),NodeArrays->NodeSolid[j],0, idx,pos,Rho_tmp,U_,alpha);
// Initialise the blue and red densities
		Rhor[idx]=alpha*Rho[idx];
		Rhob[idx]=(1- alpha) *Rho[idx];
		RhoN[idx]=-2;
	}
	delete [] pos;
	delete [] U_;
}
void D2Q9ColourFluid::InitColourFluidDomainBc(InitLBM& ini){
	double alpha=0;
	double* pos =new double[2];
	double* U_=new double[2];
	double Rho_tmp;
	int idx=0;
	for (int j=0;j<NodeArrays->NodeGlobalCorner.size();j++)
	{
// Set Index
		idx=NodeArrays->NodeGlobalCorner[j].Get_index();
// Get position
		pos[0]=NodeArrays->NodeGlobalCorner[j].get_x();
		pos[1]=NodeArrays->NodeGlobalCorner[j].get_y();
// Get initialise values from the user
		ini.IniDomainTwoPhases(parallel->getRank(),NodeArrays->NodeGlobalCorner[j],0, idx,pos,Rho_tmp,U_,alpha);
// Initialise the blue and red densities
		Rhor[idx]=alpha*Rho[idx];
		Rhob[idx]=(1- alpha) *Rho[idx];
		RhoN[idx]=(Rhor[idx]-Rhob[idx])/Rho[idx];
	}
	int count=0;
	for (int j=0;j<NodeArrays->NodeSpecialWall.size();j++)
	{
// Set Index
		idx=NodeArrays->NodeSpecialWall[j].Get_index();
// Get position
		pos[0]=NodeArrays->NodeSpecialWall[j].get_x();
		pos[1]=NodeArrays->NodeSpecialWall[j].get_y();
// Get initialise values from the user
		ini.IniDomainTwoPhases(parallel->getRank(),NodeArrays->NodeSpecialWall[j],0, idx,pos,Rho_tmp,U_,alpha);
// Set contact angle if needed
		if(PtrParameters->Get_ContactAngleType()==ContactAngleEnum::UserTeta)
			ini.IniContactAngle(parallel->getRank(),NodeArrays->NodeSpecialWall[j],0, idx,pos,teta[count]);
// Initialise the blue and red densities
		Rhor[idx]=alpha*Rho[idx];
		Rhob[idx]=(1- alpha) *Rho[idx];
		RhoN[idx]=(Rhor[idx]-Rhob[idx])/Rho[idx];
		count++;
	}
	for (int j=0;j<NodeArrays->NodeVelocity.size();j++)
	{
// Set Index
		idx=NodeArrays->NodeVelocity[j].Get_index();
// Get position
		pos[0]=NodeArrays->NodeVelocity[j].get_x();
		pos[1]=NodeArrays->NodeVelocity[j].get_y();
// Get initialise values from the user
		ini.IniDomainTwoPhases(parallel->getRank(),NodeArrays->NodeVelocity[j],0, idx,pos,Rho_tmp,U_,alpha);
// Initialise the blue and red densities
		Rhor[idx]=alpha*Rho[idx];
		Rhob[idx]=(1- alpha) *Rho[idx];
		RhoN[idx]=(Rhor[idx]-Rhob[idx])/Rho[idx];
	}

	for (int j=0;j<NodeArrays->NodePressure.size();j++)
	{
// Set Index
		idx=NodeArrays->NodePressure[j].Get_index();
// Get position
		pos[0]=NodeArrays->NodePressure[j].get_x();
		pos[1]=NodeArrays->NodePressure[j].get_y();
// Get initialise values from the user
		ini.IniDomainTwoPhases(parallel->getRank(),NodeArrays->NodePressure[j],0, idx,pos,Rho_tmp,U_,alpha);
// Initialise the blue and red densities
		Rhor[idx]=alpha*Rho[idx];
		Rhob[idx]=(1- alpha) *Rho[idx];
		RhoN[idx]=(Rhor[idx]-Rhob[idx])/Rho[idx];
	}
	for (int j=0;j<NodeArrays->NodeSymmetry.size();j++)
	{
// Set Index
		idx=NodeArrays->NodeSymmetry[j].Get_index();
// Get position
		pos[0]=NodeArrays->NodeSymmetry[j].get_x();
		pos[1]=NodeArrays->NodeSymmetry[j].get_y();
// Get initialise values from the user
		ini.IniDomainTwoPhases(parallel->getRank(),NodeArrays->NodeSymmetry[j],0, idx,pos,Rho_tmp,U_,alpha);
// Initialise the blue and red densities
		Rhor[idx]=alpha*Rho[idx];
		Rhob[idx]=(1- alpha) *Rho[idx];
		RhoN[idx]=(Rhor[idx]-Rhob[idx])/Rho[idx];
	}
	for (int j=0;j<NodeArrays->NodePeriodic.size();j++)
	{
// Set Index
		idx=NodeArrays->NodePeriodic[j].Get_index();
// Get position
		pos[0]=NodeArrays->NodePeriodic[j].get_x();
		pos[1]=NodeArrays->NodePeriodic[j].get_y();
// Get initialise values from the user
		ini.IniDomainTwoPhases(parallel->getRank(),NodeArrays->NodePeriodic[j],0, idx,pos,Rho_tmp,U_,alpha);
// Initialise the blue and red densities
		Rhor[idx]=alpha*Rho[idx];
		Rhob[idx]=(1- alpha) *Rho[idx];
		RhoN[idx]=(Rhor[idx]-Rhob[idx])/Rho[idx];
	}
	delete [] pos;
	delete [] U_;
}
void D2Q9ColourFluid::InitColourFluidWall(InitLBM& ini){
	double alpha=0;
	double* pos =new double[2];
	double* U_=new double[2];
	double Rho_tmp;
	int idx=0;
	int count=NodeArrays->NodeSpecialWall.size();

	for (int j=0;j<NodeArrays->CornerConcave.size();j++)
	{
// Set Index
		idx=NodeArrays->NodeCorner[NodeArrays->CornerConcave[j]].Get_index();
// Get position
		pos[0]=NodeArrays->NodeCorner[NodeArrays->CornerConcave[j]].get_x();
		pos[1]=NodeArrays->NodeCorner[NodeArrays->CornerConcave[j]].get_y();
// Get initialise value from the user
		ini.IniDomainTwoPhases(parallel->getRank(),NodeArrays->NodeCorner[NodeArrays->CornerConcave[j]],0, idx,pos,Rho_tmp,U_,alpha);
// Set contact angle if needed
		if(PtrParameters->Get_ContactAngleType()==ContactAngleEnum::UserTeta)
			ini.IniContactAngle(parallel->getRank(),NodeArrays->NodeCorner[NodeArrays->CornerConcave[j]],0, idx,pos,teta[count]);
// Initialise the blue and red densities
		Rhor[idx]=alpha*Rho[idx];
		Rhob[idx]=(1- alpha) *Rho[idx];
		RhoN[idx]=(Rhor[idx]-Rhob[idx])/Rho[idx];
		count++;
	}

	for (int j=0;j<NodeArrays->NodeWall.size();j++)
	{
// Set Index
		idx=NodeArrays->NodeWall[j].Get_index();
// Get position
		pos[0]=NodeArrays->NodeWall[j].get_x();
		pos[1]=NodeArrays->NodeWall[j].get_y();
// Get initialise values from the user
		ini.IniDomainTwoPhases(parallel->getRank(),NodeArrays->NodeWall[j],0, idx,pos,Rho_tmp,U_,alpha);
// Set contact angle if needed
		if(PtrParameters->Get_ContactAngleType()==ContactAngleEnum::UserTeta)
			ini.IniContactAngle(parallel->getRank(),NodeArrays->NodeWall[j],0, idx,pos,teta[count]);
// Initialise the blue and red densities
		Rhor[idx]=alpha*Rho[idx];
		Rhob[idx]=(1- alpha) *Rho[idx];
		RhoN[idx]=(Rhor[idx]-Rhob[idx])/Rho[idx];
		count++;
	}

	for (int j=0;j<NodeArrays->CornerConvex.size();j++)
	{
// Set Index
		idx=NodeArrays->NodeCorner[NodeArrays->CornerConvex[j]].Get_index();
// Get position
		pos[0]=NodeArrays->NodeCorner[NodeArrays->CornerConvex[j]].get_x();
		pos[1]=NodeArrays->NodeCorner[NodeArrays->CornerConvex[j]].get_y();
// Get initialise value from the user
		ini.IniDomainTwoPhases(parallel->getRank(),NodeArrays->NodeCorner[NodeArrays->CornerConvex[j]],0, idx,pos,Rho_tmp,U_,alpha);
// Set contact angle if needed
		if(PtrParameters->Get_ContactAngleType()==ContactAngleEnum::UserTeta)
			ini.IniContactAngle(parallel->getRank(),NodeArrays->NodeCorner[NodeArrays->CornerConvex[j]],0, idx,pos,teta[count]);
// Initialise the blue and red densities
		Rhor[idx]=alpha*Rho[idx];
		Rhob[idx]=(1- alpha) *Rho[idx];
		RhoN[idx]=(Rhor[idx]-Rhob[idx])/Rho[idx];
		count++;
	}
	CAngle.InitContactAngle(NodeArrays,PtrParameters,Opposite,teta);
	delete [] pos;
	delete [] U_;
}
void D2Q9ColourFluid::InitColourFluidInterior(InitLBM& ini){
	double alpha=0;
	double* pos =new double[2];
	double* U_=new double[2];
	double Rho_tmp;
	int idx=0;

	for (int j=0;j<NodeArrays->NodeInterior.size();j++)
	{
// Set Index
		idx=NodeArrays->NodeInterior[j].Get_index();
// Get position
		pos[0]=NodeArrays->NodeInterior[j].get_x();
		pos[1]=NodeArrays->NodeInterior[j].get_y();
// Get initialise value from the user
		ini.IniDomainTwoPhases(parallel->getRank(),NodeArrays->NodeInterior[j],0,idx,pos,Rho_tmp,U_,alpha);
// Initialise the blue and red densities
		Rhor[idx]=alpha*Rho[idx];
		Rhob[idx]=(1- alpha) *Rho[idx];
		RhoN[idx]=(Rhor[idx]-Rhob[idx])/Rho[idx];
	}

	for (int j=0;j<NodeArrays->NodeGhost.size();j++)
	{
// Set Index
		idx=NodeArrays->NodeGhost[j].Get_index();
// Get position
		pos[0]=NodeArrays->NodeGhost[j].get_x();
		pos[1]=NodeArrays->NodeGhost[j].get_y();
// Get initialise values from the user
		ini.IniDomainTwoPhases(parallel->getRank(),NodeArrays->NodeGhost[j],0, idx,pos,Rho_tmp,U_,alpha);
// Initialise the blue and red densities
		Rhor[idx]=alpha*Rho[idx];
		Rhob[idx]=(1- alpha) *Rho[idx];
		RhoN[idx]=(Rhor[idx]-Rhob[idx])/Rho[idx];
	}

	delete [] pos;
	delete [] U_;
}
void D2Q9ColourFluid::InitDistColourFluidAllDomain(){
	InitDistColourFluidDomainBc();
	InitDistColourFluidWall();
	InitDistColourFluidInterior();

// Init Solid
	double alpha=0;
	double* pos =new double[2];
	double* U_=new double[2];
	double Rho_tmp;
	int idx=0;

	for (int j=0;j<NodeArrays->NodeSolid.size();j++)
	{
// Set Index
		idx=NodeArrays->NodeSolid[j].Get_index();
// Get initialise values from the user
		alpha=(RhoN[idx]+1.0)*0.5;
// Initialise the blue and red densities
		Rhor[idx]=alpha*Rho[idx];
		Rhob[idx]=(1- alpha) *Rho[idx];
		RhoN[idx]=-2;
		for (int i=0;i<nbvelo;i++)
		{
			f[0]->f[i][idx]=0;
			f[1]->f[i][idx]=0;
		}
	}
	delete [] pos;
	delete [] U_;
}
void D2Q9ColourFluid::InitDistColourFluidDomainBc(){
	double alpha=0;
	double* pos =new double[2];
	double* U_=new double[2];
	double Rho_tmp;
	int idx=0;
	for (int j=0;j<NodeArrays->NodeGlobalCorner.size();j++)
	{
// Set Index
		idx=NodeArrays->NodeGlobalCorner[j].Get_index();
// Get initialise values from the user
		alpha=(RhoN[idx]+1.0)*0.5;
// Initialise the blue and red densities
		Rhor[idx]=alpha*Rho[idx];
		Rhob[idx]=(1- alpha) *Rho[idx];
		RhoN[idx]=(Rhor[idx]-Rhob[idx])/Rho[idx];
		for (int i=0;i<nbvelo;i++)
		{
			f[0]->f[i][idx]=CollideLowOrder::CollideEquillibrium(Rhor[idx], U[0][idx], U[1][idx], &Ei[i][0], omega[i]);
			f[1]->f[i][idx]=CollideLowOrder::CollideEquillibrium(Rhob[idx], U[0][idx], U[1][idx], &Ei[i][0], omega[i]);
		}
	}
	for (int j=0;j<NodeArrays->NodeSpecialWall.size();j++)
	{
// Set Index
		idx=NodeArrays->NodeSpecialWall[j].Get_index();
// Get initialise values from the user
		alpha=(RhoN[idx]+1.0)*0.5;
// Initialise the blue and red densities
		Rhor[idx]=alpha*Rho[idx];
		Rhob[idx]=(1- alpha) *Rho[idx];
		RhoN[idx]=(Rhor[idx]-Rhob[idx])/Rho[idx];
		for (int i=0;i<nbvelo;i++)
		{
			f[0]->f[i][idx]=CollideLowOrder::CollideEquillibrium(Rhor[idx], U[0][idx], U[1][idx], &Ei[i][0], omega[i]);
			f[1]->f[i][idx]=CollideLowOrder::CollideEquillibrium(Rhob[idx], U[0][idx], U[1][idx], &Ei[i][0], omega[i]);
		}
	}
	for (int j=0;j<NodeArrays->NodeVelocity.size();j++)
	{
// Set Index
		idx=NodeArrays->NodeVelocity[j].Get_index();
// Get initialise values from the user
		alpha=(RhoN[idx]+1.0)*0.5;
// Initialise the blue and red densities
		Rhor[idx]=alpha*Rho[idx];
		Rhob[idx]=(1- alpha) *Rho[idx];
		RhoN[idx]=(Rhor[idx]-Rhob[idx])/Rho[idx];
		for (int i=0;i<nbvelo;i++)
		{
			f[0]->f[i][idx]=CollideLowOrder::CollideEquillibrium(Rhor[idx], U[0][idx], U[1][idx], &Ei[i][0], omega[i]);
			f[1]->f[i][idx]=CollideLowOrder::CollideEquillibrium(Rhob[idx], U[0][idx], U[1][idx], &Ei[i][0], omega[i]);
		}
	}

	for (int j=0;j<NodeArrays->NodePressure.size();j++)
	{
// Set Index
		idx=NodeArrays->NodePressure[j].Get_index();
// Get initialise values from the user
		alpha=(RhoN[idx]+1.0)*0.5;
// Initialise the blue and red densities
		Rhor[idx]=alpha*Rho[idx];
		Rhob[idx]=(1- alpha) *Rho[idx];
		RhoN[idx]=(Rhor[idx]-Rhob[idx])/Rho[idx];
		for (int i=0;i<nbvelo;i++)
		{
			f[0]->f[i][idx]=CollideLowOrder::CollideEquillibrium(Rhor[idx], U[0][idx], U[1][idx], &Ei[i][0], omega[i]);
			f[1]->f[i][idx]=CollideLowOrder::CollideEquillibrium(Rhob[idx], U[0][idx], U[1][idx], &Ei[i][0], omega[i]);
		}
	}
	for (int j=0;j<NodeArrays->NodeSymmetry.size();j++)
	{
// Set Index
		idx=NodeArrays->NodeSymmetry[j].Get_index();
// Get initialise values from the user
		alpha=(RhoN[idx]+1.0)*0.5;
// Initialise the blue and red densities
		Rhor[idx]=alpha*Rho[idx];
		Rhob[idx]=(1- alpha) *Rho[idx];
		RhoN[idx]=(Rhor[idx]-Rhob[idx])/Rho[idx];
		for (int i=0;i<nbvelo;i++)
		{
			f[0]->f[i][idx]=CollideLowOrder::CollideEquillibrium(Rhor[idx], U[0][idx], U[1][idx], &Ei[i][0], omega[i]);
			f[1]->f[i][idx]=CollideLowOrder::CollideEquillibrium(Rhob[idx], U[0][idx], U[1][idx], &Ei[i][0], omega[i]);
		}
	}
	for (int j=0;j<NodeArrays->NodePeriodic.size();j++)
	{
// Set Index
		idx=NodeArrays->NodePeriodic[j].Get_index();
// Get initialise values from the user
		alpha=(RhoN[idx]+1.0)*0.5;
// Initialise the blue and red densities
		Rhor[idx]=alpha*Rho[idx];
		Rhob[idx]=(1- alpha) *Rho[idx];
		RhoN[idx]=(Rhor[idx]-Rhob[idx])/Rho[idx];
		for (int i=0;i<nbvelo;i++)
		{
			f[0]->f[i][idx]=CollideLowOrder::CollideEquillibrium(Rhor[idx], U[0][idx], U[1][idx], &Ei[i][0], omega[i]);
			f[1]->f[i][idx]=CollideLowOrder::CollideEquillibrium(Rhob[idx], U[0][idx], U[1][idx], &Ei[i][0], omega[i]);
		}
	}
	delete [] pos;
	delete [] U_;
}
void D2Q9ColourFluid::InitDistColourFluidWall(){
	double alpha=0;
	double* pos =new double[2];
	double* U_=new double[2];
	double Rho_tmp;
	int idx=0;
	int count=0;

	for (int j=0;j<NodeArrays->CornerConcave.size();j++)
	{
// Set Index
		idx=NodeArrays->NodeCorner[NodeArrays->CornerConcave[j]].Get_index();
// Get initialise value from the user
		alpha=(RhoN[idx]+1.0)*0.5;// Initialise the blue and red densities
		Rhor[idx]=alpha*Rho[idx];
		Rhob[idx]=(1- alpha) *Rho[idx];
		RhoN[idx]=(Rhor[idx]-Rhob[idx])/Rho[idx];
		for (int i=0;i<nbvelo;i++)
		{
			f[0]->f[i][idx]=CollideLowOrder::CollideEquillibrium(Rhor[idx], U[0][idx], U[1][idx], &Ei[i][0], omega[i]);
			f[1]->f[i][idx]=CollideLowOrder::CollideEquillibrium(Rhob[idx], U[0][idx], U[1][idx], &Ei[i][0], omega[i]);
		}
		count++;
	}

	for (int j=0;j<NodeArrays->NodeWall.size();j++)
	{
// Set Index
		idx=NodeArrays->NodeWall[j].Get_index();
// Get initialise values from the user
		alpha=(RhoN[idx]+1.0)*0.5;
// Initialise the blue and red densities
		Rhor[idx]=alpha*Rho[idx];
		Rhob[idx]=(1- alpha) *Rho[idx];
		RhoN[idx]=(Rhor[idx]-Rhob[idx])/Rho[idx];
		for (int i=0;i<nbvelo;i++)
		{
			f[0]->f[i][idx]=CollideLowOrder::CollideEquillibrium(Rhor[idx], U[0][idx], U[1][idx], &Ei[i][0], omega[i]);
			f[1]->f[i][idx]=CollideLowOrder::CollideEquillibrium(Rhob[idx], U[0][idx], U[1][idx], &Ei[i][0], omega[i]);
		}
		count++;
	}

	for (int j=0;j<NodeArrays->CornerConvex.size();j++)
	{
// Set Index
		idx=NodeArrays->NodeCorner[NodeArrays->CornerConvex[j]].Get_index();
// Get initialise value from the user
		alpha=(RhoN[idx]+1.0)*0.5;
// Initialise the blue and red densities
		Rhor[idx]=alpha*Rho[idx];
		Rhob[idx]=(1- alpha) *Rho[idx];
		RhoN[idx]=(Rhor[idx]-Rhob[idx])/Rho[idx];
		for (int i=0;i<nbvelo;i++)
		{
			f[0]->f[i][idx]=CollideLowOrder::CollideEquillibrium(Rhor[idx], U[0][idx], U[1][idx], &Ei[i][0], omega[i]);
			f[1]->f[i][idx]=CollideLowOrder::CollideEquillibrium(Rhob[idx], U[0][idx], U[1][idx], &Ei[i][0], omega[i]);
		}
		count++;
	}
	delete [] pos;
	delete [] U_;
}
void D2Q9ColourFluid::InitDistColourFluidInterior(){
	double alpha=0;
	double* pos =new double[2];
	double* U_=new double[2];
	double Rho_tmp;
	int idx=0;

	for (int j=0;j<NodeArrays->NodeInterior.size();j++)
	{
// Set Index
		idx=NodeArrays->NodeInterior[j].Get_index();
// Get initialise value from the user
		alpha=(RhoN[idx]+1.0)*0.5;
// Initialise the blue and red densities
		Rhor[idx]=alpha*Rho[idx];
		Rhob[idx]=(1- alpha) *Rho[idx];
		RhoN[idx]=(Rhor[idx]-Rhob[idx])/Rho[idx];
		for (int i=0;i<nbvelo;i++)
		{
			f[0]->f[i][idx]=CollideLowOrder::CollideEquillibrium(Rhor[idx], U[0][idx], U[1][idx], &Ei[i][0], omega[i]);
			f[1]->f[i][idx]=CollideLowOrder::CollideEquillibrium(Rhob[idx], U[0][idx], U[1][idx], &Ei[i][0], omega[i]);
		}
	}

	for (int j=0;j<NodeArrays->NodeGhost.size();j++)
	{
// Set Index
		idx=NodeArrays->NodeGhost[j].Get_index();
// Get initialise values from the user
		alpha=(RhoN[idx]+1.0)*0.5;
// Initialise the blue and red densities
		Rhor[idx]=alpha*Rho[idx];
		Rhob[idx]=(1- alpha) *Rho[idx];
		RhoN[idx]=(Rhor[idx]-Rhob[idx])/Rho[idx];
		for (int i=0;i<nbvelo;i++)
		{
			f[0]->f[i][idx]=CollideLowOrder::CollideEquillibrium(Rhor[idx], U[0][idx], U[1][idx], &Ei[i][0], omega[i]);
			f[1]->f[i][idx]=CollideLowOrder::CollideEquillibrium(Rhob[idx], U[0][idx], U[1][idx], &Ei[i][0], omega[i]);
		}
	}
	delete [] pos;
	delete [] U_;
}
void D2Q9ColourFluid::UpdateAllDomain(Parameters* UpdatedParam,InitLBM& ini){
	//init field
	UpdateDomainBc(UpdatedParam,ini);
	UpdateWall(UpdatedParam,ini);
	UpdateInterior(UpdatedParam,ini);
	//init distri
	InitDistColourFluidDomainBc();
	InitDistColourFluidWall();
	InitDistColourFluidInterior();
}
void D2Q9ColourFluid::UpdateDomainBc(Parameters* UpdatedParam,InitLBM& ini){
	//init field
	InitDomainBc(ini);
	InitColourFluidDomainBc(ini);
	//init distri
	InitDistColourFluidDomainBc();
}
void D2Q9ColourFluid::UpdateWall(Parameters* UpdatedParam,InitLBM& ini){
	//init field
	InitWall(ini);
	InitColourFluidWall(ini);
	//init distri
	InitDistColourFluidWall();
}
void D2Q9ColourFluid::UpdateInterior(Parameters* UpdatedParam,InitLBM& ini){
	//init field
	InitInterior(ini);
	InitColourFluidInterior(ini);
	//init distri
	InitDistColourFluidInterior();
}
void D2Q9ColourFluid::run(Parameters* UpdatedParam){

	PtrParameters=UpdatedParam;
	IniTau(PtrParameters);
	InvTau=Get_InvTau();

	//Initialise parameters for colour fluid
	beta=PtrParameters->Get_Beta();
	//A1=PtrParameters->Get_A1();
	A1=PtrParameters->Get_ATau();
	A2=PtrParameters->Get_A2();
	tension=PtrParameters->Get_SurfaceTension();

	D2Q9ColourFluid::run();
}

void D2Q9ColourFluid::run(){
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
	{
		SyncToGhost();
	//	SyncMacroVarToGhost();
		SyncMacroVarToGhost();
	//SyncVarFromSolidGhost(RhoN);
	}
	ExtrapolDensityInSolid();
	if(parallel->getSize()>1)
	{
//		SyncVarFromSolidGhost(RhoN);
		SyncVarToSolidGhost(RhoN);
	}
	ApplyBc();
	Colour_gradient();
	if(parallel->getSize()>1)
		Synchronise_Colour_gradient();

	//Extrapolate_NormalInSolid();
	CAngle.ApplyContactAngle2D(Normal);
	UpdateColourGradientOnWalls();
	//Impose_ContactAngleInSolid();

	ColourFluid_Collision();
	ApplyBc();
	UpdateMacroVariables();
	if(CalPressure)
		UpdatePressure();
	if(parallel->getSize()>1)
	{
		SyncMacroVarToGhost();
//	SyncVarToSolidGhost(RhoN);
	}
	//if(parallel->getSize()>1)
	//	SyncVarToGhost(RhoN);

	ExtrapolDensityInSolid();

	if(parallel->getSize()>1)
	{
//	SyncVarFromSolidGhost(RhoN);
		SyncVarToSolidGhost(RhoN);
//		SyncVarSolidGhost(RhoN);
	}

//	Write_Breakpoint(PtrParameters);

	Convergence::Calcul_Error(it);
	Writer->Write_Output(it);

//	Writer->Write_breakpoint(*PtrParameters);
	it++;
	if(parallel->getSize()>1)
	{

//		for (int i=1;i<NbStep+1;i++)
		while(it<NbStep+1)
		{
			Colour_gradient();
			Synchronise_Colour_gradient();
			//Extrapolate_NormalInSolid();
			CAngle.ApplyContactAngle2D(Normal);
			UpdateColourGradientOnWalls();
			//Impose_ContactAngleInSolid();

			ColourFluid_Collision();
			SyncToGhost();
			StreamD2Q9();;
			ApplyBc();
			UpdateMacroVariables();
			if(CalGradP)
				UpdatePressure();
			SyncMacroVarToGhost();
//			SyncVarToSolidGhost(RhoN);
			ExtrapolDensityInSolid();
			SyncVarToSolidGhost(RhoN);
//			SyncVarFromSolidGhost(RhoN);
			//SyncVarToGhost(RhoN);

			if(it%OutPutNStep==0)
			{
				if(CalPressure&&!CalGradP)
					UpdatePressure();
				Writer->Write_Output(it);
			}
			if(it%listing==0 )
			{
				if(CalPressure &&!CalGradP)
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
					if(CalPressure&&!CalGradP)
						UpdatePressure();
					Writer->Write_Output(it);
					it=NbStep;
				}
			}



	/*		if(it==2)
			Write_Breakpoint(PtrParameters);
			std::cout<<"****"<<std::endl;*/
			it++;
		}
	}
	else
	{
		while(it<NbStep+1)
		{
			Colour_gradient();
			//Extrapolate_NormalInSolid();
			//Impose_ContactAngleInSolid();
			CAngle.ApplyContactAngle2D(Normal);
			UpdateColourGradientOnWalls();
			ColourFluid_Collision();
			StreamD2Q9();
			ApplyBc();
			UpdateMacroVariables();
			if(CalGradP)
				UpdatePressure();
			ExtrapolDensityInSolid();
			if(it%OutPutNStep==0)
			{
				if(CalPressure&&!CalGradP)
					UpdatePressure();
				Writer->Write_Output(it);
			}
			if(it%listing==0)
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

		/*	if(it==2)
				Write_Breakpoint(PtrParameters);*/
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

//	Write_Breakpoint(PtrParameters);
//	Writer->Write_breakpoint(*PtrParameters);
}

///Calculate the colour gradient by Gunstensen formulation, density gradient or normal density gradient
void D2Q9ColourFluid::Colour_gradient(){
	std::vector<SolverEnum::PatchType> PatchType=PatchsBc->Get_PatchTypeInType();
	std::vector<int> PatchIdInType=PatchsBc->Get_PatchIdInType();
	int idx_tmp;int normal_interior=0;
	double teta;
/*	VelocityPatchBc VelPatchBc;
	std::vector<int> NodeIdx;
	for (int i=0;i<PatchsBc->Get_NumberOfPatchBc();i++)
	{
		switch(PatchType[i])
		{
		case SolverEnum::Periodic:
//			ApplyPatchPeriodic(PatchsBc->Get_PeriodicPatch()[PatchIdInType[i]]);
			break;
		case SolverEnum::Symmetry:
//			ApplyPatchSymmetry(PatchsBc->Get_SymmetryPatch()[PatchIdInType[i]]);
			break;
		case SolverEnum::Pressure:
//			ApplyPatchPressure(PatchsBc->Get_PressurePatch()[PatchIdInType[i]]);
			break;
		case SolverEnum::Velocity:
			VelPatchBc=PatchsBc->Get_VelocityPatch()[PatchIdInType[i]];
			NodeIdx=VelPatchBc.Get_NodeIndexByType();
			if(VelPatchBc.Get_extrapolationNormal())
			{
				for (int j=0;j<NodeIdx.size();j++)
				{
				// Common variables
					idx_tmp=NodeArrays->NodeVelocity[NodeIdx[j]].Get_index();
					ExtrapolationWallToSolid(G[0],NodeArrays->NodeVelocity[NodeIdx[j]].Get_connect(),NodeArrays->NodeVelocity[NodeIdx[j]].Get_BcNormal());
					ExtrapolationWallToSolid(G[1],NodeArrays->NodeVelocity[NodeIdx[j]].Get_connect(),NodeArrays->NodeVelocity[NodeIdx[j]].Get_BcNormal());
				// Calculate gradients
					(this->*PtrColourGrad)(idx_tmp,NodeArrays->NodeVelocity[NodeIdx[j]].Get_connect(),NodeArrays->NodeVelocity[NodeIdx[j]].Get_BcNormal());
					CalculNormal(idx_tmp,NodeArrays->NodeVelocity[NodeIdx[j]].Get_connect(),NodeArrays->NodeVelocity[NodeIdx[j]].Get_BcNormal());

					ExtrapolationWallToSolid(Normal[0],NodeArrays->NodeVelocity[NodeIdx[j]].Get_connect(),NodeArrays->NodeVelocity[NodeIdx[j]].Get_BcNormal());
					ExtrapolationWallToSolid(Normal[0],NodeArrays->NodeVelocity[NodeIdx[j]].Get_connect(),NodeArrays->NodeVelocity[NodeIdx[j]].Get_BcNormal());
				}
			}
			else
			{
				for (int j=0;j<NodeIdx.size();j++)
				{
				// Common variables
					idx_tmp=NodeArrays->NodeVelocity[NodeIdx[j]].Get_index();
				// Calculate gradients
					(this->*PtrColourGradBc)(idx_tmp,NodeArrays->NodeVelocity[NodeIdx[j]].Get_connect(),NodeArrays->NodeVelocity[NodeIdx[j]].Get_BcNormal());
					CalculNormal(idx_tmp,NodeArrays->NodeVelocity[NodeIdx[j]].Get_connect(),NodeArrays->NodeVelocity[NodeIdx[j]].Get_BcNormal());
				}
			}
			break;
		case SolverEnum::Wall:
			break;
		}
	}
*/
	for (int j=0;j<NodeArrays->NodeInterior.size();j++)
	{
	// Common variables
		idx_tmp=NodeArrays->NodeInterior[j].Get_index();
	// Calculate gradients
		(this->*PtrColourGrad)(idx_tmp,NodeArrays->NodeInterior[j].Get_connect(),normal_interior);
		CalculNormal(idx_tmp,NodeArrays->NodeInterior[j].Get_connect(),normal_interior);
	}

	for (int j=0;j<NodeArrays->NodeCorner.size();j++)
	{
	// Common variables
		idx_tmp=NodeArrays->NodeCorner[j].Get_index();
	// Calculate gradients
		(this->*PtrColourGradCorner)(idx_tmp,NodeArrays->NodeCorner[j].Get_connect(),NodeArrays->NodeCorner[j].Get_BcNormal());
		CalculNormal(idx_tmp,NodeArrays->NodeCorner[j].Get_connect(),NodeArrays->NodeCorner[j].Get_BcNormal());
	}

	for (int j=0;j<NodeArrays->NodeGlobalCorner.size();j++)
	{
	// Common variables
		idx_tmp=NodeArrays->NodeGlobalCorner[j].Get_index();
	// Calculate gradients
		(this->*PtrColourGradBc)(idx_tmp,NodeArrays->NodeGlobalCorner[j].Get_connect(),NodeArrays->NodeGlobalCorner[j].Get_BcNormal());
		CalculNormal(idx_tmp,NodeArrays->NodeGlobalCorner[j].Get_connect(),NodeArrays->NodeGlobalCorner[j].Get_BcNormal());
	}

	for (int j=0;j<NodeArrays->NodeVelocity.size();j++)
	{
	// Common variables
		idx_tmp=NodeArrays->NodeVelocity[j].Get_index();
	// Calculate gradients
		(this->*PtrColourGradBc)(idx_tmp,NodeArrays->NodeVelocity[j].Get_connect(),NodeArrays->NodeVelocity[j].Get_BcNormal());
		CalculNormal(idx_tmp,NodeArrays->NodeVelocity[j].Get_connect(),NodeArrays->NodeVelocity[j].Get_BcNormal());
	}

	for (int j=0;j<NodeArrays->NodePressure.size();j++)
	{
	// Common variables
		idx_tmp=NodeArrays->NodePressure[j].Get_index();
	// Calculate gradients
		(this->*PtrColourGradBc)(idx_tmp,NodeArrays->NodePressure[j].Get_connect(),NodeArrays->NodePressure[j].Get_BcNormal());
		CalculNormal(idx_tmp,NodeArrays->NodePressure[j].Get_connect(),NodeArrays->NodePressure[j].Get_BcNormal());
	}

	for (int j=0;j<NodeArrays->NodeWall.size();j++)
	{
	// Common variables
		idx_tmp=NodeArrays->NodeWall[j].Get_index();
/*		if(NodeArrays->NodeWall[j].get_x()==34&& NodeArrays->NodeWall[j].get_y()==66)
		{
			std::cout<<"RhoN (33,67): "<<RhoN[NodeArrays->NodeWall[j].Get_connect()[6]]<<std::endl;
			std::cout<<"RhoN (34,67): "<<RhoN[NodeArrays->NodeWall[j].Get_connect()[2]]<<std::endl;
			std::cout<<"RhoN (35,67): "<<RhoN[NodeArrays->NodeWall[j].Get_connect()[5]]<<std::endl;
			std::cout<<"RhoN (35,66): "<<RhoN[NodeArrays->NodeWall[j].Get_connect()[8]]<<std::endl;
		}*/
	// Calculate gradients
		(this->*PtrColourGradWall)(idx_tmp,NodeArrays->NodeWall[j].Get_connect(),NodeArrays->NodeWall[j].Get_BcNormal());
		CalculNormal(idx_tmp,NodeArrays->NodeWall[j].Get_connect(),NodeArrays->NodeWall[j].Get_BcNormal());
	}
	for (int j=0;j<NodeArrays->NodeSpecialWall.size();j++)
	{
	// Common variables
		idx_tmp=NodeArrays->NodeSpecialWall[j].Get_index();
	// Calculate gradients
		(this->*PtrColourGradCorner)(idx_tmp,NodeArrays->NodeSpecialWall[j].Get_connect(),NodeArrays->NodeSpecialWall[j].Get_BcNormal());
		CalculNormal(idx_tmp,NodeArrays->NodeSpecialWall[j].Get_connect(),NodeArrays->NodeSpecialWall[j].Get_BcNormal());
	}
	for (int j=0;j<NodeArrays->NodeSymmetry.size();j++)
	{
	// Common variables
		idx_tmp=NodeArrays->NodeSymmetry[j].Get_index();
	// Calculate gradients
		(this->*PtrColourGradBc)(idx_tmp,NodeArrays->NodeSymmetry[j].Get_connect(),NodeArrays->NodeSymmetry[j].Get_BcNormal());
		CalculNormal(idx_tmp,NodeArrays->NodeSymmetry[j].Get_connect(),NodeArrays->NodeSymmetry[j].Get_BcNormal());
	}
	for (int j=0;j<NodeArrays->NodePeriodic.size();j++)
	{
	// Common variables
		idx_tmp=NodeArrays->NodePeriodic[j].Get_index();
	// Calculate gradients
		(this->*PtrColourGrad)(idx_tmp,NodeArrays->NodePeriodic[j].Get_connect(),NodeArrays->NodePeriodic[j].Get_BcNormal());
		CalculNormal(idx_tmp,NodeArrays->NodePeriodic[j].Get_connect(),NodeArrays->NodePeriodic[j].Get_BcNormal());
	}
}
void D2Q9ColourFluid::ColourFluid_Collision()
{
//	double Ak=0.65;
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
		(this->*PtrCollisionCorner)(idx_tmp, NodeArrays->NodeCorner[j].Get_connect(),NodeArrays->NodeCorner[j].Get_BcNormal(),&fi_tmp[0]);
		(this->*PtrRecolour)(idx_tmp,&fi_tmp[0]);
	}

	for (int j=0;j<NodeArrays->NodeGlobalCorner.size();j++)
	{
	// Common variables
		idx_tmp=NodeArrays->NodeGlobalCorner[j].Get_index();
	//Model
		(this->*PtrCollisionCorner)(idx_tmp, NodeArrays->NodeGlobalCorner[j].Get_connect(),NodeArrays->NodeGlobalCorner[j].Get_BcNormal(),&fi_tmp[0]);
		(this->*PtrRecolour)(idx_tmp,&fi_tmp[0]);
	}

	for (int j=0;j<NodeArrays->NodeVelocity.size();j++)
	{
	// Common variables
		idx_tmp=NodeArrays->NodeVelocity[j].Get_index();
	//Model
		(this->*PtrCollisionBc)(idx_tmp, NodeArrays->NodeVelocity[j].Get_connect(),NodeArrays->NodeVelocity[j].Get_BcNormal(),&fi_tmp[0]);
		(this->*PtrRecolour)(idx_tmp,&fi_tmp[0]);
	}

	for (int j=0;j<NodeArrays->NodePressure.size();j++)
	{
	// Common variables
		idx_tmp=NodeArrays->NodePressure[j].Get_index();
		//Model
		(this->*PtrCollisionBc)(idx_tmp, NodeArrays->NodePressure[j].Get_connect(),NodeArrays->NodePressure[j].Get_BcNormal(),&fi_tmp[0]);
		(this->*PtrRecolour)(idx_tmp,&fi_tmp[0]);
	}

	for (int j=0;j<NodeArrays->NodeWall.size();j++)
	{
	// Common variables
		idx_tmp=NodeArrays->NodeWall[j].Get_index();
	//Model
		(this->*PtrCollisionBc)(idx_tmp, NodeArrays->NodeWall[j].Get_connect(),NodeArrays->NodeWall[j].Get_BcNormal(),&fi_tmp[0]);
		(this->*PtrRecolour)(idx_tmp,&fi_tmp[0]);
	}
	for (int j=0;j<NodeArrays->NodeSpecialWall.size();j++)
	{
	// Common variables
		idx_tmp=NodeArrays->NodeSpecialWall[j].Get_index();
	//Model
		(this->*PtrCollisionCorner)(idx_tmp, NodeArrays->NodeSpecialWall[j].Get_connect(),NodeArrays->NodeSpecialWall[j].Get_BcNormal(),&fi_tmp[0]);
		(this->*PtrRecolour)(idx_tmp,&fi_tmp[0]);
	}
	for (int j=0;j<NodeArrays->NodeSymmetry.size();j++)
	{
	// Common variables
		idx_tmp=NodeArrays->NodeSymmetry[j].Get_index();

	//Model
		(this->*PtrCollisionBc)(idx_tmp, NodeArrays->NodeSymmetry[j].Get_connect(),NodeArrays->NodeSymmetry[j].Get_BcNormal(),&fi_tmp[0]);
		(this->*PtrRecolour)(idx_tmp,&fi_tmp[0]);

	}
	for (int j=0;j<NodeArrays->NodePeriodic.size();j++)
	{
	// Common variables
		idx_tmp=NodeArrays->NodePeriodic[j].Get_index();

	//Model
		(this->*PtrCollisionBc)(idx_tmp, NodeArrays->NodePeriodic[j].Get_connect(),NodeArrays->NodePeriodic[j].Get_BcNormal(),&fi_tmp[0]);
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
/*
void D2Q9ColourFluid::CalculNormal_NoTeta(int & nodenumber, int* connect,int & normal){
	if(G_Norm[nodenumber]>0)
	{
		Normal[0][nodenumber]=G[0][nodenumber]/G_Norm[nodenumber];
		Normal[1][nodenumber]=G[1][nodenumber]/G_Norm[nodenumber];
	}

}
void D2Q9ColourFluid::CalculNormal_FixTeta(int & nodenumber, int* connect,int & normal){

	if(G_Norm[nodenumber]>0)
	{
		Normal[0][nodenumber]=G[0][nodenumber]/G_Norm[nodenumber];
		Normal[1][nodenumber]=G[1][nodenumber]/G_Norm[nodenumber];
	}

}
void D2Q9ColourFluid::CalculNormal_NoCstTeta(int & nodenumber, int* connect,int & normal){
	switch(normal)
	{
	case 1:
		Normal[0][nodenumber]=std::cos(teta[nodenumber]);
		Normal[1][nodenumber]=std::sin(teta[nodenumber]);
		break;
	case 2:
		Normal[0][nodenumber]=std::sin(teta[nodenumber]);
		Normal[1][nodenumber]=std::cos(teta[nodenumber]);
		break;
	case 3:
		Normal[0][nodenumber]=-std::cos(teta[nodenumber]);
		Normal[1][nodenumber]=std::sin(teta[nodenumber]);
		break;
	case 4:
		Normal[0][nodenumber]=std::sin(teta[nodenumber]);
		Normal[1][nodenumber]=-std::cos(teta[nodenumber]);
		break;
	default:
		std::cerr<<" Normal for Colour gradient non constant Teta is not found"<<std::endl;
	}
}
//Extrapolation of the normal of the interface at the first layer in the solid
void D2Q9ColourFluid::Extrapolate_NormalInSolid(){
	for (int j=0;j<NodeArrays->CornerConcave.size();j++)
	{
		ExtrapolationCornerConcaveToSolid(Normal[0],NodeArrays->NodeCorner[NodeArrays->CornerConcave[j]].Get_connect(),NodeArrays->NodeCorner[NodeArrays->CornerConcave[j]].Get_BcNormal());
		ExtrapolationCornerConcaveToSolid(Normal[1],NodeArrays->NodeCorner[NodeArrays->CornerConcave[j]].Get_connect(),NodeArrays->NodeCorner[NodeArrays->CornerConcave[j]].Get_BcNormal());
		Normalise(Normal[0],Normal[1],NodeArrays->NodeCorner[NodeArrays->CornerConcave[j]].Get_connect()[Opposite[NodeArrays->NodeCorner[NodeArrays->CornerConcave[j]].Get_BcNormal()]]);
	}
	for (int j=0;j<NodeArrays->NodeWall.size();j++)
	{
		ExtrapolationWallToSolid(Normal[0],NodeArrays->NodeWall[j].Get_connect(),NodeArrays->NodeWall[j].Get_BcNormal());
		ExtrapolationWallToSolid(Normal[1],NodeArrays->NodeWall[j].Get_connect(),NodeArrays->NodeWall[j].Get_BcNormal());
		Normalise(Normal[0],Normal[1],NodeArrays->NodeWall[j].Get_connect()[Opposite[NodeArrays->NodeWall[j].Get_BcNormal()]]);
	}
	for (int j=0;j<NodeArrays->CornerConvex.size();j++)
	{
		ExtrapolationCornerConvexToSolid(Normal[0],NodeArrays->NodeCorner[NodeArrays->CornerConvex[j]].Get_connect(),NodeArrays->NodeCorner[NodeArrays->CornerConvex[j]].Get_BcNormal());
		ExtrapolationCornerConvexToSolid(Normal[1],NodeArrays->NodeCorner[NodeArrays->CornerConvex[j]].Get_connect(),NodeArrays->NodeCorner[NodeArrays->CornerConvex[j]].Get_BcNormal());
		Normalise(Normal[0],Normal[1],NodeArrays->NodeCorner[NodeArrays->CornerConvex[j]].Get_connect()[Opposite[NodeArrays->NodeCorner[NodeArrays->CornerConvex[j]].Get_BcNormal()]]);
	}
}*/
void D2Q9ColourFluid::UpdateColourGradientOnWalls(){
	for (int j=0;j<NodeArrays->CornerConcave.size();j++)
	{
		UpdateColourGradient(NodeArrays->NodeCorner[NodeArrays->CornerConcave[j]].Get_connect(),NodeArrays->NodeCorner[NodeArrays->CornerConcave[j]].Get_BcNormal());
	}
	for (int j=0;j<NodeArrays->NodeWall.size();j++)
	{
		UpdateColourGradient(NodeArrays->NodeWall[j].Get_connect(),NodeArrays->NodeWall[j].Get_BcNormal());
	}
	for (int j=0;j<NodeArrays->CornerConvex.size();j++)
	{
		UpdateColourGradient(NodeArrays->NodeCorner[NodeArrays->CornerConvex[j]].Get_connect(),NodeArrays->NodeCorner[NodeArrays->CornerConvex[j]].Get_BcNormal());
	}
}
void D2Q9ColourFluid::UpdateColourGradient(int* connect, int & normal){
//Update the colour gradient
	G[0][connect[0]]=Normal[0][connect[0]]*G_Norm[connect[0]];
	G[1][connect[0]]=Normal[1][connect[0]]*G_Norm[connect[0]];
}
/*
void D2Q9ColourFluid::Impose_ContactAngleInSolid(){
	//define the contact angle in the solid
	for (int j=0;j<NodeArrays->CornerConcave.size();j++)
	{
		ContactAngleConcaveCornerInSolid(NodeArrays->NodeCorner[NodeArrays->CornerConcave[j]].Get_connect(),NodeArrays->NodeCorner[NodeArrays->CornerConcave[j]].Get_BcNormal());
	}
	for (int j=0;j<NodeArrays->NodeWall.size();j++)
	{
		ContactAngleWallInSolid(NodeArrays->NodeWall[j].Get_connect(),NodeArrays->NodeWall[j].Get_BcNormal());
	}
	for (int j=0;j<NodeArrays->CornerConvex.size();j++)
	{
		ContactAngleConvexCornerInSolid(NodeArrays->NodeCorner[NodeArrays->CornerConvex[j]].Get_connect(),NodeArrays->NodeCorner[NodeArrays->CornerConvex[j]].Get_BcNormal());
	}
	//interpolate the contact angle on the boundary from in the solid and in the interior fluid
	for (int j=0;j<NodeArrays->CornerConcave.size();j++)
	{
		ContactAngleConcaveCornerOnBc(NodeArrays->NodeCorner[NodeArrays->CornerConcave[j]].Get_connect(),NodeArrays->NodeCorner[NodeArrays->CornerConcave[j]].Get_BcNormal());
	}
	for (int j=0;j<NodeArrays->NodeWall.size();j++)
	{
		ContactAngleWallOnBc(NodeArrays->NodeWall[j].Get_connect(),NodeArrays->NodeWall[j].Get_BcNormal());
	}
	for (int j=0;j<NodeArrays->CornerConvex.size();j++)
	{
		ContactAngleConvexCornerOnBc(NodeArrays->NodeCorner[NodeArrays->CornerConvex[j]].Get_connect(),NodeArrays->NodeCorner[NodeArrays->CornerConvex[j]].Get_BcNormal());
	}
}
void D2Q9ColourFluid::Select_ContactAngle(int & normal,double & Nx, double & Ny){
	D1=std::sqrt((Nx-n1[normal][0])*(Nx-n1[normal][0])+(Ny-n1[normal][1])*(Ny-n1[normal][1]));
	D2=std::sqrt((Nx-n2[normal][0])*(Nx-n2[normal][0])+(Ny-n2[normal][1])*(Ny-n2[normal][1]));
	r=D1/(D1+D2);rMinus1=1.0-r;
	Nx=rMinus1*n1[normal][0]+r*n2[normal][0];
	Ny=rMinus1*n1[normal][1]+r*n2[normal][1];
}
void D2Q9ColourFluid::ContactAngleConcaveCornerInSolid(int* connect, int & normal){
	Select_ContactAngle(normal,Normal[0][connect[Opposite[normal]]],Normal[1][connect[Opposite[normal]]]);
	switch(normal)
	{
	case 5:
		Select_ContactAngle(IntRef(1),Normal[0][connect[3]],Normal[1][connect[3]]);
		Select_ContactAngle(IntRef(2),Normal[0][connect[4]],Normal[1][connect[4]]);

		break;
	case 6:
		Select_ContactAngle(IntRef(3),Normal[0][connect[1]],Normal[1][connect[1]]);
		Select_ContactAngle(IntRef(2),Normal[0][connect[4]],Normal[1][connect[4]]);
		break;
	case 7:
		Select_ContactAngle(IntRef(3),Normal[0][connect[1]],Normal[1][connect[1]]);
		Select_ContactAngle(IntRef(4),Normal[0][connect[2]],Normal[1][connect[2]]);
		break;
	case 8:
		Select_ContactAngle(IntRef(1),Normal[0][connect[3]],Normal[1][connect[3]]);
		Select_ContactAngle(IntRef(4),Normal[0][connect[2]],Normal[1][connect[2]]);
		break;
	}
}
void D2Q9ColourFluid::ContactAngleConvexCornerInSolid(int* connect, int & normal){
	Select_ContactAngle(normal,Normal[0][connect[Opposite[normal]]],Normal[1][connect[Opposite[normal]]]);
}
void D2Q9ColourFluid::ContactAngleWallInSolid(int* Connect, int & normal){
	Select_ContactAngle(normal,Normal[0][Connect[Opposite[normal]]],Normal[1][Connect[Opposite[normal]]]);
}


void D2Q9ColourFluid::ContactAngleConcaveCornerOnBc(int* connect, int & normal){
	Normal[0][connect[0]]=0.5*(Normal[0][connect[normal]]+Normal[0][connect[Opposite[normal]]]);
	Normal[1][connect[0]]=0.5*(Normal[1][connect[normal]]+Normal[1][connect[Opposite[normal]]]);
	// Normalise
	D_tmp=sqrt(Normal[0][connect[0]]*Normal[0][connect[0]]+Normal[1][connect[0]]*Normal[1][connect[0]]);
	if(D_tmp>0)
		{Normal[0][connect[0]]/=D_tmp;Normal[1][connect[0]]/=D_tmp;}
	else
		{Normal[0][connect[0]]=0.0; Normal[1][connect[0]]=0.0;}
	//Update the colour gradient
	G[0][connect[0]]=Normal[0][connect[0]]*G_Norm[connect[0]];
	G[1][connect[0]]=Normal[1][connect[0]]*G_Norm[connect[0]];
}
void D2Q9ColourFluid::ContactAngleConvexCornerOnBc(int* connect, int & normal){
	Normal[0][connect[0]]=0.5*(Normal[0][connect[normal]]+Normal[0][connect[Opposite[normal]]]);
	Normal[1][connect[0]]=0.5*(Normal[1][connect[normal]]+Normal[1][connect[Opposite[normal]]]);
	// Normalise
	D_tmp=sqrt(Normal[0][connect[0]]*Normal[0][connect[0]]+Normal[1][connect[0]]*Normal[1][connect[0]]);
	if(D_tmp>0)
		{Normal[0][connect[0]]/=D_tmp;Normal[1][connect[0]]/=D_tmp;}
	else
		{Normal[0][connect[0]]=0.0; Normal[1][connect[0]]=0.0;}
	//Update the colour gradient
	G[0][connect[0]]=Normal[0][connect[0]]*G_Norm[connect[0]];
	G[1][connect[0]]=Normal[1][connect[0]]*G_Norm[connect[0]];
}
void D2Q9ColourFluid::ContactAngleWallOnBc(int* connect, int & normal){
	Normal[0][connect[0]]=0.5*(Normal[0][connect[normal]]+Normal[0][connect[Opposite[normal]]]);
	Normal[1][connect[0]]=0.5*(Normal[1][connect[normal]]+Normal[1][connect[Opposite[normal]]]);
	// Normalise
	D_tmp=sqrt(Normal[0][connect[0]]*Normal[0][connect[0]]+Normal[1][connect[0]]*Normal[1][connect[0]]);
	if(D_tmp>0)
		{Normal[0][connect[0]]/=D_tmp;Normal[1][connect[0]]/=D_tmp;}
	else
		{Normal[0][connect[0]]=0.0; Normal[1][connect[0]]=0.0;}
	//Update the colour gradient
	G[0][connect[0]]=Normal[0][connect[0]]*G_Norm[connect[0]];
	G[1][connect[0]]=Normal[1][connect[0]]*G_Norm[connect[0]];
}
*/
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
	Grad(&tmp[0],&RhoN[0],connect,normal);
	G[0][nodenumber]=tmp[0];G[1][nodenumber]=tmp[1];
	G_Norm[nodenumber]=std::sqrt(G[0][nodenumber]*G[0][nodenumber]+G[1][nodenumber]*G[1][nodenumber]);
}
void D2Q9ColourFluid::Colour_gradient_DensityNormalGradBc(int & nodenumber, int* connect,int & normal){
	double tmp[2];
	GradBc(&tmp[0],&RhoN[0],connect,normal);
	G[0][nodenumber]=tmp[0];G[1][nodenumber]=tmp[1];
	G_Norm[nodenumber]=std::sqrt(G[0][nodenumber]*G[0][nodenumber]+G[1][nodenumber]*G[1][nodenumber]);
}
void D2Q9ColourFluid::Colour_gradient_DensityNormalGradCorner(int & nodenumber, int* connect,int & normal){
	double tmp[2];
	GradCorner(&tmp[0],&RhoN[0],connect,normal);
	G[0][nodenumber]=tmp[0];G[1][nodenumber]=tmp[1];
	G_Norm[nodenumber]=std::sqrt(G[0][nodenumber]*G[0][nodenumber]+G[1][nodenumber]*G[1][nodenumber]);
}
void D2Q9ColourFluid::CalculNormal(int & nodenumber, int* connect,int & normal){
	//Normalise grad(RhoN)
	if(G_Norm[nodenumber]>0)
	{
		Normal[0][nodenumber]=G[0][nodenumber]/G_Norm[nodenumber];
		Normal[1][nodenumber]=G[1][nodenumber]/G_Norm[nodenumber];
	}
}

void D2Q9ColourFluid::Recolouring_Latva(int & nodenumber, double * fi_tmp){
	if(G_Norm[nodenumber]<LimitGNorm)
	{
		if(Rhor[nodenumber]<Rhob[nodenumber])
			for(int i=0;i<nbvelo;i++)
				{f[0]->f[i][nodenumber]=0.0;f[1]->f[i][nodenumber]=fi_tmp[i];}
		else
			for(int i=0;i<nbvelo;i++)
				{f[0]->f[i][nodenumber]=fi_tmp[i];f[1]->f[i][nodenumber]=0.0;}
	}
	else
	{

		double Rhor_Rho=Rhor[nodenumber]/Rho[nodenumber];
		double Rhob_Rho=Rhob[nodenumber]/Rho[nodenumber];
		double factor=beta*Rhor[nodenumber]*Rhob[nodenumber]/Rho[nodenumber];

		f[0]->f[0][nodenumber]=Rhor_Rho*fi_tmp[0];
		f[1]->f[0][nodenumber]=Rhob_Rho*fi_tmp[0];
		for(int i=1;i<nbvelo;i++)
		{
				f[0]->f[i][nodenumber]=Rhor_Rho*fi_tmp[i]
							+factor*omega[i]*CosPhi(nodenumber,i,G_Norm[nodenumber]);
				f[1]->f[i][nodenumber]=fi_tmp[i]-f[0]->f[i][nodenumber];//Rhob_Rho*fi_tmp[i]
		}
	}
}
void D2Q9ColourFluid::Recolouring_LatvaWithEstimator(int & nodenumber, double * fi_tmp){
	if(G_Norm[nodenumber]<LimitGNorm)
	{
		if(Rhor[nodenumber]<Rhob[nodenumber])
			for(int i=0;i<nbvelo;i++)
				{f[0]->f[i][nodenumber]=0.0;f[1]->f[i][nodenumber]=fi_tmp[i];}
		else
			for(int i=0;i<nbvelo;i++)
				{f[0]->f[i][nodenumber]=fi_tmp[i];f[1]->f[i][nodenumber]=0.0;}
	}
	else
	{

		double RhoNStar=RhoN[nodenumber]*(1.0+G[0][nodenumber]*U[0][nodenumber]+G[1][nodenumber]*U[1][nodenumber]);
		double RhorStar=Rho[nodenumber]*(RhoNStar+1.0)*0.5;
		double RhobStar=Rho[nodenumber]*(1.0-RhoNStar)*0.5;
		double Rhor_Rho=RhorStar/Rho[nodenumber];
		double Rhob_Rho=RhobStar/Rho[nodenumber];
		double factor=beta*RhorStar*RhobStar/Rho[nodenumber];
		f[0]->f[0][nodenumber]=Rhor_Rho*fi_tmp[0];
		f[1]->f[0][nodenumber]=Rhob_Rho*fi_tmp[0];
		for(int i=1;i<nbvelo;i++)
		{
				f[0]->f[i][nodenumber]=Rhor_Rho*fi_tmp[i]
							+factor*omega[i]*CosPhi(nodenumber,i,G_Norm[nodenumber]);
				f[1]->f[i][nodenumber]=fi_tmp[i]-f[0]->f[i][nodenumber];//Rhob_Rho*fi_tmp[i]
		}
	}
}
/*void D2Q9ColourFluid::Recolouring_Wall(int & nodenumber, double * fi_tmp){
		f[0]->f[0][nodenumber]=Rhor[nodenumber]*fi_tmp[0]/Rho[nodenumber];
		f[1]->f[0][nodenumber]=Rhob[nodenumber]*fi_tmp[0]/Rho[nodenumber];
		if(G_Norm[nodenumber]>0)
			for(int i=1;i<nbvelo;i++)
			{
					f[0]->f[i][nodenumber]=Rhor[nodenumber]*fi_tmp[i]/Rho[nodenumber]
								+beta*omega[i]*Rhor[nodenumber]*Rhob[nodenumber]*CosPhi(nodenumber,i,G_Norm[nodenumber])/Rho[nodenumber];
					f[1]->f[i][nodenumber]=Rhob[nodenumber]*fi_tmp[i]/Rho[nodenumber]
								-beta*omega[i]*Rhor[nodenumber]*Rhob[nodenumber]*CosPhi(nodenumber,i,G_Norm[nodenumber])/Rho[nodenumber];
			}
		else
			for(int i=1;i<nbvelo;i++)
			{
				f[0]->f[i][nodenumber]=Rhor[nodenumber]*fi_tmp[i]/Rho[nodenumber];
				f[1]->f[i][nodenumber]=Rhob[nodenumber]*fi_tmp[i]/Rho[nodenumber];
			}
		if(Rhor[nodenumber]<=Rho_limiter)
			for(int i=0;i<nbvelo;i++)
				{f[0]->f[i][nodenumber]=0.0;f[1]->f[i][nodenumber]=fi_tmp[i];}
		if(Rhob[nodenumber]<=Rho_limiter)
			for(int i=0;i<nbvelo;i++)
				{f[0]->f[i][nodenumber]=fi_tmp[i];;f[1]->f[i][nodenumber]=0.0;}



}*/
double D2Q9ColourFluid::CosPhi(int nodenumber, int & direction,double & F_Norm){
	return (Normal[0][nodenumber]* Ei[direction][0]+Normal[1][nodenumber]* Ei[direction][1])/(Ei_Norm[direction]);//*G_Norm[nodenumber]///(F_Norm*std::sqrt(Ei[direction][0]*Ei[direction][0]+Ei[direction][1]*Ei[direction][1]));

//	return (G[0][nodenumber]* Ei[direction][0]+G[1][nodenumber]* Ei[direction][1])/(Ei_Norm[direction]*G_Norm[nodenumber]);//*G_Norm[nodenumber]///(F_Norm*std::sqrt(Ei[direction][0]*Ei[direction][0]+Ei[direction][1]*Ei[direction][1]));
}

double& D2Q9ColourFluid::Collision_operator_Grunau(int & i, int & nodenumber, double Ak){
	if(G_Norm[nodenumber]>0)
		D_tmp=Ak*G_Norm[nodenumber]*(((G[0][nodenumber]*Ei[i][0]+G[1][nodenumber]*Ei[i][1])/G_Norm[nodenumber])*((G[0][nodenumber]*Ei[i][0]+G[1][nodenumber]*Ei[i][1])/G_Norm[nodenumber])-3/4);
	else
		D_tmp=0;
	return D_tmp;
}
double& D2Q9ColourFluid::Collision_operator_Reis(int & i, int & nodenumber, double Ak){
	if(G_Norm[nodenumber]>0)
		D_tmp=Ak*G_Norm[nodenumber]*((omega[i]*(G[0][nodenumber]*Ei[i][0]+G[1][nodenumber]*Ei[i][1])/G_Norm[nodenumber])*((G[0][nodenumber]*Ei[i][0]+G[1][nodenumber]*Ei[i][1])/G_Norm[nodenumber])-Bi[i]);
	else
		D_tmp=0;
	return D_tmp;
}
void D2Q9ColourFluid::Collision_Grunau(int & nodenumber, int* connect,int & normal,double* fi){
//	int i0=0;
	double Fx=0,Fy=0;
//	fi[0]=f[0]->f[0][nodenumber]+f[1]->f[0][nodenumber];
//	DVec_2D_tmp[0]=0;DVec_2D_tmp[1]=0;
//	Collide_2D(i0, fi[0],Rho[nodenumber], U[0][nodenumber], U[1][nodenumber], DVec_2D_tmp[0],DVec_2D_tmp[1], Get_InvTau(Rho[nodenumber],RhoN[nodenumber]));
	for (int i=0;i<9;i++)
//	for (int i=1;i<9;i++)
	{
		//Save the mixture distribution for recolouring
		fi[i]=f[0]->f[i][nodenumber]+f[1]->f[i][nodenumber];
		Fi[i]=Collision_operator_Grunau(i, nodenumber, A1);
//		Collide_2D(i, fi[i],Rho[nodenumber], U[0][nodenumber], U[1][nodenumber], Collision_operator_Grunau(i, nodenumber, A1),DVec_2D_tmp[1], Get_InvTau(Rho[nodenumber],RhoN[nodenumber]));
	}
	Collide_2D_V2(fi,Rho[nodenumber], U[0][nodenumber], U[1][nodenumber], Fi, Fx,Fy, Get_InvTau(Rho[nodenumber],RhoN[nodenumber]),Get_Mu(Rho[nodenumber],RhoN[nodenumber]));
}
void D2Q9ColourFluid::Collision_Reis(int & nodenumber, int* connect,int & normal,double* fi){
	//int i0=0;
	double invtau=Get_InvTau(Rho[nodenumber],RhoN[nodenumber]);
	double Ak=A1*invtau;
	double Fx=0,Fy=0;
	for (int i=0;i<9;i++)
	{
		//Save the mixture distribution for recolouring
		fi[i]=f[0]->f[i][nodenumber]+f[1]->f[i][nodenumber];
		Fi[i]=Collision_operator_Reis(i, nodenumber, Ak);
//		Collide_2D(i, fi[i],Rho[nodenumber], U[0][nodenumber], U[1][nodenumber], Collision_operator_Reis(i, nodenumber, Ak),DVec_2D_tmp[1], invtau);
	}
	Collide_2D_V2(fi,Rho[nodenumber], U[0][nodenumber], U[1][nodenumber], Fi, Fx,Fy, invtau,Get_Mu(Rho[nodenumber],RhoN[nodenumber]));
}
void D2Q9ColourFluid::Collision_SurfaceForce(int & nodenumber, int* connect,int & normal,double* fi){
	if(G_Norm[nodenumber]>0)
	{SurfaceForce(nodenumber,connect,normal,F[0][nodenumber],F[1][nodenumber]);}
	else
		{F[0][nodenumber]=0;F[1][nodenumber]=0;	}
//	U[0][nodenumber]=U[0][nodenumber]+0.5*F[0][nodenumber]/Rho[nodenumber];
//	U[1][nodenumber]=U[1][nodenumber]+0.5*F[1][nodenumber]/Rho[nodenumber];
	for (int i=0;i<9;i++)
	{
		//Save the mixture distribution for recolouring
		fi[i]=f[0]->f[i][nodenumber]+f[1]->f[i][nodenumber];
		Fi[i]=0;
//		Collide_2D(i, fi[i],Rho[nodenumber], U[0][nodenumber], U[1][nodenumber], F[0][nodenumber],F[1][nodenumber], Get_InvTau(Rho[nodenumber],RhoN[nodenumber]));
	}
	Collide_2D_V2(fi,Rho[nodenumber], U[0][nodenumber], U[1][nodenumber], Fi, F[0][nodenumber],F[1][nodenumber], Get_InvTau(Rho[nodenumber],RhoN[nodenumber]),Get_Mu(Rho[nodenumber],RhoN[nodenumber]));
}

 void D2Q9ColourFluid::Collision_SurfaceForceBc(int & nodenumber, int* connect,int & normal,double* fi){
	 	if(G_Norm[nodenumber]>0)
	 	{
		SurfaceForceBc(nodenumber,connect,normal,F[0][nodenumber],F[1][nodenumber]);}
	else
		{F[0][nodenumber]=0;F[1][nodenumber]=0;}
//	U[0][nodenumber]=U[0][nodenumber]+0.5*F[0][nodenumber]/Rho[nodenumber];
//	U[1][nodenumber]=U[1][nodenumber]+0.5*F[1][nodenumber]/Rho[nodenumber];
	for (int i=0;i<9;i++)
	{
		//Save the mixture distribution for recolouring
		fi[i]=f[0]->f[i][nodenumber]+f[1]->f[i][nodenumber];
		Fi[i]=0;
//		Collide_2D(i, fi[i],Rho[nodenumber], U[0][nodenumber], U[1][nodenumber], F[0][nodenumber],F[1][nodenumber], Get_InvTau(Rho[nodenumber],RhoN[nodenumber]));
	}
	Collide_2D_V2(fi,Rho[nodenumber], U[0][nodenumber], U[1][nodenumber], Fi, F[0][nodenumber],F[1][nodenumber], Get_InvTau(Rho[nodenumber],RhoN[nodenumber]),Get_Mu(Rho[nodenumber],RhoN[nodenumber]));
}
 void D2Q9ColourFluid::Collision_SurfaceForceCorner(int & nodenumber, int* connect,int & normal,double* fi){
	if(G_Norm[nodenumber]>0)
	{	SurfaceForceCorner(nodenumber,connect,normal,F[0][nodenumber],F[1][nodenumber]);}
	else
		{F[0][nodenumber]=0;F[1][nodenumber]=0;	}
//	U[0][nodenumber]=U[0][nodenumber]+0.5*F[0][nodenumber]/Rho[nodenumber];
//	U[1][nodenumber]=U[1][nodenumber]+0.5*F[1][nodenumber]/Rho[nodenumber];
	for (int i=0;i<9;i++)
	{
		//Save the mixture distribution for recolouring
		fi[i]=f[0]->f[i][nodenumber]+f[1]->f[i][nodenumber];
		Fi[i]=0;
//		Collide_2D(i, fi[i],Rho[nodenumber], U[0][nodenumber], U[1][nodenumber], F[0][nodenumber],F[1][nodenumber], Get_InvTau(Rho[nodenumber],RhoN[nodenumber]));
	}
	Collide_2D_V2(fi,Rho[nodenumber], U[0][nodenumber], U[1][nodenumber], Fi, F[0][nodenumber],F[1][nodenumber], Get_InvTau(Rho[nodenumber],RhoN[nodenumber]),Get_Mu(Rho[nodenumber],RhoN[nodenumber]));
}
 void D2Q9ColourFluid::SurfaceForce(int & nodenumber, int* connect,int & normal,double & Fx,double & Fy){
 	Curv[nodenumber]=Curvature(nodenumber,connect,normal);
 	D_tmp=0.5*tension*Curv[nodenumber];//*G_Norm[nodenumber]; //G_Norm is to get back the density gradient and not the normalise one
 	Fx=D_tmp*G[0][nodenumber];
 	Fy=D_tmp*G[1][nodenumber];
 }
void D2Q9ColourFluid::SurfaceForceBc(int & nodenumber, int* connect,int & normal,double & Fx,double & Fy){

	Curv[nodenumber]=CurvatureBc(nodenumber,connect,normal);
	D_tmp=0.5*tension*Curv[nodenumber];
	Fx=D_tmp*G[0][nodenumber];
	Fy=D_tmp*G[1][nodenumber];
}
void D2Q9ColourFluid::SurfaceForceCorner(int & nodenumber, int* connect,int & normal,double & Fx,double & Fy){

	Curv[nodenumber]=CurvatureCorner(nodenumber,connect,normal);
	D_tmp=0.5*tension*Curv[nodenumber];
	Fx=D_tmp*G[0][nodenumber];
	Fy=D_tmp*G[1][nodenumber];
}
double D2Q9ColourFluid::Curvature(int & nodenumber, int* connect,int & normal){
	double tmp_0[2],tmp_1[2];
	Grad(&tmp_0[0],&Normal[0][0],connect,normal);//d(Nx)/dx and d(Nx)/dy
	Grad(&tmp_1[0],&Normal[1][0],connect,normal);//d(Ny)/dx and d(Ny)/dy
	return -(tmp_0[0]+tmp_1[1]);//Dot product

}
double D2Q9ColourFluid::CurvatureBc(int & nodenumber, int* connect,int & normal){
	double tmp_0[2],tmp_1[2];
	GradBc(&tmp_0[0],&Normal[0][0],connect,normal);//d(Nx)/dx and d(Nx)/dy
	GradBc(&tmp_1[0],&Normal[1][0],connect,normal);//d(Ny)/dx and d(Ny)/dy
	return -(tmp_0[0]+tmp_1[1]);//Dot product
}
double D2Q9ColourFluid::CurvatureCorner(int & nodenumber, int* connect,int & normal){
	double tmp_0[2],tmp_1[2];
	GradCorner(&tmp_0[0],&Normal[0][0],connect,normal);//d(Nx)/dx and d(Nx)/dy
	GradCorner(&tmp_1[0],&Normal[1][0],connect,normal);//d(Ny)/dx and d(Ny)/dy
	return -(tmp_0[0]+tmp_1[1]);//Dot product
}
//
void D2Q9ColourFluid::ApplyPatchVelocity(VelocityPatchBc& VelPatchBc){
	SetVelocity(VelPatchBc.Get_VelocityModel(),VelPatchBc.Get_VelocityType());
	std::vector<int> NodeIdx=VelPatchBc.Get_NodeIndexByType();
	std::vector<int> NodeIdxSpecialWalls=VelPatchBc.Get_NodeIndexByTypeSpecialWalls();
	std::vector<int> NodeIdxGlobalCorner=VelPatchBc.Get_NodeIndexByTypeGlobalCorner();
	double alpha=0;
	if(VelPatchBc.Get_extrapolationAlpha())
	{
		for (int j=0;j<NodeIdx.size();j++)
		{
			ApplyVelocity(NodeArrays->NodeVelocity[NodeIdx[j]].Get_BcNormal(),NodeArrays->NodeVelocity[NodeIdx[j]].Get_connect(),NodeArrays->NodeVelocity[NodeIdx[j]].Get_UDef(), f[0],Rhor,U[0],U[1]);
			ApplyVelocity(NodeArrays->NodeVelocity[NodeIdx[j]].Get_BcNormal(),NodeArrays->NodeVelocity[NodeIdx[j]].Get_connect(),NodeArrays->NodeVelocity[NodeIdx[j]].Get_UDef(), f[1],Rhob,U[0],U[1]);
			//Impose alpha
			ExtrapolationOnWall(RhoN,NodeArrays->NodeVelocity[NodeIdx[j]].Get_connect(),NodeArrays->NodeVelocity[NodeIdx[j]].Get_BcNormal());
			alpha=(RhoN[NodeArrays->NodeVelocity[NodeIdx[j]].Get_index()]+1.0)*0.5;
			ImposeAlpha(NodeArrays->NodeVelocity[NodeIdx[j]].Get_index(), alpha);
		}
		for (int j=0;j<NodeIdxSpecialWalls.size();j++)
		{
//			ExtrapolationOnCornerConcave(Rho,NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]].Get_connect(),NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]].Get_BcNormal());
//			ExtrapolationOnCornerConcave(Rhor,NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]].Get_connect(),NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]].Get_BcNormal());
//			ExtrapolationOnCornerConcave(Rhob,NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]].Get_connect(),NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]].Get_BcNormal());
			ExtrapolationOnCornerConcave(RhoN,NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]].Get_connect(),NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]].Get_BcNormal());
			alpha=(RhoN[NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]].Get_index()]+1.0)*0.5;
			ApplyVelocityWall(NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]],U[0][NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]].Get_index()],U[1][NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]].Get_index()],NodeArrays->TypeOfNode,f[0],Rhor,U[0],U[1]);
			ApplyVelocityWall(NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]],U[0][NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]].Get_index()],U[1][NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]].Get_index()],NodeArrays->TypeOfNode,f[1],Rhob,U[0],U[1]);
			ImposeAlpha(NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]].Get_index(), alpha);
		}
		for (int j=0;j<NodeIdxGlobalCorner.size();j++)
		{
			ApplyGlobalCorner(NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]],NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]].Get_AlphaDef()*NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]].Get_RhoDef(),NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]].Get_UDef(),NodeArrays->TypeOfNode,f[0],Rhor,U[0],U[1],Rhor[NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]].Get_index()]/Rho[NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]].Get_index()]);
			ApplyGlobalCorner(NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]],(1-NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]].Get_AlphaDef())*NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]].Get_RhoDef(),NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]].Get_UDef(),NodeArrays->TypeOfNode,f[1],Rhob,U[0],U[1],Rhob[NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]].Get_index()]/Rho[NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]].Get_index()]);
		}
	}
	else
	{
		for (int j=0;j<NodeIdx.size();j++)
		{
			ApplyVelocity(NodeArrays->NodeVelocity[NodeIdx[j]].Get_BcNormal(),NodeArrays->NodeVelocity[NodeIdx[j]].Get_connect(),NodeArrays->NodeVelocity[NodeIdx[j]].Get_UDef(), f[0],Rhor,U[0],U[1]);
			ApplyVelocity(NodeArrays->NodeVelocity[NodeIdx[j]].Get_BcNormal(),NodeArrays->NodeVelocity[NodeIdx[j]].Get_connect(),NodeArrays->NodeVelocity[NodeIdx[j]].Get_UDef(), f[1],Rhob,U[0],U[1]);
			//Impose alpha
			ImposeAlpha(NodeArrays->NodeVelocity[NodeIdx[j]].Get_index(), NodeArrays->NodeVelocity[NodeIdx[j]].Get_AlphaDef());
		}
		for (int j=0;j<NodeIdxSpecialWalls.size();j++)
		{
//			ExtrapolationOnCornerConcave(Rho,NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]].Get_connect(),NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]].Get_BcNormal());
//			ExtrapolationOnCornerConcave(Rhor,NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]].Get_connect(),NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]].Get_BcNormal());
//			ExtrapolationOnCornerConcave(Rhob,NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]].Get_connect(),NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]].Get_BcNormal());
			ExtrapolationOnCornerConcave(RhoN,NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]].Get_connect(),NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]].Get_BcNormal());
			alpha=(RhoN[NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]].Get_index()]+1.0)*0.5;
			ApplyVelocityWall(NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]],U[0][NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]].Get_index()],U[1][NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]].Get_index()],NodeArrays->TypeOfNode,f[0],Rhor,U[0],U[1]);
			ApplyVelocityWall(NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]],U[0][NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]].Get_index()],U[1][NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]].Get_index()],NodeArrays->TypeOfNode,f[1],Rhob,U[0],U[1]);
			ImposeAlpha(NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]].Get_index(), alpha);
		}
		for (int j=0;j<NodeIdxGlobalCorner.size();j++)
		{
			ApplyGlobalCorner(NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]],NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]].Get_AlphaDef()*NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]].Get_RhoDef(),NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]].Get_UDef(),NodeArrays->TypeOfNode,f[0],Rhor,U[0],U[1],Rhor[NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]].Get_index()]/Rho[NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]].Get_index()]);
			ApplyGlobalCorner(NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]],(1-NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]].Get_AlphaDef())*NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]].Get_RhoDef(),NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]].Get_UDef(),NodeArrays->TypeOfNode,f[1],Rhob,U[0],U[1],Rhob[NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]].Get_index()]/Rho[NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]].Get_index()]);
		}
	}
}
void D2Q9ColourFluid::ApplyPatchPressure(PressurePatchBc& PresPatchBc){
	SetPressure(PresPatchBc.Get_PressureModel(),PresPatchBc.Get_PressureType());
	std::vector<int> NodeIdx=PresPatchBc.Get_NodeIndexByType();
	std::vector<int> NodeIdxSpecialWalls=PresPatchBc.Get_NodeIndexByTypeSpecialWalls();
	std::vector<int> NodeIdxGlobalCorner=PresPatchBc.Get_NodeIndexByTypeGlobalCorner();
	double alpha=0;
	if(PresPatchBc.Get_extrapolationAlpha())
	{
		for (int j=0;j<NodeIdx.size();j++)
		{
			ExtrapolationOnWall(RhoN,NodeArrays->NodePressure[NodeIdx[j]].Get_connect(),NodeArrays->NodePressure[NodeIdx[j]].Get_BcNormal());
			alpha=(RhoN[NodeArrays->NodePressure[NodeIdx[j]].Get_index()]+1.0)*0.5;
			ApplyPressure(NodeArrays->NodePressure[NodeIdx[j]].Get_BcNormal(),NodeArrays->NodePressure[NodeIdx[j]].Get_connect(),alpha*NodeArrays->NodePressure[NodeIdx[j]].Get_RhoDef(), f[0],Rhor,U[0],U[1]);
			ApplyPressure(NodeArrays->NodePressure[NodeIdx[j]].Get_BcNormal(),NodeArrays->NodePressure[NodeIdx[j]].Get_connect(),(1-alpha)*NodeArrays->NodePressure[NodeIdx[j]].Get_RhoDef(), f[1],Rhob,U[0],U[1]);
		}
		for (int j=0;j<NodeIdxSpecialWalls.size();j++)
		{
			ExtrapolationOnCornerConcave(Rho,NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]].Get_connect(),NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]].Get_BcNormal());
			ExtrapolationOnCornerConcave(RhoN,NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]].Get_connect(),NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]].Get_BcNormal());
			alpha=(RhoN[NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]].Get_index()]+1.0)*0.5;
			ApplyPressureWall(NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]],alpha*Rho[NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]].Get_index()],NodeArrays->TypeOfNode,f[0],Rhor,U[0],U[1]);
			ApplyPressureWall(NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]],(1.0-alpha)*Rho[NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]].Get_index()],NodeArrays->TypeOfNode,f[1],Rhob,U[0],U[1]);
		}
		for (int j=0;j<NodeIdxGlobalCorner.size();j++)
		{
			ApplyGlobalCorner(NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]],NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]].Get_AlphaDef()*NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]].Get_RhoDef(),NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]].Get_UDef(),NodeArrays->TypeOfNode,f[0],Rhor,U[0],U[1],Rhor[NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]].Get_index()]/Rho[NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]].Get_index()]);
			ApplyGlobalCorner(NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]],(1-NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]].Get_AlphaDef())*NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]].Get_RhoDef(),NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]].Get_UDef(),NodeArrays->TypeOfNode,f[1],Rhob,U[0],U[1],Rhob[NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]].Get_index()]/Rho[NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]].Get_index()]);
/*		//Keep alpha
			Rhor[NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]].Get_index()]=0.0;
			Rhob[NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]].Get_index()]=0.0;
			for(int i=0;i<9;i++)
			{
				fi_tmp=f[0]->f[i][NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]].Get_index()]+f[1]->f[i][NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]].Get_index()];
				f[0]->f[i][NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]].Get_index()]=NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]].Get_AlphaDef()*fi_tmp;
				f[1]->f[i][NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]].Get_index()]=(1-NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]].Get_AlphaDef())*fi_tmp;
				Rhor[NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]].Get_index()]+=f[0]->f[i][NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]].Get_index()];
				Rhob[NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]].Get_index()]+=f[1]->f[i][NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]].Get_index()];
			}*/
		}
	}
	else
	{
		for (int j=0;j<NodeIdx.size();j++)
		{
			alpha=NodeArrays->NodePressure[NodeIdx[j]].Get_AlphaDef();
			ApplyPressure(NodeArrays->NodePressure[NodeIdx[j]].Get_BcNormal(),NodeArrays->NodePressure[NodeIdx[j]].Get_connect(),alpha*NodeArrays->NodePressure[NodeIdx[j]].Get_RhoDef(), f[0],Rhor,U[0],U[1]);
			ApplyPressure(NodeArrays->NodePressure[NodeIdx[j]].Get_BcNormal(),NodeArrays->NodePressure[NodeIdx[j]].Get_connect(),(1-alpha)*NodeArrays->NodePressure[NodeIdx[j]].Get_RhoDef(), f[1],Rhob,U[0],U[1]);
		}
		for (int j=0;j<NodeIdxSpecialWalls.size();j++)
		{
			ExtrapolationOnCornerConcave(Rho,NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]].Get_connect(),NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]].Get_BcNormal());
//			ExtrapolationOnCornerConcave(Rhob,NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]].Get_connect(),NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]].Get_BcNormal());
			alpha=(RhoN[NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]].Get_index()]+1.0)*0.5;
			ApplyPressureWall(NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]],alpha*Rho[NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]].Get_index()],NodeArrays->TypeOfNode,f[0],Rhor,U[0],U[1]);
			ApplyPressureWall(NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]],(1.0-alpha)*Rho[NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]].Get_index()],NodeArrays->TypeOfNode,f[1],Rhob,U[0],U[1]);
		}
		for (int j=0;j<NodeIdxGlobalCorner.size();j++)
		{
			ApplyGlobalCorner(NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]],NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]].Get_AlphaDef()*NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]].Get_RhoDef(),NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]].Get_UDef(),NodeArrays->TypeOfNode,f[0],Rhor,U[0],U[1],Rhor[NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]].Get_index()]/Rho[NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]].Get_index()]);
			ApplyGlobalCorner(NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]],(1-NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]].Get_AlphaDef())*NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]].Get_RhoDef(),NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]].Get_UDef(),NodeArrays->TypeOfNode,f[1],Rhob,U[0],U[1],Rhob[NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]].Get_index()]/Rho[NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]].Get_index()]);
		}
	}
}
void D2Q9ColourFluid::ApplyPatchSymmetry(SymmetryPatchBc& SymPatchBc){
	SetSymmetry(SymPatchBc.Get_SymmetryType());
	std::vector<int> NodeIdx=SymPatchBc.Get_NodeIndexByType();
	std::vector<int> NodeIdxSpecialWalls=SymPatchBc.Get_NodeIndexByTypeSpecialWalls();
	std::vector<int> NodeIdxGlobalCorner=SymPatchBc.Get_NodeIndexByTypeGlobalCorner();
	double alpha=0;
	if(SymPatchBc.Get_extrapolationAlpha())
	{
		for (int j=0;j<NodeIdx.size();j++)
		{
				ApplySymmetry(NodeArrays->NodeSymmetry[NodeIdx[j]].Get_BcNormal(),NodeArrays->NodeSymmetry[NodeIdx[j]].Get_connect(),NodeArrays->NodeSymmetry[NodeIdx[j]].Get_RhoDef(),NodeArrays->NodeSymmetry[NodeIdx[j]].Get_UDef(),f[0],Rhor,U[0],U[1]);
				ApplySymmetry(NodeArrays->NodeSymmetry[NodeIdx[j]].Get_BcNormal(),NodeArrays->NodeSymmetry[NodeIdx[j]].Get_connect(),NodeArrays->NodeSymmetry[NodeIdx[j]].Get_RhoDef(),NodeArrays->NodeSymmetry[NodeIdx[j]].Get_UDef(),f[1],Rhob,U[0],U[1]);
		}
		for (int j=0;j<NodeIdxSpecialWalls.size();j++)
		{
			ExtrapolationOnCornerConcave(Rho,NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]].Get_connect(),NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]].Get_BcNormal());
			ExtrapolationOnCornerConcave(Rhor,NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]].Get_connect(),NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]].Get_BcNormal());
			ExtrapolationOnCornerConcave(Rhob,NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]].Get_connect(),NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]].Get_BcNormal());
			alpha=(RhoN[NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]].Get_index()]+1.0)*0.5;
			ApplySymmetryWall(NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]],alpha*Rho[NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]].Get_index()],U[0][NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]].Get_index()],U[1][NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]].Get_index()],NodeArrays->TypeOfNode,f[0],Rhor,U[0],U[1]);
			ApplySymmetryWall(NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]],(1.0-alpha)*Rho[NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]].Get_index()],U[0][NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]].Get_index()],U[1][NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]].Get_index()],NodeArrays->TypeOfNode,f[1],Rhob,U[0],U[1]);
		}
		for (int j=0;j<NodeIdxGlobalCorner.size();j++)
		{
			ApplyGlobalCorner(NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]],NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]].Get_AlphaDef()*NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]].Get_RhoDef(),NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]].Get_UDef(),NodeArrays->TypeOfNode,f[0],Rhor,U[0],U[1],Rhor[NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]].Get_index()]/Rho[NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]].Get_index()]);
			ApplyGlobalCorner(NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]],(1-NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]].Get_AlphaDef())*NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]].Get_RhoDef(),NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]].Get_UDef(),NodeArrays->TypeOfNode,f[1],Rhob,U[0],U[1],Rhob[NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]].Get_index()]/Rho[NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]].Get_index()]);
		}
	}
	else
	{
		for (int j=0;j<NodeIdx.size();j++)
		{
				ApplySymmetry(NodeArrays->NodeSymmetry[NodeIdx[j]].Get_BcNormal(),NodeArrays->NodeSymmetry[NodeIdx[j]].Get_connect(),NodeArrays->NodeSymmetry[NodeIdx[j]].Get_RhoDef(),NodeArrays->NodeSymmetry[NodeIdx[j]].Get_UDef(),f[0],Rhor,U[0],U[1]);
				ApplySymmetry(NodeArrays->NodeSymmetry[NodeIdx[j]].Get_BcNormal(),NodeArrays->NodeSymmetry[NodeIdx[j]].Get_connect(),NodeArrays->NodeSymmetry[NodeIdx[j]].Get_RhoDef(),NodeArrays->NodeSymmetry[NodeIdx[j]].Get_UDef(),f[1],Rhob,U[0],U[1]);
		}
		for (int j=0;j<NodeIdxSpecialWalls.size();j++)
		{
			ExtrapolationOnCornerConcave(Rho,NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]].Get_connect(),NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]].Get_BcNormal());
			ExtrapolationOnCornerConcave(Rhor,NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]].Get_connect(),NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]].Get_BcNormal());
			ExtrapolationOnCornerConcave(Rhob,NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]].Get_connect(),NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]].Get_BcNormal());
			alpha=(RhoN[NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]].Get_index()]+1.0)*0.5;
			ApplySymmetryWall(NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]],alpha*Rho[NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]].Get_index()],U[0][NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]].Get_index()],U[1][NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]].Get_index()],NodeArrays->TypeOfNode,f[0],Rhor,U[0],U[1]);
			ApplySymmetryWall(NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]],(1.0-alpha)*Rho[NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]].Get_index()],U[0][NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]].Get_index()],U[1][NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]].Get_index()],NodeArrays->TypeOfNode,f[1],Rhob,U[0],U[1]);
		}
		for (int j=0;j<NodeIdxGlobalCorner.size();j++)
		{
			ApplyGlobalCorner(NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]],NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]].Get_AlphaDef()*NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]].Get_RhoDef(),NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]].Get_UDef(),NodeArrays->TypeOfNode,f[0],Rhor,U[0],U[1],Rhor[NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]].Get_index()]/Rho[NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]].Get_index()]);
			ApplyGlobalCorner(NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]],(1-NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]].Get_AlphaDef())*NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]].Get_RhoDef(),NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]].Get_UDef(),NodeArrays->TypeOfNode,f[1],Rhob,U[0],U[1],Rhob[NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]].Get_index()]/Rho[NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]].Get_index()]);
		}
	}
}
void D2Q9ColourFluid::ApplyPatchPeriodic(PeriodicPatchBc& PerPatchBc){
	SetPeriodic(PerPatchBc.Get_PeriodicType());
	std::vector<int> NodeIdx=PerPatchBc.Get_NodeIndexByType();
	std::vector<int> NodeIdxSpecialWalls=PerPatchBc.Get_NodeIndexByTypeSpecialWalls();
	std::vector<int> NodeIdxGlobalCorner=PerPatchBc.Get_NodeIndexByTypeGlobalCorner();
	double alpha=0;
	if(PerPatchBc.Get_extrapolationAlpha())
	{
		for (int j=0;j<NodeIdx.size();j++)
		{
				ApplyPeriodic(NodeArrays->NodePeriodic[NodeIdx[j]].Get_BcNormal(),NodeArrays->NodePeriodic[NodeIdx[j]].Get_connect(),NodeArrays->NodePeriodic[NodeIdx[j]].Get_RhoDef(),NodeArrays->NodePeriodic[NodeIdx[j]].Get_UDef(),f[0],Rhor,U[0],U[1],Rhor[NodeArrays->NodePeriodic[NodeIdx[j]].Get_index()]/Rho[NodeArrays->NodePeriodic[NodeIdx[j]].Get_index()]);
				ApplyPeriodic(NodeArrays->NodePeriodic[NodeIdx[j]].Get_BcNormal(),NodeArrays->NodePeriodic[NodeIdx[j]].Get_connect(),NodeArrays->NodePeriodic[NodeIdx[j]].Get_RhoDef(),NodeArrays->NodePeriodic[NodeIdx[j]].Get_UDef(),f[1],Rhob,U[0],U[1],Rhor[NodeArrays->NodePeriodic[NodeIdx[j]].Get_index()]/Rho[NodeArrays->NodePeriodic[NodeIdx[j]].Get_index()]);
		}
		for (int j=0;j<NodeIdxSpecialWalls.size();j++)
		{
			ExtrapolationOnCornerConcave(Rho,NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]].Get_connect(),NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]].Get_BcNormal());
			ExtrapolationOnCornerConcave(Rhor,NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]].Get_connect(),NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]].Get_BcNormal());
			ExtrapolationOnCornerConcave(Rhob,NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]].Get_connect(),NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]].Get_BcNormal());
			alpha=(RhoN[NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]].Get_index()]+1.0)*0.5;
			ApplyPeriodicWall(NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]],alpha*Rho[NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]].Get_index()],U[0][NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]].Get_index()],U[1][NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]].Get_index()],NodeArrays->TypeOfNode,f[0],Rhor,U[0],U[1]);
			ApplyPeriodicWall(NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]],(1.0-alpha)*Rho[NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]].Get_index()],U[0][NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]].Get_index()],U[1][NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]].Get_index()],NodeArrays->TypeOfNode,f[1],Rhob,U[0],U[1]);
		}
		for (int j=0;j<NodeIdxGlobalCorner.size();j++)
		{
			ApplyGlobalCorner(NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]],NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]].Get_AlphaDef()*NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]].Get_RhoDef(),NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]].Get_UDef(),NodeArrays->TypeOfNode,f[0],Rhor,U[0],U[1],Rhor[NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]].Get_index()]/Rho[NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]].Get_index()]);
			ApplyGlobalCorner(NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]],(1-NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]].Get_AlphaDef())*NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]].Get_RhoDef(),NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]].Get_UDef(),NodeArrays->TypeOfNode,f[1],Rhob,U[0],U[1],Rhob[NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]].Get_index()]/Rho[NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]].Get_index()]);
		}
	}
	else
	{
		for (int j=0;j<NodeIdx.size();j++)
		{
				ApplyPeriodic(NodeArrays->NodePeriodic[NodeIdx[j]].Get_BcNormal(),NodeArrays->NodePeriodic[NodeIdx[j]].Get_connect(),NodeArrays->NodePeriodic[NodeIdx[j]].Get_RhoDef(),NodeArrays->NodePeriodic[NodeIdx[j]].Get_UDef(),f[0],Rhor,U[0],U[1],Rhor[NodeArrays->NodePeriodic[NodeIdx[j]].Get_index()]/Rho[NodeArrays->NodePeriodic[NodeIdx[j]].Get_index()]);
				ApplyPeriodic(NodeArrays->NodePeriodic[NodeIdx[j]].Get_BcNormal(),NodeArrays->NodePeriodic[NodeIdx[j]].Get_connect(),NodeArrays->NodePeriodic[NodeIdx[j]].Get_RhoDef(),NodeArrays->NodePeriodic[NodeIdx[j]].Get_UDef(),f[1],Rhob,U[0],U[1],Rhor[NodeArrays->NodePeriodic[NodeIdx[j]].Get_index()]/Rho[NodeArrays->NodePeriodic[NodeIdx[j]].Get_index()]);
		}
		for (int j=0;j<NodeIdxSpecialWalls.size();j++)
		{
			ExtrapolationOnCornerConcave(Rho,NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]].Get_connect(),NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]].Get_BcNormal());
			ExtrapolationOnCornerConcave(Rhor,NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]].Get_connect(),NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]].Get_BcNormal());
			ExtrapolationOnCornerConcave(Rhob,NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]].Get_connect(),NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]].Get_BcNormal());
			alpha=(RhoN[NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]].Get_index()]+1.0)*0.5;
			ApplyPeriodicWall(NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]],alpha*Rho[NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]].Get_index()],U[0][NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]].Get_index()],U[1][NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]].Get_index()],NodeArrays->TypeOfNode,f[0],Rhor,U[0],U[1]);
			ApplyPeriodicWall(NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]],(1.0-alpha)*Rho[NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]].Get_index()],U[0][NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]].Get_index()],U[1][NodeArrays->NodeSpecialWall[NodeIdxSpecialWalls[j]].Get_index()],NodeArrays->TypeOfNode,f[1],Rhob,U[0],U[1]);
		}
		for (int j=0;j<NodeIdxGlobalCorner.size();j++)
		{
			ApplyGlobalCorner(NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]],NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]].Get_AlphaDef()*NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]].Get_RhoDef(),NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]].Get_UDef(),NodeArrays->TypeOfNode,f[0],Rhor,U[0],U[1],Rhor[NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]].Get_index()]/Rho[NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]].Get_index()]);
			ApplyGlobalCorner(NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]],(1-NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]].Get_AlphaDef())*NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]].Get_RhoDef(),NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]].Get_UDef(),NodeArrays->TypeOfNode,f[1],Rhob,U[0],U[1],Rhob[NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]].Get_index()]/Rho[NodeArrays->NodeGlobalCorner[NodeIdxGlobalCorner[j]].Get_index()]);
		}
	}
}
/*
void D2Q9ColourFluid::ApplyPatchGlobalCorner(GlobalCornerPatchBc& PerPatchBc){
	double fi_tmp=0;
	for (int j=0;j<NodeArrays->NodeGlobalCorner.size();j++)
	{
		ApplyGlobalCorner(NodeArrays->NodeGlobalCorner[j],NodeArrays->NodeGlobalCorner[j].Get_AlphaDef()*NodeArrays->NodeGlobalCorner[j].Get_RhoDef(),NodeArrays->NodeGlobalCorner[j].Get_UDef(),NodeArrays->TypeOfNode,f[0],Rhor,U[0],U[1],Rhor[NodeArrays->NodeGlobalCorner[j].Get_index()]/Rho[NodeArrays->NodeGlobalCorner[j].Get_index()]);
		ApplyGlobalCorner(NodeArrays->NodeGlobalCorner[j],(1-NodeArrays->NodeGlobalCorner[j].Get_AlphaDef())*NodeArrays->NodeGlobalCorner[j].Get_RhoDef(),NodeArrays->NodeGlobalCorner[j].Get_UDef(),NodeArrays->TypeOfNode,f[1],Rhob,U[0],U[1],Rhob[NodeArrays->NodeGlobalCorner[j].Get_index()]/Rho[NodeArrays->NodeGlobalCorner[j].Get_index()]);
	//Keep alpha
		Rhor[NodeArrays->NodeGlobalCorner[j].Get_index()]=0.0;
		Rhob[NodeArrays->NodeGlobalCorner[j].Get_index()]=0.0;
		for(int i=0;i<9;i++)
		{
			fi_tmp=f[0]->f[i][NodeArrays->NodeGlobalCorner[j].Get_index()]+f[1]->f[i][NodeArrays->NodeGlobalCorner[j].Get_index()];
			f[0]->f[i][NodeArrays->NodeGlobalCorner[j].Get_index()]=NodeArrays->NodeGlobalCorner[j].Get_AlphaDef()*fi_tmp;
			f[1]->f[i][NodeArrays->NodeGlobalCorner[j].Get_index()]=(1-NodeArrays->NodeGlobalCorner[j].Get_AlphaDef())*fi_tmp;
			Rhor[NodeArrays->NodeGlobalCorner[j].Get_index()]+=f[0]->f[i][NodeArrays->NodeGlobalCorner[j].Get_index()];
			Rhob[NodeArrays->NodeGlobalCorner[j].Get_index()]+=f[1]->f[i][NodeArrays->NodeGlobalCorner[j].Get_index()];
		}
	}
}*/

void  D2Q9ColourFluid::ImposeAlpha(int &index, double alpha)
{
	Rhor[index]=0.0;
	Rhob[index]=0.0;
	for(int i=0;i<9;i++)
	{
		f[0]->f[i][index]=alpha*f[0]->f[i][index];//fi_tmp;
		f[1]->f[i][index]=(1-alpha)*f[1]->f[i][index];//fi_tmp;
		Rhor[index]+=f[0]->f[i][index];
		Rhob[index]+=f[1]->f[i][index];
	}
}
///Select and apply boundary conditions
void D2Q9ColourFluid::ApplyBc(){
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
	double fi_tmp;
	double alpha=0;
/*
	for (int j=0;j<NodeArrays->NodeVelocity.size();j++)
	{
		ApplyVelocity(NodeArrays->NodeVelocity[j].Get_BcNormal(),NodeArrays->NodeVelocity[j].Get_connect(),NodeArrays->NodeVelocity[j].Get_UDef(), f[0],Rhor,U[0],U[1]);
		ApplyVelocity(NodeArrays->NodeVelocity[j].Get_BcNormal(),NodeArrays->NodeVelocity[j].Get_connect(),NodeArrays->NodeVelocity[j].Get_UDef(), f[1],Rhob,U[0],U[1]);

//Impose alpha
		Rhor[NodeArrays->NodeVelocity[j].Get_index()]=0.0;
		Rhob[NodeArrays->NodeVelocity[j].Get_index()]=0.0;
		for(int i=0;i<9;i++)
		{
			//fi_tmp=f[0]->f[i][NodeArrays->NodeVelocity[j].Get_index()];//+f[1]->f[i][NodeArrays->NodeVelocity[j].Get_index()];
			f[0]->f[i][NodeArrays->NodeVelocity[j].Get_index()]=NodeArrays->NodeVelocity[j].Get_AlphaDef()*f[0]->f[i][NodeArrays->NodeVelocity[j].Get_index()];//fi_tmp;
			f[1]->f[i][NodeArrays->NodeVelocity[j].Get_index()]=(1-NodeArrays->NodeVelocity[j].Get_AlphaDef())*f[1]->f[i][NodeArrays->NodeVelocity[j].Get_index()];//fi_tmp;
			Rhor[NodeArrays->NodeVelocity[j].Get_index()]+=f[0]->f[i][NodeArrays->NodeVelocity[j].Get_index()];
			Rhob[NodeArrays->NodeVelocity[j].Get_index()]+=f[1]->f[i][NodeArrays->NodeVelocity[j].Get_index()];

		}

	}


	for (int j=0;j<NodeArrays->NodePressure.size();j++)
	{
//		ApplyPressure(NodeArrays->NodePressure[j].Get_BcNormal(),NodeArrays->NodePressure[j].Get_connect(),NodeArrays->NodePressure[j].Get_AlphaDef()*NodeArrays->NodePressure[j].Get_RhoDef(), f[0],Rhor,U[0],U[1]);
//		ApplyPressure(NodeArrays->NodePressure[j].Get_BcNormal(),NodeArrays->NodePressure[j].Get_connect(),(1-NodeArrays->NodePressure[j].Get_AlphaDef())*NodeArrays->NodePressure[j].Get_RhoDef(), f[1],Rhob,U[0],U[1]);
		ExtrapolationOnWall(RhoN,NodeArrays->NodePressure[j].Get_connect(),NodeArrays->NodePressure[j].Get_BcNormal());
		alpha=(RhoN[NodeArrays->NodePressure[j].Get_index()]+1.0)*0.5;
		ApplyPressure(NodeArrays->NodePressure[j].Get_BcNormal(),NodeArrays->NodePressure[j].Get_connect(),alpha*NodeArrays->NodePressure[j].Get_RhoDef(), f[0],Rhor,U[0],U[1]);
		ApplyPressure(NodeArrays->NodePressure[j].Get_BcNormal(),NodeArrays->NodePressure[j].Get_connect(),(1-alpha)*NodeArrays->NodePressure[j].Get_RhoDef(), f[1],Rhob,U[0],U[1]);

	}
*/
	for (int j=0;j<NodeArrays->NodeWall.size();j++)
	{
		ApplyWall(NodeArrays->NodeWall[j].Get_BcNormal(),NodeArrays->NodeWall[j].Get_connect(),f[0],Rhor,U[0],U[1]);
		ApplyWall(NodeArrays->NodeWall[j].Get_BcNormal(),NodeArrays->NodeWall[j].Get_connect(),f[1],Rhob,U[0],U[1]);

	}

/*
	for (int j=0;j<NodeArrays->NodeSpecialWall.size();j++)
	{

		ExtrapolationOnCornerConcave(Rho,NodeArrays->NodeSpecialWall[j].Get_connect(),NodeArrays->NodeSpecialWall[j].Get_BcNormal());
		ExtrapolationOnCornerConcave(Rhor,NodeArrays->NodeSpecialWall[j].Get_connect(),NodeArrays->NodeSpecialWall[j].Get_BcNormal());
		ExtrapolationOnCornerConcave(Rhob,NodeArrays->NodeSpecialWall[j].Get_connect(),NodeArrays->NodeSpecialWall[j].Get_BcNormal());
		alpha=(RhoN[NodeArrays->NodeSpecialWall[j].Get_index()]+1.0)*0.5;

		ApplySpecialWall(NodeArrays->NodeSpecialWall[j],alpha*Rho[NodeArrays->NodeSpecialWall[j].Get_index()],U[0][NodeArrays->NodeSpecialWall[j].Get_index()],U[1][NodeArrays->NodeSpecialWall[j].Get_index()],NodeArrays->TypeOfNode,f[0],Rhor,U[0],U[1]);
		ApplySpecialWall(NodeArrays->NodeSpecialWall[j],(1.0-alpha)*Rho[NodeArrays->NodeSpecialWall[j].Get_index()],U[0][NodeArrays->NodeSpecialWall[j].Get_index()],U[1][NodeArrays->NodeSpecialWall[j].Get_index()],NodeArrays->TypeOfNode,f[1],Rhob,U[0],U[1]);
	}
*/
//	std::cout<<"Apply bc "<<"f6: "<<f[1]->f[6][idxcheck]<<" f8: "<<f[1]->f[8][idxcheck];

	for (int j=0;j<NodeArrays->NodeCorner.size();j++)
	{
		ApplyCornerWall(NodeArrays->NodeCorner[j], f[0],Rhor,U[0],U[1]);
		ApplyCornerWall(NodeArrays->NodeCorner[j], f[1],Rhob,U[0],U[1]);
	}
/*
	for (int j=0;j<NodeArrays->NodeSymmetry.size();j++)
	{
			ApplySymmetry(NodeArrays->NodeSymmetry[j].Get_BcNormal(),NodeArrays->NodeSymmetry[j].Get_connect(),NodeArrays->NodeSymmetry[j].Get_RhoDef(),NodeArrays->NodeSymmetry[j].Get_UDef(),f[0],Rhor,U[0],U[1]);
			ApplySymmetry(NodeArrays->NodeSymmetry[j].Get_BcNormal(),NodeArrays->NodeSymmetry[j].Get_connect(),NodeArrays->NodeSymmetry[j].Get_RhoDef(),NodeArrays->NodeSymmetry[j].Get_UDef(),f[1],Rhob,U[0],U[1]);
	}
	for (int j=0;j<NodeArrays->NodePeriodic.size();j++)
	{
			ApplyPeriodic(NodeArrays->NodePeriodic[j].Get_BcNormal(),NodeArrays->NodePeriodic[j].Get_connect(),NodeArrays->NodePeriodic[j].Get_RhoDef(),NodeArrays->NodePeriodic[j].Get_UDef(),f[0],Rhor,U[0],U[1],Rhor[NodeArrays->NodePeriodic[j].Get_index()]/Rho[NodeArrays->NodePeriodic[j].Get_index()]);
			ApplyPeriodic(NodeArrays->NodePeriodic[j].Get_BcNormal(),NodeArrays->NodePeriodic[j].Get_connect(),NodeArrays->NodePeriodic[j].Get_RhoDef(),NodeArrays->NodePeriodic[j].Get_UDef(),f[1],Rhob,U[0],U[1],Rhor[NodeArrays->NodePeriodic[j].Get_index()]/Rho[NodeArrays->NodePeriodic[j].Get_index()]);
	}
*/
	/*
	for (int j=0;j<NodeArrays->NodeGlobalCorner.size();j++)
	{
		ApplyGlobalCorner(NodeArrays->NodeGlobalCorner[j],NodeArrays->NodeGlobalCorner[j].Get_AlphaDef()*NodeArrays->NodeGlobalCorner[j].Get_RhoDef(),NodeArrays->NodeGlobalCorner[j].Get_UDef(),NodeArrays->TypeOfNode,f[0],Rhor,U[0],U[1],Rhor[NodeArrays->NodeGlobalCorner[j].Get_index()]/Rho[NodeArrays->NodeGlobalCorner[j].Get_index()]);
		ApplyGlobalCorner(NodeArrays->NodeGlobalCorner[j],(1-NodeArrays->NodeGlobalCorner[j].Get_AlphaDef())*NodeArrays->NodeGlobalCorner[j].Get_RhoDef(),NodeArrays->NodeGlobalCorner[j].Get_UDef(),NodeArrays->TypeOfNode,f[1],Rhob,U[0],U[1],Rhob[NodeArrays->NodeGlobalCorner[j].Get_index()]/Rho[NodeArrays->NodeGlobalCorner[j].Get_index()]);
//Keep alpha
		Rhor[NodeArrays->NodeGlobalCorner[j].Get_index()]=0.0;
		Rhob[NodeArrays->NodeGlobalCorner[j].Get_index()]=0.0;
		for(int i=0;i<9;i++)
		{
			fi_tmp=f[0]->f[i][NodeArrays->NodeGlobalCorner[j].Get_index()]+f[1]->f[i][NodeArrays->NodeGlobalCorner[j].Get_index()];
			f[0]->f[i][NodeArrays->NodeGlobalCorner[j].Get_index()]=NodeArrays->NodeGlobalCorner[j].Get_AlphaDef()*fi_tmp;
			f[1]->f[i][NodeArrays->NodeGlobalCorner[j].Get_index()]=(1-NodeArrays->NodeGlobalCorner[j].Get_AlphaDef())*fi_tmp;
			Rhor[NodeArrays->NodeGlobalCorner[j].Get_index()]+=f[0]->f[i][NodeArrays->NodeGlobalCorner[j].Get_index()];
			Rhob[NodeArrays->NodeGlobalCorner[j].Get_index()]+=f[1]->f[i][NodeArrays->NodeGlobalCorner[j].Get_index()];
		}
	}
*/
	parallel->barrier();


}
void D2Q9ColourFluid::ExtrapolDensityInSolid(){
	//Extrapolation of density or normal density at the first layer in the solid
	for (int j=0;j<NodeArrays->CornerConcave.size();j++)
	{
		ExtrapolDensity.ExtrapolationCornerConcaveToSolid(RhoN,NodeArrays->NodeCorner[NodeArrays->CornerConcave[j]].Get_connect(),NodeArrays->NodeCorner[NodeArrays->CornerConcave[j]].Get_BcNormal());
	}
	for (int j=0;j<NodeArrays->NodeWall.size();j++)
	{
		ExtrapolDensity.ExtrapolationWallToSolid(RhoN,NodeArrays->NodeWall[j].Get_connect(),NodeArrays->NodeWall[j].Get_BcNormal());
	}
	for (int j=0;j<NodeArrays->CornerConvex.size();j++)
	{
		ExtrapolDensity.ExtrapolationCornerConvexToSolid(RhoN,NodeArrays->NodeCorner[NodeArrays->CornerConvex[j]].Get_connect(),NodeArrays->NodeCorner[NodeArrays->CornerConvex[j]].Get_BcNormal());
	}
	for (int j=0;j<NodeArrays->NodeSpecialWall.size();j++)
	{
		ExtrapolDensity.ExtrapolationCornerConcaveToSolid(RhoN,NodeArrays->NodeSpecialWall[j].Get_connect(),NodeArrays->NodeSpecialWall[j].Get_BcNormal());
	}
/*	for (int j=0;j<NodeArrays->NodeVelocity.size();j++)
	{
		ExtrapolDensity.ExtrapolationCornerConcaveToSolid(RhoN,NodeArrays->NodeVelocity[j].Get_connect(),NodeArrays->NodeVelocity[j].Get_BcNormal());
	}*/
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
		for (int j=0;j<NodeArrays->NodeSpecialWall.size();j++)
		{
			(this->*PtrMacro)(NodeArrays->NodeSpecialWall[j].Get_index());
		}
		for (int j=0;j<NodeArrays->NodeSymmetry.size();j++)
		{
			(this->*PtrMacro)(NodeArrays->NodeSymmetry[j].Get_index());
		}
		for (int j=0;j<NodeArrays->NodePeriodic.size();j++)
		{
			(this->*PtrMacro)(NodeArrays->NodePeriodic[j].Get_index());
		}
	//	std::cout<<"Apply bc "<<"f6: "<<f[0]->f[6][idxcheck]<<" f8: "<<f[0]->f[8][idxcheck];
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
			U[j][idx]+=(f[0]->f[k][idx] + f[1]->f[k][idx])*Ei[k][j];
		}
	}
	Rho[idx]=Rhor[idx]+Rhob[idx];
	RhoN[idx]=(Rhor[idx]-Rhob[idx])/Rho[idx];
	U[0][idx]=U[0][idx]/Rho[idx];
	U[1][idx]=U[1][idx]/Rho[idx];
//	U[0][idx]=(U[0][idx]+0.5*F[0][idx])/Rho[idx];
//	U[1][idx]=(U[1][idx]+0.5*F[1][idx])/Rho[idx];
}
void D2Q9ColourFluid::UpdatePressure(){
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
void D2Q9ColourFluid::CalculatePressure(int const &idx){
	P[idx]=IdealGazIsothermalPressure(Rho[idx]);
}
double D2Q9ColourFluid::IdealGazIsothermalPressure(double const &Rho){
	return Rho/3.0;
}

void D2Q9ColourFluid::NormalDensityExtrapolationSpacial2ndOrder(int const & idxNodeArray, int & nodenumber, int* connect,int & normal){
	switch(normal)
	{
	case 1:
		RhoN[connect[3]]=RhoN[connect[0]]+std::tan(pi/2.0-teta[0])*std::abs(ExtrapolationSpacial2ndOrder(G[0],connect[0],connect[1]));
		testVar[connect[0]]=RhoN[connect[3]];
		if(RhoN[connect[3]]<-1.0) RhoN[connect[3]]=-1;
		if(RhoN[connect[3]]>1.0) RhoN[connect[3]]=1;
		break;
	case 2:
		RhoN[connect[4]]=RhoN[connect[0]]+std::tan(pi/2.0-teta[0])*std::abs(ExtrapolationSpacial2ndOrder(G[0],connect[0],connect[2]));
		testVar[connect[0]]=RhoN[connect[0]];
		testVar[connect[2]]=RhoN[connect[4]];
		if(RhoN[connect[4]]<-1.0) RhoN[connect[4]]=-1;
		if(RhoN[connect[4]]>1.0) RhoN[connect[4]]=1;
		break;
	case 3:
		RhoN[connect[1]]=RhoN[connect[0]]+std::tan(pi/2.0-teta[0])*std::abs(ExtrapolationSpacial2ndOrder(G[0],connect[0],connect[3]));
		testVar[connect[0]]=RhoN[connect[1]];
		if(RhoN[connect[1]]<-1.0) RhoN[connect[1]]=-1;
		if(RhoN[connect[1]]>1.0) RhoN[connect[1]]=1;
		break;
	case 4:
		RhoN[connect[2]]=RhoN[connect[0]]+std::tan(pi/2.0-teta[0])*std::abs(ExtrapolationSpacial2ndOrder(G[0],connect[0],connect[4]));
		testVar[connect[0]]=RhoN[connect[2]];
		if(RhoN[connect[2]]<-1.0) RhoN[connect[2]]=-1;
		if(RhoN[connect[2]]>1.0) RhoN[connect[2]]=1;
		break;
	case 5:
		RhoN[connect[7]]=RhoN[connect[0]]+std::tan(pi/2.0-teta[0])*std::abs(ExtrapolationSpacial2ndOrder(G[0],connect[0],connect[5]));
		if(RhoN[connect[7]]<-1.0) RhoN[connect[7]]=-1;
		if(RhoN[connect[7]]>1.0) RhoN[connect[7]]=1;
		break;
	case 6:
		RhoN[connect[8]]=RhoN[connect[0]]+std::tan(pi/2.0-teta[0])*std::abs(ExtrapolationSpacial2ndOrder(G[0],connect[0],connect[6]));
		if(RhoN[connect[8]]<-1.0) RhoN[connect[8]]=-1;
		if(RhoN[connect[8]]>1.0) RhoN[connect[8]]=1;
		break;
	case 7:
		RhoN[connect[5]]=RhoN[connect[0]]+std::tan(pi/2.0-teta[0])*std::abs(ExtrapolationSpacial2ndOrder(G[0],connect[0],connect[7]));
		if(RhoN[connect[5]]<-1.0) RhoN[connect[5]]=-1;
		if(RhoN[connect[5]]>1.0) RhoN[connect[5]]=1;
		break;
	case 8:
		RhoN[connect[6]]=RhoN[connect[0]]+std::tan(pi/2.0-teta[0])*std::abs(ExtrapolationSpacial2ndOrder(G[0],connect[0],connect[8]));
		if(RhoN[connect[6]]<-1.0) RhoN[connect[6]]=-1;
		if(RhoN[connect[6]]>1.0) RhoN[connect[6]]=1;
		break;
	}

}
void D2Q9ColourFluid::NormalDensityExtrapolationWeight(int const & idxNodeArray, int & nodenumber, int* connect,int & normal){
	double InvSqrt2=std::sqrt(0.5);
	double InvSqrt2_5=std::sqrt(0.4);//0.4=2/5
	double InvSumWeightWall=1.0/(1.0+2.0*InvSqrt2);//std::sqrt(2.0)/(2+std::sqrt(2.0));
	double InvSumWeightCornerConvex=1.0/(2.0+1.0*InvSqrt2);
	double InvSumWeightCornerConcave=1.0/(1.0+2.0*InvSqrt2_5);//0.4=2/5
	double InvSumWeightCornerWall=1.0/(1.0+1.0*InvSqrt2);

	switch(normal)
	{
	case 1:
		RhoN[connect[3]]=InvSumWeightWall*(InvSqrt2*(RhoN[connect[2]]+RhoN[connect[4]])+RhoN[connect[0]]);
		testVar[connect[0]]=RhoN[connect[3]];
		break;
	case 2:
		RhoN[connect[4]]=InvSumWeightWall*(InvSqrt2*(RhoN[connect[1]]+RhoN[connect[3]])+RhoN[connect[0]]);
		testVar[connect[0]]=RhoN[connect[4]];
		break;
	case 3:
		RhoN[connect[1]]=InvSumWeightWall*(InvSqrt2*(RhoN[connect[2]]+RhoN[connect[4]])+RhoN[connect[0]]);
		testVar[connect[0]]=RhoN[connect[1]];
		break;
	case 4:
		RhoN[connect[2]]=InvSumWeightWall*(InvSqrt2*(RhoN[connect[1]]+RhoN[connect[3]])+RhoN[connect[0]]);
		testVar[connect[0]]=RhoN[connect[2]];
		break;
	case 5:
		if(NodeArrays->NodeCorner[idxNodeArray].Get_CornerType()==Convex)
			RhoN[connect[7]]=InvSumWeightCornerConvex*(RhoN[connect[3]]+RhoN[connect[4]]+InvSqrt2*RhoN[connect[0]]);
		else
		{
			RhoN[connect[3]]=InvSumWeightCornerWall*(InvSqrt2*RhoN[connect[2]]+RhoN[connect[0]]);
			RhoN[connect[4]]=InvSumWeightCornerWall*(InvSqrt2*RhoN[connect[1]]+RhoN[connect[0]]);
			RhoN[connect[7]]=InvSumWeightCornerConcave*(InvSqrt2_5*(RhoN[connect[3]]+RhoN[connect[4]])+RhoN[connect[0]]);
		}
		break;
	case 6:
		if(NodeArrays->NodeCorner[idxNodeArray].Get_CornerType()==Convex)
			RhoN[connect[8]]=InvSumWeightCornerConvex*(RhoN[connect[1]]+RhoN[connect[4]]+InvSqrt2*RhoN[connect[0]]);
		else
		{
			RhoN[connect[1]]=InvSumWeightCornerWall*(InvSqrt2*RhoN[connect[2]]+RhoN[connect[0]]);
			RhoN[connect[4]]=InvSumWeightCornerWall*(InvSqrt2*RhoN[connect[3]]+RhoN[connect[0]]);
			RhoN[connect[8]]=InvSumWeightCornerConcave*(InvSqrt2_5*(RhoN[connect[3]]+RhoN[connect[2]])+RhoN[connect[0]]);
		}
		break;
	case 7:
		if(NodeArrays->NodeCorner[idxNodeArray].Get_CornerType()==Convex)
			RhoN[connect[5]]=InvSumWeightCornerConvex*(RhoN[connect[1]]+RhoN[connect[2]]+InvSqrt2*RhoN[connect[0]]);
		else
		{
			RhoN[connect[1]]=InvSumWeightCornerWall*(InvSqrt2*RhoN[connect[4]]+RhoN[connect[0]]);
			RhoN[connect[2]]=InvSumWeightCornerWall*(InvSqrt2*RhoN[connect[3]]+RhoN[connect[0]]);
			RhoN[connect[5]]=InvSumWeightCornerConcave*(InvSqrt2_5*(RhoN[connect[3]]+RhoN[connect[4]])+RhoN[connect[0]]);
		}
		break;
	case 8:
		if(NodeArrays->NodeCorner[idxNodeArray].Get_CornerType()==Convex)
			RhoN[connect[6]]=InvSumWeightCornerConvex*(RhoN[connect[3]]+RhoN[connect[2]]+InvSqrt2*RhoN[connect[0]]);
		else
		{
			RhoN[connect[3]]=InvSumWeightCornerWall*(InvSqrt2*RhoN[connect[4]]+RhoN[connect[0]]);
			RhoN[connect[2]]=InvSumWeightCornerWall*(InvSqrt2*RhoN[connect[1]]+RhoN[connect[0]]);
			RhoN[connect[6]]=InvSumWeightCornerConcave*(InvSqrt2_5*(RhoN[connect[1]]+RhoN[connect[4]])+RhoN[connect[0]]);
		}
		break;
	}
}
void  D2Q9ColourFluid::NormalDensityNoTeta(int const & idxNodeArray, int & nodenumber, int* connect,int & normal){
	switch(normal)
	{
	testVar[connect[0]]=RhoN[connect[0]];
	case 1:
		RhoN[connect[3]]=RhoN[connect[0]];
		testVar[connect[1]]=RhoN[connect[3]];
		break;
	case 2:
		RhoN[connect[4]]=RhoN[connect[0]];
		testVar[connect[2]]=RhoN[connect[4]];
		break;
	case 3:
		RhoN[connect[1]]=RhoN[connect[0]];
		testVar[connect[3]]=RhoN[connect[1]];
		break;
	case 4:
		RhoN[connect[2]]=RhoN[connect[0]];
		testVar[connect[4]]=RhoN[connect[2]];
		break;
	case 5:
		RhoN[connect[7]]=RhoN[connect[0]];
		testVar[connect[5]]=RhoN[connect[7]];
		break;
	case 6:
		RhoN[connect[8]]=RhoN[connect[0]];
		testVar[connect[6]]=RhoN[connect[8]];
		break;
	case 7:
		RhoN[connect[5]]=RhoN[connect[0]];
		testVar[connect[7]]=RhoN[connect[5]];
		break;
	case 8:
		RhoN[connect[6]]=RhoN[connect[0]];
		testVar[connect[8]]=RhoN[connect[6]];
		break;
	}
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
void D2Q9ColourFluid::Synchronise_Colour_gradient(){
	//Put variables in buffer arrays
		for (unsigned int i=0;i<IdRNodeW.size();i++)
		{
				buf_MacroSend[0][1][i]=Normal[0][IdRNodeW[i]];
				buf_MacroSend[1][1][i]=Normal[1][IdRNodeW[i]];
		}
		for (unsigned int i=0;i<IdRNodeS.size();i++)
		{
				buf_MacroSend[0][2][i]=Normal[0][IdRNodeS[i]];
				buf_MacroSend[1][2][i]=Normal[1][IdRNodeS[i]];
		}
		for (unsigned int i=0;i<IdRNodeE.size();i++)
		{
				buf_MacroSend[0][0][i]=Normal[0][IdRNodeE[i]];
				buf_MacroSend[1][0][i]=Normal[1][IdRNodeE[i]];
		}
		for (unsigned int i=0;i<IdRNodeN.size();i++)
		{
				buf_MacroSend[0][3][i]=Normal[0][IdRNodeN[i]];
				buf_MacroSend[1][3][i]=Normal[1][IdRNodeN[i]];
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
				MultiBlock_->Send(&buf_MacroSend[0][0][0],IdRNodeE.size(),1,tag_x_r);
				MultiBlock_->Send(&buf_MacroSend[1][0][0],IdRNodeE.size(),1,tag_x_r);
		}
		if(IdGNodeW.size()>=1)
		{
				MultiBlock_->Recv(&buf_MacroRecv[0][0][0],IdGNodeW.size(),3,tag_x_r,status);
				MultiBlock_->Recv(&buf_MacroRecv[1][0][0],IdGNodeW.size(),3,tag_x_r,status);
		}

		if(IdRNodeW.size()>=1)
		{
				MultiBlock_->Send(&buf_MacroSend[0][1][0],IdRNodeW.size(),3,tag_x_l);
				MultiBlock_->Send(&buf_MacroSend[1][1][0],IdRNodeW.size(),3,tag_x_l);
		}
		if(IdGNodeE.size()>=1)
		{
				MultiBlock_->Recv(&buf_MacroRecv[0][1][0],IdGNodeE.size(),1,tag_x_l,status);
				MultiBlock_->Recv(&buf_MacroRecv[1][1][0],IdGNodeE.size(),1,tag_x_l,status);
		}
			if(IdRNodeN.size()>=1)
		{
				MultiBlock_->Send(&buf_MacroSend[0][3][0],IdRNodeN.size(),0,tag_y_t);
				MultiBlock_->Send(&buf_MacroSend[1][3][0],IdRNodeN.size(),0,tag_y_t);
		}
		if(IdGNodeS.size()>=1)
		{
				MultiBlock_->Recv(&buf_MacroRecv[0][3][0],IdGNodeS.size(),2,tag_y_t,status);
				MultiBlock_->Recv(&buf_MacroRecv[1][3][0],IdGNodeS.size(),2,tag_y_t,status);
		}
		if(IdRNodeS.size()>=1)
		{
				MultiBlock_->Send(&buf_MacroSend[0][2][0],IdRNodeS.size(),2,tag_y_b);
				MultiBlock_->Send(&buf_MacroSend[1][2][0],IdRNodeS.size(),2,tag_y_b);
		}
		if(IdGNodeN.size()>=1)
		{
				MultiBlock_->Recv(&buf_MacroRecv[0][2][0],IdGNodeN.size(),0,tag_y_b,status);
				MultiBlock_->Recv(&buf_MacroRecv[1][2][0],IdGNodeN.size(),0,tag_y_b,status);
		}
	//Set variables from buffer to real variables
		for (unsigned int i=0;i<IdGNodeE.size();i++)
		{
			Normal[0][IdGNodeE[i]]=buf_MacroRecv[0][1][i];
			Normal[1][IdGNodeE[i]]=buf_MacroRecv[1][1][i];
		}
		for (unsigned int i=0;i<IdGNodeN.size();i++)
		{
			Normal[0][IdGNodeN[i]]=buf_MacroRecv[0][2][i];
			Normal[1][IdGNodeN[i]]=buf_MacroRecv[1][2][i];
		}
		for (unsigned int i=0;i<IdGNodeW.size();i++)
		{
			Normal[0][IdGNodeW[i]]=buf_MacroRecv[0][0][i];
			Normal[1][IdGNodeW[i]]=buf_MacroRecv[1][0][i];
		}

		for (unsigned int i=0;i<IdGNodeS.size();i++)
		{
			Normal[0][IdGNodeS[i]]=buf_MacroRecv[0][3][i];
			Normal[1][IdGNodeS[i]]=buf_MacroRecv[1][3][i];
		}

	//Corners of the local domain
		if(IdRNodeSE.size()>=1)
		{
				MultiBlock_->Send(&Normal[0][IdRNodeSE[0]],1,5,tag_x_l);
				MultiBlock_->Send(&Normal[1][IdRNodeSE[0]],1,5,tag_x_l);
		}
		if(IdGNodeNW.size()>=1)
		{
				MultiBlock_->Recv(&Normal[0][IdGNodeNW[0]],1,6,tag_x_l,status);
				MultiBlock_->Recv(&Normal[1][IdGNodeNW[0]],1,6,tag_x_l,status);
		}
		if(IdRNodeSW.size()>=1)
		{
				MultiBlock_->Send(&Normal[0][IdRNodeSW[0]],1,7,tag_x_l);
				MultiBlock_->Send(&Normal[1][IdRNodeSW[0]],1,7,tag_x_l);
		}
		if(IdGNodeNE.size()>=1)
		{
				MultiBlock_->Recv(&Normal[0][IdGNodeNE[0]],1,4,tag_x_l,status);
				MultiBlock_->Recv(&Normal[1][IdGNodeNE[0]],1,4,tag_x_l,status);
		}
		if(IdRNodeNE.size()>=1)
		{
				MultiBlock_->Send(&Normal[0][IdRNodeNE[0]],1,4,tag_x_l);
				MultiBlock_->Send(&Normal[1][IdRNodeNE[0]],1,4,tag_x_l);
		}
		if(IdGNodeSW.size()>=1)
		{
				MultiBlock_->Recv(&Normal[0][IdGNodeSW[0]],1,7,tag_x_l,status);
				MultiBlock_->Recv(&Normal[1][IdGNodeSW[0]],1,7,tag_x_l,status);
		}
		if(IdRNodeNW.size()>=1)
		{
				MultiBlock_->Send(&Normal[0][IdRNodeNW[0]],1,6,tag_d_bl);
				MultiBlock_->Send(&Normal[1][IdRNodeNW[0]],1,6,tag_d_bl);
		}
		if(IdGNodeSE.size()>=1)
		{
				MultiBlock_->Recv(&Normal[0][IdGNodeSE[0]],1,5,tag_d_bl,status);
				MultiBlock_->Recv(&Normal[1][IdGNodeSE[0]],1,5,tag_d_bl,status);
		}
}

/*
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
*/
