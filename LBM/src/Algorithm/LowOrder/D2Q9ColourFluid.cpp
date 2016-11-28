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
	PtrCalNormal=0;
	PtrCalNormalWall=0;
	PtrCollisionWall=0;
	PtrExtrapolDensity=0;
	Curv=0;
	D_tmp=0;
	PtrD_tmp=0;
	PtrRecolour=0;
	PtrMacro=0;
	beta=0.7;
	Rho_limiter=1.e-10;
	A1=0;
	A2=0;
	F=0;
	G=0;
	Normal=0;
	G_Norm=0;
	tension=0;
	teta=0;


}

D2Q9ColourFluid::~D2Q9ColourFluid() {
	delete [] Rhor;
	delete [] Rhob;
}

D2Q9ColourFluid::D2Q9ColourFluid(MultiBlock* MultiBlock__,ParallelManager* parallel__,WriterManager* Writer__, Parameters* Parameters__,InitLBM& ini){
	Rho_limiter=1.e-10;
	//Initialise the common variables for the two phase models.
	InitD2Q9TwoPhases(MultiBlock__,parallel__,Writer__, Parameters__, ini);
//Set Pointers On Functions for selecting the right model dynamically
	Set_PointersOnFunctions();
//Initialise variables of the colour fluid model.
	InitColourFluid(ini);
//Initialise boundary conditions.
	InitD2Q9Bc(Dic, Parameters__);
//Initialise communication between processors
	IniComVariables();
//Set the variables names and the variable pointers for output in solution
	Solution2D::Set_output();
//Set the variables names and the variable pointers for breakpoints in solution
	Solution2D::Set_breakpoint();
//Set_Convergence
	Set_Convergence();
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
// Select the method for calculating the contact angle
	Set_NormalWall();
}

void D2Q9ColourFluid::Select_Colour_Operator(ColourOperatorType OperatorType_){
	switch(OperatorType_)
	{
	case ::SurfaceForce:
		PtrCollision=&D2Q9ColourFluid::Collision_SurfaceForce;
		PtrCollisionWall=&D2Q9ColourFluid::Collision_SurfaceForce;
		break;
	case ::Grunau:
		PtrCollision=&D2Q9ColourFluid::Collision_Grunau;
		PtrCollisionWall=PtrCollision;
		break;
	case ::Reis:
		PtrCollision=&D2Q9ColourFluid::Collision_Reis;
		PtrCollisionWall=PtrCollision;
		break;
	default:
		std::cerr<<" Colour operator not found."<<std::endl;
		break;
	}
}

void D2Q9ColourFluid::Set_NormalWall(){

	switch(PtrParameters->Get_ContactAngleType())
	{
	case ::NoTeta:
		PtrCalNormal=&D2Q9ColourFluid::CalculNormal;
		PtrCalNormalWall=&D2Q9ColourFluid::CalculNormal_NoTeta;
		PtrColourGradWall=PtrColourGrad;
		PtrExtrapolDensity=&D2Q9ColourFluid::NormalDensityNoTeta;
		break;
	case ::FixTeta:
		PtrCalNormalWall=&D2Q9ColourFluid::CalculNormal_FixTeta;
		PtrCalNormal=&D2Q9ColourFluid::CalculNormal;
		PtrColourGradWall=PtrColourGrad;
		teta=new double[1];
		PtrExtrapolDensity=&D2Q9ColourFluid::NormalDensityExtrapolationWeight;
		//PtrExtrapolDensity=&D2Q9ColourFluid::NormalDensityExtrapolationSpacial2ndOrder;//&D2Q9ColourFluid::NormalDensityExtrapolationWeight;
		break;
	case ::NonCstTeta:
		PtrCalNormalWall=&D2Q9ColourFluid::CalculNormal_NoCstTeta;
		PtrCalNormal=&D2Q9ColourFluid::CalculNormal;
		Dic->AddVar(Scalar,"Teta",true,true,false,teta);
		PtrExtrapolDensity=&D2Q9ColourFluid::NormalDensityExtrapolationSpacial2ndOrder;
		break;
	default:
		std::cerr<<" Contact angle type not found."<<std::endl;
		break;
	}
}
void D2Q9ColourFluid::Set_Collide(){
	PtrDicCollide=Dic;


if(PtrParameters->Get_ColourOperatorType()== ::SurfaceForce && PtrParameters->Get_UserForceType()== ::LocalForce)
{
	std::cout<<" Warming: The local user force will be ignored. The colour model with the surface force is incompatible with a local user force"<<std::endl;
}
	if(PtrParameters->Get_ViscosityType()==ConstViscosity)
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
		case DensityGrad:
			PtrColourGrad =&D2Q9ColourFluid::Colour_gradient_DensityGrad;
			PtrColourGradWall =&D2Q9ColourFluid::Colour_gradient_DensityGrad;
			PtrColourGradBc =&D2Q9ColourFluid::Colour_gradient_DensityGradBc;
			PtrColourGradCorner =&D2Q9ColourFluid::Colour_gradient_DensityGradCorner;
			break;
		case DensityNormalGrad:
			PtrColourGrad =&D2Q9ColourFluid::Colour_gradient_DensityNormalGrad;
			PtrColourGradWall =&D2Q9ColourFluid::Colour_gradient_DensityNormalGrad;
			PtrColourGradBc =&D2Q9ColourFluid::Colour_gradient_DensityNormalGradBc;
			PtrColourGradCorner =&D2Q9ColourFluid::Colour_gradient_DensityNormalGradCorner;
			break;
		}
	// Force teta value
/*		if(PtrParameters->Get_ContactAngleType()==FixTeta)
			PtrColourGradWall =&D2Q9ColourFluid::Colour_gradient_FixTeta;
		else if(PtrParameters->Get_ContactAngleType()==NonCstTeta)
			PtrColourGradWall =&D2Q9ColourFluid::Colour_gradient_NoCstTeta;
*/
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
	PtrCalNormal=&D2Q9ColourFluid::CalculNormal;
}
void D2Q9ColourFluid::InitColourFluid(InitLBM& ini){
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
	Bi[0]=-4/27;
	for(int i=1;i<5;i++) Bi[i]=2/27;
	for(int i=5;i<9;i++) Bi[i]=5/108;
	G=new double*[2];
	F=new double*[2];
	Normal=new double*[2];
	Dic->AddSync("Rho",Rho);
	Dic->AddVar(Vector,"ColourGrad",false,true,false,G[0],G[1]);
	Dic->AddVar(Scalar,"ColourGrad_Norm",true,true,false,G_Norm);
	Dic->AddVar(Vector,"Normal",false,true,false,Normal[0],Normal[1]);
	Dic->AddVar(Vector,"SurfaceForce",false,true,false,F[0],F[1]);
	Dic->AddVar(Scalar,"RhoRed",false,true,true,Rhor);
	Dic->AddVar(Scalar,"RhoBlue",false,true,true,Rhob);
	Dic->AddVar(Scalar,"RhoN",true,true,true,RhoN);
	Dic->AddVar(Scalar,"Curvature",false,false,false,Curv);
	Dic->AddVar(Scalar,"Check",false,false,false,testVar);
	for(int i=0;i<nbnodes_total;i++)
	{
		G[0][i]=0;
		G[1][i]=0;
		G_Norm[i]=0;
		Normal[0][i]=0;
		Normal[1][i]=0;
		testVar[i]=0.0;
		F[0][i]=0;
		F[1][i]=0;
	}

std::cout<<"Number of Solid first layer="<<NodeArrays->Solid1stLayer.size()<<std::endl
		<<"Number of concave corner"<<NodeArrays->CornerConcave.size()<<std::endl
		<<"Number of convex corner"<<NodeArrays->CornerConvex.size()<<std::endl;

	InitColourFluidAllDomain(ini);
}
void D2Q9ColourFluid::InitColourFluidAllDomain(InitLBM& ini){
	InitColourFluidDomainBc(ini);
	InitColourFluidWall(ini);
	InitColourFluidInterior(ini);

// Init Solid
	double alpha=0;
	double* pos =new double[2];
	double* U_=new double[2];
	int idx=0;

	for (int j=0;j<NodeArrays->NodeSolid.size();j++)
	{
// Set Index
		idx=NodeArrays->NodeSolid[j].Get_index();
// Get position
		pos[0]=NodeArrays->NodeSolid[j].get_x();
		pos[1]=NodeArrays->NodeSolid[j].get_y();
// Get initialise values from the user
		ini.IniDomainTwoPhases(parallel->getRank(),NodeArrays->NodeSolid[j],0, idx,pos,Rho[idx],U_,alpha);
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
void D2Q9ColourFluid::InitColourFluidDomainBc(InitLBM& ini){
	double alpha=0;
	double* pos =new double[2];
	double* U_=new double[2];
	int idx=0;
	for (int j=0;j<NodeArrays->NodeGlobalCorner.size();j++)
	{
// Set Index
		idx=NodeArrays->NodeGlobalCorner[j].Get_index();
// Get position
		pos[0]=NodeArrays->NodeGlobalCorner[j].get_x();
		pos[1]=NodeArrays->NodeGlobalCorner[j].get_y();
// Get initialise values from the user
		ini.IniDomainTwoPhases(parallel->getRank(),NodeArrays->NodeGlobalCorner[j],0, idx,pos,Rho[idx],U_,alpha);
// Initialise the blue and red densities
		Rhor[idx]=alpha*Rho[idx];
		Rhob[idx]=(1- alpha) *Rho[idx];
		RhoN[idx]=(Rhor[idx]-Rhob[idx])/Rho[idx];
		for (int i=0;i<nbvelo;i++)
		{
			f[0]->f[i][idx]=CollideLowOrder::EquiDistriFunct2D(Rhor[idx], U[0][idx], U[1][idx], &Ei[i][0], omega[i]);
			f[1]->f[i][idx]=CollideLowOrder::EquiDistriFunct2D(Rhob[idx], U[0][idx], U[1][idx], &Ei[i][0], omega[i]);
		}
	}
	for (int j=0;j<NodeArrays->NodeSpecialWall.size();j++)
	{
// Set Index
		idx=NodeArrays->NodeSpecialWall[j].Get_index();
// Get position
		pos[0]=NodeArrays->NodeSpecialWall[j].get_x();
		pos[1]=NodeArrays->NodeSpecialWall[j].get_y();
// Get initialise values from the user
		ini.IniDomainTwoPhases(parallel->getRank(),NodeArrays->NodeSpecialWall[j],0, idx,pos,Rho[idx],U_,alpha);
// Initialise the blue and red densities
		Rhor[idx]=alpha*Rho[idx];
		Rhob[idx]=(1- alpha) *Rho[idx];
		RhoN[idx]=(Rhor[idx]-Rhob[idx])/Rho[idx];
		for (int i=0;i<nbvelo;i++)
		{
			f[0]->f[i][idx]=CollideLowOrder::EquiDistriFunct2D(Rhor[idx], U[0][idx], U[1][idx], &Ei[i][0], omega[i]);
			f[1]->f[i][idx]=CollideLowOrder::EquiDistriFunct2D(Rhob[idx], U[0][idx], U[1][idx], &Ei[i][0], omega[i]);
		}
	}
	for (int j=0;j<NodeArrays->NodeVelocity.size();j++)
	{
// Set Index
		idx=NodeArrays->NodeVelocity[j].Get_index();
// Get position
		pos[0]=NodeArrays->NodeVelocity[j].get_x();
		pos[1]=NodeArrays->NodeVelocity[j].get_y();
// Get initialise values from the user
		ini.IniDomainTwoPhases(parallel->getRank(),NodeArrays->NodeVelocity[j],0, idx,pos,Rho[idx],U_,alpha);
// Initialise the blue and red densities
		Rhor[idx]=alpha*Rho[idx];
		Rhob[idx]=(1- alpha) *Rho[idx];
		RhoN[idx]=(Rhor[idx]-Rhob[idx])/Rho[idx];
		for (int i=0;i<nbvelo;i++)
		{
			f[0]->f[i][idx]=CollideLowOrder::EquiDistriFunct2D(Rhor[idx], U[0][idx], U[1][idx], &Ei[i][0], omega[i]);
			f[1]->f[i][idx]=CollideLowOrder::EquiDistriFunct2D(Rhob[idx], U[0][idx], U[1][idx], &Ei[i][0], omega[i]);
		}
	}

	for (int j=0;j<NodeArrays->NodePressure.size();j++)
	{
// Set Index
		idx=NodeArrays->NodePressure[j].Get_index();
// Get position
		pos[0]=NodeArrays->NodePressure[j].get_x();
		pos[1]=NodeArrays->NodePressure[j].get_y();
// Get initialise values from the user
		ini.IniDomainTwoPhases(parallel->getRank(),NodeArrays->NodePressure[j],0, idx,pos,Rho[idx],U_,alpha);
// Initialise the blue and red densities
		Rhor[idx]=alpha*Rho[idx];
		Rhob[idx]=(1- alpha) *Rho[idx];
		RhoN[idx]=(Rhor[idx]-Rhob[idx])/Rho[idx];
		for (int i=0;i<nbvelo;i++)
		{
			f[0]->f[i][idx]=CollideLowOrder::EquiDistriFunct2D(Rhor[idx], U[0][idx], U[1][idx], &Ei[i][0], omega[i]);
			f[1]->f[i][idx]=CollideLowOrder::EquiDistriFunct2D(Rhob[idx], U[0][idx], U[1][idx], &Ei[i][0], omega[i]);
		}
	}
	for (int j=0;j<NodeArrays->NodeSymmetry.size();j++)
	{
// Set Index
		idx=NodeArrays->NodeSymmetry[j].Get_index();
// Get position
		pos[0]=NodeArrays->NodeSymmetry[j].get_x();
		pos[1]=NodeArrays->NodeSymmetry[j].get_y();
// Get initialise values from the user
		ini.IniDomainTwoPhases(parallel->getRank(),NodeArrays->NodeSymmetry[j],0, idx,pos,Rho[idx],U_,alpha);
// Initialise the blue and red densities
		Rhor[idx]=alpha*Rho[idx];
		Rhob[idx]=(1- alpha) *Rho[idx];
		RhoN[idx]=(Rhor[idx]-Rhob[idx])/Rho[idx];
		for (int i=0;i<nbvelo;i++)
		{
			f[0]->f[i][idx]=CollideLowOrder::EquiDistriFunct2D(Rhor[idx], U[0][idx], U[1][idx], &Ei[i][0], omega[i]);
			f[1]->f[i][idx]=CollideLowOrder::EquiDistriFunct2D(Rhob[idx], U[0][idx], U[1][idx], &Ei[i][0], omega[i]);
		}
	}
	for (int j=0;j<NodeArrays->NodePeriodic.size();j++)
	{
// Set Index
		idx=NodeArrays->NodePeriodic[j].Get_index();
// Get position
		pos[0]=NodeArrays->NodePeriodic[j].get_x();
		pos[1]=NodeArrays->NodePeriodic[j].get_y();
// Get initialise values from the user
		ini.IniDomainTwoPhases(parallel->getRank(),NodeArrays->NodePeriodic[j],0, idx,pos,Rho[idx],U_,alpha);
// Initialise the blue and red densities
		Rhor[idx]=alpha*Rho[idx];
		Rhob[idx]=(1- alpha) *Rho[idx];
		RhoN[idx]=(Rhor[idx]-Rhob[idx])/Rho[idx];
		for (int i=0;i<nbvelo;i++)
		{
			f[0]->f[i][idx]=CollideLowOrder::EquiDistriFunct2D(Rhor[idx], U[0][idx], U[1][idx], &Ei[i][0], omega[i]);
			f[1]->f[i][idx]=CollideLowOrder::EquiDistriFunct2D(Rhob[idx], U[0][idx], U[1][idx], &Ei[i][0], omega[i]);
		}
	}
	delete [] pos;
	delete [] U_;
}
void D2Q9ColourFluid::InitColourFluidWall(InitLBM& ini){
	double alpha=0;
	double* pos =new double[2];
	double* U_=new double[2];
	int idx=0;

	for (int j=0;j<NodeArrays->NodeCorner.size();j++)
	{
// Set Index
		idx=NodeArrays->NodeCorner[j].Get_index();
// Get position
		pos[0]=NodeArrays->NodeCorner[j].get_x();
		pos[1]=NodeArrays->NodeCorner[j].get_y();
// Get initialise value from the user
		ini.IniDomainTwoPhases(parallel->getRank(),NodeArrays->NodeCorner[j],0, idx,pos,Rho[idx],U_,alpha);
// Set contact angle if needed
		if(PtrParameters->Get_ContactAngleType()==::NonCstTeta)
			ini.IniContactAngle(parallel->getRank(),NodeArrays->NodeWall[j],0, idx,pos,teta[idx]);
// Initialise the blue and red densities
		Rhor[idx]=alpha*Rho[idx];
		Rhob[idx]=(1- alpha) *Rho[idx];
		RhoN[idx]=(Rhor[idx]-Rhob[idx])/Rho[idx];
		for (int i=0;i<nbvelo;i++)
		{
			f[0]->f[i][idx]=CollideLowOrder::EquiDistriFunct2D(Rhor[idx], U[0][idx], U[1][idx], &Ei[i][0], omega[i]);
			f[1]->f[i][idx]=CollideLowOrder::EquiDistriFunct2D(Rhob[idx], U[0][idx], U[1][idx], &Ei[i][0], omega[i]);
		}
	if(PtrParameters->Get_ContactAngleType()==::FixTeta)
		{teta[0]=PtrParameters->Get_ContactAngle();costeta=std::cos(teta[0]);sinteta=std::sin(teta[0]);}
	}

	for (int j=0;j<NodeArrays->NodeWall.size();j++)
	{
// Set Index
		idx=NodeArrays->NodeWall[j].Get_index();
// Get position
		pos[0]=NodeArrays->NodeWall[j].get_x();
		pos[1]=NodeArrays->NodeWall[j].get_y();
// Get initialise values from the user
		ini.IniDomainTwoPhases(parallel->getRank(),NodeArrays->NodeWall[j],0, idx,pos,Rho[idx],U_,alpha);
// Set contact angle if needed
		if(PtrParameters->Get_ContactAngleType()==::NonCstTeta)
			ini.IniContactAngle(parallel->getRank(),NodeArrays->NodeWall[j],0, idx,pos,teta[idx]);
// Initialise the blue and red densities
		Rhor[idx]=alpha*Rho[idx];
		Rhob[idx]=(1- alpha) *Rho[idx];
		RhoN[idx]=(Rhor[idx]-Rhob[idx])/Rho[idx];
		for (int i=0;i<nbvelo;i++)
		{
			f[0]->f[i][idx]=CollideLowOrder::EquiDistriFunct2D(Rhor[idx], U[0][idx], U[1][idx], &Ei[i][0], omega[i]);
			f[1]->f[i][idx]=CollideLowOrder::EquiDistriFunct2D(Rhob[idx], U[0][idx], U[1][idx], &Ei[i][0], omega[i]);
		}
	}

	if(PtrParameters->Get_ContactAngleType()==::FixTeta)
	{
		teta[0]=PtrParameters->Get_ContactAngle();
		costeta=std::cos(teta[0]);
		sinteta=std::sin(teta[0]);
		if (std::abs(costeta)<epsilon) {
			costeta=0.0;
			sinteta=copysign(1.0,sinteta);
		}
		if (std::abs(costeta)>1.0-epsilon) {
			costeta=copysign(1.0,costeta);
			sinteta=0.0;
		}
		n1[0][0]=0.0;
		n2[0][0]=n1[0][0];
		n1[0][1]=0.0;
		n2[0][1]=n1[0][1];
		n1[1][0]=costeta;
		n2[1][0]=n1[1][0];
		n1[1][1]=sinteta;
		n2[1][1]=-n1[1][1];
		n1[2][0]=sinteta;
		n2[2][0]=-n1[2][0];
		n1[2][1]=costeta;
		n2[2][1]=n1[2][1];
		n1[3][0]=-costeta;
		n2[3][0]=n1[3][0];
		n1[3][1]=sinteta;
		n2[3][1]=-n1[3][1];
		n1[4][0]=sinteta;
		n2[4][0]=-n1[4][0];
		n1[4][1]=-costeta;
		n2[4][1]=n1[4][1];
		double InvSqrt2=1.0/std::sqrt(2.0);
		n1[5][0]=InvSqrt2*(costeta-sinteta);
		n1[5][1]=InvSqrt2*(costeta+sinteta);
		n2[5][0]=InvSqrt2*(costeta+sinteta);
		n2[5][1]=InvSqrt2*(costeta-sinteta);

		n1[6][0]=InvSqrt2*(-costeta+sinteta);
		n1[6][1]=InvSqrt2*(costeta+sinteta);
		n2[6][0]=InvSqrt2*(-costeta-sinteta);
		n2[6][1]=InvSqrt2*(costeta-sinteta);

		n1[7][0]=InvSqrt2*(-costeta-sinteta);
		n1[7][1]=InvSqrt2*(-costeta+sinteta);
		n2[7][0]=InvSqrt2*(-costeta+sinteta);
		n2[7][1]=InvSqrt2*(-costeta-sinteta);

		n1[8][0]=InvSqrt2*(costeta+sinteta);
		n1[8][1]=InvSqrt2*(-costeta+sinteta);
		n2[8][0]=InvSqrt2*(costeta-sinteta);
		n2[8][1]=InvSqrt2*(-costeta-sinteta);
	}
	delete [] pos;
	delete [] U_;
}
void D2Q9ColourFluid::InitColourFluidInterior(InitLBM& ini){
	double alpha=0;
	double* pos =new double[2];
	double* U_=new double[2];
	int idx=0;

	for (int j=0;j<NodeArrays->NodeInterior.size();j++)
	{
// Set Index
		idx=NodeArrays->NodeInterior[j].Get_index();
// Get position
		pos[0]=NodeArrays->NodeInterior[j].get_x();
		pos[1]=NodeArrays->NodeInterior[j].get_y();
// Get initialise value from the user
		ini.IniDomainTwoPhases(parallel->getRank(),NodeArrays->NodeInterior[j],0,idx,pos,Rho[idx],U_,alpha);
// Initialise the blue and red densities
		Rhor[idx]=alpha*Rho[idx];
		Rhob[idx]=(1- alpha) *Rho[idx];
		RhoN[idx]=(Rhor[idx]-Rhob[idx])/Rho[idx];
		for (int i=0;i<nbvelo;i++)
		{
			f[0]->f[i][idx]=CollideLowOrder::EquiDistriFunct2D(Rhor[idx], U[0][idx], U[1][idx], &Ei[i][0], omega[i]);
			f[1]->f[i][idx]=CollideLowOrder::EquiDistriFunct2D(Rhob[idx], U[0][idx], U[1][idx], &Ei[i][0], omega[i]);
		}
	}

	for (int j=0;j<NodeArrays->NodeGhost.size();j++)
	{
// Set Index
		idx=NodeArrays->NodeGhost[j].Get_index();
// Get position
		pos[0]=NodeArrays->NodeGhost[j].get_x();
		pos[1]=NodeArrays->NodeGhost[j].get_y();
// Get initialise values from the user
		ini.IniDomainTwoPhases(parallel->getRank(),NodeArrays->NodeGhost[j],0, idx,pos,Rho[idx],U_,alpha);
// Initialise the blue and red densities
		Rhor[idx]=alpha*Rho[idx];
		Rhob[idx]=(1- alpha) *Rho[idx];
		RhoN[idx]=(Rhor[idx]-Rhob[idx])/Rho[idx];
		for (int i=0;i<nbvelo;i++)
		{
			f[0]->f[i][idx]=CollideLowOrder::EquiDistriFunct2D(Rhor[idx], U[0][idx], U[1][idx], &Ei[i][0], omega[i]);
			f[1]->f[i][idx]=CollideLowOrder::EquiDistriFunct2D(Rhob[idx], U[0][idx], U[1][idx], &Ei[i][0], omega[i]);
		}
	}

	delete [] pos;
	delete [] U_;
}
void D2Q9ColourFluid::UpdateAllDomain(Parameters* UpdatedParam,InitLBM& ini){
	UpdateDomainBc(UpdatedParam,ini);
	UpdateWall(UpdatedParam,ini);
	UpdateInterior(UpdatedParam,ini);
}
void D2Q9ColourFluid::UpdateDomainBc(Parameters* UpdatedParam,InitLBM& ini){
	InitDomainBc(ini);
	InitColourFluidDomainBc(ini);
}
void D2Q9ColourFluid::UpdateWall(Parameters* UpdatedParam,InitLBM& ini){
	InitWall(ini);
	InitColourFluidWall(ini);
}
void D2Q9ColourFluid::UpdateInterior(Parameters* UpdatedParam,InitLBM& ini){
	InitInterior(ini);
	InitColourFluidInterior(ini);
}
void D2Q9ColourFluid::run(Parameters* UpdatedParam){

	PtrParameters=UpdatedParam;
	IniTau(PtrParameters);
	InvTau=Get_InvTau();

	//Initialise parameters for colour fluid
	beta=PtrParameters->Get_Beta();
	A1=PtrParameters->Get_A1();
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
		SyncToGhost();
	ApplyBc();
	Colour_gradient();
	Synchronise_Colour_gradient();
	Extrapolate_NormalInSolid();
	Impose_ContactAngleInSolid();
	ColourFluid_Collision();
	UpdateMacroVariables();
	if(parallel->getSize()>1)
		SyncMacroVarToGhost();
	Extrapol_Density_Corner();
	Writer->Write_Output(it);
	it++;
	if(parallel->getSize()>1)
	{

//		for (int i=1;i<NbStep+1;i++)
		while(it<NbStep+1)
		{
			Colour_gradient();
			Synchronise_Colour_gradient();
			Extrapolate_NormalInSolid();
			Impose_ContactAngleInSolid();
			ColourFluid_Collision();
			SyncToGhost();
			StreamD2Q9();;
			ApplyBc();
			UpdateMacroVariables();
			SyncMacroVarToGhost();
			Extrapol_Density_Corner();
			if(it%OutPutNStep==0)
			{
				Writer->Write_Output(it);
			}
			if(it%listing==0 )
			{
				Convergence::Calcul_Error();
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
					Writer->Write_Output(it);
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
			Colour_gradient();
			Extrapolate_NormalInSolid();
			Impose_ContactAngleInSolid();
			ColourFluid_Collision();
			StreamD2Q9();
			ApplyBc();
			UpdateMacroVariables();
			Extrapol_Density_Corner();
			if(it%OutPutNStep==0)
			{
				Writer->Write_Output(it);
			}
			if(it%listing==0)
			{
				Convergence::Calcul_Error();
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
					Writer->Write_Output(it);
					it=NbStep;
				}
			}
			it++;
		}
	}

	Writer->Write_breakpoint(*PtrParameters);
}
/*
double D2Q9ColourFluid::Calcul_Error(){
	double error=0;
	for (int j=0;j<NodeArrays->NodeInterior.size();j++)
	{
		error+=Error_RhoN(NodeArrays->NodeInterior[j].Get_index());
		//error+=std::abs(RhoN[NodeArrays->NodeInterior[j].Get_index()]-RhoN[NodeArrays->NodeInterior[j].Get_index()])/(std::abs(RhoN[NodeArrays->NodeInterior[j].Get_index()])+1e-15);
	}
	for (int j=0;j<NodeArrays->NodeCorner.size();j++)
	{
		error+=Error_RhoN(NodeArrays->NodeCorner[j].Get_index());
	}
	for (int j=0;j<NodeArrays->NodeGlobalCorner.size();j++)
	{
		error+=Error_RhoN(NodeArrays->NodeGlobalCorner[j].Get_index());
	}
	for (int j=0;j<NodeArrays->NodeVelocity.size();j++)
	{
		error+=Error_RhoN(NodeArrays->NodeVelocity[j].Get_index());
	}
	for (int j=0;j<NodeArrays->NodePressure.size();j++)
	{
		error+=Error_RhoN(NodeArrays->NodePressure[j].Get_index());
	}
	for (int j=0;j<NodeArrays->NodeWall.size();j++)
	{
		error+=Error_RhoN(NodeArrays->NodeWall[j].Get_index());
	}
	for (int j=0;j<NodeArrays->NodeSymmetry.size();j++)
	{
		error+=Error_RhoN(NodeArrays->NodeSymmetry[j].Get_index());
	}
	for (int j=0;j<NodeArrays->NodePeriodic.size();j++)
	{
		error+=Error_RhoN(NodeArrays->NodePeriodic[j].Get_index());
	}
	return error;
}
*/
///Calculate the colour gradient by Gunstensen formulation, density gradient or normal density gradient
void D2Q9ColourFluid::Colour_gradient(){
	int idx_tmp;int normal_interior=0;
	double teta;
	for (int j=0;j<NodeArrays->NodeInterior.size();j++)
	{
	// Common variables
		idx_tmp=NodeArrays->NodeInterior[j].Get_index();
	// Calculate gradients
		(this->*PtrColourGrad)(idx_tmp,NodeArrays->NodeInterior[j].Get_connect(),normal_interior);
		(this->*PtrCalNormal)(idx_tmp,NodeArrays->NodeInterior[j].Get_connect(),normal_interior);
	}

	for (int j=0;j<NodeArrays->NodeCorner.size();j++)
	{
	// Common variables
		idx_tmp=NodeArrays->NodeCorner[j].Get_index();
	// Calculate gradients
		(this->*PtrColourGradCorner)(idx_tmp,NodeArrays->NodeCorner[j].Get_connect(),NodeArrays->NodeCorner[j].Get_BcNormal());
		(this->*PtrCalNormal)(idx_tmp,NodeArrays->NodeCorner[j].Get_connect(),NodeArrays->NodeCorner[j].Get_BcNormal());
	}

	for (int j=0;j<NodeArrays->NodeGlobalCorner.size();j++)
	{
	// Common variables
		idx_tmp=NodeArrays->NodeGlobalCorner[j].Get_index();
	// Calculate gradients
		(this->*PtrColourGradBc)(idx_tmp,NodeArrays->NodeGlobalCorner[j].Get_connect(),NodeArrays->NodeGlobalCorner[j].Get_BcNormal());
		(this->*PtrCalNormal)(idx_tmp,NodeArrays->NodeGlobalCorner[j].Get_connect(),NodeArrays->NodeGlobalCorner[j].Get_BcNormal());
	}

	for (int j=0;j<NodeArrays->NodeVelocity.size();j++)
	{
	// Common variables
		idx_tmp=NodeArrays->NodeVelocity[j].Get_index();
	// Calculate gradients
		(this->*PtrColourGradBc)(idx_tmp,NodeArrays->NodeVelocity[j].Get_connect(),NodeArrays->NodeVelocity[j].Get_BcNormal());
		(this->*PtrCalNormal)(idx_tmp,NodeArrays->NodeVelocity[j].Get_connect(),NodeArrays->NodeVelocity[j].Get_BcNormal());
	}

	for (int j=0;j<NodeArrays->NodePressure.size();j++)
	{
	// Common variables
		idx_tmp=NodeArrays->NodePressure[j].Get_index();
	// Calculate gradients
		(this->*PtrColourGradBc)(idx_tmp,NodeArrays->NodePressure[j].Get_connect(),NodeArrays->NodePressure[j].Get_BcNormal());
		(this->*PtrCalNormal)(idx_tmp,NodeArrays->NodePressure[j].Get_connect(),NodeArrays->NodePressure[j].Get_BcNormal());
	}

	for (int j=0;j<NodeArrays->NodeWall.size();j++)
	{
	// Common variables
		idx_tmp=NodeArrays->NodeWall[j].Get_index();
	// Calculate gradients
		(this->*PtrColourGradWall)(idx_tmp,NodeArrays->NodeWall[j].Get_connect(),NodeArrays->NodeWall[j].Get_BcNormal());
		(this->*PtrCalNormalWall)(idx_tmp,NodeArrays->NodeWall[j].Get_connect(),NodeArrays->NodeWall[j].Get_BcNormal());
	}

	for (int j=0;j<NodeArrays->NodeSymmetry.size();j++)
	{
	// Common variables
		idx_tmp=NodeArrays->NodeSymmetry[j].Get_index();
	// Calculate gradients
		(this->*PtrColourGradBc)(idx_tmp,NodeArrays->NodeSymmetry[j].Get_connect(),NodeArrays->NodeSymmetry[j].Get_BcNormal());
		(this->*PtrCalNormal)(idx_tmp,NodeArrays->NodeSymmetry[j].Get_connect(),NodeArrays->NodeSymmetry[j].Get_BcNormal());
	}
	for (int j=0;j<NodeArrays->NodePeriodic.size();j++)
	{
	// Common variables
		idx_tmp=NodeArrays->NodePeriodic[j].Get_index();
	// Calculate gradients
		(this->*PtrColourGrad)(idx_tmp,NodeArrays->NodePeriodic[j].Get_connect(),NodeArrays->NodePeriodic[j].Get_BcNormal());
		(this->*PtrCalNormal)(idx_tmp,NodeArrays->NodePeriodic[j].Get_connect(),NodeArrays->NodePeriodic[j].Get_BcNormal());
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
//		G_Norm[idx_tmp]=0;
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
//		G_Norm[idx_tmp]=0;
//		G_Norm=std::sqrt(G[idx_tmp][0]*G[idx_tmp][0]+G[idx_tmp][1]*G[idx_tmp][1]);
	//Model
		(this->*PtrCollisionWall)(idx_tmp, NodeArrays->NodeWall[j].Get_connect(),NodeArrays->NodeWall[j].Get_BcNormal(),&fi_tmp[0]);
		//(this->*PtrCalNormalWall)(idx_tmp,NodeArrays->NodeWall[j].Get_BcNormal());
//		(this->*PtrCalNormalWall)(idx_tmp,NodeArrays->NodeWall[j].Get_connect(),NodeArrays->NodeWall[j].Get_BcNormal());

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
	for (int j=0;j<NodeArrays->NodePeriodic.size();j++)
	{
	// Common variables
		idx_tmp=NodeArrays->NodePeriodic[j].Get_index();

	//Model
		(this->*PtrCollision)(idx_tmp, NodeArrays->NodePeriodic[j].Get_connect(),NodeArrays->NodePeriodic[j].Get_BcNormal(),&fi_tmp[0]);
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

void D2Q9ColourFluid::CalculNormal_NoTeta(int & nodenumber, int* connect,int & normal){
	if(G_Norm[nodenumber]>0)
	{
		Normal[0][nodenumber]=G[0][nodenumber]/G_Norm[nodenumber];
		Normal[1][nodenumber]=G[1][nodenumber]/G_Norm[nodenumber];
	}
	/*if(G_Norm[nodenumber]>0)
	{
		switch(normal)
		{
		case 1:

			Normal[0][nodenumber]=0.0;
			Normal[1][nodenumber]=copysign(1.0,G[1][nodenumber]);
			break;
		case 2:
			Normal[0][nodenumber]=copysign(1.0,G[0][nodenumber]);
			Normal[1][nodenumber]=0.0;
			Normal[0][connect[4]]=Normal[0][nodenumber];
			Normal[1][connect[4]]=Normal[1][nodenumber];
			break;
		case 3:
			Normal[0][nodenumber]=0.0;
			Normal[1][nodenumber]=copysign(1.0,G[1][nodenumber]);
			break;
		case 4:
			Normal[0][nodenumber]=copysign(1.0,G[0][nodenumber]);
			Normal[1][nodenumber]=0.0;
			break;
		default:
			std::cerr<<" Normal for Colour gradient Fix Teta is not found"<<std::endl;
		}
	}
	else
	{
		Normal[0][nodenumber]=0;
		Normal[1][nodenumber]=0;
	}*/
}
void D2Q9ColourFluid::CalculNormal_FixTeta(int & nodenumber, int* connect,int & normal){
/*	switch(normal)
	{
	case 1:
		G[0][nodenumber]=-cos(teta[0])*G_Norm[nodenumber];
		G[1][nodenumber]=copysign(sin(teta[0])*G_Norm[nodenumber],G[1][nodenumber]);
		break;
	case 2:
		G[0][nodenumber]=copysign(sin(teta[0])*G_Norm[nodenumber],G[0][nodenumber]);
		G[1][nodenumber]=cos(teta[0])*G_Norm[nodenumber];
		G[0][connect[2]]=G[0][nodenumber];//copysign(sin(teta[0])*G_Norm[connect[4]],G[0][connect[4]]);
		G[1][connect[2]]=G[1][nodenumber];//-cos(teta[0])*G_Norm[connect[2]];
		if(G_Norm[nodenumber]>0)
		{
			Normal[0][connect[2]]=G[0][connect[2]]/G_Norm[nodenumber];
			Normal[1][connect[2]]=G[1][connect[2]]/G_Norm[nodenumber];
		}
		break;
	case 3:
		G[0][nodenumber]=cos(teta[0])*G_Norm[nodenumber];
		G[1][nodenumber]=copysign(sin(teta[0])*G_Norm[nodenumber],G[1][nodenumber]);
		break;
	case 4:
		G[0][nodenumber]=copysign(sin(teta[0])*G_Norm[nodenumber],G[0][nodenumber]);
		G[1][nodenumber]=cos(teta[0])*G_Norm[nodenumber];
		break;
	default:
		std::cerr<<" Normal for Colour gradient Fix Teta is not found"<<std::endl;
	}*/
	if(G_Norm[nodenumber]>0)
	{
		Normal[0][nodenumber]=G[0][nodenumber]/G_Norm[nodenumber];
		Normal[1][nodenumber]=G[1][nodenumber]/G_Norm[nodenumber];
	}
/*	if(G_Norm[nodenumber]>0.0)
	{
		Normal[0][nodenumber]*=copysign(G_Norm[nodenumber],G[0][nodenumber]);
		Normal[1][nodenumber]*=-G_Norm[nodenumber];//copysign(G_Norm[nodenumber],G[1][nodenumber]);
	}
	else
	{
		Normal[0][nodenumber]=0;
		Normal[1][nodenumber]=0;
	}*/
//	G[0][nodenumber]=G_Norm[nodenumber]*Normal[0][nodenumber];
//	G[1][nodenumber]=-G_Norm[nodenumber]*Normal[1][nodenumber];
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
}
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

	if(Rhor[nodenumber]<=Rho_limiter)
	{
		for(int i=0;i<nbvelo;i++)
			{f[0]->f[i][nodenumber]=0.0;f[1]->f[i][nodenumber]=fi_tmp[i];}
	}
	else if(Rhob[nodenumber]<=Rho_limiter)
	{
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
		if(G_Norm[nodenumber]>0)
			for(int i=1;i<nbvelo;i++)
			{
					f[0]->f[i][nodenumber]=Rhor_Rho*fi_tmp[i]
								+factor*omega[i]*CosPhi(nodenumber,i,G_Norm[nodenumber]);
					f[1]->f[i][nodenumber]=fi_tmp[i]-f[0]->f[i][nodenumber];//Rhob_Rho*fi_tmp[i]
								//-factor*omega[i]*CosPhi(nodenumber,i,G_Norm[nodenumber]);
			}
		else
			for(int i=1;i<nbvelo;i++)
			{
				f[0]->f[i][nodenumber]=Rhor_Rho*fi_tmp[i];
				f[1]->f[i][nodenumber]=fi_tmp[i]-f[0]->f[i][nodenumber];//=Rhob_Rho*fi_tmp[i];
			}
	}
}
void D2Q9ColourFluid::Recolouring_Wall(int & nodenumber, double * fi_tmp){
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



}
double D2Q9ColourFluid::CosPhi(int nodenumber, int & direction,double & F_Norm){
		return (G[0][nodenumber]* Ei[direction][0]+G[1][nodenumber]* Ei[direction][1])/(Ei_Norm[direction]*G_Norm[nodenumber]);//*G_Norm[nodenumber]///(F_Norm*std::sqrt(Ei[direction][0]*Ei[direction][0]+Ei[direction][1]*Ei[direction][1]));
}

double& D2Q9ColourFluid::Collision_operator_Grunau(int & i, int & nodenumber, double Ak){
	if(G_Norm[nodenumber]>0)
		D_tmp=Ak*0.5*G_Norm[nodenumber]*(((G[0][nodenumber]*Ei[i][0]+G[1][nodenumber]*Ei[i][1])/G_Norm[nodenumber])*((G[0][nodenumber]*Ei[i][0]+G[1][nodenumber]*Ei[i][1])/G_Norm[nodenumber])-3/4);
	return D_tmp;
}
double& D2Q9ColourFluid::Collision_operator_Reis(int & i, int & nodenumber, double Ak){
	if(G_Norm[nodenumber]>0)
		D_tmp=Ak*0.5*G_Norm[nodenumber]*((omega[i]*(G[0][nodenumber]*Ei[i][0]+G[1][nodenumber]*Ei[i][1])/G_Norm[nodenumber])*((G[0][nodenumber]*Ei[i][0]+G[1][nodenumber]*Ei[i][1])/G_Norm[nodenumber])-Bi[i]);
	return D_tmp;
}
void D2Q9ColourFluid::Collision_Grunau(int & nodenumber, int* connect,int & normal,double* fi){
	int i0=0;
	fi[0]=f[0]->f[0][nodenumber]+f[1]->f[0][nodenumber];
	DVec_2D_tmp[0]=0;DVec_2D_tmp[1]=0;
	Collide_2D(i0, fi[0],Rho[nodenumber], U[0][nodenumber], U[1][nodenumber], DVec_2D_tmp[0],DVec_2D_tmp[1], Get_InvTau(Rho[nodenumber],RhoN[nodenumber]));
	for (int i=1;i<9;i++)
	{
		//Save the mixture distribution for recolouring
		fi[i]=f[0]->f[i][nodenumber]+f[1]->f[i][nodenumber];
		Collide_2D(i, fi[i],Rho[nodenumber], U[0][nodenumber], U[1][nodenumber], Collision_operator_Grunau(i, nodenumber, A1),DVec_2D_tmp[1], Get_InvTau(Rho[nodenumber],RhoN[nodenumber]));
	}
}
void D2Q9ColourFluid::Collision_Reis(int & nodenumber, int* connect,int & normal,double* fi){
	int i0=0;
	fi[0]=f[0]->f[0][nodenumber]+f[1]->f[0][nodenumber];
	DVec_2D_tmp[0]=0;DVec_2D_tmp[1]=0;
	Collide_2D(i0, fi[0],Rho[nodenumber], U[0][nodenumber], U[1][nodenumber], DVec_2D_tmp[0],DVec_2D_tmp[1], Get_InvTau(Rho[nodenumber],RhoN[nodenumber]));
	for (int i=1;i<9;i++)
	{
		//Save the mixture distribution for recolouring
		fi[i]=f[0]->f[i][nodenumber]+f[1]->f[i][nodenumber];
		Collide_2D(i, fi[i],Rho[nodenumber], U[0][nodenumber], U[1][nodenumber], Collision_operator_Reis(i, nodenumber, A1),DVec_2D_tmp[1], Get_InvTau(Rho[nodenumber],RhoN[nodenumber]));
	}
}
void D2Q9ColourFluid::Collision_SurfaceForce(int & nodenumber, int* connect,int & normal,double* fi){
	if(G_Norm[nodenumber]>0)
	{
		SurfaceForce(nodenumber,connect,normal,F[0][nodenumber],F[1][nodenumber]);
	}
	else
		{F[0][nodenumber]=0;F[1][nodenumber]=0;}
	U[0][nodenumber]=U[0][nodenumber]+0.5*F[0][nodenumber]/Rho[nodenumber];
	U[1][nodenumber]=U[1][nodenumber]+0.5*F[1][nodenumber]/Rho[nodenumber];
	for (int i=0;i<9;i++)
	{
		//Save the mixture distribution for recolouring
		fi[i]=f[0]->f[i][nodenumber]+f[1]->f[i][nodenumber];
		Collide_2D(i, fi[i],Rho[nodenumber], U[0][nodenumber], U[1][nodenumber], F[0][nodenumber],F[1][nodenumber], Get_InvTau(Rho[nodenumber],RhoN[nodenumber]));
	}
}

void D2Q9ColourFluid::SurfaceForce(int & nodenumber, int* connect,int & normal,double & Fx,double & Fy){
	Curv[nodenumber]=Curvature(nodenumber,connect,normal);
	D_tmp=0.5*tension*Curv[nodenumber];//*G_Norm[nodenumber]; //G_Norm is to get back the density gradient and not the normalise one
	Fx=D_tmp*G[0][nodenumber];
	Fy=D_tmp*G[1][nodenumber];
}


 void D2Q9ColourFluid::Collision_SurfaceForceWall(int & nodenumber, int* connect,int & normal,double* fi){
	if(G_Norm[nodenumber]>0)
	{
		SurfaceForceWall(nodenumber,connect,normal,F[0][nodenumber],F[1][nodenumber]);
	}
	else
		{F[0][nodenumber]=0;F[1][nodenumber]=0;}
	U[0][nodenumber]=U[0][nodenumber]+0.5*F[0][nodenumber]/Rho[nodenumber];
	U[1][nodenumber]=U[1][nodenumber]+0.5*F[1][nodenumber]/Rho[nodenumber];
	for (int i=0;i<9;i++)
	{
		//Save the mixture distribution for recolouring
		fi[i]=f[0]->f[i][nodenumber]+f[1]->f[i][nodenumber];
		//fi[i]=fi[i]+Collide_2D_BodyForce_Non_Constant_Tau(i, U[0][nodenumber], U[1][nodenumber], F[0][nodenumber],F[1][nodenumber], Get_InvTau(Rho[nodenumber],RhoN[nodenumber]));
		Collide_2D(i, fi[i],Rho[nodenumber], U[0][nodenumber], U[1][nodenumber], F[0][nodenumber],F[1][nodenumber], Get_InvTau(Rho[nodenumber],RhoN[nodenumber]));
	}
}

void D2Q9ColourFluid::SurfaceForceWall(int & nodenumber, int* connect,int & normal,double & Fx,double & Fy){

	Curv[nodenumber]=CurvatureWall(nodenumber,connect,normal);
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
double D2Q9ColourFluid::CurvatureWall(int & nodenumber, int* connect,int & normal){
	double tmp_0[2],tmp_1[2];
	GradBc(&tmp_0[0],&Normal[0][0],connect,normal);//d(Nx)/dx and d(Nx)/dy
	GradBc(&tmp_1[0],&Normal[1][0],connect,normal);//d(Ny)/dx and d(Ny)/dy
	return -(tmp_0[0]+tmp_1[1]);//Dot product

}

///Select and apply boundary conditions
void D2Q9ColourFluid::ApplyBc(){

	double fi_tmp;
	for (int j=0;j<NodeArrays->NodeVelocity.size();j++)
	{
		ApplyVelocity(NodeArrays->NodeVelocity[j].Get_BcNormal(),NodeArrays->NodeVelocity[j].Get_connect(),NodeArrays->NodeVelocity[j].Get_UDef(), f[0],Rhor,U[0],U[1]);
		ApplyVelocity(NodeArrays->NodeVelocity[j].Get_BcNormal(),NodeArrays->NodeVelocity[j].Get_connect(),NodeArrays->NodeVelocity[j].Get_UDef(), f[1],Rhob,U[0],U[1]);

//Impose alpha
		Rhor[NodeArrays->NodeVelocity[j].Get_index()]=0.0;
		Rhob[NodeArrays->NodeVelocity[j].Get_index()]=0.0;
		for(int i=0;i<9;i++)
		{
			fi_tmp=f[0]->f[i][NodeArrays->NodeVelocity[j].Get_index()]+f[1]->f[i][NodeArrays->NodeVelocity[j].Get_index()];
			f[0]->f[i][NodeArrays->NodeVelocity[j].Get_index()]=NodeArrays->NodeVelocity[j].Get_AlphaDef()*fi_tmp;
			f[1]->f[i][NodeArrays->NodeVelocity[j].Get_index()]=(1-NodeArrays->NodeVelocity[j].Get_AlphaDef())*fi_tmp;
			Rhor[NodeArrays->NodeVelocity[j].Get_index()]+=f[0]->f[i][NodeArrays->NodeVelocity[j].Get_index()];
			Rhob[NodeArrays->NodeVelocity[j].Get_index()]+=f[1]->f[i][NodeArrays->NodeVelocity[j].Get_index()];

		}

	}

	double alpha=0;
	for (int j=0;j<NodeArrays->NodePressure.size();j++)
	{
//		ApplyPressure(NodeArrays->NodePressure[j].Get_BcNormal(),NodeArrays->NodePressure[j].Get_connect(),NodeArrays->NodePressure[j].Get_AlphaDef()*NodeArrays->NodePressure[j].Get_RhoDef(), f[0],Rhor,U[0],U[1]);
//		ApplyPressure(NodeArrays->NodePressure[j].Get_BcNormal(),NodeArrays->NodePressure[j].Get_connect(),(1-NodeArrays->NodePressure[j].Get_AlphaDef())*NodeArrays->NodePressure[j].Get_RhoDef(), f[1],Rhob,U[0],U[1]);
		ExtrapolationOnWall(RhoN,NodeArrays->NodePressure[j].Get_connect(),NodeArrays->NodePressure[j].Get_BcNormal());
		alpha=(RhoN[NodeArrays->NodePressure[j].Get_index()]+1.0)*0.5;
		ApplyPressure(NodeArrays->NodePressure[j].Get_BcNormal(),NodeArrays->NodePressure[j].Get_connect(),alpha*NodeArrays->NodePressure[j].Get_RhoDef(), f[0],Rhor,U[0],U[1]);
		ApplyPressure(NodeArrays->NodePressure[j].Get_BcNormal(),NodeArrays->NodePressure[j].Get_connect(),(1-alpha)*NodeArrays->NodePressure[j].Get_RhoDef(), f[1],Rhob,U[0],U[1]);

	}
	for (int j=0;j<NodeArrays->NodeWall.size();j++)
	{
		ApplyWall(NodeArrays->NodeWall[j].Get_BcNormal(),NodeArrays->NodeWall[j].Get_connect(),f[0],Rhor,U[0],U[1]);
		ApplyWall(NodeArrays->NodeWall[j].Get_BcNormal(),NodeArrays->NodeWall[j].Get_connect(),f[1],Rhob,U[0],U[1]);
	}
	for (int j=0;j<NodeArrays->NodeSpecialWall.size();j++)
	{
		ExtrapolationCornerConcaveToSolid(Rhor,NodeArrays->NodeSpecialWall[j].Get_connect(),NodeArrays->NodeSpecialWall[j].Get_BcNormal());
		ExtrapolationCornerConcaveToSolid(Rhob,NodeArrays->NodeSpecialWall[j].Get_connect(),NodeArrays->NodeSpecialWall[j].Get_BcNormal());
		ExtrapolationCornerConcaveToSolid(U[0],NodeArrays->NodeSpecialWall[j].Get_connect(),NodeArrays->NodeSpecialWall[j].Get_BcNormal());
		ExtrapolationCornerConcaveToSolid(U[1],NodeArrays->NodeSpecialWall[j].Get_connect(),NodeArrays->NodeSpecialWall[j].Get_BcNormal());

//		ApplyGlobalCorner(NodeArrays->NodeGlobalCorner[j],NodeArrays->NodeGlobalCorner[j].Get_AlphaDef()*NodeArrays->NodeGlobalCorner[j].Get_RhoDef(),NodeArrays->NodeGlobalCorner[j].Get_UDef(),NodeArrays->TypeOfNode,f[0],Rhor,U[0],U[1]);

		ApplySpecialWall(NodeArrays->NodeSpecialWall[j],Rhor[NodeArrays->NodeSpecialWall[j].Get_index()],U[0][NodeArrays->NodeSpecialWall[j].Get_index()],U[1][NodeArrays->NodeSpecialWall[j].Get_index()],NodeArrays->TypeOfNode,f[0],Rhor,U[0],U[1]);
		ApplySpecialWall(NodeArrays->NodeSpecialWall[j],Rhob[NodeArrays->NodeSpecialWall[j].Get_index()],U[0][NodeArrays->NodeSpecialWall[j].Get_index()],U[1][NodeArrays->NodeSpecialWall[j].Get_index()],NodeArrays->TypeOfNode,f[1],Rhob,U[0],U[1]);

//		ApplySpecialWall(NodeArrays->NodeSpecialWall[j],NodeArrays->TypeOfNode,f[1],Rhob,U[0],U[1]);
//		ApplySpecialWall(NodeArrays->NodeSpecialWall[j],f[0],NodeArrays->TypeOfNode,Rhor,U[0],U[1]);
//		ApplySpecialWall(NodeArrays->NodeSpecialWall[j],f[1],NodeArrays->TypeOfNode,Rhob,U[0],U[1]);
	}
	for (int j=0;j<NodeArrays->NodeCorner.size();j++)
	{
		ApplyCornerWall(NodeArrays->NodeCorner[j], f[0],Rhor,U[0],U[1]);
		ApplyCornerWall(NodeArrays->NodeCorner[j], f[1],Rhob,U[0],U[1]);
	}
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
	parallel->barrier();

}
void D2Q9ColourFluid::Extrapol_Density_Corner(){
	//Extrapolation of density or normal density at the first layer in the solid
			for (int j=0;j<NodeArrays->NodeWall.size();j++)
			{
				(this->*PtrExtrapolDensity)(j,NodeArrays->NodeWall[j].Get_index(),NodeArrays->NodeWall[j].Get_connect(),NodeArrays->NodeWall[j].Get_BcNormal());
			}
			for (int j=0;j<NodeArrays->NodeCorner.size();j++)
			{
				(this->*PtrExtrapolDensity)(j,NodeArrays->NodeCorner[j].Get_index(),NodeArrays->NodeCorner[j].Get_connect(),NodeArrays->NodeCorner[j].Get_BcNormal());
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
		for (int j=0;j<NodeArrays->NodePeriodic.size();j++)
		{
			(this->*PtrMacro)(NodeArrays->NodePeriodic[j].Get_index());
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
