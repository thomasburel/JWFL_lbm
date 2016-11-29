/*
 ============================================================================
 Name        : LBM.c
 Author      : thomas
 Version     :
 Copyright   : Your copyright notice
 Description : LBM Code Developing by JWFL at University of Strathclyde
 ============================================================================
*/

#include <cmath>
#include <iostream>
#include <stdlib.h>

#include "Core/Simulation.h"

using namespace std;

 
int main(int argc, char *argv[]) {

	int idx=0;
/*	if (argv[1]==0)
	{
		std::cout<<"No argument. Index set to 1"<<std::endl;
		idx=1;
	}
	else
	{
		idx=atoi(argv[1]);
	}*/

	double pi=atan(1.0)*4.0 ;
stringstream FileExportStream;
// fluid 2 is the continuous fluid	and fluid 1 the droplet
//Domain size
	double L=1289;
	double H=370;

	double Diameter=H;
	double Ca=0;

	double U2_ref=0.01;
	double Re=1;

	double Rho1_ref=1;
	double Rho2_ref=Rho1_ref;

	double nu_2=1.0/6.0;
	double lambda=1;
	double nu_1=lambda*nu_2;

	double tau1=(6.0*nu_1+1.0)*0.5;
	double tau2=(6.0*nu_2+1.0)*0.5;
	double La=0;
	double sigma=0.0001;//La*(Rho1_ref*nu_1)*2/Diameter;//0.001;//0.01;//0.001;


	double Mach =U2_ref*std::sqrt(3);
	double Kn=(Mach/Re)*std::sqrt(pi/2.0);
//Pressure drop
	double deltaP=0.000001;
// Pressure inlet
	double Pmax=1+deltaP;
// Pressure outlet
	double Pmin=1-deltaP;

/// Create object for the simulation.
	Simulation simu;
	Parameters Param;

/// Set Simulation Arguments
	Param.Set_Arguments(&argc,&argv,false);

/// Set Solver type
// Dimension
	Param.Set_Dimension(SolverEnum::D2);
// Scheme
	Param.Set_Scheme(SolverEnum::Q9);

// Set Domain size
	Param.Set_Domain_Size((int)L,(int)H); //Cells

// Set User Parameters
	//U2_ref=0.01;Pmax=1;Pmin=1;
	double contactangle=120*pi/180.0;
	Param.Set_ContactAngleType(FixTeta);//NoTeta, FixTeta or NonCstTeta
	if(Param.Get_ContactAngleType()==NoTeta)
		contactangle=pi/2.0;
	Param.Set_ContactAngle(contactangle);

	Param.Set_UserDroplets(contactangle,sigma,Diameter,Re,Ca);
	Param.Set_UserParameters(U2_ref,H,L,Pmax,Pmin);

// Set User Force type
	Param.Set_UserForceType(None);
// Set delta x for output
	Param.Set_deltax(1);


/// Set Boundary condition type for the boundaries of the domain
/// Boundary condition accepted: Wall, Pressure, Velocity and Symmetry
/// Order Bottom, Right, Top, Left, (Front, Back for 3D)
	Param.Set_BcType(Symmetry,Pressure,Symmetry,Pressure);


/// Set Pressure Type
	Param.Set_PressureType(FixP);//FixP,zeroPGrad1st
/// Set Global Corner type
	Param.Set_CornerPressureType(ExtrapolCP);//FixCP,ExtrapolCP
/// Wall boundary condition type (Implemented BounceBack and Diffuse)
	Param.Set_WallType(BounceBack);//BounceBack,HeZouWall
	Param.Set_VelocityModel(HeZouV);
/// Set Periodic boundary condition to add a pressure drop term
	Param.Set_PeriodicType(PressureForce);//Simple,PressureForce
	Param.Set_PressureDrop(deltaP);
/// Number of maximum timestep
	Param.Set_NbStep(500000);
/// Interval for output
	Param.Set_OutPutNSteps(10000);// interval
///Display information during the calculation every N iteration
	Param.Set_listing(500);
	Param.Set_ErrorMax(1e-11);
	Param.Set_ErrorVariable(SolverEnum::VelocityX);

///Selection of variables to export
	Param.Set_VariablesOutput(true,true);// export Rho,U

/// Define the Output filename
/*	FileExportStream<<"Droplet_shear_"<< fixed << setprecision(0)<<H<<"x"<<L
			<<setprecision(2)<<"Re_"<<Re<<"_Ca_"<<Ca<<"_Conf_"<<confinement<<"_lambda_"<<viscosity_ratio
			<< scientific<<setprecision(3)<<"sigma_"<<sigma;*/
	FileExportStream.str("");
	FileExportStream<<"Debug_poiseuille_serpentine";
	Param.Set_OutputFileName(FileExportStream.str());

	// Multiphase model (SinglePhase or ColourFluid)
	Param.Set_Model(SolverEnum::SinglePhase);
	Param.Set_ViscosityType(HarmonicViscosity);//ConstViscosity,HarmonicViscosity

	//Gradient definition
	Param.Set_GradientType(LBMStencil); //FD or LBMStencil
	Param.Set_ExtrapolationType(WeightDistanceExtrapol);//NoExtrapol,TailorExtrapol,WeightDistanceExtrapol

/// Singlephase Parameters
	Param.Set_Tau(tau1);
	Param.Set_Rho(Rho1_ref);

/// Multiphase Parameters
	// Normal density output
	Param.Set_NormalDensityOutput(true);
	//Density of each fluid
	Param.Set_Rho_1(Rho1_ref);
	Param.Set_Rho_2(Rho2_ref);
	//Relaxation time for each fluid
	Param.Set_Tau_1(tau1);
	Param.Set_Tau_2(tau2);
	//Surface tension
	Param.Set_SurfaceTension(sigma);
	//Colour fluid Parameters
	Param.Set_A1(0.00001);
	Param.Set_A2(Param.Get_A1());
	Param.Set_Beta(0.7);// Between 0 and 1
	Param.Set_ColourGradType(DensityNormalGrad);//Gunstensen or DensityGrad or DensityNormalGrad
	Param.Set_RecolouringType(LatvaKokkoRothman);
	Param.Set_ColourOperatorType(SurfaceForce);//Grunau or Reis or SurfaceForce





/// Initialise the simulation with parameters
	simu.InitSimu(Param, true);
	simu.barrier();

		if(simu.Is_MainProcessor())
			cout<<"Reynolds: "<<Re<<" Surface tension is: "<<sigma<< " Tau 2: "<<tau2<<" U 1: "<<U2_ref<<" Number of cell in x direction: "<<L<<" Number of cell in y direction: "<<H<<std::endl;

//	simu.UpdateAllDomain(Param);
/// Run the simulation
	simu.RunSimu(Param);

//	}
}

