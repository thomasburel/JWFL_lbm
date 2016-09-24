/*
 ============================================================================
 Name        : LBM.c
 Author      : thomas
 Version     :
 Copyright   : Your copyright notice
 Description : LBM Code Developing by JWFL at University of Strathclyde
 ============================================================================
*/

#include <math.h> 
//#include "mpi.h"
#include <iostream>
//#include "Mesh/SingleBlock/Block2D.h"
//#include "Parrallelism/MpiManager.h"
//#include "Core/Parameters.h"
#include "Core/Simulation.h"
//#include <boost/archive/tmpdir.hpp>
//#include <boost/archive/xml_oarchive.hpp>
//#include "Core/DefSerialization.h"
using namespace std;

 
int main(int argc, char *argv[]) {


	double start,end;





	//H=re*nu/Umax-1;
	//H=round(H);
	//L=H;
	//re=(H+1)*Umax/nu;
	//re=re/100;
	//Pmax=1.0004;
	//Umax=0.005;
	//re=(H+1)*Umax/nu;

	//Umax=0.002;//re*nu/(H+1);
	//********Knudsen formulation*****
/*	double Kn=0.01; //Knudsen number
	double pi=atan(1)*4 ;
	double dt=1.0/(L)/sqrt(3.0);
	double omega = dt /(Kn*sqrt(2.0/pi)+0.50*dt);
	tau=1.0/omega;
	nu=(2.0*tau-1.0)/6.0 ;
//	nu=Kn*sqrt(2.0/pi);
	nu=0.01667;
	tau=(6.0*nu+1.0)*0.5;
	re=Umax*(H/2+1)/nu;
	//Umax=re*nu/(H+1);
	Umax=0.003;*/

// fluid 2 is the continuous fluid	and fluid 1 the droplet
//Domain size
	double H=50;
	double L=100;
/*
// Global parameter
	double Ca=0.2;
	double Re=0.01;
	double diameter=50;
	double R=diameter/2;

//Fluid viscosity
	double Nu_2=0.8;
	double viscosity_ratio=1;

	double Nu_1=Nu_2*viscosity_ratio;
	double tau1=(6.0*Nu_1+1.0)*0.5;
	double tau2=(6.0*Nu_2+1.0)*0.5;

	double shear_rate=Re*Nu_2/(R*R);

//Reference values
	double Rho1_ref=1.0;
	double Rho2_ref=1.0;
	double U1_ref=0.00001;
	double U2_ref=shear_rate*H/2;


//Surface tension
	double sigma=(Re/Ca)*Nu_2*Nu_2*Rho2_ref/R;

// Pressure inlet
	double Pmax=1;
// Pressure outlet
	double Pmin=1;



	cout<<"Reynolds: "<<Re<<" Ca is: "<<Ca<<" Surface tension is: "<<sigma<< " Tau 1: "<<tau1<< " Nu 1: "<<Nu_1<<" U 1: "<<U2_ref<<" Number of cell in x direction: "<<L<<" Number of cell in y direction: "<<H<<endl;
//	if(std::abs(Re-)>1e-3)
		std::cout<<"Reynolds imposed: "<<Re<<" Reynolds calculated: "<< shear_rate*R*R/(Nu_2)<<std::endl;

		std::cout<<"Capilary number imposed: "<<Ca<<" Capilary number calculated: "<< shear_rate*R*Nu_2/(sigma*Rho2_ref)<<std::endl;
*/
/// Create object for the simulation.
	Simulation simu;
	Parameters Param;

/// Set Simulation Arguments
	Param.Set_Arguments(&argc,&argv,false);

/// Set Solver type
// Dimension
	Param.Set_Dimension(D2);
// Scheme
	Param.Set_Scheme(Q9);

// Set Domain size
	Param.Set_Domain_Size((int)L,(int)H); //Cells

// Set Reynolds
	double re=20;
	double u=0.0505051;
	double nu=u*H/re;
	double tau1=(6.0*nu+1.0)*0.5;
	double Rho1_ref=1.0;

// Set User Parameters
	Param.Set_UserParameters(u,H,L,1.0,1.0);
//	Param.Set_UserParameters(U2_ref,H,L,Pmax,Pmin);
//	Param.Set_TwoPhaseUserParameters(Re,Ca,diameter, sigma);
// Set User Force type
	Param.Set_UserForceType(None);
// Set delta x for output
	Param.Set_deltax(1);


/// Set Boundary condition type for the boundaries of the domain
/// Boundary condition accepted: Wall, Pressure, Velocity and Symmetry
/// Order Bottom, Right, Top, Left, (Front, Back for 3D)
//	Param.Set_BcType(Velocity,Periodic,Velocity,Periodic);
//	Param.Set_BcType(Periodic,Periodic,Periodic,Periodic);
	Param.Set_BcType(Wall,Pressure,Wall,Velocity);

//	Param.Set_BcType(Wall,Pressure,Wall,Pressure);
/// Set Pressure Type
	Param.Set_PressureType(FixP);
/// Set Global Corner type
	Param.Set_CornerPressureType(FixCP);
/// Wall boundary condition type (Implemented BounceBack and Diffuse)
	Param.Set_WallType(BounceBack);

/// Number of maximum timestep
	Param.Set_NbStep(50000);
/// Interval for output
	Param.Set_OutPutNSteps(5000);// interval
///Display information during the calculation every N iteration
	Param.Set_listing(500);

///Selection of variables to export
	Param.Set_VariablesOutput(true,true);// export Rho,U

/// Define the Output filename
	Param.Set_OutputFileName("Poiseuille");

	// Multiphase model (SinglePhase or ColourFluid)
	Param.Set_Model(SinglePhase);

	//Gradient definition
	Param.Set_GradientType(LBMStencil); //FD or LBMStencil

/// Singlephase Parameters
	Param.Set_Tau(tau1);
	Param.Set_Rho(Rho1_ref);
/*
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
	Param.Set_ColourOperatorType(SurfaceForce);//Grunau or SurfaceForce

*/

/// Initialise the simulation with parameters
	simu.InitSimu(Param, true);
	simu.barrier();

/// Run the simulation
	simu.RunSimu();
}

