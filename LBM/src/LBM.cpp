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

	double Umax=0.132013201320132;//0.01;
	double H=299;
	double L=299;
	double Pmax=1.00001;
	double Pmin=0.99999;
	//double Umaxtmp,Htmp,Ltmp,Pmaxtmp,Pmintmp;

	double tau=0.6; // Relaxation time
	double re=10; // Reynold Number
	double ma=0.01; // Mach number in lattice
	//double Kn=0.001; //Knudsen number
	double nu;
	nu=(2.0*tau-1.0)/6.0 ;
	Umax=ma/std::sqrt(3.0);
	//Umax=re*nu/(H+1);
	H=427;//599
	L=532;//820
	H=854;
	L=1064;
	H=100;
	L=100;
	//H=re*nu/Umax-1;
	//H=round(H);
	//L=H;
	//re=(H+1)*Umax/nu;
	//re=re/100;
	//Pmax=1.0004;
	//Umax=0.005;
	//re=(H+1)*Umax/nu;

	Umax=0.002;//re*nu/(H+1);
	//********Knudsen formulation*****
	double Kn=0.01; //Knudsen number
	double pi=atan(1)*4 ;
	double dt=1.0/(L)/sqrt(3.0);
	double omega = dt /(Kn*sqrt(2.0/pi)+0.50*dt);
	tau=1.0/omega;
	nu=(2.0*tau-1.0)/6.0 ;
//	nu=Kn*sqrt(2.0/pi);
//	tau=(6.0*nu+1.0)*0.5;
	re=Umax*(H/2+1)/nu;
	//Umax=re*nu/(H+1);
	cout<<"Reynolds: "<<re<< " Tau: "<<tau<< " Nu: "<<nu<<" U: "<<Umax<<" Number of cell in x direction: "<<L<<" Number of cell in y direction: "<<H<<endl;




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

// Set User Parameters
	Param.Set_UserParameters(Umax,H,L,Pmax,Pmin);
// Set delta x for output
	Param.Set_deltax(1/L);

// Relaxation time for single phase
	Param.Set_Tau(tau);

/// Set Boundary condition type for the boundaries of the domain
/// Boundary condition accepted: Wall, Pressure, Velocity and Symmetry
/// Order Bottom, Right, Top, Left, (Front, Back for 3D)
	Param.Set_BcType(Wall,Wall,Wall,Wall);
/// Wall boundary condition type (Implemented BounceBack and Diffuse)
	Param.Set_WallType(BounceBack);

/// Number of maximum timestep
	Param.Set_NbStep(10);
/// Interval for output
	Param.Set_OutPutNSteps(1);// interval
///Display information during the calculation every N iteration
	Param.Set_listing(500);

///Selection of variables to export
	Param.Set_VariablesOutput(true,true);// export Rho,U

/// Define the Output filename
	Param.Set_OutputFileName("Two_colour");

	// Multiphase model (SinglePhase or ColourFluid)
	Param.Set_Model(ColourFluid);


/// Singlephase Parameters
	Param.Set_Rho(1.0);
	Param.Set_Tau(tau);

/// Multiphase Parameters
	//Density of each fluid
	Param.Set_Rho_1(1.001);
	Param.Set_Rho_2(1);
	//Relaxation time for each fluid
	Param.Set_Tau_1(tau);
	Param.Set_Tau_2(tau);
	//Surface tension
	Param.Set_SurfaceTension(0.01);
		//Colour fluid Parameters
		Param.Set_A1(0.01);
		Param.Set_A2(Param.Get_A1());
		Param.Set_Beta(0.7);
		Param.Set_ColourGradType(Gunstensen);
		Param.Set_RecolouringType(LatvaKokkoRothman);
		Param.Set_ColourOperatorType(Grunau);


/// Initialise the simulation with parameters
	simu.InitSimu(Param, true);
	simu.barrier();
/// Run the simulation
	simu.RunSimu();

}

