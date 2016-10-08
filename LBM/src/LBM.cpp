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
//#include "mpi.h"
#include <iostream>
#include <stdlib.h>
//#include "Mesh/SingleBlock/Block2D.h"
//#include "Parrallelism/MpiManager.h"
//#include "Core/Parameters.h"
#include "Core/Simulation.h"
//#include <boost/archive/tmpdir.hpp>
//#include <boost/archive/xml_oarchive.hpp>
//#include "Core/DefSerialization.h"
using namespace std;

 
int main(int argc, char *argv[]) {

	int idx=0;
	if (argv[1]==0)
	{
		std::cout<<"No argument. Index set to 1"<<std::endl;
		idx=1;
	}
	else
	{
		idx=atoi(argv[1]);
	}

	double H_vec[]={400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,600,600,600,600,600,600,600,600,600,600,600,600};
	double Lambda_vec[]={1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,11.5,11.5,11.5,11.5,11.5,11.5,11.5,11.5,11.5,11.5,11.5,11.5};
	double Confinement_vec[]={0.74,0.74,0.74,0.74,0.59,0.59,0.59,0.59,0.34,0.34,0.34,0.34,0.15,0.15,0.15,0.15,0.74,0.74,0.74,0.74,0.59,0.59,0.59,0.59,0.34,0.34,0.34,0.34,0.15,0.15,0.15,0.15,0.74,0.74,0.74,0.74,0.59,0.59,0.59,0.59,0.34,0.34,0.34,0.34};
	double Ca_vec[]={0.4,0.43,0.47,0.5,0.4,0.43,0.47,0.5,0.4,0.43,0.47,0.5,0.4,0.43,0.47,0.5,0.4,0.47,0.53,0.6,0.5,0.8,1.1,1.5,2,2.6,3.4,4,12,13.5,15,16,0.4,0.47,0.53,0.6,30,33,36,40,70,73,76,80};
	double Rho_2_vec[]={1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
	double Sigma_vec[]={0.0003,0.0003,0.0003,0.0003,0.0003,0.0003,0.0003,0.0003,0.001,0.001,0.001,0.001,0.002,0.002,0.002,0.002,0.0003,0.0003,0.0003,0.0003,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.00005,0.00005,0.00005,0.00005,0.0005,0.0005,0.0005,0.0005,0.00001,0.00001,0.00001,0.00001,0.00001,0.00001,0.00001,0.00001};
	double Mu_2_vec[]={0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3};

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
	double pi=atan(1)*4 ;
stringstream FileExportStream;
// fluid 2 is the continuous fluid	and fluid 1 the droplet
//Domain size
	double H=150;
	double L=3*H;

// Global parameter
	double Ca=0.2;
	double Re=0.1;
	//Surface tension
	double sigma=0.001;//(Re/Ca)*Nu_2*Nu_2*Rho2_ref/R;
	double confinement=0.5;

	double diameter=confinement*H;
	double R=diameter/2;

	//Reference values
	double Rho2_ref=1.0;
	double Rho1_ref=Rho2_ref-3.0*sigma/R;

//Fluid viscosity
	double Nu_2=std::sqrt(Ca*sigma*R*Rho2_ref/Re);
	double viscosity_ratio=1;

	double Nu_1=Nu_2*viscosity_ratio;
	double tau1=(6.0*Nu_1+1.0)*0.5;
	double tau2=(6.0*Nu_2+1.0)*0.5;

	double shear_rate=Re*Nu_2/(R*R);



	double U1_ref=shear_rate*R;
	double U2_ref=shear_rate*H/2;

/*

	U2_ref=0.01;
	shear_rate=2*U2_ref/H;
	U1_ref=shear_rate*R;
	Nu_2=shear_rate*R*R/Re;
	Rho2_ref=1.0;
	sigma=Nu_2*shear_rate*R/(Ca*Rho2_ref);
	Rho1_ref=Rho2_ref-3.0*sigma/R;
	Nu_1=Nu_2*viscosity_ratio;
	tau1=(6.0*Nu_1+1.0)*0.5;
	tau2=(6.0*Nu_2+1.0)*0.5;*/

	//Domain size
		H=H_vec[idx];
		L=3*H;

	// Global parameter
		Ca=Ca_vec[idx];
		//Surface tension
		sigma=Sigma_vec[idx];//(Re/Ca)*Nu_2*Nu_2*Rho2_ref/R;
		confinement=Confinement_vec[idx];
		viscosity_ratio=Lambda_vec[idx];

		diameter=confinement*H;
		R=diameter/2;

		//Reference values
		Rho2_ref=Rho_2_vec[idx];
		double mu_2=Mu_2_vec[idx];
		Rho1_ref=Rho2_ref;//-3.0*sigma/R;
		double mu_1=mu_2*viscosity_ratio;

	//Fluid viscosity
		Nu_2=mu_2*Rho2_ref;
		Nu_1=mu_1*Rho1_ref;

		tau1=(6.0*Nu_1+1.0)*0.5;
		tau2=(6.0*Nu_2+1.0)*0.5;

		shear_rate=Ca*sigma/(mu_2*R);
		U1_ref=shear_rate*R;
		U2_ref=shear_rate*H/2;

		Re=shear_rate*R*R/Nu_2;

	double Mach =U2_ref*std::sqrt(3);
	double Kn=(Mach/Re)*std::sqrt(pi/2.0);
// Pressure inlet
	double Pmax=1;
// Pressure outlet
	double Pmin=1;



//	cout<<"Reynolds: "<<Re<<" Ca is: "<<Ca<<" Surface tension is: "<<sigma<< " Tau 2: "<<tau2<< " Nu 2: "<<Nu_2<<" U 1: "<<U2_ref<<" Number of cell in x direction: "<<L<<" Number of cell in y direction: "<<H<<std::endl;


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
/*	double re=20;
	double u=0.0505051;
	double nu=u*H/re;
	double tau1=(6.0*nu+1.0)*0.5;
	double Rho1_ref=1.0;

// Set User Parameters
	Param.Set_UserParameters(u,H,L,1.0,1.0);*/
	Param.Set_UserParameters(U2_ref,H,L,Pmax,Pmin);
	Param.Set_TwoPhaseUserParameters(Re,Ca,diameter, sigma);
// Set User Force type
	Param.Set_UserForceType(None);
// Set delta x for output
	Param.Set_deltax(1/H);


/// Set Boundary condition type for the boundaries of the domain
/// Boundary condition accepted: Wall, Pressure, Velocity and Symmetry
/// Order Bottom, Right, Top, Left, (Front, Back for 3D)
	Param.Set_BcType(Velocity,Periodic,Velocity,Periodic);
//	Param.Set_BcType(Periodic,Periodic,Periodic,Periodic);
//	Param.Set_BcType(Wall,Wall,Wall,Wall);
//	Param.Set_BcType(Wall,Pressure,Wall,Pressure);

/// Set Pressure Type
	Param.Set_PressureType(FixP);
/// Set Global Corner type
	Param.Set_CornerPressureType(FixCP);
/// Wall boundary condition type (Implemented BounceBack and Diffuse)
	Param.Set_WallType(BounceBack);
	Param.Set_VelocityModel(Ladd);

/// Number of maximum timestep
	Param.Set_NbStep(500000);
/// Interval for output
	Param.Set_OutPutNSteps(20000);// interval
///Display information during the calculation every N iteration
	Param.Set_listing(2000);
	Param.Set_ErrorMax(1e-10);

///Selection of variables to export
	Param.Set_VariablesOutput(true,true);// export Rho,U

/// Define the Output filename
/*	FileExportStream<<"Droplet_shear_"<< fixed << setprecision(0)<<H<<"x"<<L
			<<setprecision(2)<<"Re_"<<Re<<"_Ca_"<<Ca<<"_Conf_"<<confinement<<"_lambda_"<<viscosity_ratio
			<< scientific<<setprecision(3)<<"sigma_"<<sigma;*/
	FileExportStream.str("MeshIni");
	Param.Set_OutputFileName(FileExportStream.str());

	// Multiphase model (SinglePhase or ColourFluid)
	Param.Set_Model(ColourFluid);
	Param.Set_ViscosityType(ConstViscosity);

	//Gradient definition
	Param.Set_GradientType(LBMStencil); //FD or LBMStencil

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
	Param.Set_ColourOperatorType(SurfaceForce);//Grunau or SurfaceForce



/// Initialise the simulation with parameters
	simu.InitSimu(Param, true);
	simu.barrier();

	//double myUref[] = {0.001,0.0002,0.0004,0.0006,0.0008,0.001,0.002,0.004,0.006,0.008,0.01,0.02,0.04};//{0.00001,0.0001,0.001,0.01};
	//double myUref[] = {0.001,0.002,0.003,0.004,0.005,0.006,0.007,0.008,0.009,0.01};
	//double myUref[] = {0.008,0.009,0.01};
/*	double myUref[] = {0.3,0.25,0.2,0.15,0.1,0.05,0.001};
	std::vector<double> vectloop (myUref, myUref + sizeof(myUref) / sizeof(double) );


	for(int i=0;i<vectloop.size();i++)
	{
		Ca=vectloop[i];
		sigma=0.001;
		Nu_2=std::sqrt(R*sigma*Ca/(Rho2_ref*Re));
		shear_rate=Re*Nu_2/(R*R);
		U2_ref=shear_rate*H/2;
		Rho1_ref=Rho2_ref;//-3.0*sigma/R;
		Nu_1=Nu_2*viscosity_ratio;
		tau1=(6.0*Nu_1+1.0)*0.5;
		tau2=(6.0*Nu_2+1.0)*0.5;*/
		FileExportStream.str("");
		FileExportStream<<"Droplet_shear_"<< fixed << setprecision(0)<<H<<"x"<<L
				<< scientific<<setprecision(2)<<"Re_"<<Re<<"_Ca_"<<Ca<<"_Conf_"<<confinement<<"_lambda_"<<viscosity_ratio
			<< scientific<<setprecision(3)<<"sigma_"<<sigma;
		Param.Set_OutputFileName(FileExportStream.str());
		Param.Set_UserParameters(U2_ref,H,L,Pmax,Pmin);
		Param.Set_TwoPhaseUserParameters(Re,Ca,diameter, sigma);
		Param.Set_SurfaceTension(sigma);
		Param.Set_Rho_1(Rho1_ref);
		Param.Set_Rho_1(Rho2_ref);
		Param.Set_Tau_1(tau1);
		Param.Set_Tau_2(tau2);

		if(simu.Is_MainProcessor())
			cout<<"Reynolds: "<<Re<<" Ca is: "<<Ca<<" Surface tension is: "<<sigma<< " Tau 2: "<<tau2<< " Nu 2: "<<Nu_2<<" U 1: "<<U2_ref<<" Number of cell in x direction: "<<L<<" Number of cell in y direction: "<<H<<std::endl;



//	simu.UpdateAllDomain(Param);
/// Run the simulation
	simu.RunSimu(Param);

//	}
}

