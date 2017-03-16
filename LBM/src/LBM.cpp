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

/// Create object for the simulation.
	Simulation simu;
	Parameters Param;
	double pi=atan(1.0)*4.0 ;
	stringstream FileExportStream;
/// Set Simulation Arguments
	Param.Set_Arguments(&argc,&argv,false);

/// Set Solver type
// Dimension
	Param.Set_Dimension(SolverEnum::D2);
// Scheme
	Param.Set_Scheme(SolverEnum::Q9);



// fluid 2 is the continuous fluid	and fluid 1 the droplet
// Set Domain size
	double L=600;
	double H=195;
	Param.Set_Domain_Size((int)L,(int)H); //Cells
//Set Lattice unity
	double deltaTLattice=1;
	double deltaXLattice=1;
// Set lattice size and lattice time step (default value is 1 for both)
	 Param.Set_deltaT(deltaTLattice);Param.Set_deltaX(deltaXLattice);

	 double Diameter=H;
	double Ca=0;

	double U2_ref=0.01;
	double Re=1;

	double Rho1_ref=1;
	double Rho2_ref=Rho1_ref;


	double lambda=1.0/43.0;
	double nu_1=4.91*1.e-2/lambda;
	double nu_2=lambda*nu_1;

	double tau1=0.5+nu_1*deltaTLattice*3.0/(deltaXLattice*deltaXLattice);
	double tau2=0.5+nu_2*deltaTLattice*3.0/(deltaXLattice*deltaXLattice);
	double La=0;
	double sigma=1.96*1.e-3;//La*(Rho1_ref*nu_1)*2/Diameter;//0.001;//0.01;//0.001;


	double Mach =U2_ref*std::sqrt(3);
	double Kn=(Mach/Re)*std::sqrt(pi/2.0);
//Pressure drop
	double deltaP=0.0;
	double refPressure=1.0/3.0;
// Pressure inlet
	double Pmax=refPressure+deltaP;
// Pressure outlet
	double Pmin=refPressure-deltaP;






// Set User Parameters
	//U2_ref=0.01;Pmax=1;Pmin=1;
//Contact angle parameters
	double contactangle=30*pi/180.0;
	Param.Set_ContactAngleType(ContactAngleEnum::FixTeta);//NoTeta, FixTeta or NonCstTeta
	Param.Set_ContactAngleModel(ContactAngleEnum::Interpol);//Standard or Interpol
	Param.Set_SwitchSelectTeta(ContactAngleEnum::Linear);//Binary or Linear
	Param.Set_NormalExtrapolType(ContactAngleEnum::WeightDistanceExtrapol);//NoExtrapol,TailorExtrapol,or WeightDistanceExtrapol
	Param.Set_NormalInterpolType(ContactAngleEnum::LinearLeastSquareInterpol);//NoInterpol,LinearInterpol,LinearLeastSquareInterpol
	Param.Set_NumberOfInterpolNodeInSolid(3);
	Param.Set_NumberOfInterpolNodeInFluid(3);
	if(Param.Get_ContactAngleType()==ContactAngleEnum::NoTeta)
		contactangle=pi/2.0;
	Param.Set_ContactAngle(contactangle);

	Param.Set_UserDroplets(contactangle,sigma,Diameter,Re,Ca);
	Param.Set_UserParameters(U2_ref,H,L,Pmax,Pmin);

// Set User Force type
	Param.Set_UserForceType(None);
// Set delta x for output
	Param.Set_deltax(1);


/// Set Boundary condition type for the boundaries of the domain
/// Boundary condition accepted: Wall, Pressure, Velocity, Periodic and Symmetry
/// Order Bottom, Right, Top, Left, (Front, Back for 3D)
	Param.Set_BcType(Symmetry,Pressure,Symmetry,Velocity);


/// Set Pressure Type
	Param.Set_PressureType(FixP);//FixP,zeroPGrad1st
/// Set Global Corner type
	Param.Set_CornerPressureType(ExtrapolCP);//FixCP,ExtrapolCP
/// Wall boundary condition type (Implemented BounceBack and Diffuse)
	Param.Set_WallType(BounceBack);//BounceBack,HeZouWall
	Param.Set_VelocityModel(HeZouV);
/// Set Periodic boundary condition to add a pressure drop term
	Param.Set_PeriodicType(Simple);//Simple,PressureForce
	Param.Set_PressureDrop(deltaP);
/// Number of maximum timestep
	Param.Set_NbStep(10000);
/// Interval for output
	Param.Set_OutPutNSteps(100);// interval
///Display information during the calculation every N iteration
	Param.Set_listing(100);
	Param.Set_ErrorMax(1e-10);
	Param.Set_ErrorVariable(SolverEnum::VelocityX);

///Selection of variables to export
	Param.Set_VariablesOutput(true,true);// export Rho,U

	// Multiphase model (SinglePhase or ColourFluid)
	Param.Set_Model(SolverEnum::SinglePhase);
	Param.Set_ViscosityType(HarmonicViscosity);//ConstViscosity,HarmonicViscosity
	if(Param.Get_ViscosityType()==ConstViscosity)
		lambda=1;

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
	Param.Set_ColourGradLimiter(0.001);
	Param.Set_A1(0.00001);
	Param.Set_A2(Param.Get_A1());
	Param.Set_Beta(0.7);// Between 0 and 1
	Param.Set_ColourGradType(ColourFluidEnum::DensityNormalGrad);//Gunstensen or DensityGrad or DensityNormalGrad
	Param.Set_RecolouringType(ColourFluidEnum::LatvaKokkoRothman);
	Param.Set_ColourOperatorType(ColourFluidEnum::SurfaceForce);//Grunau or Reis or SurfaceForce

	/// Define the Output filename
	/*	FileExportStream<<"Droplet_shear_"<< fixed << setprecision(0)<<H<<"x"<<L
				<<setprecision(2)<<"Re_"<<Re<<"_Ca_"<<Ca<<"_Conf_"<<confinement<<"_lambda_"<<viscosity_ratio
				<< scientific<<setprecision(3)<<"sigma_"<<sigma;*/
		FileExportStream.str("");
/*		FileExportStream<<"5ContractionsSymmetry_"<< fixed << setprecision(0)<<L<<"x"<<H;
//		FileExportStream<<"LinearLeastSquareInterpol_linearswitch_teta170_beta0.7_5fluids_3Solid";
		if(Param.Get_ContactAngleModel()==ContactAngleEnum::Interpol)
		{
			if(Param.Get_NormalInterpolType()==ContactAngleEnum::LinearLeastSquareInterpol)
			{
				FileExportStream<<"LinearLeastSquareInterpol_"
						<<Param.Get_NumberOfInterpolNodeInFluid()<<"-Fluids_"
						<<Param.Get_NumberOfInterpolNodeInSolid()<<"-Solids_";
			}
			if(Param.Get_NormalInterpolType()==ContactAngleEnum::LinearInterpol)
					FileExportStream<<"LinearInterpol_";
		}
		else
			FileExportStream<<"Standard_";
		if(Param.Get_SwitchSelectTeta()==ContactAngleEnum::Linear)
			FileExportStream<<"LinearSwitch_";
		if(Param.Get_SwitchSelectTeta()==ContactAngleEnum::Binary)
			FileExportStream<<"BinarySwitch_";
		FileExportStream<<setprecision(2)<<lambda<<"-lambda_"
				<<Param.Get_Beta()<<"-beta_"<<Param.Get_ContactAngle()*180.0/pi<<"-teta_"
						<< scientific<<setprecision(3)<<"sigma_"<<sigma<<"_ColourGradLimiter"<<Param.Get_ColourGradLimiter();
*/
		FileExportStream.str("Contraction1162_test");
		Param.Set_OutputFileName(FileExportStream.str());

//	Param.Set_OutputFileName("Testread");

//	Param.Add_VariableToInit("Contraction1162_50000.cgns",SolverEnum::Density);
//	Param.Add_VariableToInit("Contraction1162_50000.cgns",SolverEnum::VelocityX);
//	Param.Add_VariableToInit("Contraction1162_50000.cgns",SolverEnum::VelocityY);
//	Param.Add_VariableToInit("Contraction1162_testdeltaT_50000.cgns",SolverEnum::RhoN);

/// Initialise the simulation with parameters
	simu.InitSimu(Param, true);
/*	double *d_=0;
	std::string variablename("Density");
	std::string filename("Debug_poiseuille_serpentine_105000.cgns");
	simu.Test_ReadData(d_,variablename,filename);
	if(d_==0)
		std::cout<<"Error reading data"<<std::endl;
	else
		if(simu.Is_MainProcessor())
		for (int i=0;i<20;i++)
			std::cout<<i<<" "<<d_[i]<<std::endl;
*/
	simu.barrier();

		if(simu.Is_MainProcessor())
			cout<<" Surface tension is: "<<sigma<< " Tau 1: "<<tau1<<" Tau 2: "<<tau2<<" Number of cell in x direction: "<<L<<" Number of cell in y direction: "<<H<<std::endl;

//	simu.UpdateAllDomain(Param);
/// Run the simulation
	simu.RunSimu(Param);

//	}
}

