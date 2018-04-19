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
//Domain size
	double L=243;//118;//57;//431;
	double H=73;//37;//18;//131;
// Set Domain size
	Param.Set_Domain_Size((int)L,(int)H); //Cells
//Set Lattice unity
	double deltaTLattice=1;
	double deltaXLattice=1;
// Set lattice size and lattice time step (default value is 1 for both)
	 Param.Set_deltaT(deltaTLattice);Param.Set_deltaX(deltaXLattice);

	double Diameter=H;
//	double Ca=0;

	double U2_ref=0.01;


	double Rho1_ref=1;
	double Rho2_ref=Rho1_ref;

	double nu_1=1.0/(6.0*6.0);
	double lambda=6.0;///5;
	double nu_2=nu_1*lambda;

	double tau1=0.5+nu_1*deltaTLattice*3.0/(deltaXLattice*deltaXLattice);
	double tau2=0.5+nu_2*deltaTLattice*3.0/(deltaXLattice*deltaXLattice);
	double La=0;
	double sigma=0.0001;//La*(Rho1_ref*nu_1)*2/Diameter;//0.001;//0.01;//0.001;


//	double Mach =U2_ref*std::sqrt(3);
//	double Kn=(Mach/Re)*std::sqrt(pi/2.0);
//Pressure drop
	double deltaP=0.0001;
// Pressure inlet
	double Pmax=1+deltaP;
// Pressure outlet
	double Pmin=1-deltaP;



/// Set User Parameters

//Contact angle parameters
	double contactangle=90*pi/180.0;//30
	Param.Set_ContactAngleType(ContactAngleEnum::UniformTeta);//NoTeta, UniformTeta or UserTeta
	Param.Set_ContactAngleModel(ContactAngleEnum::Interpol);//Standard or Interpol
	Param.Set_SwitchSelectTeta(ContactAngleEnum::Linear);//Binary or Linear
	Param.Set_NormalExtrapolType(ContactAngleEnum::WeightDistanceExtrapol);//NoExtrapol,TailorExtrapol,or WeightDistanceExtrapol
	Param.Set_NormalInterpolType(ContactAngleEnum::LinearLeastSquareInterpol);//NoInterpol,LinearInterpol,LinearLeastSquareInterpol
	Param.Set_NumberOfInterpolNodeInSolid(1);
	Param.Set_NumberOfInterpolNodeInFluid(3);
	if(Param.Get_ContactAngleType()==ContactAngleEnum::NoTeta)
		contactangle=pi/2.0;
	Param.Set_ContactAngle(contactangle);

//	Param.Set_UserDroplets(contactangle,sigma,Diameter,Re,Ca);
	Param.Set_UserParameters(U2_ref,H,L,Pmax,Pmin);

// Set User Force type
	Param.Set_UserForceType(ModelEnum::None);//None,LocalForce or BodyForce
// Set delta x for conversion from lattice unit to physical unit
	Param.Set_deltax(1.0);//1774*1e-6/L


/// Set Boundary condition type for the boundaries of the domain
/// Boundary condition accepted: Wall, Pressure, Velocity, Periodic and Symmetry
/// Order Bottom, Right, Top, Left, (Front, Back for 3D)
	Param.Set_BcType(Symmetry,Pressure,Symmetry,Velocity);

	Param.Set_ModelOfFluid(SolverEnum::Compressible);
/// Set Pressure Type
	Param.Set_PressureType(FixP);//FixP,zeroPGrad1st
/// Set Global Corner type
	Param.Set_CornerPressureType(ExtrapolCP);//FixCP,ExtrapolCP
/// Wall boundary condition type (Implemented BounceBack and Diffuse)
	Param.Set_WallType(HalfWayBounceBack);//BounceBack,HeZouWall
	Param.Set_CollisionOnWalls(true);
	Param.Set_VelocityModel(HeZouV);
/// Set Periodic boundary condition to add a pressure drop term
	Param.Set_PeriodicType(Simple);//Simple,PressureForce
	Param.Set_PressureDrop(deltaP);
/// Number of maximum timestep
	Param.Set_NbStep(5000000);
/// Interval for output
	Param.Set_OutPutNSteps(40000);// interval
///Display information during the calculation every N iteration
	Param.Set_listing(20000);
	Param.Set_ErrorMax(5e-18*H*L);
	Param.Set_ErrorVariable(SolverEnum::RhoN);
	Param.SetErrorLpCase(false);
///Selection of variables to export
	Param.Set_VariablesOutput(true,true,false);// export Rho,U,


// Multiphase model (SinglePhase or ColourFluid)
	Param.Set_Model(SolverEnum::ColourFluid);
	if(tau1==tau2 || Param.Get_Model()==SolverEnum::SinglePhase)
		Param.Set_ViscosityType(ConstViscosity);//ConstViscosity,HarmonicViscosity
	else
		Param.Set_ViscosityType(HarmonicViscosity);//ConstViscosity,HarmonicViscosity
	if(Param.Get_ViscosityType()==ConstViscosity)
		lambda=1;

//Gradient definition
	Param.Set_GradientType(ModelEnum::LBMStencil); //FD or LBMStencil
	Param.Set_ExtrapolationType(ModelEnum::WeightDistanceExtrapol);//NoExtrapol,TailorExtrapol,WeightDistanceExtrapol

/// Singlephase Parameters
	Param.Set_Tau(tau1);
	Param.Set_Rho(Rho1_ref);


/// Multiphase Parameters
	// outputs
	//Colour gradient
	Param.Set_ColourGradientOutput(true);
	//Density of Red fluid
	Param.Set_BlueDensityOutput(false);
	//Density of blue fluid
	Param.Set_RedDensityOutput(false);
	//Normal of the interface
	Param.Set_NormalOutput(true);
	//RhoN
	Param.Set_NormalDensityOutput(true);

	//Density of each fluid
	Param.Set_Rho_1(Rho1_ref);
	Param.Set_Rho_2(Rho2_ref);
	//Relaxation time for each fluid
	Param.Set_Tau_1(tau1);
	Param.Set_Tau_2(tau2);
	//Surface tension
	Param.Set_SurfaceTension(sigma);
	Param.Set_ATau(Param.Convert_SigmaToATau());
	//Colour fluid Parameters
	Param.Set_ColourGradLimiter(0.00000001);
	Param.Set_Beta(0.95);// Between 0 and 1
	Param.Set_ColourGradType(ColourFluidEnum::DensityNormalGrad);//Gunstensen or DensityGrad or DensityNormalGrad
	Param.Set_RecolouringType(ColourFluidEnum::LatvaKokkoRothmanEstimator);//LatvaKokkoRothman or LatvaKokkoRothmanEstimator
	Param.Set_ColourOperatorType(ColourFluidEnum::SurfaceForce);//Grunau or Reis or SurfaceForce

// Porous media
	Param.PorousMediaCase(false);
	Param.CalculatePorosity(true);
	Param.CalculateProductionRate(true);
	Param.CalculatePermeability(true);
	Param.CalculateDarcyPermeability(true);
	Param.Set_SectionMedia(H);
	Param.Set_LenghtMedia(L);
	Param.Set_PositionMedia(0,L,0,H);
	Param.Set_HeleShawBodyForce(PorousMediaEnum::no);
	Param.Set_Depth(24.54*1e-06/Param.Get_deltax());

/// Define the Output filename
	FileExportStream.str("");
//	FileExportStream<<"Serpentine_Two_Phase_90Degree_Init_SPhase"<< fixed << setprecision(0)<<H<<"x"<<L;
	FileExportStream<<"Serpentine_Two_Phase_90Degree"<< fixed << setprecision(0)<<H<<"x"<<L;
	Param.Set_OutputFileName(FileExportStream.str());
/*
	Param.Add_VariableToInit("../../../SinglePhase/run/Serpentine_Single_Phase_73x243_15000.cgns",SolverEnum::Density);
	Param.Add_VariableToInit("../../../SinglePhase/run/Serpentine_Single_Phase_73x243_15000.cgns",SolverEnum::VelocityX);
	Param.Add_VariableToInit("../../../SinglePhase/run/Serpentine_Single_Phase_73x243_15000.cgns",SolverEnum::VelocityY);
*/

/// Initialise the simulation with parameters

	simu.InitSimu(Param, true);
	if(simu.Is_MainProcessor())
		simu.Save_Parameters(Param,"iniParamInclined.xml");
	double Re=1;

	if(simu.Is_MainProcessor())
		cout<<"Reynolds: "<<Re<<" Surface tension is: "<<sigma<< " Tau 1: "<<tau1<< " Tau 2: "<<tau2<<" U max: "<<U2_ref<<" Number of cell in x direction: "<<L<<" Number of cell in y direction: "<<H<<std::endl;


	simu.barrier();

	if(simu.Is_MainProcessor())
		cout<<"Tau 1: "<<tau1<<" Number of cell in x direction: "<<L<<" Number of cell in y direction: "<<H<<std::endl;
//Setup for simulation loop	
	double HeighChannel=19;
	double mySigma[] = {1e-02,1e-03,1e-04,1e-05};

	double Ca=0;
	std::vector<double> vectloop (mySigma, mySigma + sizeof(mySigma) / sizeof(double) );
//Check first value
	sigma=vectloop[0];
	U2_ref=Re*nu_1/HeighChannel;
	Ca=U2_ref*nu_1/sigma;
	if(simu.Is_MainProcessor())
		cout<<"U: "<<U2_ref<<" Ca: "<<Ca<<endl;
	Param.Set_UserParameters(U2_ref,H,L,Pmax,Pmin);
	Param.Set_UserDroplets(contactangle,sigma,Diameter,Re,Ca);
	// reinit the field
	 //  	simu.UpdateAllDomain(Param);
	for(int i=0;i<vectloop.size();i++)
	{
		
	//Surface tension
		sigma=vectloop[i];
		Param.Set_SurfaceTension(sigma);
		Param.Set_ATau(Param.Convert_SigmaToATau());

		U2_ref=Re*nu_1/HeighChannel;
		Ca=U2_ref*nu_1/sigma;
		FileExportStream.str("");
		FileExportStream<<"Serpentine_Half-Way_"<< fixed << setprecision(0)<<H<<"x"<<L<<"_Hchannel_"<<HeighChannel
						<< scientific<<setprecision(2)<<"_Re_"<<Re<<"_Ca_"<<Ca<<"_Sigma_"<<sigma;
		Param.Set_OutputFileName(FileExportStream.str());
		Param.Set_UserParameters(U2_ref,H,L,Pmax,Pmin);
		Param.Set_UserDroplets(contactangle,sigma,Diameter,Re,Ca);
		if(simu.Is_MainProcessor())
			cout<<"U 1: "<<U2_ref<<" Ca: "<<Ca<<
			" Number of cell in x direction: "<<L<<" Number of cell in y direction: "<<H<<std::endl;

// reinit the field
   		simu.UpdateAllDomain(Param);
   	//	simu.UpdateDomainBc(Param);
/// Run the simulation
		simu.RunSimu(Param);
	}
	simu.FinalizeSimu();

	if(simu.Is_MainProcessor())
		cout<<" END Simulation"<<endl;

}

