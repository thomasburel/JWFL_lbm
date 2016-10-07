/*
 * Parameters.h
 *
 *  Created on: 5 May 2015
 *      Author: thomas
 */

#ifndef CORE_PARAMETERS_H_
#define CORE_PARAMETERS_H_

#include "../User/UserParameters.h"
#include "../Mesh/SingleBlock/Node2D.h"
#include <iomanip>
#include <iostream>
#include <fstream>
#include <cmath>
//#include <string>
#include <boost/serialization/string.hpp>
#include <boost/serialization/nvp.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/list.hpp>
#include <boost/serialization/version.hpp>
#include <boost/serialization/assume_abstract.hpp>
#include <boost/archive/tmpdir.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>


/*
 *
 */


using namespace std;

//General enumeration
enum parralleltype {Serial,Mpi,Openmp,Hybrid};
enum OutputFormat {CGNSFormat,TecplotFormat};
enum DensityExport{Density,NormalDensity,BothDensity};

//Solver enumeration
enum dimension {D2,D3};
enum schemetype {Q9,Q16};
enum modeltype {SinglePhase, ColourFluid};

//model enumeration
enum FluidType{Newtonian};
enum UserForceType{None,LocalForce,BodyForce};
enum GradientType{FD,LBMStencil};

//Boundary conditions enumeration
enum WallType {BounceBack, HalfWayBounceBack, Diffuse, Specular,HeZouWall,HeZouWallVel};
enum PressureModel{HeZouP};
enum PressureType{FixP,zeroPGrad1st};
enum VelocityModel{HeZouV,Ladd};
enum VelocityType{FixV,zeroVGrad1st};
enum CornerModel{HoChan};
enum CornerPressureType{FixCP,ExtrapolCP};

//Two phases enumeration
enum TetaType{NoTeta, FixTeta, NonCstTeta};
enum ViscosityType{ConstViscosity,LinearViscosity,HarmonicViscosity};
//Colour fluid enumeration
enum ColourGradType{Gunstensen,DensityGrad,DensityNormalGrad};
enum RecolouringType{LatvaKokkoRothman};
enum ColourOperatorType {Grunau,Reis,SurfaceForce};

class ColourFluid{
private:
	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive &ar, const unsigned int version)
	{
		ar & BOOST_SERIALIZATION_NVP(beta)
		 & BOOST_SERIALIZATION_NVP(A1)
		 & BOOST_SERIALIZATION_NVP(A2)
		 & BOOST_SERIALIZATION_NVP(ColourGrad)
		 & BOOST_SERIALIZATION_NVP(Recolouring)
		 & BOOST_SERIALIZATION_NVP(ColourOperator);
	}
public:
	void Set_Beta(double beta_=0.7){beta=beta_;};
	double Get_Beta() const{return beta;};
	void Set_A1(double A_=0){A2=A_;};
	double Get_A1() const{return A2;};
	void Set_A2(double A_=0){A1=A_;};
	double Get_A2() const{return A1;};
	void Set_ColourGradType(ColourGradType type){ColourGrad=type;};
	ColourGradType Get_ColourGradType(){return ColourGrad;};
	void Set_RecolouringType(RecolouringType type){Recolouring=type;};
	RecolouringType Get_RecolouringType(){return Recolouring;};
	void Set_ColourOperatorType(ColourOperatorType type){ColourOperator=type;};
	ColourOperatorType Get_ColourOperatorType(){return ColourOperator;};
protected:
	double beta;
	double A1,A2;
	ColourGradType ColourGrad;
	RecolouringType Recolouring;
	ColourOperatorType ColourOperator;

};

class MultiphaseModels: public ColourFluid{
private:
	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive &ar, const unsigned int version)
	{
		ar & BOOST_SERIALIZATION_NVP(tension)
		 & BOOST_SERIALIZATION_NVP(NormalDensityOutput)
		 & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ColourFluid);
	}
public:
	void Set_SurfaceTension(double tension_=0){tension=tension_;};
	double Get_SurfaceTension() const{return tension;};
	void Set_NormalDensityOutput(bool NormalDensityOutput_=0){NormalDensityOutput=NormalDensityOutput_;};
	bool Get_NormalDensityOutput() const{return NormalDensityOutput;};
protected:
	bool NormalDensityOutput;
	double tension;
};
class MeshParameters {
private:
	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive &ar, const unsigned int version)
	{
		ar & BOOST_SERIALIZATION_NVP(nx)
		   & BOOST_SERIALIZATION_NVP(ny)
		   & BOOST_SERIALIZATION_NVP(nz)
		   & BOOST_SERIALIZATION_NVP(dimension_)
		   & BOOST_SERIALIZATION_NVP(MeshFile)
		   & BOOST_SERIALIZATION_NVP(NbNodes);
	}
public:
	void Set_Dimension(dimension dim_);
	dimension Get_Dimension() const;
	void set_MeshParameters();
	void Set_Domain_Size(int nx_=1,int ny_=1,int nz_=1);
	int Get_Nx() const;
	int Get_Ny() const;
	int Get_Nz() const;
	int Get_NbNodes() const;
	void Set_MeshFile(std::string MeshFile_);
	std::string Get_MeshFile() const;
protected:
	int nx,ny,nz;
	dimension dimension_; //2D or 3D
	std::string MeshFile;
	int NbNodes;
};

class ModelsParameters:public MultiphaseModels {
public:
//Single phase
	void Set_Rho(double Rho_=1){Rho_1=Rho_;};
	double Get_Rho() const{return Rho_1;};
	void Set_Tau(double Tau_=1){Tau_1=Tau_;};
	double Get_Tau() const{return Tau_1;};
//multiphase
	void Set_Rho_1(double Rho){Rho_1=Rho;};
	double Get_Rho_1() const {return Rho_1;};
	void Set_Rho_2(double Rho){Rho_2=Rho;};
	double Get_Rho_2() const {return Rho_2;};
	void Set_Tau_1(double Tau){Tau_1=Tau;};
	double Get_Tau_1() const {return Tau_1;};
	void Set_Tau_2(double Tau){Tau_2=Tau;};
	double Get_Tau_2() const {return Tau_2;};
private:
	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive &ar, const unsigned int version)
	{
		ar & BOOST_SERIALIZATION_NVP(Rho_1)
		   & BOOST_SERIALIZATION_NVP(Rho_2)
		   & BOOST_SERIALIZATION_NVP(Tau_1)
		   & BOOST_SERIALIZATION_NVP(Tau_2)
		   & BOOST_SERIALIZATION_BASE_OBJECT_NVP(MultiphaseModels);
	}

protected:
	double Rho_1, Rho_2;
	double Tau_1, Tau_2;
};

class PhysicalParameters {
public:
	void Set_deltax(double dx_);
	double Get_deltax() const ;
	void Set_deltat(double dt_){dt=dt_;};
	double Get_deltat() const {return dt;};
private:
	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive &ar, const unsigned int version)
	{
		ar & BOOST_SERIALIZATION_NVP(nu0)
		   & BOOST_SERIALIZATION_NVP(l0)
		   & BOOST_SERIALIZATION_NVP(t0)
		   & BOOST_SERIALIZATION_NVP(Re)
		   & BOOST_SERIALIZATION_NVP(l)
		   & BOOST_SERIALIZATION_NVP(t)
		   & BOOST_SERIALIZATION_NVP(dx)
		   & BOOST_SERIALIZATION_NVP(dt);
	}

protected:
	double nu0, l0,t0;//physical reference values
	double Re,l,t;// non dimension reference values
	double dx,dt; // convert from lattice unit to non dimension value

};

class IniParameters {
private:
	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive &ar, const unsigned int version)
	{
		ar & BOOST_SERIALIZATION_NVP(U_ini)
		   & BOOST_SERIALIZATION_NVP(V_ini)
		   & BOOST_SERIALIZATION_NVP(W_ini);
	}

protected:
	double U_ini,V_ini,W_ini;

};

class BoundaryParameters {
private:
	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive &ar, const unsigned int version)
	{
		ar & BOOST_SERIALIZATION_NVP(WallTypeParam)
		   & BOOST_SERIALIZATION_NVP(NbGlobalBcType);
		if (NbGlobalBcType==4)
		{
			ar & boost::serialization::make_nvp("GlobalBcType_0",GlobalBcType[0] );
			ar & boost::serialization::make_nvp("GlobalBcType_1",GlobalBcType[1] );
			ar & boost::serialization::make_nvp("GlobalBcType_2",GlobalBcType[2] );
			ar & boost::serialization::make_nvp("GlobalBcType_3",GlobalBcType[3] );
		}
		else
		{
			ar & boost::serialization::make_nvp("GlobalBcType_0",GlobalBcType[0] );
			ar & boost::serialization::make_nvp("GlobalBcType_1",GlobalBcType[1] );
			ar & boost::serialization::make_nvp("GlobalBcType_2",GlobalBcType[2] );
			ar & boost::serialization::make_nvp("GlobalBcType_3",GlobalBcType[3] );
			ar & boost::serialization::make_nvp("GlobalBcType_4",GlobalBcType[4] );
			ar & boost::serialization::make_nvp("GlobalBcType_5",GlobalBcType[5] );
		}
		ar & BOOST_SERIALIZATION_NVP(SymmetryTypeParam)
		   & BOOST_SERIALIZATION_NVP(PressureModelParam)
		   & BOOST_SERIALIZATION_NVP(PressureTypeParam)
		   & BOOST_SERIALIZATION_NVP(VelocityModelParam)
		   & BOOST_SERIALIZATION_NVP(VelocityTypeParam)
		   & BOOST_SERIALIZATION_NVP(CornerModelParam)
		   & BOOST_SERIALIZATION_NVP(CornerPressureTypeParam);
	}
public:
	NodeType Get_GlobalBcType(int side) const;
	void Set_WallType(WallType WallType_);
	WallType Get_WallType() const;
	void Set_SymmetryType(SymmetryType SymmetryType_);
	SymmetryType Get_SymmetryType() const;
	void Set_PressureModel(PressureModel PressureModelParam_){PressureModelParam=PressureModelParam_;};
	PressureModel Get_PressureModel() const {return PressureModelParam;};
	void Set_PressureType(PressureType PressureType_){PressureTypeParam=PressureType_;};
	PressureType Get_PressureType() const {return PressureTypeParam;};
	void Set_VelocityModel(VelocityModel VelocityModel_){VelocityModelParam=VelocityModel_;};
	VelocityModel Get_VelocityModel() const {return VelocityModelParam;};
	void Set_VelocityType(VelocityType VelocityType_){VelocityTypeParam=VelocityType_;};
	VelocityType Get_VelocityType() const {return VelocityTypeParam;};
	void Set_CornerModel(CornerModel CornerModel_){CornerModelParam=CornerModel_;};
	CornerModel Get_CornerModel() const {return CornerModelParam;};
	void Set_CornerPressureType(CornerPressureType CornerPressureType_){CornerPressureTypeParam=CornerPressureType_;};
	CornerPressureType Get_CornerPressureType() const {return CornerPressureTypeParam;};
protected:
	WallType WallTypeParam;
	SymmetryType SymmetryTypeParam;
	int NbGlobalBcType;
	NodeType * GlobalBcType;
	PressureModel PressureModelParam;
	PressureType PressureTypeParam;
	VelocityModel VelocityModelParam;
	VelocityType VelocityTypeParam;
	CornerModel CornerModelParam;
	CornerPressureType CornerPressureTypeParam;
};

class SolverParameters {
private:
	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive &ar, const unsigned int version)
	{
		ar & BOOST_SERIALIZATION_NVP(scheme)
		   & BOOST_SERIALIZATION_NVP(model)
		   & BOOST_SERIALIZATION_NVP(fluid)
		   & BOOST_SERIALIZATION_NVP(Gradient)
		   & BOOST_SERIALIZATION_NVP(UserForce)
		   & BOOST_SERIALIZATION_NVP(NbVelocities)
		   & BOOST_SERIALIZATION_NVP(ErrorMax);


	}

public:
	void set_SolverParameters();
	void Set_Scheme(schemetype scheme_);
	schemetype Get_Scheme() const;
	void Set_Model(modeltype model_);
	modeltype Get_Model() const;
	void Set_FluidType(FluidType fluidtype_){fluid=fluidtype_;};
	FluidType Get_FluidType() const{return fluid;};
	void Set_UserForceType(UserForceType UserForceType_){UserForce=UserForceType_;};
	UserForceType Get_UserForceType() const{return UserForce;};
	void Set_ErrorMax(double ErrorMax_){ErrorMax=ErrorMax_;};
	double Get_ErrorMax(){return ErrorMax;};
	void Set_GradientType(GradientType GradientType_){Gradient=GradientType_;};
	GradientType Get_GradientType() const{return Gradient;};
	void Set_ContactAngle(double teta_){teta=teta_;};
	double Get_ContactAngle(){return teta;};
	void Set_ContactAngleType(TetaType tetaType_){tetaType=tetaType_;};
	TetaType Get_ContactAngleType(){return tetaType;};
	void Set_ViscosityType(ViscosityType viscosityType_){viscosityType=viscosityType_;};
	ViscosityType Get_ViscosityType(){return viscosityType;};
	int Get_NbVelocities() const;
	double Convert_TauToNu(double Tau_){return (Tau_-0.5)*cs2*deltaT;};
	double Convert_MutoNu(double mu_,double rho){return mu_/rho;};
	double Convert_NutoMu(double nu_,double rho){return nu_*rho;};
	double Convert_NuToTau(double nu_){return nu_/(cs2*deltaT)+0.5;};

	double Get_cs() const {return cs;};
	double Get_cs2() const {return cs2;};
	double Get_deltaT() const {return deltaT;};
protected:
	schemetype scheme;
	modeltype model;
	FluidType fluid;
	GradientType Gradient;
	UserForceType UserForce;
	int NbVelocities;
	double cs,cs2;
	double deltaT;
	double ErrorMax;

	TetaType tetaType;
	double teta;
	ViscosityType viscosityType;



};

class CalculationParameters {
private:
	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive &ar, const unsigned int version)
	{
		ar & BOOST_SERIALIZATION_NVP(parrallel)
		   & BOOST_SERIALIZATION_NVP(NbSteps)
		   & BOOST_SERIALIZATION_NVP(IntervalOutput)
		   & BOOST_SERIALIZATION_NVP(IntervalListing)
		   & BOOST_SERIALIZATION_NVP(NbVariablesOutput)
		   & BOOST_SERIALIZATION_NVP(OutputFileName)
		   & BOOST_SERIALIZATION_NVP(density)
		   & BOOST_SERIALIZATION_NVP(velocity);
		/*  ar & BOOST_SERIALIZATION_NVP(argc)
		   & BOOST_SERIALIZATION_NVP(argv);*/
	}
public:
	void Set_OutputFormat(OutputFormat OutputFormat_) {Format=OutputFormat_;};
	OutputFormat Get_OutputFormat(){return Format;};
	std::string Get_OutputFileName() const;
	void Set_OutputFileName(std::string string_);
	void Set_Parallel(parralleltype parrallel_);
	parralleltype Get_Parallel() const;
	std::string* Get_PtrVariablesOutput() const;
	int Get_NbVariablesOutput() const;
	void Set_NbStep(int NbStep_);
	int Get_NbStep() const;
	void Set_OutPutNSteps(int NbStep_);
	int Get_OutPutNSteps() const;
	void Set_listing(int IntervalListning);
	int Get_listing() const;
	bool Get_output_density() const {return density;};
	bool Get_output_velocity() const {return velocity;};
protected:
	void Set_VariablesOutput(int nbvar, std::string * strinput);
protected:
	parralleltype parrallel;
	int *argc;
	char ***argv;
	OutputFormat Format;
	int NbSteps, IntervalOutput,IntervalListing;
	std::string *VariablesOutput;
	int NbVariablesOutput;
	std::string OutputFileName;
	bool density,velocity;
};
class RestartParameters {
private:
	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive &ar, const unsigned int version)
	{
		ar & BOOST_SERIALIZATION_NVP(restart);
		ar & BOOST_SERIALIZATION_NVP(RestartFile);
	}
public:
	void set_RestartParameters();
protected:
	bool restart;
	std::string RestartFile;
};
class Parameters :
		public UserParameters,
		public MeshParameters,
		public ModelsParameters,
		public PhysicalParameters,
		public IniParameters,
		public BoundaryParameters,
		public SolverParameters,
		public CalculationParameters,
		public RestartParameters
{
private:
	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive &ar, const unsigned int version)
	{
		ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(UserParameters);
		ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(MeshParameters);
		ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ModelsParameters);
		ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(PhysicalParameters);
		ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(IniParameters);
		ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(BoundaryParameters);
		ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(SolverParameters);
		ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(CalculationParameters);
		ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(RestartParameters);

		ar & BOOST_SERIALIZATION_NVP(verbous);
	}
public:
	Parameters();
	virtual ~Parameters();

	void Change_Dimension(dimension dim_);
	void Set_BcType(NodeType Bc0,NodeType Bc1,NodeType Bc2,NodeType Bc3,NodeType Bc4=Wall,NodeType Bc5=Wall);

	void Set_Arguments(int *argc_,char ***argv_,bool verbous=false);

	int* Get_Argc() const;
	char*** Get_Argv() const;
	bool Get_Verbous() const;

	void Set_VariablesOutput(bool Rho, bool U){density=Rho;velocity=U;};

	//void Save_Parameters();
	//void Load_Parameters();

private:
	bool verbous;

};

#endif /* CORE_PARAMETERS_H_ */
