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

enum dimension {D2,D3};

enum schemetype {Q9,Q16};

enum modeltype {SinglePhase, ColourFluid};

enum parralleltype {Serial,Mpi,Openmp,Hybrid};

enum WallType {BounceBack, HalfWayBounceBack, Diffuse, Specular};

enum OutputFormat {CGNSFormat,TecplotFormat};


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

class ModelsParameters {
public:
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
		   & BOOST_SERIALIZATION_NVP(Tau_2);
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
		   //for (int i=0;i<NbGlobalBcType;i++)
			//   ar & BOOST_SERIALIZATION_NVP(GlobalBcType[i]);
	}
public:
	NodeType Get_GlobalBcType(int side) const;
	void Set_WallType(WallType WallType_);
	WallType Get_WallType() const;
	void Set_SymmetryType(SymmetryType SymmetryType_);
	SymmetryType Get_SymmetryType() const;

protected:
	WallType WallTypeParam;
	SymmetryType SymmetryTypeParam;
	int NbGlobalBcType;
	NodeType * GlobalBcType;
};

class SolverParameters {
private:
	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive &ar, const unsigned int version)
	{
		ar & BOOST_SERIALIZATION_NVP(scheme)
		   & BOOST_SERIALIZATION_NVP(NbVelocities)
		   & BOOST_SERIALIZATION_NVP(Tau);
	}

public:
	void set_SolverParameters();
	void Set_Scheme(schemetype scheme_);
	schemetype Get_Scheme() const;
	void Set_Model(modeltype model_);
	modeltype Get_Model() const;
	int Get_NbVelocities() const;
	void Set_Tau(double Tau_=1);
	double Get_Tau() const;

protected:
	schemetype scheme;
	modeltype model;
	int NbVelocities;
	double Tau;
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

	void Set_VariablesOutput(bool Rho, bool U);

	//void Save_Parameters();
	//void Load_Parameters();

private:
	bool verbous;

};

#endif /* CORE_PARAMETERS_H_ */
