/*
 * Simulation.h
 *
 *  Created on: 17 Apr 2015
 *      Author: thomas
 */

#ifndef CORE_SIMULATION_H_
#define CORE_SIMULATION_H_


#include "../Mesh/MultiBlock.h"
#include "../Algorithm/LowOrder.h"
#include "Solver.h"
#include <string>
#include "../Parallelism/ParallelManager.h"
#include "../IO/Writers.h"
#include "Parameters.h"
#include "InitLBM.h"
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <fstream>
#include <ios>

using namespace std;
//



class Simulation {
public:
	Simulation();
	virtual ~Simulation();
	void InitSimu(Parameters &Parameters_, bool create_mesh);
	InitLBM& InitSimu();
	void Create_Mesh();
	void Create_Mesh(SolverEnum::dimension dim, int Nx, int Ny, int Nz);
	void Import_Mesh(string MeshFile);
	void RunSimu();
	void RunSimu(Parameters &UpdatedParam);

	void UpdateAllDomain(Parameters &UpdatedParam);
	void UpdateDomainBc(Parameters &UpdatedParam);
	void UpdateWall(Parameters &UpdatedParam);
	void UpdateInterior(Parameters &UpdatedParam);

	void FinalizeSimu();

	double get_time();
	bool Is_MainProcessor(){return parallel->isMainProcessor();};

	void Save_Simulation();

	void barrier();
private:
	void SelectSolver();
	void SelectOutputFormat();
private:
	double time;
	//Generic objects
	InitLBM ini;
	Solver* Solver_;
	MultiBlock* MultiBlock_;
	ParallelManager* parallel;
	WriterManager* Writer;
	Parameters* PtrParameters;

	//real object created => add here new writer, solver, block....
	CGNS* WriterCGNS;
	TecplotIO* WriterTecplot;
	MultiBlock2D* MultiBlock2D_;
	D2Q9* SolverD2Q9;
	D2Q9TwoPhases* SolverD2Q9TwoPhases;
/*private:
	friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
       ar & time & ini & Solver_ &  MultiBlock_ & parrallel & Writer & PtrParameters;
       ar & time & ini & PtrParameters;
    }*/
private:
	void Save_Parameters(Parameters & object);
	void Save_Parameters(Parameters & object,std::string &filename);
	void Load_Parameters(Parameters & object);
	void Load_Parameters(Parameters & object,std::string &filename);
	void Save_InitLBM(InitLBM & object);
	void Save_InitLBM(InitLBM & object,std::string &filename);
	void Load_InitLBM(InitLBM & object);
	void Load_InitLBM(InitLBM & object,std::string &filename);
	void Save_Writer();
	void Save_Writer(std::string &filename,ios_base::openmode mode = ios_base::out);
	void Load_Writer();
	void Load_Writer(std::string &filename);
	void Save_MultiBlock();
	void Save_MultiBlock(std::string &filename,ios_base::openmode mode = ios_base::out);
	void Load_MultiBlock();
	void Load_MultiBlock(std::string &filename);
	void Save_Solver(std::string &filename,ios_base::openmode mode);


};


#endif /* CORE_SIMULATION_H_ */
