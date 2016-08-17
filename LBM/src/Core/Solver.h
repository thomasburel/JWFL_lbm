/*
 * Solver.h
 *
 *  Created on: 17 Apr 2015
 *      Author: thomas
 */

#ifndef CORE_SOLVER_H_
#define CORE_SOLVER_H_

#include "Parameters.h"
#include "Solution.h"
#include "../Algorithm/LowOrder/CollideLowOrder.h"
#include "../Algorithm/LowOrder/StreamLowOrder.h"
#include "InitLBM.h"
#include "../Algorithm/LowOrder/Gradients.h"

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/list.hpp>
#include <boost/serialization/assume_abstract.hpp>

class Solver {
public:
	Solver();
	virtual ~Solver();
	virtual void get_time()=0;
	virtual void run()=0;
protected:
	void saveParameters();
protected:
	Parameters *PtrParameters;
	//double InvTau;
	int nbvelo;
	int nbnode;
//	double **Ei;//Ei[dimension][nb_velocity]
//	double *omega;//omega[nb_velocity]
private:
	friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
       ar & nbvelo & nbnode & PtrParameters;
    }
};

class SolverSinglePhase: public Solver {
public:
	SolverSinglePhase();
	virtual ~SolverSinglePhase();
	virtual void get_time()=0;
	virtual void run()=0;
	void set_f_ini();
	double** get_f_ini();
	void set_f_name();
	std::string* get_f_name();
protected:
	DistriFunct* f;
	double* ftmp; //tmp array for streaming (one dimensional to reduce memory consumption)
private:
	double** f_ini; //save pointers for breakpoints
	std::string* f_name;
private:
	friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
       ar & boost::serialization::base_object<Solver>(*this) & f;
       for (int i=0;i<nbvelo+3;i++)
    	   ar & f_name[i];
    }
};

class SolverTwoPhases: public Solver, public Gradients {
public:
	SolverTwoPhases();
	virtual ~SolverTwoPhases();
	virtual void get_time()=0;
	virtual void run()=0;
	void set_f_ini();
	double** get_f_ini();
	void set_f_name();
	std::string* get_f_name();
protected:
	DistriFunct** f;
	double* ftmp; //tmp array for streaming (one dimensional to reduce memory consumption)
private:
	double** f_ini; //save pointers for breakpoints
	std::string* f_name;
private:
	friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
       ar & boost::serialization::base_object<Solver>(*this);
    	ar & f[0]& f[1];
       for (int i=0;i<nbvelo+3;i++)
    	   ar & f_name[i];
    }
};

class SolverSinglePhaseLowOrder: public SolverSinglePhase, public StreamLowOrder, public CollideLowOrder {
public:
	SolverSinglePhaseLowOrder();
	virtual ~SolverSinglePhaseLowOrder();
	virtual void get_time()=0;
	virtual void run()=0;
private:
	friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
       ar & boost::serialization::base_object<SolverSinglePhase>(*this);
       ar & boost::serialization::base_object<StreamLowOrder>(*this);
       ar & boost::serialization::base_object<CollideLowOrder>(*this);
    }
};
class SolverTwoPhasesLowOrder: public SolverTwoPhases, public StreamLowOrder, public CollideLowOrder {
public:
	SolverTwoPhasesLowOrder();
	virtual ~SolverTwoPhasesLowOrder();
	virtual void get_time()=0;
	virtual void run()=0;
private:
	friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
       ar & boost::serialization::base_object<SolverTwoPhases>(*this);
       ar & boost::serialization::base_object<StreamLowOrder>(*this);
      // ar & boost::serialization::base_object<CollideLowOrder>(*this);
    }
};
class SolverSinglePhaseLowOrder2D: public SolverSinglePhaseLowOrder, public Solution2D {
public:
	SolverSinglePhaseLowOrder2D();
	virtual ~SolverSinglePhaseLowOrder2D();
	virtual void get_time();
	void Set_Solver(MultiBlock* PtrMultiBlock_,ParallelManager* PtrParallel_,WriterManager* PtrWriter_, Parameters* PtrParameters_);
	virtual void run()=0;

protected:
	std::vector<int>::iterator itBc;
	std::vector<int> IdNodeN,IdNodeE,IdNodeS,IdNodeW;//Block connections
	std::vector<int> IdNodeSW,IdNodeSE,IdNodeNW,IdNodeNE;//Block connections
	//Mark real and ghost nodes for communication
	std::vector<int> IdRNodeN,IdRNodeE,IdRNodeS,IdRNodeW,IdGNodeN,IdGNodeE,IdGNodeS,IdGNodeW;
	std::vector<int> IdRNodeSW,IdRNodeSE,IdRNodeNW,IdRNodeNE,IdGNodeSW,IdGNodeSE,IdGNodeNW,IdGNodeNE;
private:
	friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
       ar & boost::serialization::base_object<SolverSinglePhaseLowOrder>(*this);
       ar & boost::serialization::base_object<Solution2D>(*this);
       ar & IdNodeN & IdNodeE & IdNodeS & IdNodeW;
       ar & IdNodeSW & IdNodeSE & IdNodeNW & IdNodeNE;
    }
};
class SolverTwoPhasesLowOrder2D: public SolverTwoPhasesLowOrder, public Solution2D {
public:
	SolverTwoPhasesLowOrder2D();
	virtual ~SolverTwoPhasesLowOrder2D();
	virtual void get_time();
	void Set_Solver(MultiBlock* PtrMultiBlock_,ParallelManager* PtrParallel_,WriterManager* PtrWriter_, Parameters* PtrParameters_);
	virtual void run()=0;

protected:
	std::vector<int>::iterator itBc;
	std::vector<int> IdNodeN,IdNodeE,IdNodeS,IdNodeW;//Block connections
	std::vector<int> IdNodeSW,IdNodeSE,IdNodeNW,IdNodeNE;//Block connections
	//Mark real and ghost nodes for communication
	std::vector<int> IdRNodeN,IdRNodeE,IdRNodeS,IdRNodeW,IdGNodeN,IdGNodeE,IdGNodeS,IdGNodeW;
	std::vector<int> IdRNodeSW,IdRNodeSE,IdRNodeNW,IdRNodeNE,IdGNodeSW,IdGNodeSE,IdGNodeNW,IdGNodeNE;
private:
	friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
       ar & boost::serialization::base_object<SolverTwoPhasesLowOrder>(*this);
       ar & boost::serialization::base_object<Solution2D>(*this);
       ar & IdNodeN & IdNodeE & IdNodeS & IdNodeW;
       ar & IdNodeSW & IdNodeSE & IdNodeNW & IdNodeNE;
    }
};
#endif /* CORE_SOLVER_H_ */
