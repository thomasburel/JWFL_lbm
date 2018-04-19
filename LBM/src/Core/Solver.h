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
#include "Tau.h"
#include "Convergence.h"
#include "../Algorithm/LowOrder/CollideLowOrder.h"
#include "../Algorithm/LowOrder/StreamLowOrder.h"
#include "InitLBM.h"
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/list.hpp>
#include <boost/serialization/assume_abstract.hpp>
#include <boost/math/special_functions/fpclassify.hpp> // isnan
#include "../Algorithm/Tools/Gradients.h"
#include "../Algorithm/Tools/Extrapolation.h"

class Solver: public Tau, public Convergence, public Extrapolation{
public:
	Solver();
	virtual ~Solver();
	virtual void get_time()=0;
	virtual void run()=0;
	virtual void run(Parameters *UpdatedParam)=0;
	virtual void UpdateAllDomainFromFile(Parameters* UpdatedParam,InitLBM& ini)=0;
	virtual void UpdateAllDomain(Parameters* UpdatedParam,InitLBM& ini)=0;
	virtual void UpdateDomainBc(Parameters* UpdatedParam,InitLBM& ini)=0;
	virtual void UpdateWall(Parameters* UpdatedParam,InitLBM& ini)=0;
	virtual void UpdateInterior(Parameters* UpdatedParam,InitLBM& ini)=0;
	virtual void SetUserForce(Parameters* UpdatedParam)=0;
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
	virtual void run(Parameters *UpdatedParam)=0;
	virtual void UpdateAllDomainFromFile(Parameters* UpdatedParam,InitLBM& ini)=0;
	virtual void UpdateAllDomain(Parameters* UpdatedParam,InitLBM& ini)=0;
	virtual void UpdateDomainBc(Parameters* UpdatedParam,InitLBM& ini)=0;
	virtual void UpdateWall(Parameters* UpdatedParam,InitLBM& ini)=0;
	virtual void UpdateInterior(Parameters* UpdatedParam,InitLBM& ini)=0;
	virtual void SetUserForce(Parameters* UpdatedParam)=0;
	void set_f_ini();
	double** get_f_ini();
	void set_f_name();
	std::string* get_f_name();
protected:
	DistriFunct* f;
	double* ftmp; //tmp array for streaming (one dimensional to reduce memory consumption)

private:
	friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
       ar & boost::serialization::base_object<Solver>(*this) & f;
    }
};

class SolverTwoPhases: public Solver, public Gradients {
public:
	SolverTwoPhases();
	virtual ~SolverTwoPhases();
	virtual void get_time()=0;
	virtual void run()=0;
	virtual void run(Parameters *UpdatedParam)=0;
	virtual void UpdateAllDomainFromFile(Parameters* UpdatedParam,InitLBM& ini)=0;
	virtual void UpdateAllDomain(Parameters* UpdatedParam,InitLBM& ini)=0;
	virtual void UpdateDomainBc(Parameters* UpdatedParam,InitLBM& ini)=0;
	virtual void UpdateWall(Parameters* UpdatedParam,InitLBM& ini)=0;
	virtual void UpdateInterior(Parameters* UpdatedParam,InitLBM& ini)=0;
	virtual void SetUserForce(Parameters* UpdatedParam)=0;

protected:
	DistriFunct** f;
	double* ftmp; //tmp array for streaming (one dimensional to reduce memory consumption)
private:
	friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
       ar & boost::serialization::base_object<Solver>(*this);
    	ar & f[0]& f[1];

    }
};

class SolverSinglePhaseLowOrder: public SolverSinglePhase, public StreamLowOrder, public CollideLowOrder {
public:
	SolverSinglePhaseLowOrder();
	virtual ~SolverSinglePhaseLowOrder();
	virtual void get_time()=0;
	virtual void run()=0;
	virtual void run(Parameters *UpdatedParam)=0;
	virtual void UpdateAllDomainFromFile(Parameters* UpdatedParam,InitLBM& ini)=0;
	virtual void UpdateAllDomain(Parameters* UpdatedParam,InitLBM& ini)=0;
	virtual void UpdateDomainBc(Parameters* UpdatedParam,InitLBM& ini)=0;
	virtual void UpdateWall(Parameters* UpdatedParam,InitLBM& ini)=0;
	virtual void UpdateInterior(Parameters* UpdatedParam,InitLBM& ini)=0;
	virtual void SetUserForce(Parameters* UpdatedParam){Set_UserForce(UpdatedParam);};
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
	virtual void run(Parameters *UpdatedParam)=0;
	virtual void UpdateAllDomainFromFile(Parameters* UpdatedParam,InitLBM& ini)=0;
	virtual void UpdateAllDomain(Parameters* UpdatedParam,InitLBM& ini)=0;
	virtual void UpdateDomainBc(Parameters* UpdatedParam,InitLBM& ini)=0;
	virtual void UpdateWall(Parameters* UpdatedParam,InitLBM& ini)=0;
	virtual void UpdateInterior(Parameters* UpdatedParam,InitLBM& ini)=0;
	virtual void SetUserForce(Parameters* UpdatedParam){Set_UserForce(UpdatedParam);};
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
	void Add_OneDistributionToDictionary();
	void Updated_OneDistributionToDictionary();
	void Write_Breakpoint(Parameters *Param);
	void Set_Solver(MultiBlock* PtrMultiBlock_,ParallelManager* PtrParallel_,WriterManager* PtrWriter_, Parameters* PtrParameters_);
	virtual void run()=0;
	virtual void run(Parameters *UpdatedParam)=0;
	virtual void UpdateAllDomainFromFile(Parameters* UpdatedParam,InitLBM& ini)=0;
	virtual void UpdateAllDomain(Parameters* UpdatedParam,InitLBM& ini)=0;
	virtual void UpdateDomainBc(Parameters* UpdatedParam,InitLBM& ini)=0;
	virtual void UpdateWall(Parameters* UpdatedParam,InitLBM& ini)=0;
	virtual void UpdateInterior(Parameters* UpdatedParam,InitLBM& ini)=0;
protected:
	std::vector<int>::iterator itBc;
	std::vector<int> IdNodeN,IdNodeE,IdNodeS,IdNodeW;//Block connections
	std::vector<int> IdNodeSW,IdNodeSE,IdNodeNW,IdNodeNE;//Block connections
	//Mark real and ghost nodes for communication
	std::vector<int> IdRNodeN,IdRNodeE,IdRNodeS,IdRNodeW,IdGNodeN,IdGNodeE,IdGNodeS,IdGNodeW;
	std::vector<int> IdRNodeSW,IdRNodeSE,IdRNodeNW,IdRNodeNE,IdGNodeSW,IdGNodeSE,IdGNodeNW,IdGNodeNE;
	std::vector<int> SolidIdRNodeN,SolidIdRNodeE,SolidIdRNodeS,SolidIdRNodeW,SolidIdGNodeN,SolidIdGNodeE,SolidIdGNodeS,SolidIdGNodeW;
	std::vector<int> SolidIdRNodeSW,SolidIdRNodeSE,SolidIdRNodeNW,SolidIdRNodeNE,SolidIdGNodeSW,SolidIdGNodeSE,SolidIdGNodeNW,SolidIdGNodeNE;

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
	void Add_TwoDistributionsToDictionary();
	void Updated_TwoDistributionsToDictionary();
	void Write_Breakpoint(Parameters *Param);
	void Set_Solver(MultiBlock* PtrMultiBlock_,ParallelManager* PtrParallel_,WriterManager* PtrWriter_, Parameters* PtrParameters_);
	virtual void run()=0;
	virtual void run(Parameters *UpdatedParam)=0;
	virtual void UpdateAllDomainFromFile(Parameters* UpdatedParam,InitLBM& ini)=0;
	virtual void UpdateAllDomain(Parameters* UpdatedParam,InitLBM& ini)=0;
	virtual void UpdateDomainBc(Parameters* UpdatedParam,InitLBM& ini)=0;
	virtual void UpdateWall(Parameters* UpdatedParam,InitLBM& ini)=0;
	virtual void UpdateInterior(Parameters* UpdatedParam,InitLBM& ini)=0;
protected:
	std::vector<int>::iterator itBc;
	std::vector<int> IdNodeN,IdNodeE,IdNodeS,IdNodeW;//Block connections
	std::vector<int> IdNodeSW,IdNodeSE,IdNodeNW,IdNodeNE;//Block connections
	//Mark real and ghost nodes for communication
	std::vector<int> IdRNodeN,IdRNodeE,IdRNodeS,IdRNodeW,IdGNodeN,IdGNodeE,IdGNodeS,IdGNodeW;
	std::vector<int> IdRNodeSW,IdRNodeSE,IdRNodeNW,IdRNodeNE,IdGNodeSW,IdGNodeSE,IdGNodeNW,IdGNodeNE;
	std::vector<int> SolidIdRNodeN,SolidIdRNodeE,SolidIdRNodeS,SolidIdRNodeW,SolidIdGNodeN,SolidIdGNodeE,SolidIdGNodeS,SolidIdGNodeW;
	std::vector<int> SolidIdRNodeSW,SolidIdRNodeSE,SolidIdRNodeNW,SolidIdRNodeNE,SolidIdGNodeSW,SolidIdGNodeSE,SolidIdGNodeNW,SolidIdGNodeNE;
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
