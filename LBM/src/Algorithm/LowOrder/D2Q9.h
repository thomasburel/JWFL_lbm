/*
 * D2Q9.h
 *
 *  Created on: 9 Jun 2015
 *      Author: thomas
 */

#ifndef ALGORITHM_LOWORDER_D2Q9_H_
#define ALGORITHM_LOWORDER_D2Q9_H_
#include "D2Q9CommonVar.h"
#include "../../Core/Parameters.h"
#include "../../Core/Solver.h"
#include "../../Core/GlobalDef.h"
#include "Boundaries/D2Q9Bc.h"
#include <iostream>
#include <cmath>

class D2Q9: public SolverSinglePhaseLowOrder2D, public D2Q9Bc, protected D2Q9CommonVar {
public:
	D2Q9();
	D2Q9(MultiBlock* MultiBlock__,ParallelManager* parallel__,WriterManager* Writer__, Parameters* Parameters__,InitLBM& ini);
	virtual ~D2Q9();
	virtual void init(InitLBM& ini);
	virtual void run();
	virtual void run(Parameters* UpdatedParam);
	virtual void UpdateAllDomainFromFile(Parameters* UpdatedParam,InitLBM& ini);
	virtual void UpdateAllDomain(Parameters* UpdatedParam,InitLBM& ini);
	virtual void UpdateDomainBc(Parameters* UpdatedParam,InitLBM& ini);
	virtual void UpdateWall(Parameters* UpdatedParam,InitLBM& ini);
	virtual void UpdateInterior(Parameters* UpdatedParam,InitLBM& ini);


private:
	//! Initialise D2Q9.
    /*!
      \param ini : initialisation class (generic initialised methods).
    */
	void InitD2Q9(MultiBlock* MultiBlock__,ParallelManager* parallel__,WriterManager* Writer__, Parameters* Parameters_ ,InitLBM& ini);
	void InitialiseFromFile();
	void InitAllDomain(InitLBM& ini);
	void InitDomainBc(InitLBM& ini);
	void InitWall(InitLBM& ini);
	void InitInterior(InitLBM& ini);
	//! Initialise the distributions.
	void InitDistAllDomain();
	void InitDistDomainBc();
	void InitDistWall();
	void InitDistInterior();

	void Set_Collide();
	void CollideD2Q9();
	void CollideD2Q9_NoBodyForce();
	void CollideD2Q9_WithBodyForce();
	void CollideD2Q9_ExcludeWalls_NoBodyForce();
	void CollideD2Q9_ExcludeWalls_WithBodyForce();
	void StreamD2Q9();

	void Set_PointersOnFunctions();

// Macroscopic calculation
	void Set_Macro();
	//void UpdateMacroVariables();
	void UpdateMacroVariables();
	//void UpdateMacroVariables without body force();
	void UpdateMacroVariables_NoBodyForce();
	//void UpdateMacroVariables with body force();
	void UpdateMacroVariables_WithBodyForce();
	/// Calculate \f$\rho\f$ and \f$\vec{U}\f$ in the local domain
	void MacroVariables(int& idx);
	/// Calculate \f$\rho\f$ and \f$\vec{U}\f$ in the local domain with including the external force
	void MacroVariablesWithForce(int& idx);

	//Calculate Pressure
	void UpdatePressure();
	void CalculatePressure(int const &idx);
	double IdealGazIsothermalPressure(double const &Rho);

private:
	//Streaming by type of node
	void SelectStream(int & nodenumber, unsigned int& direction);
	void TmptoDistri(unsigned int& direction);
	void InteriorStream(int & nodenumber, unsigned int& direction);
	void GhostStream(int & nodenumber, unsigned int& direction);
	void WallStream(int & nodenumber, unsigned int& direction);
	void CornerStream(int & nodenumber, unsigned int& direction);
	void PeriodicStream(int & nodenumber, unsigned int& direction);
	void VelocityStream(int & nodenumber, unsigned int& direction);
	void PressureStream(int & nodenumber, unsigned int& direction);
	void SolidStream(int & nodenumber, unsigned int& direction);

	//Set streaming for Boundaries conditions (Stream only what we are able to stream)
	void Set_BcType(int & nodenumber);
	void Set_GhostType(int & nodenumber);
	void Set_WallType(int & nodenumber);
	void Set_CornerType(int & nodenumber);
	void Set_VelocityType(int & nodenumber);
	void Set_PressureType(int & nodenumber);
	void Set_BcType();
	void Set_GhostType(NodeGhost2D& Node);
	void Set_WallType(NodeWall2D& Node);
	void Set_CornerType(NodeCorner2D& Node);
	void Set_VelocityType(NodeVelocity2D& Node);
	void Set_PressureType(NodePressure2D& Node);
	void Set_SymmetryType(NodeSymmetry2D& NodeIn);
	void Set_PeriodicType(NodePeriodic2D& NodeIn);

	//Set Orientation of Boundaries

	void StreamingOrientation(int & nodenumber, bool GhostStreaming[9]);
	void StreamingOrientation(NodeGhost2D& Node, bool GhostStreaming[9]);
	void StreamingOrientation(NodeWall2D& Node, bool WallStreaming[9]);
	void StreamingOrientation(NodeCorner2D& Node, bool CornerStreaming[9]);
	void StreamingOrientation(NodeVelocity2D& Node, bool VelocityStreaming[9]);
	void StreamingOrientation(NodePressure2D& Node, bool PressureStreaming[9]);
	void StreamingOrientation(NodeSymmetry2D& Node, bool SymmetryStreaming[9]);
	void StreamingOrientation(NodePeriodic2D& Node, bool PeriodicStreaming[9]);

	//Apply boundary conditions
	void ApplyBc();
	void ApplyPatchPressure(PressurePatchBc& PresPatchBc);
	void ApplyPatchVelocity(VelocityPatchBc& VelPatchBc);
	void ApplyPatchSymmetry(SymmetryPatchBc& SymPatchBc);
	void ApplyPatchPeriodic(PeriodicPatchBc& PerPatchBc);
	//Apply conditions before streaming
	void ApplyPreStreamBc();
	void ApplyPreStreamBc_on();
	void ApplyPreStreamBc_off(){};	
	void ApplyPatchPressurePreStream(PressurePatchBc& PresPatchBc);
	void ApplyPatchVelocityPreStream(VelocityPatchBc& VelPatchBc);
	void ApplyPatchSymmetryPreStream(SymmetryPatchBc& SymPatchBc);
	void ApplyPatchPeriodicPreStream(PeriodicPatchBc& PerPatchBc);
	unsigned int& Connect (int &NodeNumber,unsigned int& direction);
	//void Set_Connect(int &NodeNumber,unsigned int& direction);
	double Cal_RhoCorner(int &normalBc, int &nodenumber);
	double Cal_RhoCorner(NodeCorner2D& Node);

	//Communications
	void IniComVariables();
	void SyncFromGhost();
	void SyncToGhost();
	void GhostNodesSyncFromGhost();
	void GhostNodesSyncToGhost();
	void CornerNodesSyncFromGhost();
	void CornerNodesSyncToGhost();
	void SyncMacroVarToGhost();

private:
	double **F;///< Surface Force and Colour gradient/density gradient
	double* tmp;// variable to copy tmp to distribution function
//	double Ei[9][2]; //Velocity in the distribution function
//	double omega[9];//Weight in the distribution function
//	unsigned int Opposite[9]; //opposite direction in the distribution function

	int** DiagConnect; //not use
	unsigned int tmpreturn; //return for connection between nodes
	int intTmpReturn;
	double doubleTmpReturn;

	int Nd_variables_sync;//number of variable has to be synchronise
	std::vector<double*> SyncVar;
	double ***buf_send, ***buf_recv; //buffers to send and receive
	int *size_buf; // size of buffers

	int Nd_MacroVariables_sync;//number of variable has to be synchronise
	double ***buf_MacroSend, ***buf_MacroRecv; //buffers to send and receive
	int *size_MacroBuf; // size of buffers

///Simplify notation for pointer on a member function of D2Q9ColourFluid class for macroscopic variables calculation methods

	typedef void(D2Q9::*Macro)();
	typedef void(D2Q9::*Collision)();
	typedef void(D2Q9::*PreStream)();
	Macro PtrMacro;///< Macroscopic pointer
	Collision PtrCollision;///< Collision pointer
	PreStream PtrPreStream;///< Pre-streaming pointer for applying BC in two times as for half-way bounce-back
private:
	friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
       ar & boost::serialization::base_object<SolverSinglePhaseLowOrder2D>(*this);
       ar & Nd_variables_sync;
    }

};






#endif /* ALGORITHM_LOWORDER_D2Q9_H_ */
