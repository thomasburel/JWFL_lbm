/*
 * ============================================================================
 * Convergence.h
 *
 *  Created on: 30 Sep 2016
 *      Author: Thomas Burel
 *      Version     :
 *      Copyright   : Your copyright notice
 *      Description : LBM Code Developing by JWFL at University of Strathclyde
 * ============================================================================
 */

#ifndef SRC_CORE_CONVERGENCE_H_
#define SRC_CORE_CONVERGENCE_H_
#include "../Mesh/MultiBlock.h"
#include "Dictionary.h"
#include "../Mesh/SingleBlock/Patch/PatchBc.h"
#include "Viscosity.h"
enum TypeConverge{none,GlobalConvergence,FieldConvergence};
enum TypeConvergeScalar{SacalarGlobal, ScalarField};
enum TypeConvergeVector{VectorGlobal,VectorField};
class Convergence {
public:
	Convergence();
	virtual ~Convergence();
	void Set_Convergence();
	void Calcul_Error();
	double Get_Error(){return Error;};

private:
	void Calcul_Error_ScalarField();
	void Calcul_Error_VectorField();
	void NoCalcul_Error(){};
	double Calcul_Error_ScalarFieldInOneProc();
	void Set_ConvergencePatchBc();
	void Calcul_PorousMediaConvergence();
	double Calcul_EOR();
	double Calcul_Permeability_SinglePhase();

	void Avg_ScalarNodeArray(std::vector<int>& NodeArray,std::vector<int>& NodeArraySpecialWall,std::vector<int>& NodeArrayGlobalCorner,double *&Var1,double &sum);
	void Avg_VectorNodeArray(std::vector<int>& NodeArray,std::vector<int>& NodeArraySpecialWall,std::vector<int>& NodeArrayGlobalCorner,double *&Var1,double *&Var2,double &sum,double &sum2);
	void Sum_ScalarNodeArray(std::vector<int>& NodeArray,std::vector<int>& NodeArraySpecialWall,std::vector<int>& NodeArrayGlobalCorner,double *&Var1,double &sum);
	void Sum_VectorNodeArray(std::vector<int>& NodeArray,std::vector<int>& NodeArraySpecialWall,std::vector<int>& NodeArrayGlobalCorner,double *&Var1,double *&Var2,double &sum1,double &sum2);


protected:
	MultiBlock *PtrMultiBlockConv;
	Dictionary *PtrDicConv;
	Parameters *PtrParmConv;
	PatchBc *PtrPatchBcConv;
	Viscosity *PtrViscosityConv;

private:
	double Sum_Current;
	double Error;///< Save error and it used inside the sum if needed
	double Error_sum, Error_avg;
	double Error_tmp;

	std::vector<int> InletPatchId,OutletPatchId;
	double *RhoNEOR;
	std::vector<int> *NodeId,*NodeIdSpeWall,*NodeIdGloCorner;
	double avg;


	int NbNodes;

	double *Scalar_CurrentTime;///< Pointer on the scalar of the Current Time
	double **Vector_CurrentTime;///< Pointer on the scalar of the Current Time
	//Save previous time step
	double *Scalar_last;///< Save previous time step for scalars
	double **Vector_last;///< Save previous time step for vectors

};

#endif /* SRC_CORE_CONVERGENCE_H_ */
