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
	double Calcul_Error_ScalarField();



protected:
	MultiBlock *PtrMultiBlockConv;
	Dictionary *PtrDicConv;

private:
	double Sum_Current;
	double Error;///< Save error and it used inside the sum if needed
	double Error_sum, Error_avg;
	double Error_tmp;


	int NbNodes;

	double *Scalar_CurrentTime;///< Pointer on the scalar of the Current Time
	double **Vector_CurrentTime;///< Pointer on the scalar of the Current Time
	//Save previous time step
	double *Scalar_last;///< Save previous time step for scalars
	double **Vector_last;///< Save previous time step for vectors

};

#endif /* SRC_CORE_CONVERGENCE_H_ */
