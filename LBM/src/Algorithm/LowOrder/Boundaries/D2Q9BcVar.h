/*
 * ============================================================================
 * D2Q9BcVar.h
 *
 *  Created on: 29 Aug 2016
 *      Author: Thomas Burel
 *      Version     :
 *      Copyright   : Your copyright notice
 *      Description : LBM Code Developing by JWFL at University of Strathclyde
 * ============================================================================
 */

#ifndef SRC_ALGORITHM_LOWORDER_BOUNDARIES_D2Q9BCVAR_H_
#define SRC_ALGORITHM_LOWORDER_BOUNDARIES_D2Q9BCVAR_H_

#include "../../Tools/Extrapolation.h"
#include "../../../Mesh/SingleBlock/Node2D.h"
#include "../../../Core/Parameters.h"
#include "../../../Core/GlobalDef.h"
#include <map>
#include "../../../Core/Dictionary.h"
#include "../../../Mesh/SingleBlock/NodeArrays.h"

class D2Q9BcVar{
public:
	D2Q9BcVar();
	virtual ~D2Q9BcVar();

protected:
	double ru1,ru2;
	double q23,q16,q13;
	short int OppositeBc[9]; //opposite direction in the distribution function
	double omegaBc[9]; //Weight of the distribution function. It is need for some boundary conditions
	double **EiBc;
	//double EiBc[9][2];
	double rhodiff;
	double SumWeightS,SumWeightE,SumWeightN,SumWeightW;
	double SumWeightConcaveSE,SumWeightConcaveNE,SumWeightConcaveNW,SumWeightConcaveSW;
	double SumWeightConvexSE,SumWeightConvexNE,SumWeightConvexNW,SumWeightConvexSW;
	int *BounceBackWallConnect[9];///Keep the rule of bounce back to avoid a "if" condition
	int DiffusekWallConnect[9][3];///Keep the rule of diffuse to avoid a "if" condition
	double invCs2;
	
};

#endif /* SRC_ALGORITHM_LOWORDER_BOUNDARIES_D2Q9BCVAR_H_ */
