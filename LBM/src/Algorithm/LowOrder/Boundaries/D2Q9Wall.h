/*
 * ============================================================================
 * D2Q9Wall.h
 *
 *  Created on: 29 Aug 2016
 *      Author: Thomas Burel
 *      Version     :
 *      Copyright   : Your copyright notice
 *      Description : LBM Code Developing by JWFL at University of Strathclyde
 * ============================================================================
 */

#ifndef SRC_ALGORITHM_LOWORDER_BOUNDARIES_D2Q9WALL_H_
#define SRC_ALGORITHM_LOWORDER_BOUNDARIES_D2Q9WALL_H_

#include "D2Q9BcVar.h"

class D2Q9Wall: public D2Q9BcVar {
public:
	D2Q9Wall();
	virtual ~D2Q9Wall();

private:
	double rhodiff;
	double SumWeightS,SumWeightE,SumWeightN,SumWeightW;
	double SumWeightConcaveSE,SumWeightConcaveNE,SumWeightConcaveNW,SumWeightConcaveSW;
	double SumWeightConvexSE,SumWeightConvexNE,SumWeightConvexNW,SumWeightConvexSW;
};

#endif /* SRC_ALGORITHM_LOWORDER_BOUNDARIES_D2Q9WALL_H_ */
