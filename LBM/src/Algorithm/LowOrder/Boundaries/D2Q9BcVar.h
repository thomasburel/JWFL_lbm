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
#include "../../../Mesh/SingleBlock/Node2D.h"
class D2Q9BcVar {
public:
	D2Q9BcVar();
	virtual ~D2Q9BcVar();

protected:
	double ru1,ru2;
	double q23,q16,q13;
	short int Opposite[9]; //opposite direction in the distribution function
};

#endif /* SRC_ALGORITHM_LOWORDER_BOUNDARIES_D2Q9BCVAR_H_ */
