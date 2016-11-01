/*
 * ============================================================================
 * D2Q9CommonVar.h
 *
 *  Created on: 6 Sep 2016
 *      Author: Thomas Burel
 *      Version     :
 *      Copyright   : Your copyright notice
 *      Description : LBM Code Developing by JWFL at University of Strathclyde
 * ============================================================================
 */

#ifndef SRC_ALGORITHM_LOWORDER_D2Q9COMMONVAR_H_
#define SRC_ALGORITHM_LOWORDER_D2Q9COMMONVAR_H_
#include <cmath>

class D2Q9CommonVar {
public:
	D2Q9CommonVar();
	virtual ~D2Q9CommonVar();

	short int *Opposite; //< Opposite direction in the distribution function
	double *omega; //< Weight of the distribution function. It is need for some boundary conditions
	double **Ei; //< Velocity in the distribution space
	double *Ei_Norm;//< Velocity Magnitude in the distribution space (only need it in few models)
	double pi;
};

#endif /* SRC_ALGORITHM_LOWORDER_D2Q9COMMONVAR_H_ */
