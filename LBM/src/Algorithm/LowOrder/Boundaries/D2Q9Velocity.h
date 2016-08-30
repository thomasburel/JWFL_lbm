/*
 * ============================================================================
 * D2Q9Velocity.h
 *
 *  Created on: 29 Aug 2016
 *      Author: Thomas Burel
 *      Version     :
 *      Copyright   : Your copyright notice
 *      Description : LBM Code Developing by JWFL at University of Strathclyde
 * ============================================================================
 */

#ifndef SRC_ALGORITHM_LOWORDER_BOUNDARIES_D2Q9VELOCITY_H_
#define SRC_ALGORITHM_LOWORDER_BOUNDARIES_D2Q9VELOCITY_H_

#include "D2Q9BcVar.h"

class D2Q9Velocity: public D2Q9BcVar {
public:
	D2Q9Velocity();
	virtual ~D2Q9Velocity();

	void BC_HeZou_U(int & NormalBc, double* fi, double U, double V);

private:
	void FUNC_HeZou_U (double & a,double & b,double & c,double & d,double & e,double & f,double & g,double & h,double & i,double & U,double & V);


};

#endif /* SRC_ALGORITHM_LOWORDER_BOUNDARIES_D2Q9VELOCITY_H_ */
