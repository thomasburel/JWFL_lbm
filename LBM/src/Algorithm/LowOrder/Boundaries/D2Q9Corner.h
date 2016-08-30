/*
 * ============================================================================
 * D2Q9Corner.h
 *
 *  Created on: 29 Aug 2016
 *      Author: Thomas Burel
 *      Version     :
 *      Copyright   : Your copyright notice
 *      Description : LBM Code Developing by JWFL at University of Strathclyde
 * ============================================================================
 */

#ifndef SRC_ALGORITHM_LOWORDER_BOUNDARIES_D2Q9CORNER_H_
#define SRC_ALGORITHM_LOWORDER_BOUNDARIES_D2Q9CORNER_H_

#include "D2Q9BcVar.h"

class D2Q9Corner: public D2Q9BcVar {
public:
	D2Q9Corner();
	virtual ~D2Q9Corner();

	void BC_corner(int & NormalCorner, double* fi, double Rho, double U, double V);
	void BC_corner_no_vel_concave(int & NormalCorner, double* fi, double Rho);
	void BC_corner_no_vel_convex(int & NormalCorner, double* fi);

private:
	void FUNC_corner (double & a,double & b,double & c,double & d,double & e,double & f,double & g,double & h,double & i,double & U,double & V,double & Rho);
	void FUNC_corner_no_vel_concave (double & a,double & b,double & c,double & d,double & e,double & f,double & g,double & h,double & i,double & Rho);
	void FUNC_corner_no_vel_convex (double & a,double & b,double & c,double & d,double & e,double & f);

private:
		double InvRho,InvU,InvV;
};

#endif /* SRC_ALGORITHM_LOWORDER_BOUNDARIES_D2Q9CORNER_H_ */
