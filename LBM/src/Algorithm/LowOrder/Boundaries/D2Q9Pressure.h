/*
 * ============================================================================
 * D2Q9Pressure.h
 *
 *  Created on: 29 Aug 2016
 *      Author: Thomas Burel
 *      Version     :
 *      Copyright   : Your copyright notice
 *      Description : LBM Code Developing by JWFL at University of Strathclyde
 * ============================================================================
 */

#ifndef SRC_ALGORITHM_LOWORDER_BOUNDARIES_D2Q9PRESSURE_H_
#define SRC_ALGORITHM_LOWORDER_BOUNDARIES_D2Q9PRESSURE_H_
#include "D2Q9BcVar.h"

class D2Q9Pressure: public D2Q9BcVar {
public:
	D2Q9Pressure();
	virtual ~D2Q9Pressure();


	void ApplyPressure(NodePressure2D& Node, double* fi, double Rho, double & U, double & V);

	void BC_HeZou_P(int & NormalBc, double* fi, double Rho, double & U, double & V);

private:
	void FUNC_HeZou_P (double & a,double & b,double & c,double & d,double & e,double & f,double & g,double & h,double & i,double & V,double & Rho);


// Pointers on function
///Simplify notation for pointer on a member function of D2Q9Pressure class for Pressure model used
	typedef void(D2Q9Pressure::*PressureMethod)(int & NormalBc, double* fi, double Rho, double & U, double & V);
	typedef double(D2Q9Pressure::*CalculRho)(int & NormalBc, double* fi, double Rho, double & U, double & V);
//Define name for pointers on functions
	PressureMethod PtrPressureMethod;
	CalculRho PtrCalculRho;
};

#endif /* SRC_ALGORITHM_LOWORDER_BOUNDARIES_D2Q9PRESSURE_H_ */
