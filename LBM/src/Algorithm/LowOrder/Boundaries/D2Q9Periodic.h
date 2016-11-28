/*
 * ============================================================================
 * D2Q9Periodic.h
 *
 *  Created on: 2 Sep 2016
 *      Author: Thomas Burel
 *      Version     :
 *      Copyright   : Your copyright notice
 *      Description : LBM Code Developing by JWFL at University of Strathclyde
 * ============================================================================
 */

#ifndef SRC_ALGORITHM_LOWORDER_BOUNDARIES_D2Q9PERIODIC_H_
#define SRC_ALGORITHM_LOWORDER_BOUNDARIES_D2Q9PERIODIC_H_

#include "D2Q9BcVar.h"

class D2Q9Periodic: public D2Q9BcVar {
public:
	D2Q9Periodic();
	virtual ~D2Q9Periodic();

	void Set_Periodic(Parameters *Param);
	void ApplyPeriodic(int const &BcNormal,int const *Connect, double const &Rho_def, double weightDensity, double const *UDef, DistriFunct* f_in, double *Rho, double *U, double *V);

private:
	void ApplyPeriodicBc(int const &BcNormal,int const *Connect, double const &Rho_def, double & weightDensity, double const *UDef, double *LocalForce, DistriFunct* f_in);
	void ApplyPeriodicBc_PressureForce(int const &BcNormal,int const *Connect, double const &Rho_def, double & weightDensity, double const *UDef, double *LocalForce, DistriFunct* f_in);

// Pointers on function
///Simplify notation for pointer on a member function of D2Q9Periodic class for Periodic model used
	typedef void(D2Q9Periodic::*PeriodicMethod)(int const &BcNormal,int const *Connect, double const &Rho_def, double & weightDensity, double const *UDef, double *LocalForce, DistriFunct* f_in);
//Define name for pointers on functions
	PeriodicMethod PtrPeriodicMethod;
	double *LocalForce;
	double PressureDrop,PressureDropAdjust;
};

#endif /* SRC_ALGORITHM_LOWORDER_BOUNDARIES_D2Q9PERIODIC_H_ */
