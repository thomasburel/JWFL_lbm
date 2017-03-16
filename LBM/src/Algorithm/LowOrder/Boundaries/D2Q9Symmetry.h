/*
 * ============================================================================
 * D2Q9Symmetry.h
 *
 *  Created on: 2 Sep 2016
 *      Author: Thomas Burel
 *      Version     :
 *      Copyright   : Your copyright notice
 *      Description : LBM Code Developing by JWFL at University of Strathclyde
 * ============================================================================
 */

#ifndef SRC_ALGORITHM_LOWORDER_BOUNDARIES_D2Q9SYMMETRY_H_
#define SRC_ALGORITHM_LOWORDER_BOUNDARIES_D2Q9SYMMETRY_H_

#include "D2Q9BcVar.h"

class D2Q9Symmetry: public D2Q9BcVar {
public:
	D2Q9Symmetry();
	virtual ~D2Q9Symmetry();

	void Set_Symmetry(Parameters *Param, double ** &Ei);
	void SetSymmetry(SymmetryType SymmetryType_);
	void SetSymmetryType(SymmetryType SymmetryType_);
	void ApplySymmetry(int const &BcNormal,int const *Connect, double const &Rho_def, double const *UDef, DistriFunct* f_in, double *Rho, double *U, double *V);
private:
	void ApplySymmetryOnNode(int const &BcNormal,int const *Connect, double const &Rho_def, double const *UDef, double *LocalForce, DistriFunct* f_in);


// Pointers on function
///Simplify notation for pointer on a member function of D2Q9Symmetry class for Symmetry model used
	typedef void(D2Q9Symmetry::*SymmetryMethod)(int const &BcNormal,int const *Connect, double const &Rho_def, double const *UDef, double *LocalForce, DistriFunct* f_in);
//Define name for pointers on functions
	SymmetryMethod PtrSymmetryMethod;
	double *LocalForce;
};

#endif /* SRC_ALGORITHM_LOWORDER_BOUNDARIES_D2Q9SYMMETRY_H_ */
