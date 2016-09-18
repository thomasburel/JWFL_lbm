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

	void SetVelocity(Parameters *Param);
	void ApplyVelocity(int const &BcNormal,int const *Connect, double const *UDef, DistriFunct * & f_in, double *Rho, double *U, double *V);
	inline void Get_CalculU(int const &BcNormal,int const *Connect, double const *UDef, double *U, double *V, double & Ureturn,double & Vreturn){(this->*PtrCalculU)(BcNormal,Connect,UDef,U,V); Ureturn=U[Connect[0]];Vreturn=V[Connect[0]];};

private:
//Function for calculating the velocity
	void FixU(int const &BcNormal,int const *Connect, double const *UDef, double *U, double *V);
	void NoGradU_1stOrder(int const &BcNormal,int const *Connect, double const *UDef, double *U, double *V);

	void BC_HeZou_U(int const &BcNormal,int const *Connect, double const *UDef, DistriFunct * & f_in, double & U, double & V);
	void FUNC_HeZou_U (double & a,double & b,double & c,double & d,double & e,double & f,double & g,double & h,double & i,double & U,double & V);

// Pointers on function
///Simplify notation for pointer on a member function of D2Q9Velocity class for Velocity model used
	typedef void(D2Q9Velocity::*VelocityMethod)(int const &BcNormal,int const *Connect, double const *UDef, DistriFunct * & f_in, double & U, double & V);
	typedef void(D2Q9Velocity::*CalculU)(int const &BcNormal,int const *Connect, double const *UDef, double *U, double *V);

//Define name for pointers on functions
	VelocityMethod PtrVelocityMethod;
	CalculU PtrCalculU;

	double feq,Rho;
	double InvU,InvV,rho,ru,rv,fequi;
};

#endif /* SRC_ALGORITHM_LOWORDER_BOUNDARIES_D2Q9VELOCITY_H_ */
