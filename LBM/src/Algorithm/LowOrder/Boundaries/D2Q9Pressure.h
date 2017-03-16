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
#include <cmath>
class D2Q9Pressure: public D2Q9BcVar {
public:
	D2Q9Pressure();
	virtual ~D2Q9Pressure();

	void Set_PressureBcs(Parameters *Param, double ** &Ei);
	void SetPressure(PressureModel PressureModel_,PressureType PressureType_);
	void SetPressureModel(PressureModel PressureModel_);
	void SetPressureType(PressureType PressureType_);
	void ApplyPressure(int const &BcNormal,int const *Connect, double const &Rho_def, DistriFunct* f_in, double *Rho, double *U, double *V);
	inline void Get_CalculRho(int const &BcNormal,int const *Connect, double const &Rho_def, double *Rho, double & Rhoreturn){(this->*PtrCalculRho)(BcNormal,Connect,Rho_def,Rho);Rhoreturn=Rho[Connect[0]];};
private:
//Methods for pressure models (generic)
	//He Zou with pressure formulation
	void BC_HeZou_P(int const &BcNormal,int const *Connect, double const &Rho_def, DistriFunct* f_in, double Rho, double & U, double & V);


//Functions for pressure models	(Calculations)
	void FUNC_HeZou_P (double & a,double & b,double & c,double & d,double & e,double & f,double & g,double & h,double & i,double & V,double & Rho);


//Function for calculating the pressure and setting in the global variable
	void FixRho(int const &BcNormal,int const *Connect, double const &Rho_def, double *Rho);
	void NoGradRho_1stOrder(int const &BcNormal,int const *Connect, double const &Rho_def, double *Rho);

//Private variables
	double 	InvRho,feq,U;

// Pointers on function
///Simplify notation for pointer on a member function of D2Q9Pressure class for Pressure model used
	typedef void(D2Q9Pressure::*PressureMethod)(int const &BcNormal,int const *Connect, double const &Rho_def, DistriFunct* f_in, double Rho, double & U, double & V);
	typedef void(D2Q9Pressure::*CalculRho)(int const &BcNormal,int const *Connect, double const &Rho_def, double *Rho);
//Define name for pointers on functions
	PressureMethod PtrPressureMethod;
	CalculRho PtrCalculRho;
};

#endif /* SRC_ALGORITHM_LOWORDER_BOUNDARIES_D2Q9PRESSURE_H_ */
