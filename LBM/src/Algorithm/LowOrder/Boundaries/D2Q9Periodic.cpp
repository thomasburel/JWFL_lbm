/*
 * ============================================================================
 * D2Q9Periodic.cpp
 *
 *  Created on: 2 Sep 2016
 *      Author: Thomas Burel
 *      Version     :
 *      Copyright   : Your copyright notice
 *      Description : LBM Code Developing by JWFL at University of Strathclyde
 * ============================================================================
 */

#include "D2Q9Periodic.h"

D2Q9Periodic::D2Q9Periodic() {
	PtrPeriodicMethod=0;
	LocalForce=0;
}

D2Q9Periodic::~D2Q9Periodic() {
	// TODO Auto-generated destructor stub
}

void D2Q9Periodic::Set_Periodic(Parameters *Param){
	PtrPeriodicMethod=&D2Q9Periodic::ApplyPeriodicBc;
}
void D2Q9Periodic::ApplyPeriodic(int const &BcNormal,int const *Connect, double const &Rho_def, double const *UDef, DistriFunct* f_in, double *Rho, double *U, double *V){
	(this->*PtrPeriodicMethod)(BcNormal,Connect,Rho_def,UDef,LocalForce,f_in);
}

void D2Q9Periodic::ApplyPeriodicBc(int const &BcNormal,int const *Connect, double const &Rho_def, double const *UDef, double *LocalForce, DistriFunct* f_in){

}
