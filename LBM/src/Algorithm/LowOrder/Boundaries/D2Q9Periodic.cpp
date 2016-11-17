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
	PressureDrop=0;
	PressureDropAdjust=0;
}

D2Q9Periodic::~D2Q9Periodic() {
	// TODO Auto-generated destructor stub
}

void D2Q9Periodic::Set_Periodic(Parameters *Param){
	PtrPeriodicMethod=&D2Q9Periodic::ApplyPeriodicBc;
}
void D2Q9Periodic::ApplyPeriodic(int const &BcNormal,int const *Connect, double const &Rho_def, double const *UDef, DistriFunct* f_in, double *Rho, double *U, double *V, double weightDensity){
	(this->*PtrPeriodicMethod)(BcNormal,Connect,Rho_def,weightDensity,UDef,LocalForce,f_in);
}

void D2Q9Periodic::ApplyPeriodicBc(int const &BcNormal,int const *Connect, double const &Rho_def, double & weightDensity, double const *UDef, double *LocalForce, DistriFunct* f_in){

}
void D2Q9Periodic::ApplyPeriodicBc_PressureForce(int const &BcNormal,int const *Connect, double const &Rho_def, double & weightDensity, double const *UDef, double *LocalForce, DistriFunct* f_in){
	PressureDropAdjust=PressureDrop*weightDensity;
	switch(BcNormal)
			{
			case 2:
				f_in->f[0][Connect[0]]+=PressureDropAdjust*omegaBc[0];
				f_in->f[2][Connect[0]]+=PressureDropAdjust*omegaBc[2];
				f_in->f[5][Connect[0]]+=PressureDropAdjust*omegaBc[5];
				f_in->f[6][Connect[0]]+=PressureDropAdjust*omegaBc[6];
				break;
			case 4:
				f_in->f[0][Connect[0]]-=PressureDropAdjust*omegaBc[0];
				f_in->f[4][Connect[0]]-=PressureDropAdjust*omegaBc[4];
				f_in->f[7][Connect[0]]-=PressureDropAdjust*omegaBc[7];
				f_in->f[8][Connect[0]]-=PressureDropAdjust*omegaBc[8];
				break;
			case 1:
				f_in->f[0][Connect[0]]+=PressureDropAdjust*omegaBc[0];
				f_in->f[1][Connect[0]]+=PressureDropAdjust*omegaBc[1];
				f_in->f[5][Connect[0]]+=PressureDropAdjust*omegaBc[5];
				f_in->f[8][Connect[0]]+=PressureDropAdjust*omegaBc[8];
				break;
			case 3:
				f_in->f[0][Connect[0]]-=PressureDropAdjust*omegaBc[0];
				f_in->f[3][Connect[0]]-=PressureDropAdjust*omegaBc[3];
				f_in->f[6][Connect[0]]-=PressureDropAdjust*omegaBc[6];
				f_in->f[7][Connect[0]]-=PressureDropAdjust*omegaBc[7];
				break;
			default:
				std::cerr<<"Direction PeriodicBc-PressureForce not found"<<std::endl;
				break;
			}
}
