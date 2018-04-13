/*
 * UserInitLBM.cpp
 *
 *  Created on: 12 Jun 2015
 *      Author: Thomas Burel
 *
 *  This file is used to initialise the domain for all models (2D, 3D and single or multiphases).
 *  For single phase, Rho is the density, U[0] is Ux, U[1] is Uy and U[2] is Uz
 *  For Multiphase, Rho is the mass fraction (0 <= Rho <= 1), U[0] is Ux, U[1] is Uy and U[2] is Uz
 *  	The mass Fraction (Alpha) is defined as Rho= Rho_1 * Alpha +Rho_2 * (1-Alpha)
 *
 *  Notes: element number is not used at this time and set to 0.
 */

#include "UserInit.h"

UserInit::UserInit() {
	Block_=0;

}

UserInit::~UserInit() {
	// TODO Auto-generated destructor stub
}

void UserInit::UserBc(Parameters& PtrParameters, int elem, int nodenumber, double* pos,double& Rho, double* U, double& alpha){

	double Umax,H,L,Pmax,Pmin;
	double teta,sigma,Diameter,Re,Ca;
	PtrParameters.Get_UserParameters(Umax,H,L,Pmax,Pmin);
	PtrParameters.Get_UserDroplets(teta,sigma,Diameter,Re,Ca);
	double HeighChannel=19.0;
	double first_wall=3;
	double visco=0.1;
	U[0]=0.0;
	U[1]=0.0;
	Rho=1.0;
	alpha=0.0;
	if(pos[0]<=5)
		alpha=1.0;
	if(pos[0]<=1)
		U[0]=Umax*(1.0-pow(2.0*(pos[1]-HeighChannel/2.0-first_wall)/HeighChannel,2.0));
	else
		U[0]=0.0;
}

void UserInit::UserIc (Parameters& PtrParameters, int elem, int nodenumber, double* pos ,double& Rho, double* U, double& alpha){
	double Umax,H,L,Pmax,Pmin;
	double teta,sigma,Diameter,Re,Ca;
	PtrParameters.Get_UserParameters(Umax,H,L,Pmax,Pmin);
	PtrParameters.Get_UserDroplets(teta,sigma,Diameter,Re,Ca);


	U[0]=0.0;
	U[1]=0.0;

	Rho=1.0;
	alpha=0.0;

	double R=8;

	if(pos[0]<=5)
	{

		alpha=1.0;
	}

}

