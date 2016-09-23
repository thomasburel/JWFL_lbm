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
	PtrParameters.Get_UserParameters(Umax,H,L,Pmax,Pmin);
	if(pos[0<=0])
		U[0]=Umax*4.0*(pos[1]/H)*(1-(pos[1]/H));
	else
		U[0]=0.0;
	U[1]=0;
	Rho=Pmax-pos[0]*(Pmax-Pmin)/L;
	/*
//*********** Poiseuille ini***************
	double Ubot=-Umax;//-0.01;
	double Utop=Umax;//0.01;
	double Vleft=0;//-0.01;
	double Vright=0;//0.01;
//Left side of the domain
	if(pos[0]<=0)
	// Global Corner at the bottom left side
 		if(pos[1]<=0)
  		{
  			U[0]=Ubot;
  			U[1]=Vleft;
  			Rho=Pmax;
  			alpha=0;
  		}
  		else
 	// Global Corner at the Top left side
 		if(pos[1]>=H)
  		{
  			U[0]=Utop;
  			U[1]=Vleft;
  			Rho=Pmax;
  			alpha=0;
  		}
  	//Left side of the domain and excluding the two Global Corners
  		else
  		{
   			U[0]=(Utop-Ubot)/H*pos[1]+Ubot;//Umax*4.0*(pos[1]/H)*(1-(pos[1]/H));
  			U[1]=Vleft;
  			Rho=Pmax;
  			alpha=0;
 		}
	else
//Right side of the domain
 	if(pos[0]>=L)
	// Global Corner at the bottom right side
  		if(pos[1]<=0)
  		{
  			U[0]=Ubot;
  			U[1]=Vright;
  			Rho=Pmin;
  			alpha=0;
  		}
  		else
  	// Global Corner at the Top right side
  		if(pos[1]>=H)
  		{
  			U[0]=Utop;
  			U[1]=Vright;
  			Rho=Pmin;
  			alpha=0;
  		}
  	//Right side of the domain and excluding the two Global Corner
  		else
  		{
    		U[0]=(Utop-Ubot)/H*pos[1]+Ubot;//Umax*4.0*(pos[1]/H)*(1-(pos[1]/H));
  			U[1]=Vright;
  			Rho=Pmin;
  			alpha=0;
  		}
  	else
//Bottom side of the domain
 	if(pos[1]<=0)
  	{
  		U[0]=Ubot;
  		U[1]=(Vright-Vleft)/L*pos[0]+Vleft;
  		Rho=Pmax-pos[0]*(Pmax-Pmin)/L;//Pmax-pos[0]*(Pmax-Pmin)/L;
  		alpha=0;
  	}
	//Top side of the domain
  	else
  	{
  		U[0]=Utop;
  		U[1]=(Vright-Vleft)/L*pos[0]+Vleft;
  		Rho=Pmax-pos[0]*(Pmax-Pmin)/L;//Pmax-pos[0]*(Pmax-Pmin)/L;
  		alpha=0;
  	}
*/
}

void UserInit::UserIc (Parameters& PtrParameters, int elem, int nodenumber, double* pos ,double& Rho, double* U, double& alpha){
	double Umax,H,L,Pmax,Pmin;
	double Re,Ca,diameter, sigma;
	PtrParameters.Get_UserParameters(Umax,H,L,Pmax,Pmin);
	PtrParameters.Get_TwoPhaseUserParameters(Re,Ca,diameter, sigma);


		double Ubot=0;//-0.01;
		double Utop=0;//0.01;
		double Vleft=0;//-0.01;
		double Vright=0;//0.01;
		U[0]=(Utop-Ubot)/H*pos[1]+Ubot;//Umax*4.0*(pos[1]/H)*(1-(pos[1]/H));
		U[1]=(Vright-Vleft)/L*pos[0]+Vleft;
		Rho=Pmax-pos[0]*(Pmax-Pmin)/L;

		U[0]=0;
		U[1]=0;
		Rho=1;
		alpha=0;

		double R=diameter/2;
		if(pow(pos[0]-(L)/2,2.0)+pow(pos[1]-(H)/2,2.0)<=R*R+0.000001)
		{
			alpha=1;
			Rho=1;
		}
}

