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

		//if (pos[1]>0.1 && pos[1]<H-0.1)
		//U[0]=Umax;//*4.0*(pos[1]/H)*(1-(pos[1]/H));
	//else
		//U[0]=0.0;
		U[0]=0.0;//002;
	U[1]=0.0;
	//Rho=Pmax-pos[0]*(Pmax-Pmin)/L;
	//Rho=(Pmin+Pmax)/2.0;
	//Rho=Pmax;
	//Rho=Pmin;
	Rho=Pmin;

/*	if (pos[0]==L)
		Rho=Pmin;
	if (pos[0]==0)
		Rho=Pmax;
*/
//	Rho=1;
/*	if(pos[0]<=3)
		Rho=1.1;
	else
		Rho=0.9;*/
	//Rho=Pmax;

	//Rho=0;
	//if (pos[0]==L)
	//	Rho=Pmax;
	//Rho=1.0;//Pmax-pos[0]*(Pmax-Pmin)/H;
	/*double coef=0.5;
	if(pos[0]<L/2)
	{
		if(pos[1]>H/2)
			Rho=Pmax+coef*(Pmax-Pmin);
		else
			Rho=Pmax-coef*(Pmax-Pmin);

		if(pos[1]==0 && pos[0]==0) std::cout<<"Pressure inlet downstream"<<Rho<<std::endl;
		if(pos[1]==H && pos[0]==0) std::cout<<"Pressure inlet upstream"<<Rho<<std::endl;
	}
	else
		Rho=Pmin;*/

//*********** Poiseuille ini***************
//Left side of the domain
/*	if(pos[0]<=0)
	// Global Corner at the bottom left side
 		if(pos[1]<=0)
  		{
  			U[0]=0;
  			U[1]=0;
  			Rho=Pmax;
  		}
  		else
 	// Global Corner at the Top left side
 		if(pos[1]>=H)
  		{
  			U[0]=0;
  			U[1]=0;
  			Rho=Pmax;
  		}
  	//Left side of the domain and excluding the two Global Corners
  		else
  		{
   			U[0]=Umax*4.0*(pos[1]/H)*(1-(pos[1]/H));;
  			U[1]=0;
  			Rho=Pmax;
 		}
	else
//Right side of the domain
 	if(pos[0]>=L)
	// Global Corner at the bottom right side
  		if(pos[1]<=0)
  		{
  			U[0]=0;
  			U[1]=0;
  			Rho=Pmin;
  		}
  		else
  	// Global Corner at the Top right side
  		if(pos[1]>=H)
  		{
  			U[0]=0;
  			U[1]=0;
  			Rho=Pmin;
  		}
  	//Right side of the domain and excluding the two Global Corner
  		else
  		{
    		U[0]=Umax*4.0*(pos[1]/H)*(1-(pos[1]/H));;
  			U[1]=0;
  			Rho=Pmin;
  		}
  	else
//Bottom side of the domain
 	if(pos[1]<=0)
  	{
  		U[0]=0;
  		U[1]=0;
  		Rho=Pmax-pos[0]*(Pmax-Pmin)/L
  	}
	//Top side of the domain
  	else
  	{
  		U[0]=0;
  		U[1]=0;
  		Rho=Pmax-pos[0]*(Pmax-Pmin)/L
  	}
*/

//*********** Driven Cavity***************
/*	if(pos[1]>0.1 )
			U[0]=Umax;
		else
			U[0]=0;
		U[1]=0.0;
		Rho=Pmin;*/

	Rho=PtrParameters.Get_Rho_2();
	alpha=0;
}

void UserInit::UserIc (Parameters& PtrParameters, int elem, int nodenumber, double* pos ,double& Rho, double* U, double& alpha){
	double Umax,H,L,Pmax,Pmin;
	PtrParameters.Get_UserParameters(Umax,H,L,Pmax,Pmin);

	//U[0]=Umax;//*4.0*(pos[1]/H)*(1-(pos[1]/H));///pos[0];
	//U[0]=0.0;
	U[0]=0.00;//02;
	U[1]=0.0;
	//Rho=(Pmax+Pmin)/2.0;
	Rho=Pmin;
	int dx=15;
	if((pos[0]>=L/2-dx&&pos[0]<=L/2+dx)&&(pos[1]>=H/2-dx&&pos[1]<=H/2+dx))
	{
		alpha=1;
		Rho=PtrParameters.Get_Rho_1();
//		std::cout<<" Node number fluid 1 are: "<<nodenumber<<" x: "<<pos[0]<<" y: "<<pos[1]<<std::endl;
	}
	else
	{
		alpha=0;
		Rho=PtrParameters.Get_Rho_2();
	}

//	Rho=1;
	//Rho=Pmax-pos[0]*(Pmax-Pmin)/L;

/*	if(pos[0]<=3)
		Rho=1.1;
	else
		Rho=0.9;*/
	//Rho=Pmax;

	//U[0]=Umax;
/*	U[0]=0.0;
	U[1]=0.0;
	Rho=1;*/
	//Rho=1;
	//U[0]=Umax*4.0*(pos[1]/H)*(1-(pos[1]/H));
}

