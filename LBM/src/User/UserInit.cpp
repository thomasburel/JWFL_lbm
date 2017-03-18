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

	U[0]=0.0;//Umax*pos[1]/H;
	U[1]=0.0;
	//Rho=Pmin;
	Rho=1.0;//(Pmax-pos[0]*(Pmax-Pmin)/L);
	if(pos[0]<=0)
		Rho=Pmax;
	else
		Rho=Pmin;
	alpha=1.0;

	double R=30;//Diameter/2.0;
	if(pow(pos[0]-50,2.0)+pow(pos[1]-cos(teta)*(R+0.0001),2.0)<=R*R+1e-8)
	{
		alpha=0.0;
	}
//	if(pow(pos[0]-(L)/2.0,2.0)+pow(pos[1]+cos(teta)*R,2.0)<=R*R)
/*	if(pow(pos[0]-(L-10)/2.0,2.0)+pow(pos[1]-H/4,2.0)<=R*R)
	{
		alpha=1.0;
		Rho=Pmin;//+sigma*3.0/R;
	}*/
	if(pos[0]<=10)
	{
		U[0]=0.0;//1.e-4;//*(1.0-pow(2.0*(pos[1]-H/2.0)/H,2.0));
		U[1]=0.0;
		Rho=Pmax;
		alpha=0.0;
	}
	else
	{
		U[0]=0.0;
		U[1]=0.0;
		Rho=Pmax;//-pos[0]*(Pmax-Pmin)/L;
		alpha=1.0;
	}
	U[0]=0.0;
	if(pos[0]<=100)
	U[0]=1.e-2;
/*	if(pow(pos[0]-(L)/2.0,2.0)+pow(pos[1],2.0)<=R*R+1e-8)
	{
		alpha=1.0;
		Rho=Pmin+sigma*3.0/R;
	}
	if(pos[1]>=H-1)
		U[0]=Umax;*/
/*	double hin=H/2.0;
	double g=0.00000001;
	double U1=g*H*H/(8*(2.0*PtrParameters.Get_Tau_1()-1.0)/6.0);
	double U2=g*hin*hin/(8*(2.0*PtrParameters.Get_Tau_2()-1.0)/6.0);

	if(pos[1]<=hin/2.0 ||pos[1]>=H-hin/2.0)
	{
		U[0]=U1*(1.0-pow(2.0*(pos[1]-H/2.0)/H,2.0));//(pos[1]/H)*(1.0-(pos[1]/H));
		U[1]=0.0;
		Rho=Pmin;
		alpha=1.0;
		//std::cout<<U[0]<<std::endl;
	}
	else
	{
		U[0]=U2*(1.0-pow(2.0*(pos[1]-H/2.0)/hin,2.0))+U1*(1.0-pow(hin/H,2.0));
		U[1]=0.0;
		Rho=Pmin;
		alpha=0.0;
	}*/
}

void UserInit::UserIc (Parameters& PtrParameters, int elem, int nodenumber, double* pos ,double& Rho, double* U, double& alpha){
	double Umax,H,L,Pmax,Pmin;
	double teta,sigma,Diameter,Re,Ca;
	PtrParameters.Get_UserParameters(Umax,H,L,Pmax,Pmin);
	PtrParameters.Get_UserDroplets(teta,sigma,Diameter,Re,Ca);
	/*

	U[0]=0.0;//Umax*pos[1]/H;
		U[1]=0.0;
		//Rho=Pmin;
		Rho=1.0;//Pmax-pos[0]*(Pmax-Pmin)/L;
		alpha=1.0;

		double R=30;//Diameter/2.0;
		if(pow(pos[0]-50,2.0)+pow(pos[1]-cos(teta)*(R+0.0001),2.0)<=R*R+1e-8)
		{
			alpha=0.0;
		}*/

//		if(pow(pos[0]-(L)/2.0,2.0)+pow(pos[1]+cos(teta)*R,2.0)<=R*R)
	//	if(pow(pos[0]-(L-10)/2.0,2.0)+pow(pos[1]-H/4.0,2.0)<=R*R+1e-8)
/*		if(pow(pos[0]-(L)/2.0,2.0)+pow(pos[1],2.0)<=R*R+1e-8)
		{
			alpha=1.0;
			Rho=Pmin+sigma*3.0/R;
		}
*/
	if(pos[0]<=10)
	{
		U[0]=0.0;
		U[1]=0.0;
		Rho=Pmax;
		alpha=0.0;
	}
	else
	{
		U[0]=0.0;
		U[1]=0.0;
		Rho=Pmax;//-pos[0]*(Pmax-Pmin)/L;
		alpha=1.0;
	}

	U[0]=0.0;//1.e-4;
	if(pos[0]<=100)
	U[0]=1.e-2;
/*	double R=10;//Diameter/2.0;
	double xc,yc;
	double delatx,deltay;
	delatx=(L/10);
	deltay=(H/10);
	xc=round(pos[0]/delatx)*delatx;
	yc=round(pos[1]/deltay)*deltay;
	if(pow(pos[0]-xc,2.0)+pow(pos[1]-yc,2.0)<=R*R+1e-8)
	{
		alpha=1.0;
	}*/

/*	double hin=H/2.0;
	double g=0.00000001;
	double U1=g*H*H/(8*(2.0*PtrParameters.Get_Tau_1()-1.0)/6.0);
	double U2=g*hin*hin/(8*(2.0*PtrParameters.Get_Tau_2()-1.0)/6.0);
//	std::cout<<"U1: "<<U1<<" U2: "<<U2<<" hin: "<<hin<< " H: "<<H<<std::endl;
	if(pos[1]<=hin/2.0 ||pos[1]>=H-hin/2.0)
	{
		U[0]=U1*(1.0-pow(2.0*(pos[1]-H/2.0)/H,2.0));
		U[1]=0.0;
		Rho=Pmin;
		alpha=1.0;
		//std::cout<<"in: "<<U[0]<<std::endl;
	}
	else
	{
		U[0]=U2*(1.0-pow(2.0*(pos[1]-H/2.0)/hin,2.0))+U1*(1.0-pow(hin/H,2.0));
		U[1]=0.0;
		Rho=Pmin;
		alpha=0.0;
	}*/
}

