/*
 * ============================================================================
 * D2Q9Pressure.cpp
 *
 *  Created on: 29 Aug 2016
 *      Author: Thomas Burel
 *      Version     :
 *      Copyright   : Your copyright notice
 *      Description : LBM Code Developing by JWFL at University of Strathclyde
 * ============================================================================
 */

#include "D2Q9Pressure.h"

D2Q9Pressure::D2Q9Pressure() {
	InvRho=0;
	feq=0;
	U=0;;
	PtrPressureMethod=0;
	PtrCalculRho=0;
}

D2Q9Pressure::~D2Q9Pressure() {
	// TODO Auto-generated destructor stub
}

void D2Q9Pressure::Set_PressureBcs(Parameters *Param, double ** &Ei){
//Setup the pressure assumptions and the model
	SetPressure(Param->Get_PressureModel(),Param->Get_PressureType());
	EiBc=Ei;
}
void D2Q9Pressure::SetPressure(PressureModel PressureModel_,PressureType PressureType_){
	//Setup the Pressure model
	SetPressureModel(PressureModel_);
	//Setup the pressure assumptions
	SetPressureType(PressureType_);
}
void D2Q9Pressure::SetPressureModel(PressureModel PressureModel_){
	//Setup the Pressure model
		switch(PressureModel_)
		{
		case HeZouP:
			PtrPressureMethod=&D2Q9Pressure::BC_HeZou_P;
			break;
		default:
			std::cerr<<"Pressure model has not been found."<<std::endl;
		}
}
void D2Q9Pressure::SetPressureType(PressureType PressureType_){
	//Setup the pressure assumptions
		switch(PressureType_)
		{
		case FixP:
			PtrCalculRho=&D2Q9Pressure::FixRho;
			break;
		case zeroPGrad1st:
			PtrCalculRho=&D2Q9Pressure::NoGradRho_1stOrder;
			break;
		default:
			std::cerr<<"Pressure Type has not been found."<<std::endl;
		}
}
void D2Q9Pressure::ApplyPressure(int const &BcNormal,int const *Connect, double const &Rho_def, DistriFunct* f_in, double *Rho, double *U, double *V){
	(this->*PtrCalculRho)(BcNormal,Connect, Rho_def,Rho);
	(this->*PtrPressureMethod)(BcNormal,Connect, Rho_def, f_in, Rho[Connect[0]], U[Connect[0]], V[Connect[0]]);
}
void D2Q9Pressure::FixRho(int const &BcNormal,int const *Connect, double const &Rho_def, double *Rho){
	Rho[Connect[0]]=Rho_def;
}
void D2Q9Pressure::NoGradRho_1stOrder(int const &BcNormal,int const *Connect, double const &Rho_def, double *Rho){
	Rho[Connect[0]]=Rho[Connect[BcNormal]];
}
/// He Zou Boundary conditions exclude corners (Pressure imposed)
void D2Q9Pressure::FUNC_HeZou_P (double & a,double & b,double & c,double & d,double & e,double & f,double & g,double & h,double & i,double & V,double & Rho){
		U=copysign(1.0,Rho)*(1.0-((a+b+c)+2.0*(d+e+f))/std::abs(Rho));
		ru1=Rho*U;
		ru2=Rho*V;
		feq=(b-c)*0.5;
		g=d+ q23*ru1;
		h=e+ru2*0.5+ru1*q16-feq;
		i=f-ru2*0.5+ru1*q16+feq;
}
void D2Q9Pressure::BC_HeZou_P(int const &BcNormal,int const *Connect, double const &Rho_def, DistriFunct* f_in, double Rho, double & U, double & V){
    if(Rho==0)
    {
    	for (int i=0;i<9;i++) f_in->f[i][Connect[0]]=0;
    }
    else
		switch (BcNormal)
		  {
		  case 4: //Top
			  InvRho=-Rho;
			  FUNC_HeZou_P(f_in->f[0][Connect[0]],f_in->f[1][Connect[0]],f_in->f[3][Connect[0]],f_in->f[2][Connect[0]],f_in->f[6][Connect[0]],f_in->f[5][Connect[0]],f_in->f[4][Connect[0]],f_in->f[8][Connect[0]],f_in->f[7][Connect[0]],U,InvRho);
			  break;
		  case 2: //Bot
			  FUNC_HeZou_P(f_in->f[0][Connect[0]],f_in->f[1][Connect[0]],f_in->f[3][Connect[0]],f_in->f[4][Connect[0]],f_in->f[7][Connect[0]],f_in->f[8][Connect[0]],f_in->f[2][Connect[0]],f_in->f[5][Connect[0]],f_in->f[6][Connect[0]],U,Rho);
			  break;
		   case 1: //West
			  FUNC_HeZou_P(f_in->f[0][Connect[0]],f_in->f[2][Connect[0]],f_in->f[4][Connect[0]],f_in->f[3][Connect[0]],f_in->f[7][Connect[0]],f_in->f[6][Connect[0]],f_in->f[1][Connect[0]],f_in->f[5][Connect[0]],f_in->f[8][Connect[0]],V,Rho);
			  break;
		  case 3: //East
			  InvRho=-Rho;
			  FUNC_HeZou_P(f_in->f[0][Connect[0]],f_in->f[2][Connect[0]],f_in->f[4][Connect[0]],f_in->f[1][Connect[0]],f_in->f[8][Connect[0]],f_in->f[5][Connect[0]],f_in->f[3][Connect[0]],f_in->f[6][Connect[0]],f_in->f[7][Connect[0]],V,InvRho);
			  break;
		  default :
			  std::cout<<" Problem in the direction of HeZou unknown. Direction is: "<<BcNormal<<std::endl;
		  }
}
