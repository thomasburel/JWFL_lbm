/*
 * ============================================================================
 * D2Q9Velocity.cpp
 *
 *  Created on: 29 Aug 2016
 *      Author: Thomas Burel
 *      Version     :
 *      Copyright   : Your copyright notice
 *      Description : LBM Code Developing by JWFL at University of Strathclyde
 * ============================================================================
 */

#include "D2Q9Velocity.h"

D2Q9Velocity::D2Q9Velocity() {
	PtrVelocityMethod=0;
	PtrCalculU=0;
	feq=0;
	Rho=0;
	InvU=0;
	InvV=0;
	rho=0;
	ru=0;
	rv=0;
	fequi=0;
}

D2Q9Velocity::~D2Q9Velocity() {
	// TODO Auto-generated destructor stub
}

void D2Q9Velocity::SetVelocity(Dictionary *PtrDic,NodeArrays2D* NodeArrays, Parameters *Param, double ** &Ei){
//Setup the velocity assumptions
	SetVelocityModel(Param->Get_VelocityModel());
//Setup the velocity model
	SetVelocityType(Param->Get_VelocityType());

	EiBc=Ei;
}
void D2Q9Velocity::SetVelocityModel(VelocityModel VelocityModel_){
	//Setup the velocity model
		switch(VelocityModel_)
		{
		case HeZouV:
			PtrVelocityMethod=&D2Q9Velocity::BC_HeZou_U;
			break;
		case Ladd:
			PtrVelocityMethod=&D2Q9Velocity::BC_Ladd;
			break;
		default:
			std::cerr<<"Velocity model has not been found."<<std::endl;
		}
}
void D2Q9Velocity::SetVelocityType(VelocityType VelocityType_){
	//Setup the velocity assumptions
		switch(VelocityType_)
		{
		case FixV:
			PtrCalculU=&D2Q9Velocity::FixU;
			break;
		case zeroVGrad1st:
			PtrCalculU=&D2Q9Velocity::NoGradU_1stOrder;
			break;
		default:
			std::cerr<<"Velocity Type has not been found."<<std::endl;
		}
}
void D2Q9Velocity::SetVelocity(VelocityModel VelocityModel_,VelocityType VelocityType_){
	//Setup the velocity assumptions
		SetVelocityModel(VelocityModel_);
	//Setup the velocity model
		SetVelocityType(VelocityType_);
}
void D2Q9Velocity::ApplyVelocity(int const &BcNormal,int const *Connect, double const *UDef, DistriFunct * & f_in, double *Rho, double *U, double *V){
	(this->*PtrCalculU)(BcNormal,Connect, UDef,U,V);
	(this->*PtrVelocityMethod)(BcNormal,Connect, f_in,Rho[Connect[0]], U[Connect[0]], V[Connect[0]]);
}
void D2Q9Velocity::FixU(int const &BcNormal,int const *Connect, double const *UDef, double *U, double *V){
	U[Connect[0]]=UDef[0];
	V[Connect[0]]=UDef[1];
}
void D2Q9Velocity::NoGradU_1stOrder(int const &BcNormal,int const *Connect, double const *UDef, double *U, double *V){
	switch(BcNormal)
	{
		case 1:
			U[Connect[0]]=U[Connect[1]];
			V[Connect[0]]=UDef[1];
			break;
		case 2:
			U[Connect[0]]=UDef[0];
			V[Connect[0]]=V[Connect[2]];
			break;
		case 3:
			U[Connect[0]]=U[Connect[3]];
			V[Connect[0]]=UDef[1];
			break;
		case 4:
			U[Connect[0]]=UDef[0];
			V[Connect[0]]=V[Connect[4]];
			break;
		default:
			std::cerr<<"Bc normal choice for Zero Grad U boundary condition not detected."<<std::endl;
	}
}
///He Zou Boundary conditions exclude corners (Velocity imposed)
void D2Q9Velocity::FUNC_HeZou_U (double & a,double & b,double & c,double & d,double & e,double & f,double & g,double & h,double & i,double & U,double & V){

	Rho=(a+b+c+2.0*(d+e+f))/(1.0-U);
	ru1=Rho*U;
	ru2=Rho*V;
	feq=(b-c)*0.5;
	g=d+ q23*ru1;
	h=e+ru2*0.5+ru1*q16-feq;
	i=f-ru2*0.5+ru1*q16+feq;
}
//     BC He Zou with U impose
void D2Q9Velocity::BC_HeZou_U(int const &BcNormal,int const *Connect, DistriFunct * & f_in,double & Rho, double & U, double & V){

      switch (BcNormal)
      {
      case 4: //Top
    	  InvV=-V;
          FUNC_HeZou_U(f_in->f[0][Connect[0]],f_in->f[1][Connect[0]],f_in->f[3][Connect[0]],f_in->f[2][Connect[0]],f_in->f[6][Connect[0]],f_in->f[5][Connect[0]],f_in->f[4][Connect[0]],f_in->f[8][Connect[0]],f_in->f[7][Connect[0]],InvV,U);
          break;
      case 2: //Bot
    	  FUNC_HeZou_U(f_in->f[0][Connect[0]],f_in->f[1][Connect[0]],f_in->f[3][Connect[0]],f_in->f[4][Connect[0]],f_in->f[7][Connect[0]],f_in->f[8][Connect[0]],f_in->f[2][Connect[0]],f_in->f[5][Connect[0]],f_in->f[6][Connect[0]],V,U);
    	  break;
      case 1: //West
    	  FUNC_HeZou_U(f_in->f[0][Connect[0]],f_in->f[2][Connect[0]],f_in->f[4][Connect[0]],f_in->f[3][Connect[0]],f_in->f[7][Connect[0]],f_in->f[6][Connect[0]],f_in->f[1][Connect[0]],f_in->f[5][Connect[0]],f_in->f[8][Connect[0]],U,V);
    	  break;
      case 3://East
    	  InvU=-U;
    	  FUNC_HeZou_U(f_in->f[0][Connect[0]],f_in->f[2][Connect[0]],f_in->f[4][Connect[0]],f_in->f[1][Connect[0]],f_in->f[8][Connect[0]],f_in->f[5][Connect[0]],f_in->f[3][Connect[0]],f_in->f[6][Connect[0]],f_in->f[7][Connect[0]],InvU,V);

    	  break;
      default :
          std::cout<<" Problem in the direction of HeZou unknown. Direction is: "<<BcNormal<<std::endl;
          break;
      }
}
//     BC moving wall with half way bounce-back
void D2Q9Velocity::BC_Ladd(int const &BcNormal,int const *Connect, DistriFunct * & f_in,double & Rho, double & U, double & V){

	switch(BcNormal)
			{
			case 2:
				f_in->f[2][Connect[0]]=f_in->f[OppositeBc[2]][Connect[0]]-6*omegaBc[OppositeBc[2]]*Rho*(EiBc[OppositeBc[2]][1]*V);
				f_in->f[5][Connect[0]]=f_in->f[OppositeBc[5]][Connect[0]]-6*omegaBc[OppositeBc[5]]*Rho*(EiBc[OppositeBc[5]][0]*U+EiBc[OppositeBc[5]][1]*V);
				f_in->f[6][Connect[0]]=f_in->f[OppositeBc[6]][Connect[0]]-6*omegaBc[OppositeBc[6]]*Rho*(EiBc[OppositeBc[6]][0]*U+EiBc[OppositeBc[6]][1]*V);
				break;
			case 4:
				f_in->f[4][Connect[0]]=f_in->f[OppositeBc[4]][Connect[0]]-6*omegaBc[OppositeBc[4]]*Rho*(EiBc[OppositeBc[4]][1]*V);
				f_in->f[7][Connect[0]]=f_in->f[OppositeBc[7]][Connect[0]]-6*omegaBc[OppositeBc[7]]*Rho*(EiBc[OppositeBc[7]][0]*U+EiBc[OppositeBc[7]][1]*V);
				f_in->f[8][Connect[0]]=f_in->f[OppositeBc[8]][Connect[0]]-6*omegaBc[OppositeBc[8]]*Rho*(EiBc[OppositeBc[8]][0]*U+EiBc[OppositeBc[8]][1]*V);
				break;
			case 1:
				f_in->f[1][Connect[0]]=f_in->f[OppositeBc[1]][Connect[0]]-6*omegaBc[OppositeBc[1]]*Rho*(EiBc[OppositeBc[1]][0]*U);
				f_in->f[5][Connect[0]]=f_in->f[OppositeBc[5]][Connect[0]]-6*omegaBc[OppositeBc[5]]*Rho*(EiBc[OppositeBc[5]][0]*U+EiBc[OppositeBc[5]][1]*V);
				f_in->f[8][Connect[0]]=f_in->f[OppositeBc[8]][Connect[0]]-6*omegaBc[OppositeBc[8]]*Rho*(EiBc[OppositeBc[8]][0]*U+EiBc[OppositeBc[8]][1]*V);
				break;
			case 3:
				f_in->f[3][Connect[0]]=f_in->f[OppositeBc[3]][Connect[0]]-6*omegaBc[OppositeBc[3]]*Rho*(EiBc[OppositeBc[3]][0]*U);
				f_in->f[6][Connect[0]]=f_in->f[OppositeBc[6]][Connect[0]]-6*omegaBc[OppositeBc[6]]*Rho*(EiBc[OppositeBc[6]][0]*U+EiBc[OppositeBc[6]][1]*V);
				f_in->f[7][Connect[0]]=f_in->f[OppositeBc[7]][Connect[0]]-6*omegaBc[OppositeBc[7]]*Rho*(EiBc[OppositeBc[7]][0]*U+EiBc[OppositeBc[7]][1]*V);
				break;
			default:
				std::cerr<<"Direction moving wall half-way bounce back not found. Node Index is: "<<Connect[0]<<" and direction is: "<<BcNormal<<std::endl;
				break;
			}
}
