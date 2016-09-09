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
void D2Q9Velocity::SetVelocity(Parameters *Param){
//Setup the velocity assumptions
	switch(Param->Get_VelocityType())
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
//Setup the velocity model
	switch(Param->Get_VelocityModel())
	{
	case HeZouV:
		PtrVelocityMethod=&D2Q9Velocity::BC_HeZou_U;
		break;
	default:
		std::cerr<<"Velocity model has not been found."<<std::endl;
	}
}
void D2Q9Velocity::ApplyVelocity(int const &BcNormal,int const *Connect, double const *UDef, DistriFunct * & f_in, double *Rho, double *U, double *V){
	(this->*PtrCalculU)(BcNormal,Connect, UDef,U,V);
	(this->*PtrVelocityMethod)(BcNormal,Connect, UDef, f_in, U[Connect[0]], V[Connect[0]]);
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
void D2Q9Velocity::BC_HeZou_U(int const &BcNormal,int const *Connect, double const *UDef, DistriFunct * & f_in, double & U, double & V){

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
