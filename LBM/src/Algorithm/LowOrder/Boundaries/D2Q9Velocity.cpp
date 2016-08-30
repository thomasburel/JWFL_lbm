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

}

D2Q9Velocity::~D2Q9Velocity() {
	// TODO Auto-generated destructor stub
}

///He Zou Boundary conditions exclude corners (Velocity imposed)
void D2Q9Velocity::FUNC_HeZou_U (double & a,double & b,double & c,double & d,double & e,double & f,double & g,double & h,double & i,double & U,double & V){
	double feq,Rho;
	Rho=(a+b+c+2.0*(d+e+f))/(1.0-U);
	ru1=Rho*U;
	ru2=Rho*V;
	feq=(b-c)*0.5;
	g=d+ q23*ru1;
	h=e+ru2*0.5+ru1*q16-feq;
	i=f-ru2*0.5+ru1*q16+feq;
}
//     BC He Zou with U impose
void D2Q9Velocity::BC_HeZou_U(int & NormalBc, double* fi, double U, double V){

	double InvU,InvV,rho,ru,rv,fequi;
      switch (NormalBc)
      {
      case 4: //Top
    	  InvV=-V;
          FUNC_HeZou_U(fi[0],fi[1],fi[3],fi[2],fi[6],fi[5],fi[4],fi[8],fi[7],InvV,U);
          break;
      case 2: //Bot
    	  FUNC_HeZou_U(fi[0],fi[1],fi[3],fi[4],fi[7],fi[8],fi[2],fi[5],fi[6],V,U);
    	  break;
      case 1: //West
    	  FUNC_HeZou_U(fi[0],fi[2],fi[4],fi[3],fi[7],fi[6],fi[1],fi[5],fi[8],U,V);
    	  break;
      case 3://East
    	  InvU=-U;
    	  FUNC_HeZou_U(fi[0],fi[2],fi[4],fi[1],fi[8],fi[5],fi[3],fi[6],fi[7],InvU,V);

    	  break;
      default :
          std::cout<<" Problem in the direction of HeZou unknown. Direction is: "<<NormalBc<<std::endl;//,node_HeZou%parameter%ID
          break;
      }
}
