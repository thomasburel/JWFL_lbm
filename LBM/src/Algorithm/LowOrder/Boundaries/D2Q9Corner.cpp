/*
 * ============================================================================
 * D2Q9Corner.cpp
 *
 *  Created on: 29 Aug 2016
 *      Author: Thomas Burel
 *      Version     :
 *      Copyright   : Your copyright notice
 *      Description : LBM Code Developing by JWFL at University of Strathclyde
 * ============================================================================
 */

#include "D2Q9Corner.h"

D2Q9Corner::D2Q9Corner() {
	InvRho=0;
	InvU=0;
	InvV=0;
}

D2Q9Corner::~D2Q9Corner() {
	// TODO Auto-generated destructor stub
}

/// Corner treat by Chih-Fung Ho, Cheng Chang, Kuen-Hau Lin and Chao-An Lin
/// Consistent Boundary Conditions for 2D and 3D Lattice Boltzmann Simulations
void D2Q9Corner::FUNC_corner (double & a,double & b,double & c,double & d,double & e,double & f,double & g,double & h,double & i,double & U,double & V,double & Rho){
	ru1=Rho*U;
	ru2=Rho*V;
	b=g + q23*ru1;
	c=h + q23*ru2;
	d=i + ru1*q16 + ru2*q16;
	e=(Rho-a-ru1)*0.5-ru2*q13-(g + h +i);
	f=e + q16*ru1 - q16*ru2  ;
}
void D2Q9Corner::FUNC_corner_no_vel_concave (double & a,double & b,double & c,double & d,double & e,double & f,double & g,double & h,double & i,double & Rho){
	b=g;
	c=h;
	d=i;
	e=(Rho-a)*0.5-(g + h +i);
	f=e;
}
void D2Q9Corner::FUNC_corner_no_vel_convex (double & a,double & b,double & c,double & d,double & e,double & f){
	f=a-b+c-d+e;
}
void D2Q9Corner::BC_corner(int & NormalCorner, double* fi, double Rho, double U, double V){

      switch (NormalCorner)
      {

      case 8: //Top West
    	  InvV=-V;
    	  FUNC_corner(fi[0],fi[1],fi[4],fi[8],fi[7],fi[5],fi[3],fi[2],fi[6],U,InvV,Rho);//node_corner%node%fi(1),node_corner%node%fi(2),node_corner%node%fi(5), &
   	  break;
      case 5: //Bottom West
    	  FUNC_corner(fi[0],fi[1],fi[2],fi[5],fi[6],fi[8],fi[3],fi[4],fi[7],U,V,Rho);
    	  break;
      case 7 : //Top East
    	  InvU=-U;
    	  InvV=-V;
    	  FUNC_corner(fi[0],fi[3],fi[4],fi[7],fi[8],fi[6],fi[1],fi[2],fi[5],InvU,InvV,Rho);
    	  break;
      case 6: //Bottom East
    	  InvU=-U;
    	  FUNC_corner(fi[0],fi[3],fi[2],fi[6],fi[5],fi[7],fi[1],fi[4],fi[8],InvU,V,Rho);
    	  break;
      default :
          std::cout<<" Problem in the direction of HeZou Corner unknown. Direction is: "<<NormalCorner<<std::endl;//,node_HeZou%parameter%ID
      }
}
void D2Q9Corner::BC_corner_no_vel_concave(int & NormalCorner, double* fi, double Rho){
    switch (NormalCorner)
      {

      case 8: //Top West

    	  FUNC_corner_no_vel_concave(fi[0],fi[1],fi[4],fi[8],fi[7],fi[5],fi[3],fi[2],fi[6],Rho);//node_corner%node%fi(1),node_corner%node%fi(2),node_corner%node%fi(5), &
   	  break;
      case 5: //Bottom West
    	  FUNC_corner_no_vel_concave(fi[0],fi[1],fi[2],fi[5],fi[6],fi[8],fi[3],fi[4],fi[7],Rho);
    	  break;
      case 7 : //Top East

    	  FUNC_corner_no_vel_concave(fi[0],fi[3],fi[4],fi[7],fi[8],fi[6],fi[1],fi[2],fi[5],Rho);
    	  break;
      case 6: //Bottom East

    	  FUNC_corner_no_vel_concave(fi[0],fi[3],fi[2],fi[6],fi[5],fi[7],fi[1],fi[4],fi[8],Rho);
    	  break;
      default :
          std::cout<<" Problem in the direction of concave Corner unknown. Direction is: "<<NormalCorner<<std::endl;//,node_HeZou%parameter%ID
      }
}
void D2Q9Corner::BC_corner_no_vel_convex(int & NormalCorner, double* fi){
    switch (NormalCorner)
      {

      case 6: //Top West
    	  FUNC_corner_no_vel_convex(fi[1],fi[3],fi[5],fi[7],fi[8],fi[6]);
   	  break;
      case 7: //Bottom West
    	  FUNC_corner_no_vel_convex(fi[1],fi[3],fi[5],fi[6],fi[8],fi[7]);
    	  break;
      case 5 : //Top East
    	  FUNC_corner_no_vel_convex(fi[3],fi[1],fi[7],fi[8],fi[6],fi[5]);
    	  break;
      case 8: //Bottom East
    	  FUNC_corner_no_vel_convex(fi[3],fi[1],fi[6],fi[5],fi[7],fi[8]);
    	  break;
      default :
          std::cout<<" Problem in the direction of convex corner unknown. Direction is: "<<NormalCorner<<std::endl;//,node_HeZou%parameter%ID
      }
}
