/*
 * Boundaries.cpp
 *
 *  Created on: 17 Jul 2015
 *      Author: thomas
 */

#include "Boundaries.h"


Boundaries::Boundaries() {
	ru1=0;
	ru2=0;
	q23=2.0/3.0;
	q16=1.0/6.0;
	q13=1.0/3.0;
	rhodiff=0;
	SumWeightS=0;
	SumWeightE=0;
	SumWeightN=0;
	SumWeightW=0;
	SumWeightConcaveSE=0;
	SumWeightConcaveNE=0;
	SumWeightConcaveNW=0;
	SumWeightConcaveSW=0;
	SumWeightConvexSE=0;
	SumWeightConvexNE=0;
	SumWeightConvexNW=0;
	SumWeightConvexSW=0;
	/*MPI_Comm_rank(MPI_COMM_WORLD,&rank);

		char buffer[50]; // make sure it's big enough
		snprintf(buffer, sizeof(buffer), "Hezou_%d.txt", rank);
		myFlux.open(buffer);*/
}

Boundaries::~Boundaries() {
	// TODO Auto-generated destructor stub
}

///He Zou Boundary conditions exclude corners (Velocity imposed)
void Boundaries::FUNC_HeZou_U (double & a,double & b,double & c,double & d,double & e,double & f,double & g,double & h,double & i,double & U,double & V){
	double feq,Rho;
	Rho=(a+b+c+2.0*(d+e+f))/(1.0-U);
	ru1=Rho*U;
	ru2=Rho*V;
	feq=(b-c)*0.5;
	g=d+ q23*ru1;
	h=e+ru2*0.5+ru1*q16-feq;
	i=f-ru2*0.5+ru1*q16+feq;
	//Rhotmp=a+b+c+d+e+f+g+h+i;
	//std::cout<<"Rho Def:"<<Rho<<" Rho Cal: "<<Rhotmp<<" Ux: "<<U<<" Uy: "<<V <<std::endl;
}

/// He Zou Boundary conditions exclude corners (Pressure imposed)
void Boundaries::FUNC_HeZou_P (double & a,double & b,double & c,double & d,double & e,double & f,double & g,double & h,double & i,double & V,double & Rho){
	double feq,U;
	//std::cout<<copysign(1.0,Rho)<<std::endl;
	//U=copysign(1.0,Rho)*(1.0-(a+b+c+2.0*(d+e+f))/abs(Rho));
	U=copysign(1.0,Rho)*(1.0-((a+b+c)+2.0*(d+e+f))/std::abs(Rho));
	//Rho=abs(Rho);
	ru1=Rho*U;
	ru2=Rho*V;
	feq=(b-c)*0.5;
	g=d+ q23*ru1;
	h=e+ru2*0.5+ru1*q16-feq;
	i=f-ru2*0.5+ru1*q16+feq;
	/*Rhotmp=a+b+c+d+e+f+g+h+i;
	myFlux.precision(15);
	myFlux<<" Rho Def:\t"<<Rho<<"\tRho Cal:\t"<<Rhotmp<<"\tUx:\t"<<U<<"\tUy:\t"<<V <<std::endl;
	myFlux<<a<<"\t"<<b<<"\t"<<c<<"\t"<<d<<"\t"<<e<<"\t"<<f<<"\t"<<g<<"\t"<<h<<"\t"<<i<<std::endl;*/

}

/// Corner treat by Chih-Fung Ho, Cheng Chang, Kuen-Hau Lin and Chao-An Lin
/// Consistent Boundary Conditions for 2D and 3D Lattice Boltzmann Simulations
void Boundaries::FUNC_corner (double & a,double & b,double & c,double & d,double & e,double & f,double & g,double & h,double & i,double & U,double & V,double & Rho){
	ru1=Rho*U;
	ru2=Rho*V;
	b=g + q23*ru1;
	c=h + q23*ru2;
	d=i + ru1*q16 + ru2*q16;
	e=(Rho-a-ru1)*0.5-ru2*q13-(g + h +i);
	f=e + q16*ru1 - q16*ru2  ;
}
void Boundaries::FUNC_corner_no_vel_concave (double & a,double & b,double & c,double & d,double & e,double & f,double & g,double & h,double & i,double & Rho){
	b=g;
	c=h;
	d=i;
	e=(Rho-a)*0.5-(g + h +i);
	f=e;
}
void Boundaries::FUNC_corner_no_vel_convex (double & a,double & b,double & c,double & d,double & e,double & f){
	f=a-b+c-d+e;
}
/*void Boundaries::FUNC_corner_outer (double & a,double & b,double & c,double & d,double & e,double & f,double & g,double & h,double & i,double & U,double & V,double & Rho){
	ru1=Rho*U;
	ru2=Rho*V;
	b=g + q23*ru1;
	c=h + q23*ru2;
	d=i + ru1*q16 + ru2*q16;
	e=(Rho-a-ru1)*0.5-ru2*q13-(g + h +i);
	f=e + q16*ru1 - q16*ru2  ;
}*/

//     BC He Zou with U impose
void Boundaries::BC_HeZou_U(int & NormalBc, double* fi, double U, double V){

	double InvU,InvV,rho,ru,rv,fequi;
//	std::cout<<" Ux is: "<< U<<" Uy is: "<<V <<std::endl;
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

void Boundaries::BC_HeZou_P(int & NormalBc, double* fi, double Rho, double & U, double & V){


	double InvRho;

      switch (NormalBc)
      {

      case 4: //Top
    	  InvRho=-Rho;
    	  FUNC_HeZou_P(fi[0],fi[1],fi[3],fi[2],fi[6],fi[5],fi[4],fi[8],fi[7],U,InvRho);
    	  break;
      case 2: //Bot
    	  FUNC_HeZou_P(fi[0],fi[1],fi[3],fi[4],fi[7],fi[8],fi[2],fi[5],fi[6],U,Rho);
    	  break;
       case 1: //West
          FUNC_HeZou_P(fi[0],fi[2],fi[4],fi[3],fi[7],fi[6],fi[1],fi[5],fi[8],V,Rho);
          break;
      case 3: //East
    	  InvRho=-Rho;
    	  FUNC_HeZou_P(fi[0],fi[2],fi[4],fi[1],fi[8],fi[5],fi[3],fi[6],fi[7],V,InvRho);
          break;

      default :
          std::cout<<" Problem in the direction of HeZou unknown. Direction is: "<<NormalBc<<std::endl;//,node_HeZou%parameter%ID
      }
}

void Boundaries::BC_corner(int & NormalCorner, double* fi, double Rho, double U, double V){


	double InvRho,InvU,InvV;

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
void Boundaries::BC_corner_no_vel_concave(int & NormalCorner, double* fi, double Rho){
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
void Boundaries::BC_corner_no_vel_convex(int & NormalCorner, double* fi){
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
