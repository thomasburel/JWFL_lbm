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
	PtrCornerMethod=0;
	PtrCornerWallMethod=0;
	PtrCalculRhoCornerWall=0;

	doubleTmpReturn=0;
	direction1=0;
	direction2=0;
}

D2Q9Corner::~D2Q9Corner() {
	// TODO Auto-generated destructor stub
}
void D2Q9Corner::Set_Corner(Parameters *Param){
	PtrCornerMethod=&D2Q9Corner::ApplyHoChan;
	switch(Param->Get_WallType())
	{
	case Diffuse:
		PtrCornerWallMethod=&D2Q9Corner::ApplyDiffuseWall;
		break;
	case BounceBack:
		PtrCornerWallMethod=&D2Q9Corner::ApplyBounceBack;
		break;
	case HeZouWall:
		PtrCornerWallMethod=&D2Q9Corner::ApplyHoChanNoVel;
		break;
	default:
		std::cerr<<"Wall node type not found."<<std::endl;
	}
	switch(Param->Get_CornerPressureType())
	{
	case FixCP:
		PtrCalculRhoCornerWall=&D2Q9Corner::FixRho;
		break;
	case ExtrapolCP:
		PtrCalculRhoCornerWall=&D2Q9Corner::ExtrapolationAvgRho;
		break;
	default:
		std::cerr<<"Corner pressure type not found."<<std::endl;
	}

}
void D2Q9Corner::ApplyCorner(NodeCorner2D& Node, DistriFunct* f_in,double const & RhoDef,double const & UDef,double const & VDef, double *Rho, double *U, double *V){
	Rho[Node.Get_index()]=RhoDef;
	U[Node.Get_index()]=UDef;
	V[Node.Get_index()]=VDef;
	(this->*PtrCornerMethod)(Node, f_in, Rho[Node.Get_index()], U[Node.Get_index()], V[Node.Get_index()]);
}
void D2Q9Corner::ApplyCornerWall(NodeCorner2D& Node, DistriFunct* f_in, double *Rho, double *U, double *V){
	(this->*PtrCornerWallMethod)(Node, f_in, (this->*PtrCalculRhoCornerWall)(Node, Rho), U[Node.Get_index()], V[Node.Get_index()]);
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
void D2Q9Corner::FUNC_corner_no_vel (double & a,double & b,double & c,double & d,double & e,double & f,double & g,double & h,double & i,double & Rho){
	b=g;
	c=h;
	d=i;
	e=(Rho-a)*0.5-(g + h +i);
	f=e;
}

void D2Q9Corner::ApplyHoChan(NodeCorner2D& Node, DistriFunct* f_in, double Rho, double U, double V){

      switch (Node.Get_BcNormal())
      {

      case 8: //Top West
    	  InvV=-V;
    	  FUNC_corner(f_in->f[0][Node.Get_index()],f_in->f[1][Node.Get_index()],f_in->f[4][Node.Get_index()],f_in->f[8][Node.Get_index()],f_in->f[7][Node.Get_index()],f_in->f[5][Node.Get_index()],f_in->f[3][Node.Get_index()],f_in->f[2][Node.Get_index()],f_in->f[6][Node.Get_index()],U,InvV,Rho);
    	  break;
      case 5: //Bottom West
    	  FUNC_corner(f_in->f[0][Node.Get_index()],f_in->f[1][Node.Get_index()],f_in->f[2][Node.Get_index()],f_in->f[5][Node.Get_index()],f_in->f[6][Node.Get_index()],f_in->f[8][Node.Get_index()],f_in->f[3][Node.Get_index()],f_in->f[4][Node.Get_index()],f_in->f[7][Node.Get_index()],U,V,Rho);
    	  break;
      case 7 : //Top East
    	  InvU=-U;
    	  InvV=-V;
    	  FUNC_corner(f_in->f[0][Node.Get_index()],f_in->f[3][Node.Get_index()],f_in->f[4][Node.Get_index()],f_in->f[7][Node.Get_index()],f_in->f[8][Node.Get_index()],f_in->f[6][Node.Get_index()],f_in->f[1][Node.Get_index()],f_in->f[2][Node.Get_index()],f_in->f[5][Node.Get_index()],InvU,InvV,Rho);
    	  break;
      case 6: //Bottom East
    	  InvU=-U;
    	  FUNC_corner(f_in->f[0][Node.Get_index()],f_in->f[3][Node.Get_index()],f_in->f[2][Node.Get_index()],f_in->f[6][Node.Get_index()],f_in->f[5][Node.Get_index()],f_in->f[7][Node.Get_index()],f_in->f[1][Node.Get_index()],f_in->f[4][Node.Get_index()],f_in->f[8][Node.Get_index()],InvU,V,Rho);
    	  break;

      default :
          std::cout<<" Problem in the direction of HeZou Corner unknown. Direction is: "<<Node.Get_BcNormal()<<std::endl;
      }
}
void D2Q9Corner::ApplyHoChanNoVel(NodeCorner2D& Node, DistriFunct* f_in, double Rho, double U, double V){
    switch (Node.Get_BcNormal())
      {

      case 8: //Top West

    	  FUNC_corner_no_vel(f_in->f[0][Node.Get_index()],f_in->f[1][Node.Get_index()],f_in->f[4][Node.Get_index()],f_in->f[8][Node.Get_index()],f_in->f[7][Node.Get_index()],f_in->f[5][Node.Get_index()],f_in->f[3][Node.Get_index()],f_in->f[2][Node.Get_index()],f_in->f[6][Node.Get_index()],Rho);
   	  break;
      case 5: //Bottom West
    	  FUNC_corner_no_vel(f_in->f[0][Node.Get_index()],f_in->f[1][Node.Get_index()],f_in->f[2][Node.Get_index()],f_in->f[5][Node.Get_index()],f_in->f[6][Node.Get_index()],f_in->f[8][Node.Get_index()],f_in->f[3][Node.Get_index()],f_in->f[4][Node.Get_index()],f_in->f[7][Node.Get_index()],Rho);
    	  break;
      case 7 : //Top East

    	  FUNC_corner_no_vel(f_in->f[0][Node.Get_index()],f_in->f[3][Node.Get_index()],f_in->f[4][Node.Get_index()],f_in->f[7][Node.Get_index()],f_in->f[8][Node.Get_index()],f_in->f[6][Node.Get_index()],f_in->f[1][Node.Get_index()],f_in->f[2][Node.Get_index()],f_in->f[5][Node.Get_index()],Rho);
    	  break;
      case 6: //Bottom East

    	  FUNC_corner_no_vel(f_in->f[0][Node.Get_index()],f_in->f[3][Node.Get_index()],f_in->f[2][Node.Get_index()],f_in->f[6][Node.Get_index()],f_in->f[5][Node.Get_index()],f_in->f[7][Node.Get_index()],f_in->f[1][Node.Get_index()],f_in->f[4][Node.Get_index()],f_in->f[8][Node.Get_index()],Rho);
    	  break;
      default :
          std::cout<<"Corner direction is unknown. Direction is: "<<Node.Get_BcNormal()<<std::endl;
      }
}
void D2Q9Corner::ApplyBounceBack(NodeCorner2D& Node, DistriFunct* f_in, double Rho, double U, double V){
	double tmp=0;
	switch(Node.Get_BcNormal())
			{
			case 5:
				f_in->f[5][Node.Get_index()]=f_in->f[OppositeBc[5]][Node.Get_index()];
				f_in->f[1][Node.Get_index()]=f_in->f[OppositeBc[1]][Node.Get_index()];
				f_in->f[2][Node.Get_index()]=f_in->f[OppositeBc[2]][Node.Get_index()];
				if (Node.stream()[3]==false)
				{

					double sumfi=f_in->f[0][Node.Get_index()]+f_in->f[1][Node.Get_index()]+f_in->f[2][Node.Get_index()]+f_in->f[3][Node.Get_index()]+f_in->f[4][Node.Get_index()]+f_in->f[5][Node.Get_index()]+f_in->f[7][Node.Get_index()];
					f_in->f[6][Node.Get_index()]=(Rho-sumfi)/2;
					f_in->f[8][Node.Get_index()]=(Rho-sumfi)/2;
				}
				else
				{
					tmp=f_in->f[6][Node.Get_index()]+f_in->f[8][Node.Get_index()];
					f_in->f[6][Node.Get_index()]=tmp*0.5;
					f_in->f[8][Node.Get_index()]=tmp*0.5;
				}
				break;
			case 6:
				f_in->f[6][Node.Get_index()]=f_in->f[OppositeBc[6]][Node.Get_index()];
				f_in->f[2][Node.Get_index()]=f_in->f[OppositeBc[2]][Node.Get_index()];
				f_in->f[3][Node.Get_index()]=f_in->f[OppositeBc[3]][Node.Get_index()];
				if (Node.stream()[1]==false)
				{
					double sumfi=f_in->f[0][Node.Get_index()]+f_in->f[1][Node.Get_index()]+f_in->f[2][Node.Get_index()]+f_in->f[3][Node.Get_index()]+f_in->f[4][Node.Get_index()]+f_in->f[6][Node.Get_index()]+f_in->f[8][Node.Get_index()];
					f_in->f[5][Node.Get_index()]=(Rho-sumfi)/2;
					f_in->f[7][Node.Get_index()]=(Rho-sumfi)/2;
				}
				else
				{
					tmp=f_in->f[7][Node.Get_index()]+f_in->f[5][Node.Get_index()];
					f_in->f[5][Node.Get_index()]=tmp*0.5;
					f_in->f[7][Node.Get_index()]=tmp*0.5;
				}
				break;
			case 7:
				f_in->f[7][Node.Get_index()]=f_in->f[OppositeBc[7]][Node.Get_index()];
				f_in->f[3][Node.Get_index()]=f_in->f[OppositeBc[3]][Node.Get_index()];
				f_in->f[4][Node.Get_index()]=f_in->f[OppositeBc[4]][Node.Get_index()];
				if (Node.stream()[1]==false)
				{
					double sumfi=f_in->f[0][Node.Get_index()]+f_in->f[1][Node.Get_index()]+f_in->f[2][Node.Get_index()]+f_in->f[3][Node.Get_index()]+f_in->f[4][Node.Get_index()]+f_in->f[5][Node.Get_index()]+f_in->f[7][Node.Get_index()];
					f_in->f[6][Node.Get_index()]=(Rho-sumfi)/2;
					f_in->f[8][Node.Get_index()]=(Rho-sumfi)/2;
				}
				else
				{
					tmp=f_in->f[6][Node.Get_index()]+f_in->f[8][Node.Get_index()];
					f_in->f[6][Node.Get_index()]=tmp*0.5;
					f_in->f[8][Node.Get_index()]=tmp*0.5;
				}
				break;
			case 8:
				f_in->f[8][Node.Get_index()]=f_in->f[OppositeBc[8]][Node.Get_index()];
				f_in->f[1][Node.Get_index()]=f_in->f[OppositeBc[1]][Node.Get_index()];
				f_in->f[4][Node.Get_index()]=f_in->f[OppositeBc[4]][Node.Get_index()];
				if (Node.stream()[2]==false)
				{

					double sumfi=f_in->f[0][Node.Get_index()]+f_in->f[1][Node.Get_index()]+f_in->f[2][Node.Get_index()]+f_in->f[3][Node.Get_index()]+f_in->f[4][Node.Get_index()]+f_in->f[6][Node.Get_index()]+f_in->f[8][Node.Get_index()];
					f_in->f[5][Node.Get_index()]=(Rho-sumfi)/2;
					f_in->f[7][Node.Get_index()]=(Rho-sumfi)/2;
				}
				else
				{
					tmp=f_in->f[5][Node.Get_index()]+f_in->f[7][Node.Get_index()];
					f_in->f[5][Node.Get_index()]=tmp*0.5;
					f_in->f[7][Node.Get_index()]=tmp*0.5;
				}
				break;
			default:
				std::cerr<<"Direction corner bounce back not found. x:"<<Node.get_x()<<" y: "<<Node.get_y()<<std::endl;
				break;
			}
}
//Get the density from the global variable
double D2Q9Corner::FixRho(NodeCorner2D& Node, double *Rho){
	return Rho[Node.Get_index()];
}
//Get the density by using the two direct neighbours
double D2Q9Corner::ExtrapolationAvgRho(NodeCorner2D& Node, double *Rho){
	switch(Node.Get_BcNormal())
	{
	case 5:
		direction1=1;
		direction2=2;
		doubleTmpReturn=(Rho[Node.Get_connect()[direction1]]+Rho[Node.Get_connect()[direction2]])*0.5;
		break;
	case 6:
		direction1=2;
		direction2=3;
		doubleTmpReturn=(Rho[Node.Get_connect()[direction1]]+Rho[Node.Get_connect()[direction2]])*0.5;
		break;
	case 7:
		direction1=3;
		direction2=4;
		doubleTmpReturn=(Rho[Node.Get_connect()[direction1]]+Rho[Node.Get_connect()[direction2]])*0.5;
		break;
	case 8:
		direction1=1;
		direction2=4;
		doubleTmpReturn=(Rho[Node.Get_connect()[direction1]]+Rho[Node.Get_connect()[direction2]])*0.5;
		break;
	}

	return doubleTmpReturn;
}
void D2Q9Corner::ApplyDiffuseWall(NodeCorner2D& Node, DistriFunct* f_in, double Rho, double U, double V){
	double tmp=0;
	switch(Node.Get_BcNormal())
			{
			case 5:
					rhodiff=(f_in->f[3][Node.Get_index()]+f_in->f[4][Node.Get_index()]+f_in->f[7][Node.Get_index()])/SumWeightConvexNE;
					f_in->f[5][Node.Get_index()]=omegaBc[5]*rhodiff;
					f_in->f[1][Node.Get_index()]=omegaBc[1]*rhodiff;
					f_in->f[2][Node.Get_index()]=omegaBc[2]*rhodiff;
					f_in->f[6][Node.Get_index()]=omegaBc[6]*rhodiff;
					f_in->f[8][Node.Get_index()]=omegaBc[8]*rhodiff;
				break;
			case 6:

					rhodiff=(f_in->f[1][Node.Get_index()]+f_in->f[4][Node.Get_index()]+f_in->f[8][Node.Get_index()])/SumWeightConvexNW;
					f_in->f[6][Node.Get_index()]=omegaBc[6]*rhodiff;
					f_in->f[2][Node.Get_index()]=omegaBc[2]*rhodiff;
					f_in->f[3][Node.Get_index()]=omegaBc[3]*rhodiff;
					f_in->f[5][Node.Get_index()]=omegaBc[5]*rhodiff;
					f_in->f[7][Node.Get_index()]=omegaBc[7]*rhodiff;
				break;
			case 7:

					rhodiff=(f_in->f[1][Node.Get_index()]+f_in->f[2][Node.Get_index()]+f_in->f[5][Node.Get_index()])/SumWeightConvexSW;
					f_in->f[7][Node.Get_index()]=omegaBc[7]*rhodiff;
					f_in->f[3][Node.Get_index()]=omegaBc[3]*rhodiff;
					f_in->f[4][Node.Get_index()]=omegaBc[4]*rhodiff;
					f_in->f[6][Node.Get_index()]=omegaBc[6]*rhodiff;
					f_in->f[8][Node.Get_index()]=omegaBc[8]*rhodiff;
				break;
			case 8:

					rhodiff=(f_in->f[2][Node.Get_index()]+f_in->f[3][Node.Get_index()]+f_in->f[6][Node.Get_index()])/SumWeightConvexSE;
					f_in->f[8][Node.Get_index()]=omegaBc[8]*rhodiff;
					f_in->f[1][Node.Get_index()]=omegaBc[1]*rhodiff;
					f_in->f[4][Node.Get_index()]=omegaBc[4]*rhodiff;
					f_in->f[5][Node.Get_index()]=omegaBc[5]*rhodiff;
					f_in->f[7][Node.Get_index()]=omegaBc[7]*rhodiff;
				break;
			default:
				std::cerr<<"Direction: "<< Node.Get_BcNormal()<<" (Corner diffuse boundary conditions) not found"<<std::endl;
				break;
			}
}
