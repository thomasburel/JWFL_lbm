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
void D2Q9Corner::Set_Corner(Parameters *Param, double ** &Ei){
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
		Extrapol.SelectExtrapolationType(NoExtrapol);
		PtrCalculRhoCornerWall=&D2Q9Corner::FixRho;
		break;
	case ExtrapolCP:
		Extrapol.SelectExtrapolationType(WeightDistanceExtrapol);
		PtrCalculRhoCornerWall=&D2Q9Corner::ExtrapolationAvgRho;
		break;
	default:
		std::cerr<<"Corner pressure type not found."<<std::endl;
	}
	EiBc=Ei;
}
void D2Q9Corner::ApplyCorner(NodeCorner2D& Node, DistriFunct* f_in,double const & RhoDef,double const & UDef,double const & VDef, double *Rho, double *U, double *V){
	Rho[Node.Get_index()]=RhoDef;
	U[Node.Get_index()]=UDef;
	V[Node.Get_index()]=VDef;
	(this->*PtrCornerMethod)(Node.Get_BcNormal(),Node.Get_connect(),Node.Get_CornerType()==::Concave, f_in, GetRho(Node, Rho), U[Node.Get_index()], V[Node.Get_index()]);
}
void D2Q9Corner::ApplyCornerWall(NodeCorner2D& Node, DistriFunct* f_in, double *Rho, double *U, double *V){
	(this->*PtrCornerWallMethod)(Node.Get_BcNormal(),Node.Get_connect(),Node.Get_CornerType()==::Concave,f_in, GetRho(Node, Rho), U[Node.Get_index()], V[Node.Get_index()]);
}

void D2Q9Corner::ApplyCornerSpecialWall(NodeWall2D& Node, DistriFunct* f_in, double *Rho, double *U, double *V){
	(this->*PtrCornerMethod)(Node.Get_BcNormal(),Node.Get_connect(),true, f_in, Rho[Node.Get_index()], U[Node.Get_index()], V[Node.Get_index()]);
}
void D2Q9Corner::ApplyPreVelSpecialWall(NodeWall2D& Node, DistriFunct* f_in,double const & RhoDef,double const & UDef,double const & VDef){
	(this->*PtrCornerWallMethod)(Node.Get_BcNormal(),Node.Get_connect(),true,f_in, RhoDef, UDef, VDef);
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

void D2Q9Corner::ApplyHoChan(int const &BcNormal,int const *Connect, bool concave, DistriFunct* f_in, double Rho, double U, double V){

      switch (BcNormal)
      {

      case 8: //Top West
    	  InvV=-V;
    	  FUNC_corner(f_in->f[0][Connect[0]],f_in->f[1][Connect[0]],f_in->f[4][Connect[0]],f_in->f[8][Connect[0]],f_in->f[7][Connect[0]],f_in->f[5][Connect[0]],f_in->f[3][Connect[0]],f_in->f[2][Connect[0]],f_in->f[6][Connect[0]],U,InvV,Rho);
    	  break;
      case 5: //Bottom West
    	  FUNC_corner(f_in->f[0][Connect[0]],f_in->f[1][Connect[0]],f_in->f[2][Connect[0]],f_in->f[5][Connect[0]],f_in->f[6][Connect[0]],f_in->f[8][Connect[0]],f_in->f[3][Connect[0]],f_in->f[4][Connect[0]],f_in->f[7][Connect[0]],U,V,Rho);
    	  break;
      case 7 : //Top East
    	  InvU=-U;
    	  InvV=-V;
    	  FUNC_corner(f_in->f[0][Connect[0]],f_in->f[3][Connect[0]],f_in->f[4][Connect[0]],f_in->f[7][Connect[0]],f_in->f[8][Connect[0]],f_in->f[6][Connect[0]],f_in->f[1][Connect[0]],f_in->f[2][Connect[0]],f_in->f[5][Connect[0]],InvU,InvV,Rho);
    	  break;
      case 6: //Bottom East
    	  InvU=-U;
    	  FUNC_corner(f_in->f[0][Connect[0]],f_in->f[3][Connect[0]],f_in->f[2][Connect[0]],f_in->f[6][Connect[0]],f_in->f[5][Connect[0]],f_in->f[7][Connect[0]],f_in->f[1][Connect[0]],f_in->f[4][Connect[0]],f_in->f[8][Connect[0]],InvU,V,Rho);
    	  break;

      default :
          std::cout<<" Problem in the direction of HeZou Corner unknown. Direction is: "<<BcNormal<<std::endl;
      }
}
void D2Q9Corner::ApplyHoChanNoVel(int const &BcNormal,int const *Connect, bool concave, DistriFunct* f_in, double Rho, double U, double V){
    switch (BcNormal)
      {
      case 8: //Top West

    	  FUNC_corner_no_vel(f_in->f[0][Connect[0]],f_in->f[1][Connect[0]],f_in->f[4][Connect[0]],f_in->f[8][Connect[0]],f_in->f[7][Connect[0]],f_in->f[5][Connect[0]],f_in->f[3][Connect[0]],f_in->f[2][Connect[0]],f_in->f[6][Connect[0]],Rho);
   	  break;
      case 5: //Bottom West
    	  FUNC_corner_no_vel(f_in->f[0][Connect[0]],f_in->f[1][Connect[0]],f_in->f[2][Connect[0]],f_in->f[5][Connect[0]],f_in->f[6][Connect[0]],f_in->f[8][Connect[0]],f_in->f[3][Connect[0]],f_in->f[4][Connect[0]],f_in->f[7][Connect[0]],Rho);
    	  break;
      case 7 : //Top East

    	  FUNC_corner_no_vel(f_in->f[0][Connect[0]],f_in->f[3][Connect[0]],f_in->f[4][Connect[0]],f_in->f[7][Connect[0]],f_in->f[8][Connect[0]],f_in->f[6][Connect[0]],f_in->f[1][Connect[0]],f_in->f[2][Connect[0]],f_in->f[5][Connect[0]],Rho);
    	  break;
      case 6: //Bottom East

    	  FUNC_corner_no_vel(f_in->f[0][Connect[0]],f_in->f[3][Connect[0]],f_in->f[2][Connect[0]],f_in->f[6][Connect[0]],f_in->f[5][Connect[0]],f_in->f[7][Connect[0]],f_in->f[1][Connect[0]],f_in->f[4][Connect[0]],f_in->f[8][Connect[0]],Rho);
    	  break;
      default :
          std::cout<<"Corner direction is unknown. Direction is: "<<BcNormal<<std::endl;
      }
}
void D2Q9Corner::ApplyBounceBack(int const &BcNormal,int const *Connect, bool concave, DistriFunct* f_in, double Rho, double U, double V){
	double tmp=0;
	switch(BcNormal)
			{
			case 5:
				f_in->f[5][Connect[0]]=f_in->f[OppositeBc[5]][Connect[0]];
				f_in->f[1][Connect[0]]=f_in->f[OppositeBc[1]][Connect[0]];
				f_in->f[2][Connect[0]]=f_in->f[OppositeBc[2]][Connect[0]];
				if (concave)
				{

					double sumfi=f_in->f[0][Connect[0]]+f_in->f[1][Connect[0]]+f_in->f[2][Connect[0]]+f_in->f[3][Connect[0]]+f_in->f[4][Connect[0]]+f_in->f[5][Connect[0]]+f_in->f[7][Connect[0]];
					f_in->f[6][Connect[0]]=(Rho-sumfi)/2;
					f_in->f[8][Connect[0]]=(Rho-sumfi)/2;
				}
				else
				{
					tmp=f_in->f[6][Connect[0]]+f_in->f[8][Connect[0]];
					f_in->f[6][Connect[0]]=tmp*0.5;
					f_in->f[8][Connect[0]]=tmp*0.5;
				}
				break;
			case 6:
				f_in->f[6][Connect[0]]=f_in->f[OppositeBc[6]][Connect[0]];
				f_in->f[2][Connect[0]]=f_in->f[OppositeBc[2]][Connect[0]];
				f_in->f[3][Connect[0]]=f_in->f[OppositeBc[3]][Connect[0]];
				if (concave)
				{
					double sumfi=f_in->f[0][Connect[0]]+f_in->f[1][Connect[0]]+f_in->f[2][Connect[0]]+f_in->f[3][Connect[0]]+f_in->f[4][Connect[0]]+f_in->f[6][Connect[0]]+f_in->f[8][Connect[0]];
					f_in->f[5][Connect[0]]=(Rho-sumfi)/2;
					f_in->f[7][Connect[0]]=(Rho-sumfi)/2;
				}
				else
				{
					tmp=f_in->f[7][Connect[0]]+f_in->f[5][Connect[0]];
					f_in->f[5][Connect[0]]=tmp*0.5;
					f_in->f[7][Connect[0]]=tmp*0.5;
				}
				break;
			case 7:
				f_in->f[7][Connect[0]]=f_in->f[OppositeBc[7]][Connect[0]];
				f_in->f[3][Connect[0]]=f_in->f[OppositeBc[3]][Connect[0]];
				f_in->f[4][Connect[0]]=f_in->f[OppositeBc[4]][Connect[0]];
				if (concave)
				{
					double sumfi=f_in->f[0][Connect[0]]+f_in->f[1][Connect[0]]+f_in->f[2][Connect[0]]+f_in->f[3][Connect[0]]+f_in->f[4][Connect[0]]+f_in->f[5][Connect[0]]+f_in->f[7][Connect[0]];
					f_in->f[6][Connect[0]]=(Rho-sumfi)/2;
					f_in->f[8][Connect[0]]=(Rho-sumfi)/2;
				}
				else
				{
					tmp=f_in->f[6][Connect[0]]+f_in->f[8][Connect[0]];
					f_in->f[6][Connect[0]]=tmp*0.5;
					f_in->f[8][Connect[0]]=tmp*0.5;
				}
				break;
			case 8:
				f_in->f[8][Connect[0]]=f_in->f[OppositeBc[8]][Connect[0]];
				f_in->f[1][Connect[0]]=f_in->f[OppositeBc[1]][Connect[0]];
				f_in->f[4][Connect[0]]=f_in->f[OppositeBc[4]][Connect[0]];
				if (concave)
				{

					double sumfi=f_in->f[0][Connect[0]]+f_in->f[1][Connect[0]]+f_in->f[2][Connect[0]]+f_in->f[3][Connect[0]]+f_in->f[4][Connect[0]]+f_in->f[6][Connect[0]]+f_in->f[8][Connect[0]];
					f_in->f[5][Connect[0]]=(Rho-sumfi)/2;
					f_in->f[7][Connect[0]]=(Rho-sumfi)/2;
				}
				else
				{
					tmp=f_in->f[5][Connect[0]]+f_in->f[7][Connect[0]];
					f_in->f[5][Connect[0]]=tmp*0.5;
					f_in->f[7][Connect[0]]=tmp*0.5;
				}
				break;
			default:
				std::cerr<<"Direction corner bounce back not found. "<<std::endl;
				break;
			}
}

//Get the density from the global variable
double D2Q9Corner::GetRho(NodeCorner2D& Node, double *Rho){
	(this->*PtrCalculRhoCornerWall)(Node, Rho);
	return Rho[Node.Get_index()];
}
//double D2Q9Corner::FixRho(NodeCorner2D& Node, double *Rho){
void D2Q9Corner::FixRho(NodeCorner2D& Node, double *Rho){
//	Extrapol.ExtrapolationOnCornerConcave(Rho,Node.Get_connect(),Node.Get_BcNormal());
//	return Rho[Node.Get_index()];
}
//Get the density by using the two direct neighbours
void D2Q9Corner::ExtrapolationAvgRho(NodeCorner2D& Node, double *Rho){
	switch(Node.Get_BcNormal())
	{
	case 5:
		if(Node.Get_CornerType()==Concave)
			Extrapol.ExtrapolationOnCornerConcave(Rho,Node.Get_connect(),Node.Get_BcNormal());
		else
			Extrapol.ExtrapolationOnCornerConvex(Rho,Node.Get_connect(),Node.Get_BcNormal());
		break;
	case 6:
		if(Node.Get_CornerType()==Concave)
			Extrapol.ExtrapolationOnCornerConcave(Rho,Node.Get_connect(),Node.Get_BcNormal());
		else
			Extrapol.ExtrapolationOnCornerConvex(Rho,Node.Get_connect(),Node.Get_BcNormal());
		break;
	case 7:
		if(Node.Get_CornerType()==Concave)
			Extrapol.ExtrapolationOnCornerConcave(Rho,Node.Get_connect(),Node.Get_BcNormal());
		else
			Extrapol.ExtrapolationOnCornerConvex(Rho,Node.Get_connect(),Node.Get_BcNormal());
		break;
	case 8:
		if(Node.Get_CornerType()==Concave)
			Extrapol.ExtrapolationOnCornerConcave(Rho,Node.Get_connect(),Node.Get_BcNormal());
		else
			Extrapol.ExtrapolationOnCornerConvex(Rho,Node.Get_connect(),Node.Get_BcNormal());
		break;
	}
}
void D2Q9Corner::ApplyDiffuseWall(int const &BcNormal,int const *Connect, bool concave, DistriFunct* f_in, double Rho, double U, double V){
	//double tmp=0;
	switch(BcNormal)
			{
			case 5:
					rhodiff=(f_in->f[3][Connect[0]]+f_in->f[4][Connect[0]]+f_in->f[7][Connect[0]])/SumWeightConvexNE;
					f_in->f[5][Connect[0]]=omegaBc[5]*rhodiff;
					f_in->f[1][Connect[0]]=omegaBc[1]*rhodiff;
					f_in->f[2][Connect[0]]=omegaBc[2]*rhodiff;
					f_in->f[6][Connect[0]]=omegaBc[6]*rhodiff;
					f_in->f[8][Connect[0]]=omegaBc[8]*rhodiff;
				break;
			case 6:

					rhodiff=(f_in->f[1][Connect[0]]+f_in->f[4][Connect[0]]+f_in->f[8][Connect[0]])/SumWeightConvexNW;
					f_in->f[6][Connect[0]]=omegaBc[6]*rhodiff;
					f_in->f[2][Connect[0]]=omegaBc[2]*rhodiff;
					f_in->f[3][Connect[0]]=omegaBc[3]*rhodiff;
					f_in->f[5][Connect[0]]=omegaBc[5]*rhodiff;
					f_in->f[7][Connect[0]]=omegaBc[7]*rhodiff;
				break;
			case 7:

					rhodiff=(f_in->f[1][Connect[0]]+f_in->f[2][Connect[0]]+f_in->f[5][Connect[0]])/SumWeightConvexSW;
					f_in->f[7][Connect[0]]=omegaBc[7]*rhodiff;
					f_in->f[3][Connect[0]]=omegaBc[3]*rhodiff;
					f_in->f[4][Connect[0]]=omegaBc[4]*rhodiff;
					f_in->f[6][Connect[0]]=omegaBc[6]*rhodiff;
					f_in->f[8][Connect[0]]=omegaBc[8]*rhodiff;
				break;
			case 8:

					rhodiff=(f_in->f[2][Connect[0]]+f_in->f[3][Connect[0]]+f_in->f[6][Connect[0]])/SumWeightConvexSE;
					f_in->f[8][Connect[0]]=omegaBc[8]*rhodiff;
					f_in->f[1][Connect[0]]=omegaBc[1]*rhodiff;
					f_in->f[4][Connect[0]]=omegaBc[4]*rhodiff;
					f_in->f[5][Connect[0]]=omegaBc[5]*rhodiff;
					f_in->f[7][Connect[0]]=omegaBc[7]*rhodiff;
				break;
			default:
				std::cerr<<"Direction: "<< BcNormal<<" (Corner diffuse boundary conditions) not found"<<std::endl;
				break;
			}
}
