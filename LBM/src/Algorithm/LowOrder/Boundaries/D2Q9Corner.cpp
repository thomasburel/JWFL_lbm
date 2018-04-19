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
	PtrCornerSpecialWallMethod=0;
	PtrCornerWallSpecialWallMethod=0;
	PtrCalculRhoCornerWall=0;
	PtrPreStreamCornerMethod=0;
	PtrPreStreamCornerWallMethod=0;
	PtrPreStreamCornerSpecialWallMethod=0;
	doubleTmpReturn=0;
	direction1=0;
	direction2=0;
}

D2Q9Corner::~D2Q9Corner() {
	// TODO Auto-generated destructor stub
}
void D2Q9Corner::Set_Corner(Dictionary *PtrDic,NodeArrays2D* NodeArrays, Parameters *Param, double ** &Ei,unsigned int nbDistributions){
	PtrCornerMethod=&D2Q9Corner::ApplyHoChan;
	PtrCornerSpecialWallMethod=&D2Q9Corner::ApplyHoChan;
	PtrCornerWallSpecialWallMethod=&D2Q9Corner::ApplyHoChan;
	switch(Param->Get_WallType())
	{
	case Diffuse:
		PtrCornerWallMethod=&D2Q9Corner::ApplyDiffuseWall;
		//PtrCornerWallSpecialWallMethod=&D2Q9Corner::ApplyDiffuseWall;
		break;
	case BounceBack:
		PtrCornerWallMethod=&D2Q9Corner::ApplyBounceBack;
		//PtrCornerWallSpecialWallMethod=&D2Q9Corner::ApplyBounceBack;
		break;
	case HalfWayBounceBack:
		PtrCornerWallMethod=&D2Q9Corner::ApplyHalfWayBounceBack;
		//PtrCornerWallSpecialWallMethod=&D2Q9Corner::ApplyHalfWayBounceBack;
		PtrPreStreamCornerMethod=&D2Q9Corner::ApplyHalfWayBounceBackPreStream;
		PtrPreStreamCornerWallMethod=&D2Q9Corner::ApplyHalfWayBounceBackPreStream;
		PtrPreStreamCornerSpecialWallMethod=&D2Q9Corner::ApplyHalfWayBounceBackPreStream;
		SetHalfWayBounceBack(NodeArrays,nbDistributions);
		break;
	case HeZouWall:
		PtrCornerWallMethod=&D2Q9Corner::ApplyHoChanNoVel;
		//PtrCornerWallSpecialWallMethod=&D2Q9Corner::ApplyHoChanNoVel;
		break;
	default:
		std::cerr<<"Wall node type not found."<<std::endl;
	}
	switch(Param->Get_CornerPressureType())
	{
	case FixCP:
		Extrapol.SelectExtrapolationType(ModelEnum::NoExtrapol);
		PtrCalculRhoCornerWall=&D2Q9Corner::FixRho;
		break;
	case ExtrapolCP:
		Extrapol.SelectExtrapolationType(ModelEnum::WeightDistanceExtrapol);
		PtrCalculRhoCornerWall=&D2Q9Corner::ExtrapolationAvgRho;
		break;
	case LocalCP:
		PtrCalculRhoCornerWall=&D2Q9Corner::LocalRho;
		break;
	default:
		std::cerr<<"Corner pressure type not found."<<std::endl;
	}
	EiBc=Ei;
}
void D2Q9Corner::SetHalfWayBounceBack(NodeArrays2D* NodeArrays,unsigned int nbDistributions){
	//Prepare to save 3 discrete velocity
	for(unsigned int i=0;i<NodeArrays->NodeWall.size();i++)
		NodeArrays->NodeWall[i].Ini_SaveData(3*nbDistributions);
	for(unsigned int i=0;i<NodeArrays->NodeSpecialWall.size();i++)
		NodeArrays->NodeSpecialWall[i].Ini_SaveData(3*nbDistributions);
	for(unsigned int i=0;i<NodeArrays->NodeCorner.size();i++)
		NodeArrays->NodeCorner[i].Ini_SaveData(5*nbDistributions);
}
void D2Q9Corner::ApplyCorner(NodeCorner2D& Node, DistriFunct* f_in, unsigned int idxDistribution,double const & RhoDef,double const & UDef,double const & VDef, double *Rho, double *U, double *V){
	Rho[Node.Get_index()]=RhoDef;
	U[Node.Get_index()]=UDef;
	V[Node.Get_index()]=VDef;
	(this->*PtrCornerMethod)(Node,Node.Get_BcNormal(),Node.Get_connect(),Node.Get_CornerType()==::Concave, f_in, GetRho(Node, Rho,f_in), U[Node.Get_index()], V[Node.Get_index()],idxDistribution);
}
void D2Q9Corner::ApplyCornerWall(NodeCorner2D& Node, DistriFunct* f_in, unsigned int idxDistribution, double *Rho, double *U, double *V){
	(this->*PtrCornerWallMethod)(Node,Node.Get_BcNormal(),Node.Get_connect(),Node.Get_CornerType()==::Concave,f_in, GetRho(Node, Rho,f_in), U[Node.Get_index()], V[Node.Get_index()],idxDistribution);
}

void D2Q9Corner::ApplyCornerSpecialWall(NodeWall2D& Node, DistriFunct* f_in, unsigned int idxDistribution, double *Rho, double *U, double *V){
	(this->*PtrCornerSpecialWallMethod)(Node,Node.Get_BcNormal(),Node.Get_connect(),true, f_in, Rho[Node.Get_index()], U[Node.Get_index()], V[Node.Get_index()],idxDistribution);
}
void D2Q9Corner::ApplyPreVelSpecialWall(NodeWall2D& Node, DistriFunct* f_in, unsigned int idxDistribution,double const & RhoDef,double const & UDef,double const & VDef){
	(this->*PtrCornerWallSpecialWallMethod)(Node,Node.Get_BcNormal(),Node.Get_connect(),true,f_in, RhoDef, UDef, VDef,idxDistribution);
}

void D2Q9Corner::ApplyCornerPreStream(NodeCorner2D& Node, DistriFunct* f_in, unsigned int idxDistribution,double const & RhoDef,double const & UDef,double const & VDef, double *Rho, double *U, double *V){
	Rho[Node.Get_index()]=RhoDef;
	U[Node.Get_index()]=UDef;
	V[Node.Get_index()]=VDef;
	(this->*PtrPreStreamCornerMethod)(Node,Node.Get_BcNormal(),Node.Get_connect(),Node.Get_CornerType()==::Concave, f_in, GetRho(Node, Rho,f_in), U[Node.Get_index()], V[Node.Get_index()],idxDistribution);
}
void D2Q9Corner::ApplyCornerWallPreStream(NodeCorner2D& Node, DistriFunct* f_in, unsigned int idxDistribution, double *Rho, double *U, double *V){
	(this->*PtrPreStreamCornerWallMethod)(Node,Node.Get_BcNormal(),Node.Get_connect(),Node.Get_CornerType()==::Concave,f_in, GetRho(Node, Rho,f_in), U[Node.Get_index()], V[Node.Get_index()],idxDistribution);
}

void D2Q9Corner::ApplyCornerSpecialWallPreStream(NodeWall2D& Node, DistriFunct* f_in, unsigned int idxDistribution, double *Rho, double *U, double *V){
	(this->*PtrPreStreamCornerSpecialWallMethod)(Node,Node.Get_BcNormal(),Node.Get_connect(),true, f_in, Rho[Node.Get_index()], U[Node.Get_index()], V[Node.Get_index()],idxDistribution);
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
template <class T>
void D2Q9Corner::ApplyHoChan(T& Node,int const &BcNormal,int const *Connect, bool concave, DistriFunct* f_in, double Rho, double U, double V, unsigned int idxDistribution){

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
template <class T>
void D2Q9Corner::ApplyHoChanNoVel(T& Node,int const &BcNormal,int const *Connect, bool concave, DistriFunct* f_in, double Rho, double U, double V, unsigned int idxDistribution){
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
template <class T>
void D2Q9Corner::ApplyBounceBack(T& Node,int const &BcNormal,int const *Connect, bool concave, DistriFunct* f_in, double Rho, double U, double V, unsigned int idxDistribution){
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
				std::cerr<<"Direction corner bounce back not found. ID node: "<<Connect[0]<<std::endl;
				break;
			}
}
template <class T>
void D2Q9Corner::ApplyHalfWayBounceBack(T& Node,int const &BcNormal,int const *Connect, bool concave, DistriFunct* f_in, double Rho, double U, double V, unsigned int idxDistribution){
//	double tmp=0;
//	std::cout<<"corner id"<<idxDistribution<<std::endl;
	unsigned int idx0=idxDistribution*5;
	f_in->f[BcNormal][Connect[0]]=Node.Get_SaveData(idx0);//f_in->f[OppositeBc[BcNormal]][Connect[OppositeBc[BcNormal]]];
	if (concave)
	{
		idx0++;
		for (unsigned int i=0;i<4;i++){
			f_in->f[BounceBackWallConnect[BcNormal][i]][Connect[0]]= Node.Get_SaveData(idx0+i);
		}
	}

}
template <class T>
void D2Q9Corner::ApplyHalfWayBounceBackPreStream(T& Node,int const &BcNormal,int const *Connect, bool concave, DistriFunct* f_in, double Rho, double U, double V, unsigned int idxDistribution){
//	std::cout<<"Pre corner id"<<idxDistribution<<std::endl;
	unsigned int idx0=idxDistribution*5;
	Node.Set_SaveData(idx0,f_in->f[OppositeBc[BcNormal]][Connect[0]]);
	if (concave)
	{
		idx0++;
		for (unsigned int i=0;i<4;i++){
			Node.Set_SaveData(idx0+i,f_in->f[OppositeBc[BounceBackWallConnect[BcNormal][i]]][Connect[0]]);
		}
	}

}

//Get the density from the global variable
double D2Q9Corner::GetRho(NodeCorner2D& Node, double *Rho,DistriFunct* f_in){
	(this->*PtrCalculRhoCornerWall)(Node, Rho,f_in);
	return Rho[Node.Get_index()];
}
//Use the previous rho
void D2Q9Corner::FixRho(NodeCorner2D& Node, double *Rho,DistriFunct* f_in){
//	Extrapol.ExtrapolationOnCornerConcave(Rho,Node.Get_connect(),Node.Get_BcNormal());
//	return Rho[Node.Get_index()];
}
//use the known distribution to calculate rho
void D2Q9Corner::LocalRho(NodeCorner2D& Node, double *Rho,DistriFunct* f_in){
	switch(Node.Get_BcNormal())
	{
	case 5:
		Rho[Node.Get_connect()[0]]=f_in->f[0][Node.Get_connect()[0]]+2.0*(f_in->f[3][Node.Get_connect()[0]]+f_in->f[4][Node.Get_connect()[0]]+f_in->f[7][Node.Get_connect()[0]]);
		if(Node.Get_CornerType()==Convex)
			Rho[Node.Get_connect()[0]]+=f_in->f[6][Node.Get_connect()[0]]+f_in->f[8][Node.Get_connect()[0]];
		break;
	case 6:
		Rho[Node.Get_connect()[0]]=f_in->f[0][Node.Get_connect()[0]]+2.0*(f_in->f[1][Node.Get_connect()[0]]+f_in->f[4][Node.Get_connect()[0]]+f_in->f[8][Node.Get_connect()[0]]);
		if(Node.Get_CornerType()==Convex)
			Rho[Node.Get_connect()[0]]+=f_in->f[5][Node.Get_connect()[0]]+f_in->f[7][Node.Get_connect()[0]];
		break;
	case 7:
		Rho[Node.Get_connect()[0]]=f_in->f[0][Node.Get_connect()[0]]+2.0*(f_in->f[1][Node.Get_connect()[0]]+f_in->f[2][Node.Get_connect()[0]]+f_in->f[5][Node.Get_connect()[0]]);
		if(Node.Get_CornerType()==Convex)
			Rho[Node.Get_connect()[0]]+=f_in->f[6][Node.Get_connect()[0]]+f_in->f[8][Node.Get_connect()[0]];
		break;
	case 8:
		Rho[Node.Get_connect()[0]]=f_in->f[0][Node.Get_connect()[0]]+2.0*(f_in->f[2][Node.Get_connect()[0]]+f_in->f[3][Node.Get_connect()[0]]+f_in->f[6][Node.Get_connect()[0]]);
		if(Node.Get_CornerType()==Convex)
			Rho[Node.Get_connect()[0]]+=f_in->f[5][Node.Get_connect()[0]]+f_in->f[7][Node.Get_connect()[0]];
		break;
	}
}
//Get the density by using the two direct neighbours
void D2Q9Corner::ExtrapolationAvgRho(NodeCorner2D& Node, double *Rho,DistriFunct* f_in){
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
template <class T>
void D2Q9Corner::ApplyDiffuseWall(T& Node,int const &BcNormal,int const *Connect, bool concave, DistriFunct* f_in, double Rho, double U, double V, unsigned int idxDistribution){
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
