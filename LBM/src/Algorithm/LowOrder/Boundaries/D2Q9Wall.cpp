/*
 * ============================================================================
 * D2Q9Wall.cpp
 *
 *  Created on: 29 Aug 2016
 *      Author: Thomas Burel
 *      Version     :
 *      Copyright   : Your copyright notice
 *      Description : LBM Code Developing by JWFL at University of Strathclyde
 * ============================================================================
 */

#include "D2Q9Wall.h"

D2Q9Wall::D2Q9Wall() {
	feq=0;
	Rho=0;
	PtrWallMethod=0;
	PtrSpecialWallMethod=0;
}

D2Q9Wall::~D2Q9Wall() {
	// TODO Auto-generated destructor stub
}

void D2Q9Wall::SetWall(Parameters *Param){
	switch(Param->Get_WallType())
	{
	case Diffuse:
		PtrWallMethod=&D2Q9Wall::ApplyDiffuseWall;
		PtrSpecialWallMethod=&D2Q9Wall::ApplyDiffuseWallSymmetry;
		break;
	case BounceBack:
		PtrWallMethod=&D2Q9Wall::ApplyBounceBackWall;
		PtrSpecialWallMethod=&D2Q9Wall::ApplyBounceBackSymmetry;
		break;
	case HeZouWall:
		PtrWallMethod=&D2Q9Wall::ApplyHeZouWall;
		PtrSpecialWallMethod=&D2Q9Wall::ApplyBounceBackSymmetry;
		break;
	default:
		std::cerr<<"Wall node type not found."<<std::endl;
	}
}
void D2Q9Wall::ApplyWall(int const &BcNormal,int const *Connect, DistriFunct* f_in, double const *Rho, double const *U, double const *V){
	(this->*PtrWallMethod)(BcNormal,Connect,f_in);
}
void D2Q9Wall::ApplySpecialWall(NodeWall2D& Node, DistriFunct* f_in, std::map<int,NodeType> TypeOfNode_, double const *Rho, double const *U, double const *V){
	(this->*PtrSpecialWallMethod)(Node,f_in,TypeOfNode_);
}
void D2Q9Wall::ApplyDiffuseWall(int const &BcNormal,int const *Connect, DistriFunct* f_in){

	switch(BcNormal)
			{
			case 2:
				rhodiff=(f_in->f[4][Connect[0]]+f_in->f[7][Connect[0]]+f_in->f[8][Connect[0]])/SumWeightS;
				f_in->f[2][Connect[0]]=omegaBc[2]*rhodiff;
				f_in->f[5][Connect[0]]=omegaBc[5]*rhodiff;
				f_in->f[6][Connect[0]]=omegaBc[6]*rhodiff;
				break;
			case 4:
				rhodiff=(f_in->f[2][Connect[0]]+f_in->f[5][Connect[0]]+f_in->f[6][Connect[0]])/SumWeightN;
				f_in->f[4][Connect[0]]=omegaBc[4]*rhodiff;
				f_in->f[7][Connect[0]]=omegaBc[7]*rhodiff;
				f_in->f[8][Connect[0]]=omegaBc[8]*rhodiff;
				break;
			case 1:
				rhodiff=(f_in->f[3][Connect[0]]+f_in->f[6][Connect[0]]+f_in->f[7][Connect[0]])/SumWeightW;
				f_in->f[1][Connect[0]]=omegaBc[1]*rhodiff;
				f_in->f[5][Connect[0]]=omegaBc[5]*rhodiff;
				f_in->f[8][Connect[0]]=omegaBc[8]*rhodiff;

				break;
			case 3:
				rhodiff=(f_in->f[1][Connect[0]]+f_in->f[5][Connect[0]]+f_in->f[8][Connect[0]])/SumWeightE;
				f_in->f[3][Connect[0]]=omegaBc[3]*rhodiff;
				f_in->f[6][Connect[0]]=omegaBc[6]*rhodiff;
				f_in->f[7][Connect[0]]=omegaBc[7]*rhodiff;

				break;
			default:
				std::cerr<<"Direction: "<< BcNormal<<" (Wall diffuse boundary condition) not found"<<std::endl;
				break;
			}
}

///Bounceback Wall treatment
void D2Q9Wall::ApplyBounceBackWall(int const &BcNormal,int const *Connect, DistriFunct* f_in){

	switch(BcNormal)
			{
			case 2:
				f_in->f[2][Connect[0]]=f_in->f[OppositeBc[2]][Connect[0]];
				f_in->f[5][Connect[0]]=f_in->f[OppositeBc[5]][Connect[0]];
				f_in->f[6][Connect[0]]=f_in->f[OppositeBc[6]][Connect[0]];
				break;
			case 4:
				f_in->f[4][Connect[0]]=f_in->f[OppositeBc[4]][Connect[0]];
				f_in->f[7][Connect[0]]=f_in->f[OppositeBc[7]][Connect[0]];
				f_in->f[8][Connect[0]]=f_in->f[OppositeBc[8]][Connect[0]];
				break;
			case 1:
				f_in->f[1][Connect[0]]=f_in->f[OppositeBc[1]][Connect[0]];
				f_in->f[5][Connect[0]]=f_in->f[OppositeBc[5]][Connect[0]];
				f_in->f[8][Connect[0]]=f_in->f[OppositeBc[8]][Connect[0]];
				break;
			case 3:
				f_in->f[3][Connect[0]]=f_in->f[OppositeBc[3]][Connect[0]];
				f_in->f[6][Connect[0]]=f_in->f[OppositeBc[6]][Connect[0]];
				f_in->f[7][Connect[0]]=f_in->f[OppositeBc[7]][Connect[0]];
				break;
			default:
				std::cerr<<"Direction wall bounce back not found. Node Index is: "<<Connect[0]<<" and direction is: "<<BcNormal<<std::endl;
				break;
			}
}
void D2Q9Wall::ApplyHeZouWall(int const &BcNormal,int const *Connect, DistriFunct* f_in){
    switch (BcNormal)
    {
    case 4: //Top

      FUNC_HeZou_NoU(f_in->f[0][Connect[0]],f_in->f[1][Connect[0]],f_in->f[3][Connect[0]],f_in->f[2][Connect[0]],f_in->f[6][Connect[0]],f_in->f[5][Connect[0]],f_in->f[4][Connect[0]],f_in->f[8][Connect[0]],f_in->f[7][Connect[0]]);
      break;
    case 2: //Bot
  	  FUNC_HeZou_NoU(f_in->f[0][Connect[0]],f_in->f[1][Connect[0]],f_in->f[3][Connect[0]],f_in->f[4][Connect[0]],f_in->f[7][Connect[0]],f_in->f[8][Connect[0]],f_in->f[2][Connect[0]],f_in->f[5][Connect[0]],f_in->f[6][Connect[0]]);
  	  break;
    case 1: //West
  	  FUNC_HeZou_NoU(f_in->f[0][Connect[0]],f_in->f[2][Connect[0]],f_in->f[4][Connect[0]],f_in->f[3][Connect[0]],f_in->f[7][Connect[0]],f_in->f[6][Connect[0]],f_in->f[1][Connect[0]],f_in->f[5][Connect[0]],f_in->f[8][Connect[0]]);
  	  break;
    case 3://East

  	  FUNC_HeZou_NoU(f_in->f[0][Connect[0]],f_in->f[2][Connect[0]],f_in->f[4][Connect[0]],f_in->f[1][Connect[0]],f_in->f[8][Connect[0]],f_in->f[5][Connect[0]],f_in->f[3][Connect[0]],f_in->f[6][Connect[0]],f_in->f[7][Connect[0]]);

  	  break;
    default :
        std::cout<<" Problem in the direction of HeZou unknown. Direction is: "<<BcNormal<<std::endl;//,node_HeZou%parameter%ID
        break;
    }
}
///He Zou Boundary conditions(Velocity force to 0)
void D2Q9Wall::FUNC_HeZou_NoU (double & a,double & b,double & c,double & d,double & e,double & f,double & g,double & h,double & i){
	Rho=a+b+c+2.0*(d+e+f);
	feq=(b-c)*0.5;
	g=d;
	h=e-feq;
	i=f+feq;
}
void D2Q9Wall::ApplyDiffuseWallSymmetry(NodeWall2D& Node, DistriFunct* f_in, std::map<int,NodeType>& TypeOfNode_){

	switch(Node.Get_BcNormal())
			{
			case 2:
				rhodiff=(f_in->f[4][Node.Get_index()]+f_in->f[7][Node.Get_index()]+f_in->f[8][Node.Get_index()])/SumWeightS;
				f_in->f[2][Node.Get_index()]=omegaBc[2]*rhodiff;
				f_in->f[5][Node.Get_index()]=omegaBc[5]*rhodiff;
				f_in->f[6][Node.Get_index()]=omegaBc[6]*rhodiff;
				break;
			case 4:
				rhodiff=(f_in->f[2][Node.Get_index()]+f_in->f[5][Node.Get_index()]+f_in->f[6][Node.Get_index()])/SumWeightN;
				f_in->f[4][Node.Get_index()]=omegaBc[4]*rhodiff;
				f_in->f[7][Node.Get_index()]=omegaBc[7]*rhodiff;
				f_in->f[8][Node.Get_index()]=omegaBc[8]*rhodiff;
				break;
			case 1:
				//apply Symmetry
				if(TypeOfNode_[(int)Node.Get_connect(0)]==Solid)//Bottom
				{
					f_in->f[2][Node.Get_index()]=f_in->f[OppositeBc[2]][Node.Get_index()];
					f_in->f[6][Node.Get_index()]=f_in->f[7][Node.Get_index()];
				}
				else //Top
				{
					f_in->f[4][Node.Get_index()]=f_in->f[OppositeBc[4]][Node.Get_index()];
					f_in->f[7][Node.Get_index()]=f_in->f[6][Node.Get_index()];
				}
				rhodiff=(f_in->f[3][Node.Get_index()]+f_in->f[6][Node.Get_index()]+f_in->f[7][Node.Get_index()])/SumWeightW;
				f_in->f[1][Node.Get_index()]=omegaBc[1]*rhodiff;
				f_in->f[5][Node.Get_index()]=omegaBc[5]*rhodiff;
				f_in->f[8][Node.Get_index()]=omegaBc[8]*rhodiff;
				break;
			case 3:
				//apply Symmetry
				if(TypeOfNode_[(int)Node.Get_connect(0)]==Solid)//Bottom
				{
					f_in->f[2][Node.Get_index()]=f_in->f[OppositeBc[2]][Node.Get_index()];
					f_in->f[5][Node.Get_index()]=f_in->f[8][Node.Get_index()];
				}
				else //Top
				{
					f_in->f[4][Node.Get_index()]=f_in->f[OppositeBc[4]][Node.Get_index()];
					f_in->f[8][Node.Get_index()]=f_in->f[5][Node.Get_index()];
				}
				rhodiff=(f_in->f[1][Node.Get_index()]+f_in->f[5][Node.Get_index()]+f_in->f[8][Node.Get_index()])/SumWeightE;
				f_in->f[3][Node.Get_index()]=omegaBc[3]*rhodiff;
				f_in->f[6][Node.Get_index()]=omegaBc[6]*rhodiff;
				f_in->f[7][Node.Get_index()]=omegaBc[7]*rhodiff;
				break;
			default:
				std::cerr<<"Direction: "<< Node.Get_BcNormal()<<" (Wall diffuse boundary condition) not found"<<std::endl;
				break;
			}
}
void D2Q9Wall::ApplyBounceBackSymmetry(NodeWall2D & Node, DistriFunct* f_in, std::map<int,NodeType>&  TypeOfNode_){

	switch(Node.Get_BcNormal())
			{
			case 2:
				f_in->f[2][Node.Get_index()]=f_in->f[OppositeBc[2]][Node.Get_index()];
				f_in->f[5][Node.Get_index()]=f_in->f[OppositeBc[5]][Node.Get_index()];
				f_in->f[6][Node.Get_index()]=f_in->f[OppositeBc[6]][Node.Get_index()];
				break;
			case 4:
				f_in->f[4][Node.Get_index()]=f_in->f[OppositeBc[4]][Node.Get_index()];
				f_in->f[7][Node.Get_index()]=f_in->f[OppositeBc[7]][Node.Get_index()];
				f_in->f[8][Node.Get_index()]=f_in->f[OppositeBc[8]][Node.Get_index()];
				break;
			case 1:
				//apply Symmetry
				if(TypeOfNode_[(int)Node.Get_connect(0)]==Solid)//Bottom
				{
					f_in->f[2][Node.Get_index()]=f_in->f[OppositeBc[2]][Node.Get_index()];
					f_in->f[6][Node.Get_index()]=f_in->f[7][Node.Get_index()];
				}
				else //Top
				{
					f_in->f[4][Node.Get_index()]=f_in->f[OppositeBc[4]][Node.Get_index()];
					f_in->f[7][Node.Get_index()]=f_in->f[6][Node.Get_index()];
				}
				f_in->f[1][Node.Get_index()]=f_in->f[OppositeBc[1]][Node.Get_index()];
				f_in->f[5][Node.Get_index()]=f_in->f[OppositeBc[5]][Node.Get_index()];
				f_in->f[8][Node.Get_index()]=f_in->f[OppositeBc[8]][Node.Get_index()];
				break;
			case 3:
				//apply Symmetry
				if(TypeOfNode_[(int)Node.Get_connect(0)]==Solid)//Bottom
				{
					f_in->f[2][Node.Get_index()]=f_in->f[OppositeBc[2]][Node.Get_index()];
					f_in->f[5][Node.Get_index()]=f_in->f[8][Node.Get_index()];
				}
				else //Top
				{
					f_in->f[4][Node.Get_index()]=f_in->f[OppositeBc[4]][Node.Get_index()];
					f_in->f[8][Node.Get_index()]=f_in->f[5][Node.Get_index()];
				}
				f_in->f[3][Node.Get_index()]=f_in->f[OppositeBc[3]][Node.Get_index()];
				f_in->f[6][Node.Get_index()]=f_in->f[OppositeBc[6]][Node.Get_index()];
				f_in->f[7][Node.Get_index()]=f_in->f[OppositeBc[7]][Node.Get_index()];
				break;
			default:
				std::cerr<<"Direction wall bounce back not found. Index: "<<Node.Get_index()<<std::endl;
				break;
			}
}
