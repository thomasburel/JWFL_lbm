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
	//PtrSpecialWallMethod=0;
	PtrWallGlobalCornerMethod=0;
	PtrPreStreamWallMethod=0;
	PtrPreStreamWallGlobalCornerMethod=0;
}

D2Q9Wall::~D2Q9Wall() {
	// TODO Auto-generated destructor stub
}

void D2Q9Wall::SetWall(Dictionary *PtrDic,NodeArrays2D* NodeArrays, Parameters *Param, double ** &Ei,unsigned int nbDistributions){
	switch(Param->Get_WallType())
	{
	case Diffuse:
		PtrWallMethod=&D2Q9Wall::ApplyDiffuseWall;
		PtrWallGlobalCornerMethod=&D2Q9Wall::ApplyDiffuseWall;
		break;
	case BounceBack:
		PtrWallMethod=&D2Q9Wall::ApplyBounceBackWall;
		PtrWallGlobalCornerMethod=&D2Q9Wall::ApplyBounceBackWall;
		break;
	case HalfWayBounceBack:
		PtrWallMethod=&D2Q9Wall::ApplyHalfWayBounceBackWall;
		PtrWallGlobalCornerMethod=&D2Q9Wall::ApplyHalfWayBounceBackWall;
		PtrPreStreamWallMethod=&D2Q9Wall::ApplyHalfWayBounceBackWallPreStream;
		PtrPreStreamWallGlobalCornerMethod=&D2Q9Wall::ApplyHalfWayBounceBackWallPreStream;
		SetHalfWayBounceBack(NodeArrays,nbDistributions);
		break;
	case HeZouWall:
		PtrWallMethod=&D2Q9Wall::ApplyHeZouWall;
		PtrWallGlobalCornerMethod=&D2Q9Wall::ApplyHeZouWall;
		break;
	default:
		std::cerr<<"Wall node type not found."<<std::endl;
	}
	EiBc=Ei;
}
void D2Q9Wall::SetHalfWayBounceBack(NodeArrays2D* NodeArrays,unsigned int nbDistributions){
	//Prepare to save 3 discrete velocity
	for(unsigned int i=0;i<NodeArrays->NodeWall.size();i++)
		NodeArrays->NodeWall[i].Ini_SaveData(3*nbDistributions);
	for(unsigned int i=0;i<NodeArrays->NodeSpecialWall.size();i++)
		NodeArrays->NodeSpecialWall[i].Ini_SaveData(3*nbDistributions);
	for(unsigned int i=0;i<NodeArrays->NodeGlobalCorner.size();i++)
		NodeArrays->NodeGlobalCorner[i].Ini_SaveData(3*nbDistributions);
}
void D2Q9Wall::ApplyWall(NodeWall2D& Node, int const &BcNormal,int const *Connect, DistriFunct* f_in, unsigned int idxDistribution, double const *Rho, double const *U, double const *V){
	(this->*PtrWallMethod)(Node,BcNormal,Connect,f_in, idxDistribution);
}
void D2Q9Wall::ApplyWall(NodeCorner2D& Node, int const &BcNormal,int const *Connect, DistriFunct* f_in, unsigned int idxDistribution, double const *Rho, double const *U, double const *V){
	(this->*PtrWallGlobalCornerMethod)(Node,BcNormal,Connect,f_in, idxDistribution);
}
void D2Q9Wall::ApplyWallPreStream(NodeWall2D& Node, int const &BcNormal,int const *Connect, DistriFunct* f_in, unsigned int idxDistribution, double const *Rho, double const *U, double const *V){
	(this->*PtrPreStreamWallMethod)(Node,BcNormal,Connect,f_in, idxDistribution);
}
void D2Q9Wall::ApplyWallPreStream(NodeCorner2D& Node, int const &BcNormal,int const *Connect, DistriFunct* f_in, unsigned int idxDistribution, double const *Rho, double const *U, double const *V){
	(this->*PtrPreStreamWallGlobalCornerMethod)(Node,BcNormal,Connect,f_in, idxDistribution);
}
template <class T>
void D2Q9Wall::ApplyDiffuseWall(T& Node,int const &BcNormal,int const *Connect, DistriFunct* f_in, unsigned int idxDistribution){

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
template <class T>
void D2Q9Wall::ApplyBounceBackWall(T& Node,int const &BcNormal,int const *Connect, DistriFunct* f_in, unsigned int idxDistribution){
	f_in->f[BcNormal][Connect[0]]=f_in->f[OppositeBc[BcNormal]][Connect[0]];
	f_in->f[BounceBackWallConnect[BcNormal][0]][Connect[0]]=f_in->f[OppositeBc[BounceBackWallConnect[BcNormal][0]]][Connect[0]];
	f_in->f[BounceBackWallConnect[BcNormal][1]][Connect[0]]=f_in->f[OppositeBc[BounceBackWallConnect[BcNormal][1]]][Connect[0]];
}
///HalfWay Bounceback Wall treatment
template <class T>
void D2Q9Wall::ApplyHalfWayBounceBackWall(T& Node,int const &BcNormal,int const *Connect, DistriFunct* f_in, unsigned int idxDistribution){
//	std::cout<<"wall id"<<idxDistribution<<std::endl;
	unsigned int idx0=idxDistribution*3;
	unsigned int idx1=idx0+1;unsigned int idx2=idx0+2;
	f_in->f[BcNormal][Connect[0]]=Node.Get_SaveData(idx0);
	f_in->f[BounceBackWallConnect[BcNormal][0]][Connect[0]]=Node.Get_SaveData(idx1);
	f_in->f[BounceBackWallConnect[BcNormal][1]][Connect[0]]=Node.Get_SaveData(idx2);
/*	Node.Set_SaveData(idx0,f_in->f[OppositeBc[BcNormal]][Connect[0]]);
	Node.Set_SaveData(idx1,f_in->f[OppositeBc[BounceBackWallConnect[BcNormal][0]]][Connect[0]]);
	Node.Set_SaveData(idx2,f_in->f[OppositeBc[BounceBackWallConnect[BcNormal][1]]][Connect[0]]);
*/
}
///HalfWay Bounceback Wall treatment before streaming
template <class T>
void D2Q9Wall::ApplyHalfWayBounceBackWallPreStream(T& Node,int const &BcNormal,int const *Connect, DistriFunct* f_in, unsigned int idxDistribution){
//	std::cout<<"Pre wall id"<<idxDistribution<<std::endl;
	unsigned int idx0=idxDistribution*3;
	unsigned int idx1=idx0+1;unsigned int idx2=idx0+2;
	Node.Set_SaveData(idx0,f_in->f[OppositeBc[BcNormal]][Connect[0]]);
	Node.Set_SaveData(idx1,f_in->f[OppositeBc[BounceBackWallConnect[BcNormal][0]]][Connect[0]]);
	Node.Set_SaveData(idx2,f_in->f[OppositeBc[BounceBackWallConnect[BcNormal][1]]][Connect[0]]);

}
template <class T>
void D2Q9Wall::ApplyHeZouWall(T& Node,int const &BcNormal,int const *Connect, DistriFunct* f_in, unsigned int idxDistribution){
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
/*
void D2Q9Wall::ApplyDiffuseWallSymmetry(NodeWall2D& Node, DistriFunct* f_in, unsigned int idxDistribution, std::map<int,NodeType>& TypeOfNode_){

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
void D2Q9Wall::ApplyBounceBackSymmetry(NodeWall2D & Node, DistriFunct* f_in, unsigned int idxDistribution, std::map<int,NodeType>&  TypeOfNode_){
	switch(Node.Get_BcNormal())
			{
			case 2:
//apply symmetry
				if(TypeOfNode_[Node.Get_connect()[3]]==Solid)//left
				{
					f_in->f[1][Node.Get_index()]=f_in->f[3][Node.Get_index()];
					f_in->f[8][Node.Get_index()]=f_in->f[7][Node.Get_index()];
				}
				else //right
				{
					f_in->f[3][Node.Get_index()]=f_in->f[1][Node.Get_index()];
					f_in->f[7][Node.Get_index()]=f_in->f[8][Node.Get_index()];
				}
				break;
			case 4:
//apply symmetry
				if(TypeOfNode_[Node.Get_connect()[3]]==Solid)//left
				{
					f_in->f[1][Node.Get_index()]=f_in->f[3][Node.Get_index()];
					f_in->f[5][Node.Get_index()]=f_in->f[6][Node.Get_index()];
				}
				else //right
				{
					f_in->f[3][Node.Get_index()]=f_in->f[1][Node.Get_index()];
					f_in->f[6][Node.Get_index()]=f_in->f[5][Node.Get_index()];
				}
				break;
			case 1:
//apply symmetry
				if(TypeOfNode_[Node.Get_connect()[2]]==Solid)//top
				{
					f_in->f[4][Node.Get_index()]=f_in->f[2][Node.Get_index()];
					f_in->f[7][Node.Get_index()]=f_in->f[6][Node.Get_index()];
				}
				else //bottom
				{
					f_in->f[2][Node.Get_index()]=f_in->f[4][Node.Get_index()];
					f_in->f[6][Node.Get_index()]=f_in->f[7][Node.Get_index()];
				}
				break;
			case 3:
//apply symmetry
				if(TypeOfNode_[Node.Get_connect()[2]]==Solid)//top
				{
					f_in->f[4][Node.Get_index()]=f_in->f[2][Node.Get_index()];
					f_in->f[8][Node.Get_index()]=f_in->f[5][Node.Get_index()];
				}
				else //bottom
				{
					f_in->f[2][Node.Get_index()]=f_in->f[4][Node.Get_index()];
					f_in->f[7][Node.Get_index()]=f_in->f[6][Node.Get_index()];
				}
				break;
			default:
				std::cerr<<"Direction wall bounce back not found. Index: "<<Node.Get_index()<<std::endl;
				break;
			}
	f_in->f[Node.Get_BcNormal()][Node.Get_connect()[0]]=f_in->f[OppositeBc[Node.Get_BcNormal()]][Node.Get_connect()[0]];
	f_in->f[BounceBackWallConnect[Node.Get_BcNormal()][0]][Node.Get_connect()[0]]=f_in->f[OppositeBc[BounceBackWallConnect[Node.Get_BcNormal()][0]]][Node.Get_connect()[0]];
	f_in->f[BounceBackWallConnect[Node.Get_BcNormal()][1]][Node.Get_connect()[0]]=f_in->f[OppositeBc[BounceBackWallConnect[Node.Get_BcNormal()][1]]][Node.Get_connect()[0]];
}

void D2Q9Wall::ApplyHalfWayBounceBackSymmetry(NodeWall2D & Node, DistriFunct* f_in, unsigned int idxDistribution, std::map<int,NodeType>&  TypeOfNode_){
	switch(Node.Get_BcNormal())
			{
			case 2:
//apply symmetry
				if(TypeOfNode_[Node.Get_connect()[3]]==Solid)//left
				{
					f_in->f[1][Node.Get_index()]=f_in->f[3][Node.Get_index()];
					f_in->f[8][Node.Get_index()]=f_in->f[7][Node.Get_index()];
				}
				else //right
				{
					f_in->f[3][Node.Get_index()]=f_in->f[1][Node.Get_index()];
					f_in->f[7][Node.Get_index()]=f_in->f[8][Node.Get_index()];
				}
				break;
			case 4:
//apply symmetry
				if(TypeOfNode_[Node.Get_connect()[3]]==Solid)//left
				{
					f_in->f[1][Node.Get_index()]=f_in->f[3][Node.Get_index()];
					f_in->f[5][Node.Get_index()]=f_in->f[6][Node.Get_index()];
				}
				else //right
				{
					f_in->f[3][Node.Get_index()]=f_in->f[1][Node.Get_index()];
					f_in->f[6][Node.Get_index()]=f_in->f[5][Node.Get_index()];
				}
				break;
			case 1:
//apply symmetry
				if(TypeOfNode_[Node.Get_connect()[2]]==Solid)//top
				{
					f_in->f[4][Node.Get_index()]=f_in->f[2][Node.Get_index()];
					f_in->f[7][Node.Get_index()]=f_in->f[6][Node.Get_index()];
				}
				else //bottom
				{
					f_in->f[2][Node.Get_index()]=f_in->f[4][Node.Get_index()];
					f_in->f[6][Node.Get_index()]=f_in->f[7][Node.Get_index()];
				}
				break;
			case 3:
//apply symmetry
				if(TypeOfNode_[Node.Get_connect()[2]]==Solid)//top
				{
					f_in->f[4][Node.Get_index()]=f_in->f[2][Node.Get_index()];
					f_in->f[8][Node.Get_index()]=f_in->f[5][Node.Get_index()];
				}
				else //bottom
				{
					f_in->f[2][Node.Get_index()]=f_in->f[4][Node.Get_index()];
					f_in->f[7][Node.Get_index()]=f_in->f[6][Node.Get_index()];
				}
				break;
			default:
				std::cerr<<"Direction wall bounce back not found. Index: "<<Node.Get_index()<<std::endl;
				break;
			}
	f_in->f[Node.Get_BcNormal()][Node.Get_connect()[0]]=f_in->f[OppositeBc[Node.Get_BcNormal()]][Node.Get_connect()[OppositeBc[Node.Get_BcNormal()]]];
	f_in->f[BounceBackWallConnect[Node.Get_BcNormal()][0]][Node.Get_connect()[0]]=f_in->f[OppositeBc[BounceBackWallConnect[Node.Get_BcNormal()][0]]][Node.Get_connect()[OppositeBc[BounceBackWallConnect[Node.Get_BcNormal()][0]]]];
	f_in->f[BounceBackWallConnect[Node.Get_BcNormal()][1]][Node.Get_connect()[0]]=f_in->f[OppositeBc[BounceBackWallConnect[Node.Get_BcNormal()][1]]][Node.Get_connect()[OppositeBc[BounceBackWallConnect[Node.Get_BcNormal()][1]]]];
}
*/
