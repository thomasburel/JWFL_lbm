/*
 * ============================================================================
 * D2Q9GlobalCorner.cpp
 *
 *  Created on: 2 Sep 2016
 *      Author: Thomas Burel
 *      Version     :
 *      Copyright   : Your copyright notice
 *      Description : LBM Code Developing by JWFL at University of Strathclyde
 * ============================================================================
 */

#include "D2Q9GlobalCorner.h"

D2Q9GlobalCorner::D2Q9GlobalCorner() {
	// TODO Auto-generated constructor stub
	BcMethods=0;
}

D2Q9GlobalCorner::~D2Q9GlobalCorner() {
	// TODO Auto-generated destructor stub
}

void D2Q9GlobalCorner::Set_GlobalCorner(Parameters *Param,D2Q9GenericBc* D2Q9GenericBc ){
	BcMethods=D2Q9GenericBc;
}

void D2Q9GlobalCorner::ApplyGlobalCorner(NodeCorner2D& Node, std::map<int,NodeType> TypeOfNode_, DistriFunct* f_in)
{
	FunctionGlobalCorner(Node,Node.Get_RhoDef(),Node.Get_UDef(),TypeOfNode_,f_in,BcMethods->Get_Rho(),BcMethods->Get_V(),BcMethods->Get_V());
}
void D2Q9GlobalCorner::ApplyGlobalCorner(NodeCorner2D& Node, double const Rho_def, double const *UDef, std::map<int,NodeType> TypeOfNode_, DistriFunct* f_in,double * & Rho, double * &U, double * &V)
{
	FunctionGlobalCorner(Node,Rho_def,UDef,TypeOfNode_,f_in,Rho,U,V);
}
void D2Q9GlobalCorner::FunctionGlobalCorner(NodeCorner2D& Node, double const Rho_def, double const *UDef, std::map<int,NodeType> &TypeOfNode_, DistriFunct* &f_in,double * & Rho, double * &U, double * &V)
{

	NodeType NodeTypeTmp1;
	NodeType NodeTypeTmp2;
	int normal=0;
	switch(Node.Get_BcNormal())
	{
	//Bottom left corner
	case 5:
		NodeTypeTmp1=TypeOfNode_[Node.Get_connect()[2]];
		NodeTypeTmp2=TypeOfNode_[Node.Get_connect()[1]];
		if(NodeTypeTmp1==Symmetry||NodeTypeTmp2==Symmetry)
		{
			if(NodeTypeTmp1!=Symmetry)
			{
			//Apply symmetry on the bottom side
				normal=2;
				BcMethods->ApplySymmetry(normal,Node.Get_connect(),Rho_def,UDef,f_in,Rho,U,V);
			//Apply Other boundary on the left side
				normal=1;
			}
			else
			{
			//Apply symmetry on the left side
				normal=1;
				BcMethods->ApplySymmetry(normal,Node.Get_connect(),Rho_def,UDef,f_in,Rho,U,V);
			//Apply Other boundary on the bottom side
				normal=2;
			}
			switch(NodeTypeTmp2)
			{
			case Pressure:
				BcMethods->ApplyPressure(normal,Node.Get_connect(),Rho_def,f_in,Rho,U,V);
				break;
			case Velocity:
				BcMethods->ApplyVelocity(normal,Node.Get_connect(),UDef,f_in,Rho,U,V);
				break;
			case Symmetry:
				BcMethods->ApplySymmetry(normal,Node.Get_connect(),Rho_def,UDef,f_in,Rho,U,V);
				break;
			case Wall:
				BcMethods->ApplyWall(normal,Node.Get_connect(),f_in,Rho,U,V);
				break;
			case Periodic:
				BcMethods->ApplyPeriodic(normal,Node.Get_connect(),Rho_def,UDef,f_in,Rho,U,V);
				break;
			default:
				std::cerr<<"Global corner condition with one symmetry is not found"<<std::endl;
			}
		}
		else if(NodeTypeTmp1==Periodic||NodeTypeTmp2==Periodic)
		{
			if(NodeTypeTmp1!=Periodic)
			{
			//Apply symmetry on the bottom side
				normal=2;
				BcMethods->ApplyPeriodic(normal,Node.Get_connect(),Rho_def,UDef,f_in,Rho,U,V);
			//Apply Other boundary on the left side
				normal=1;
			}
			else
			{
			//Apply symmetry on the left side
				normal=1;
				BcMethods->ApplyPeriodic(normal,Node.Get_connect(),Rho_def,UDef,f_in,Rho,U,V);
			//Apply Other boundary on the bottom side
				normal=2;
			}
			switch(NodeTypeTmp2)
			{
			case Pressure:
				BcMethods->ApplyPressure(normal,Node.Get_connect(),Rho_def,f_in,Rho,U,V);
				break;
			case Velocity:
				BcMethods->ApplyVelocity(normal,Node.Get_connect(),UDef,f_in,Rho,U,V);
				break;
			case Symmetry:
				BcMethods->ApplySymmetry(normal,Node.Get_connect(),Rho_def,UDef,f_in,Rho,U,V);
				break;
			case Wall:
				BcMethods->ApplyWall(normal,Node.Get_connect(),f_in,Rho,U,V);
				break;
			case Periodic:
				BcMethods->ApplyPeriodic(normal,Node.Get_connect(),Rho_def,UDef,f_in,Rho,U,V);
				break;
			default:
				std::cerr<<"Global corner condition with one periodic is not found"<<std::endl;
			}
		}
		else
		{
			if(NodeTypeTmp1==Wall||NodeTypeTmp2==Wall)
			{
				if(NodeTypeTmp1==Wall)
					switch(NodeTypeTmp2)
					{
					case Pressure:
						BcMethods->ApplyCorner(Node,f_in,Rho_def,BcMethods->Get_U(Node.Get_connect()[2]),BcMethods->Get_V(Node.Get_connect()[2]),Rho,U,V);
						break;
					case Velocity:
						BcMethods->ApplyCorner(Node,f_in,BcMethods->Get_Rho(Node.Get_connect()[1]),UDef[0],UDef[1],Rho,U,V);
						break;
					case Wall:
						BcMethods->ApplyCornerWall(Node,f_in,Rho,U,V);
					default:
						std::cerr<<"Global corner condition with one wall is not found"<<std::endl;
					}
				else
					switch(NodeTypeTmp1)
						{
					case Pressure:
						BcMethods->ApplyCorner(Node,f_in,Rho_def,BcMethods->Get_U(Node.Get_connect()[1]),BcMethods->Get_V(Node.Get_connect()[1]),Rho,U,V);
						break;
					case Velocity:
						BcMethods->ApplyCorner(Node,f_in,BcMethods->Get_Rho(Node.Get_connect()[2]),UDef[0],UDef[1],Rho,U,V);
						break;
						default:
							std::cerr<<"Global corner condition with one wall is not found"<<std::endl;
						}
			}
			else
				BcMethods->ApplyCorner(Node,f_in,Rho_def,UDef[0],UDef[1],Rho,U,V);
		}
		break;
	case 6:
		NodeTypeTmp1=TypeOfNode_[Node.Get_connect()[2]];
		NodeTypeTmp2=TypeOfNode_[Node.Get_connect()[3]];
		if(NodeTypeTmp1==Symmetry||NodeTypeTmp2==Symmetry)
		{
			if(NodeTypeTmp1!=Symmetry)
			{
			//Apply symmetry on the bottom side
				normal=2;
				BcMethods->ApplySymmetry(normal,Node.Get_connect(),Rho_def,UDef,f_in,Rho,U,V);
			//Apply Other boundary on the right side
				normal=3;
			}
			else
			{
			//Apply symmetry on the right side
				normal=3;
				BcMethods->ApplySymmetry(normal,Node.Get_connect(),Rho_def,UDef,f_in,Rho,U,V);
			//Apply Other boundary on the bottom side
				normal=2;
			}
			switch(NodeTypeTmp2)
			{
			case Pressure:
				BcMethods->ApplyPressure(normal,Node.Get_connect(),Rho_def,f_in,Rho,U,V);
				break;
			case Velocity:
				BcMethods->ApplyVelocity(normal,Node.Get_connect(),UDef,f_in,Rho,U,V);
				break;
			case Symmetry:
				BcMethods->ApplySymmetry(normal,Node.Get_connect(),Rho_def,UDef,f_in,Rho,U,V);
				break;
			case Wall:
				BcMethods->ApplyWall(normal,Node.Get_connect(),f_in,Rho,U,V);
				break;
			case Periodic:
				BcMethods->ApplyPeriodic(normal,Node.Get_connect(),Rho_def,UDef,f_in,Rho,U,V);
				break;
			default:
				std::cerr<<"Global corner condition with one symmetry is not found"<<std::endl;
			}
		}
		else if(NodeTypeTmp1==Periodic||NodeTypeTmp2==Periodic)
		{
			if(NodeTypeTmp1!=Periodic)
			{
			//Apply symmetry on the bottom side
				normal=2;
				BcMethods->ApplyPeriodic(normal,Node.Get_connect(),Rho_def,UDef,f_in,Rho,U,V);
			//Apply Other boundary on the right side
				normal=3;
			}
			else
			{
			//Apply symmetry on the right side
				normal=3;
				BcMethods->ApplyPeriodic(normal,Node.Get_connect(),Rho_def,UDef,f_in,Rho,U,V);
			//Apply Other boundary on the bottom side
				normal=2;
			}
			switch(NodeTypeTmp2)
			{
			case Pressure:
				BcMethods->ApplyPressure(normal,Node.Get_connect(),Rho_def,f_in,Rho,U,V);
				break;
			case Velocity:
				BcMethods->ApplyVelocity(normal,Node.Get_connect(),UDef,f_in,Rho,U,V);
				break;
			case Symmetry:
				BcMethods->ApplySymmetry(normal,Node.Get_connect(),Rho_def,UDef,f_in,Rho,U,V);
				break;
			case Wall:
				BcMethods->ApplyWall(normal,Node.Get_connect(),f_in,Rho,U,V);
				break;
			case Periodic:
				BcMethods->ApplyPeriodic(normal,Node.Get_connect(),Rho_def,UDef,f_in,Rho,U,V);
				break;
			default:
				std::cerr<<"Global corner condition with one periodic is not found"<<std::endl;
			}
		}
		else
		{
			if(NodeTypeTmp1==Wall||NodeTypeTmp2==Wall)
			{
				if(NodeTypeTmp1==Wall)
					switch(NodeTypeTmp2)
					{
					case Pressure:
						BcMethods->ApplyCorner(Node,f_in,Rho_def,BcMethods->Get_U(Node.Get_connect()[2]),BcMethods->Get_V(Node.Get_connect()[2]),Rho,U,V);
						break;
					case Velocity:
						BcMethods->ApplyCorner(Node,f_in,BcMethods->Get_Rho(Node.Get_connect()[3]),UDef[0],UDef[1],Rho,U,V);
						break;
					case Wall:
						BcMethods->ApplyCornerWall(Node,f_in,Rho,U,V);
					default:
						std::cerr<<"Global corner condition with one wall is not found"<<std::endl;
					}
				else
					switch(NodeTypeTmp1)
						{
					case Pressure:
						BcMethods->ApplyCorner(Node,f_in,Rho_def,BcMethods->Get_U(Node.Get_connect()[3]),BcMethods->Get_V(Node.Get_connect()[3]),Rho,U,V);
						break;
					case Velocity:
						BcMethods->ApplyCorner(Node,f_in,BcMethods->Get_Rho(Node.Get_connect()[2]),UDef[0],UDef[1],Rho,U,V);
						break;
						default:
							std::cerr<<"Global corner condition with one wall is not found"<<std::endl;
						}
			}
			else
				BcMethods->ApplyCorner(Node,f_in,Rho_def,UDef[0],UDef[1],Rho,U,V);
		}

		break;
	case 7:
		NodeTypeTmp1=TypeOfNode_[Node.Get_connect()[4]];
		NodeTypeTmp2=TypeOfNode_[Node.Get_connect()[3]];
		if(NodeTypeTmp1==Symmetry||NodeTypeTmp2==Symmetry)
		{
			if(NodeTypeTmp1!=Symmetry)
			{
			//Apply symmetry on the top side
				normal=4;
				BcMethods->ApplySymmetry(normal,Node.Get_connect(),Rho_def,UDef,f_in,Rho,U,V);
			//Apply Other boundary on the right side
				normal=3;
			}
			else
			{
			//Apply symmetry on the right side
				normal=3;
				BcMethods->ApplySymmetry(normal,Node.Get_connect(),Rho_def,UDef,f_in,Rho,U,V);
			//Apply Other boundary on the top side
				normal=4;
			}
			switch(NodeTypeTmp2)
			{
			case Pressure:
				BcMethods->ApplyPressure(normal,Node.Get_connect(),Rho_def,f_in,Rho,U,V);
				break;
			case Velocity:
				BcMethods->ApplyVelocity(normal,Node.Get_connect(),UDef,f_in,Rho,U,V);
				break;
			case Symmetry:
				BcMethods->ApplySymmetry(normal,Node.Get_connect(),Rho_def,UDef,f_in,Rho,U,V);
				break;
			case Wall:
				BcMethods->ApplyWall(normal,Node.Get_connect(),f_in,Rho,U,V);
				break;
			case Periodic:
				BcMethods->ApplyPeriodic(normal,Node.Get_connect(),Rho_def,UDef,f_in,Rho,U,V);
				break;
			default:
				std::cerr<<"Global corner condition with one symmetry is not found"<<std::endl;
			}
		}
		else if(NodeTypeTmp1==Periodic||NodeTypeTmp2==Periodic)
		{
			if(NodeTypeTmp1!=Periodic)
			{
			//Apply symmetry on the top side
				normal=4;
				BcMethods->ApplyPeriodic(normal,Node.Get_connect(),Rho_def,UDef,f_in,Rho,U,V);
			//Apply Other boundary on the right side
				normal=3;
			}
			else
			{
			//Apply symmetry on the right side
				normal=4;
				BcMethods->ApplyPeriodic(normal,Node.Get_connect(),Rho_def,UDef,f_in,Rho,U,V);
			//Apply Other boundary on the top side
				normal=3;
			}
			switch(NodeTypeTmp2)
			{
			case Pressure:
				BcMethods->ApplyPressure(normal,Node.Get_connect(),Rho_def,f_in,Rho,U,V);
				break;
			case Velocity:
				BcMethods->ApplyVelocity(normal,Node.Get_connect(),UDef,f_in,Rho,U,V);
				break;
			case Symmetry:
				BcMethods->ApplySymmetry(normal,Node.Get_connect(),Rho_def,UDef,f_in,Rho,U,V);
				break;
			case Wall:
				BcMethods->ApplyWall(normal,Node.Get_connect(),f_in,Rho,U,V);
				break;
			case Periodic:
				BcMethods->ApplyPeriodic(normal,Node.Get_connect(),Rho_def,UDef,f_in,Rho,U,V);
				break;
			default:
				std::cerr<<"Global corner condition with one periodic is not found"<<std::endl;
			}
		}
		else
		{
			if(NodeTypeTmp1==Wall||NodeTypeTmp2==Wall)
			{
				if(NodeTypeTmp1==Wall)
					switch(NodeTypeTmp2)
					{
					case Pressure:
						BcMethods->ApplyCorner(Node,f_in,Rho_def,BcMethods->Get_U(Node.Get_connect()[4]),BcMethods->Get_V(Node.Get_connect()[4]),Rho,U,V);
						break;
					case Velocity:
						BcMethods->ApplyCorner(Node,f_in,BcMethods->Get_Rho(Node.Get_connect()[3]),UDef[0],UDef[1],Rho,U,V);
						break;
					case Wall:
						BcMethods->ApplyCornerWall(Node,f_in,Rho,U,V);
					default:
						std::cerr<<"Global corner condition with one wall is not found"<<std::endl;
					}
				else
					switch(NodeTypeTmp1)
						{
					case Pressure:
						BcMethods->ApplyCorner(Node,f_in,Rho_def,BcMethods->Get_U(Node.Get_connect()[3]),BcMethods->Get_V(Node.Get_connect()[3]),Rho,U,V);
						break;
					case Velocity:
						BcMethods->ApplyCorner(Node,f_in,BcMethods->Get_Rho(Node.Get_connect()[4]),UDef[0],UDef[1],Rho,U,V);
						break;
						default:
							std::cerr<<"Global corner condition with one wall is not found"<<std::endl;
						}
			}
			else
				BcMethods->ApplyCorner(Node,f_in,Rho_def,UDef[0],UDef[1],Rho,U,V);
		}
		break;
	case 8:
		NodeTypeTmp1=TypeOfNode_[Node.Get_connect()[4]];
		NodeTypeTmp2=TypeOfNode_[Node.Get_connect()[1]];
		if(NodeTypeTmp1==Symmetry||NodeTypeTmp2==Symmetry)
		{
			if(NodeTypeTmp1!=Symmetry)
			{
			//Apply symmetry on the top side
				normal=4;
				BcMethods->ApplySymmetry(normal,Node.Get_connect(),Rho_def,UDef,f_in,Rho,U,V);
			//Apply Other boundary on the left side
				normal=1;
			}
			else
			{
			//Apply symmetry on the left side
				normal=1;
				BcMethods->ApplySymmetry(normal,Node.Get_connect(),Rho_def,UDef,f_in,Rho,U,V);
			//Apply Other boundary on the top side
				normal=4;
			}
			switch(NodeTypeTmp2)
			{
			case Pressure:
				BcMethods->ApplyPressure(normal,Node.Get_connect(),Rho_def,f_in,Rho,U,V);
				break;
			case Velocity:
				BcMethods->ApplyVelocity(normal,Node.Get_connect(),UDef,f_in,Rho,U,V);
				break;
			case Symmetry:
				BcMethods->ApplySymmetry(normal,Node.Get_connect(),Rho_def,UDef,f_in,Rho,U,V);
				break;
			case Wall:
				BcMethods->ApplyWall(normal,Node.Get_connect(),f_in,Rho,U,V);
				break;
			case Periodic:
				BcMethods->ApplyPeriodic(normal,Node.Get_connect(),Rho_def,UDef,f_in,Rho,U,V);
				break;
			default:
				std::cerr<<"Global corner condition with one symmetry is not found"<<std::endl;
			}
		}
		else if(NodeTypeTmp1==Periodic||NodeTypeTmp2==Periodic)
		{
			if(NodeTypeTmp1!=Periodic)
			{
			//Apply symmetry on the top side
				normal=4;
				BcMethods->ApplyPeriodic(normal,Node.Get_connect(),Rho_def,UDef,f_in,Rho,U,V);
			//Apply Other boundary on the left side
				normal=1;
			}
			else
			{
			//Apply symmetry on the left side
				normal=1;
				BcMethods->ApplyPeriodic(normal,Node.Get_connect(),Rho_def,UDef,f_in,Rho,U,V);
			//Apply Other boundary on the top side
				normal=4;
			}
			switch(NodeTypeTmp2)
			{
			case Pressure:
				BcMethods->ApplyPressure(normal,Node.Get_connect(),Rho_def,f_in,Rho,U,V);
				break;
			case Velocity:
				BcMethods->ApplyVelocity(normal,Node.Get_connect(),UDef,f_in,Rho,U,V);
				break;
			case Symmetry:
				BcMethods->ApplySymmetry(normal,Node.Get_connect(),Rho_def,UDef,f_in,Rho,U,V);
				break;
			case Wall:
				BcMethods->ApplyWall(normal,Node.Get_connect(),f_in,Rho,U,V);
				break;
			case Periodic:
				BcMethods->ApplyPeriodic(normal,Node.Get_connect(),Rho_def,UDef,f_in,Rho,U,V);
				break;
			default:
				std::cerr<<"Global corner condition with one periodic is not found"<<std::endl;
			}
		}
		else
		{
			if(NodeTypeTmp1==Wall||NodeTypeTmp2==Wall)
			{
				if(NodeTypeTmp1==Wall)
					switch(NodeTypeTmp2)
					{
					case Pressure:
						BcMethods->ApplyCorner(Node,f_in,Rho_def,BcMethods->Get_U(Node.Get_connect()[4]),BcMethods->Get_V(Node.Get_connect()[4]),Rho,U,V);
						break;
					case Velocity:
						BcMethods->ApplyCorner(Node,f_in,BcMethods->Get_Rho(Node.Get_connect()[1]),UDef[0],UDef[1],Rho,U,V);
						break;
					case Wall:
						BcMethods->ApplyCornerWall(Node,f_in,Rho,U,V);
					default:
						std::cerr<<"Global corner condition with one wall is not found"<<std::endl;
					}
				else
					switch(NodeTypeTmp1)
						{
					case Pressure:
						BcMethods->ApplyCorner(Node,f_in,Rho_def,BcMethods->Get_U(Node.Get_connect()[1]),BcMethods->Get_V(Node.Get_connect()[1]),Rho,U,V);
						break;
					case Velocity:
						BcMethods->ApplyCorner(Node,f_in,BcMethods->Get_Rho(Node.Get_connect()[4]),UDef[0],UDef[1],Rho,U,V);
						break;
						default:
							std::cerr<<"Global corner condition with one wall is not found"<<std::endl;
						}
			}
			else
				BcMethods->ApplyCorner(Node,f_in,Rho_def,UDef[0],UDef[1],Rho,U,V);
		}
		break;
	default:
		std::cerr<<"Direction: "<< Node.Get_BcNormal()<<" not found for Global Corner"<<std::endl;
	}


}


