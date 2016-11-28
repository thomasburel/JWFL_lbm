/*
 * ============================================================================
 * D2Q9SpecialWall.cpp
 *
 *  Created on: 2 Sep 2016
 *      Author: Thomas Burel
 *      Version     :
 *      Copyright   : Your copyright notice
 *      Description : LBM Code Developing by JWFL at University of Strathclyde
 * ============================================================================
 */

#include "D2Q9SpecialWall.h"

D2Q9SpecialWall::D2Q9SpecialWall() {
	// TODO Auto-generated constructor stub
	BcMethods=0;
}

D2Q9SpecialWall::~D2Q9SpecialWall() {
	// TODO Auto-generated destructor stub
}

void D2Q9SpecialWall::Set_SpecialWall(Parameters *Param,D2Q9GenericBc* D2Q9GenericBc ){
	BcMethods=D2Q9GenericBc;
}

/*void D2Q9SpecialWall::ApplySpecialWall(NodeWall2D& Node, std::map<int,NodeType> TypeOfNode_, DistriFunct* f_in)
{
	FunctionSpecialWall(Node,Node.Get_RhoDef(),Node.Get_UDef(),TypeOfNode_,f_in,BcMethods->Get_Rho(),BcMethods->Get_U(),BcMethods->Get_V());
}*/
void D2Q9SpecialWall::ApplySpecialWall(NodeWall2D& Node, double const Rho_def, double const UDef, double const VDef, std::map<int,NodeType> TypeOfNode_, DistriFunct* f_in,double * & Rho, double * &U, double * &V)
{
	U_tmp[0]=UDef;U_tmp[1]=VDef;
	FunctionSpecialWall(Node,Rho_def,&U_tmp[0],TypeOfNode_,f_in,Rho,U,V);
}
void D2Q9SpecialWall::FunctionSpecialWall(NodeWall2D& Node, double const Rho_def, double const *UDef, std::map<int,NodeType> &TypeOfNode_, DistriFunct* &f_in,double * & Rho, double * &U, double * &V)
{

	NodeType NodeTypeTmp1;
	NodeType NodeTypeTmp2;
	int normal=0;
	switch(Node.Get_BcNormal())
	{
	//Bottom left corner
	case 5:
		NodeTypeTmp1=TypeOfNode_[Node.Get_connect()[1]];
		NodeTypeTmp2=TypeOfNode_[Node.Get_connect()[2]];
		if(NodeTypeTmp1==Symmetry||NodeTypeTmp2==Symmetry)
		{
			if(NodeTypeTmp1!=Symmetry)
			{
			//Apply symmetry on the bottom side
				normal=1;
				BcMethods->ApplySymmetry(normal,Node.Get_connect(),Rho_def,UDef,f_in,Rho,U,V);
			//Apply Other boundary on the left side
				normal=2;
				BcMethods->ApplyWall(normal,Node.Get_connect(),f_in,Rho,U,V);
			}
			else
			{
			//Apply symmetry on the left side
				normal=2;
				BcMethods->ApplySymmetry(normal,Node.Get_connect(),Rho_def,UDef,f_in,Rho,U,V);
			//Apply Other boundary on the bottom side
				normal=1;
				BcMethods->ApplyWall(normal,Node.Get_connect(),f_in,Rho,U,V);
			}
		}
		else if(NodeTypeTmp1==Periodic||NodeTypeTmp2==Periodic)
		{
			if(NodeTypeTmp1!=Periodic)
			{
			//Apply symmetry on the bottom side
				normal=1;
				BcMethods->ApplyPeriodic(normal,Node.Get_connect(),Rho_def,UDef,f_in,Rho,U,V);
			//Apply Other boundary on the left side
				normal=2;
				BcMethods->ApplyWall(normal,Node.Get_connect(),f_in,Rho,U,V);
			}
			else
			{
			//Apply symmetry on the left side
				normal=2;
				BcMethods->ApplyPeriodic(normal,Node.Get_connect(),Rho_def,UDef,f_in,Rho,U,V);
			//Apply Other boundary on the bottom side
				normal=1;
				BcMethods->ApplyWall(normal,Node.Get_connect(),f_in,Rho,U,V);
			}
		}
		else
			if(NodeTypeTmp1==Wall)
			{
				switch(NodeTypeTmp2)
				{
				case Pressure:
					BcMethods->ApplyPreVelSpecialWall(Node,f_in,Rho_def,BcMethods->Get_U(Node.Get_connect()[1]),BcMethods->Get_V(Node.Get_connect()[1]),Rho,U,V);
					break;
				case Velocity:
					BcMethods->ApplyPreVelSpecialWall(Node,f_in,BcMethods->Get_Rho(Node.Get_connect()[2]),UDef[0],UDef[1],Rho,U,V);
					break;
				case Wall:
					BcMethods->ApplyCornerSpecialWall(Node,f_in,Rho,U,V);
					break;
				default:
					std::cerr<<"Special Wall condition with one wall is not found. Direction of the special wall: "<<Node.Get_BcNormal()<<" Type of node 1: "<<NodeTypeTmp1<<" Type of node 2: "<<NodeTypeTmp2<<std::endl;
				}
			}
			else
				switch(NodeTypeTmp1)
				{
				case Pressure:
					BcMethods->ApplyPreVelSpecialWall(Node,f_in,Rho_def,BcMethods->Get_U(Node.Get_connect()[2]),BcMethods->Get_V(Node.Get_connect()[2]),Rho,U,V);
					break;
				case Velocity:
					BcMethods->ApplyPreVelSpecialWall(Node,f_in,BcMethods->Get_Rho(Node.Get_connect()[1]),UDef[0],UDef[1],Rho,U,V);
					break;
				default:
					std::cerr<<"Special Wall condition with one wall is not found. Direction of the special wall: "<<Node.Get_BcNormal()<<" Type of node 1: "<<NodeTypeTmp1<<" Type of node 2: "<<NodeTypeTmp2<<std::endl;
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
				BcMethods->ApplyWall(normal,Node.Get_connect(),f_in,Rho,U,V);
			}
			else
			{
			//Apply symmetry on the right side
				normal=3;
				BcMethods->ApplySymmetry(normal,Node.Get_connect(),Rho_def,UDef,f_in,Rho,U,V);
			//Apply Other boundary on the bottom side
				normal=2;
				BcMethods->ApplyWall(normal,Node.Get_connect(),f_in,Rho,U,V);
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
				BcMethods->ApplyWall(normal,Node.Get_connect(),f_in,Rho,U,V);
			}
			else
			{
			//Apply symmetry on the right side
				normal=3;
				BcMethods->ApplyPeriodic(normal,Node.Get_connect(),Rho_def,UDef,f_in,Rho,U,V);
			//Apply Other boundary on the bottom side
				normal=2;
				BcMethods->ApplyWall(normal,Node.Get_connect(),f_in,Rho,U,V);
			}

		}
		else
			if(NodeTypeTmp1==Wall)
				switch(NodeTypeTmp2)
				{
				case Pressure:
					BcMethods->ApplyPreVelSpecialWall(Node,f_in,Rho_def,BcMethods->Get_U(Node.Get_connect()[2]),BcMethods->Get_V(Node.Get_connect()[2]),Rho,U,V);
					break;
				case Velocity:
					BcMethods->ApplyPreVelSpecialWall(Node,f_in,BcMethods->Get_Rho(Node.Get_connect()[3]),UDef[0],UDef[1],Rho,U,V);
					break;
				case Wall:
					BcMethods->ApplyCornerSpecialWall(Node,f_in,Rho,U,V);
					break;
				default:
					std::cerr<<"Special Wall condition with one wall is not found. Direction of the special wall: "<<Node.Get_BcNormal()<<" Type of node 1: "<<NodeTypeTmp1<<" Type of node 2: "<<NodeTypeTmp2<<std::endl;
				}
			else
				switch(NodeTypeTmp1)
				{
				case Pressure:
					BcMethods->ApplyPreVelSpecialWall(Node,f_in,Rho_def,BcMethods->Get_U(Node.Get_connect()[3]),BcMethods->Get_V(Node.Get_connect()[3]),Rho,U,V);
					break;
				case Velocity:
					BcMethods->ApplyPreVelSpecialWall(Node,f_in,BcMethods->Get_Rho(Node.Get_connect()[2]),UDef[0],UDef[1],Rho,U,V);
					break;
				default:
					std::cerr<<"Special Wall condition with one wall is not found. Direction of the special wall: "<<Node.Get_BcNormal()<<" Type of node 1: "<<NodeTypeTmp1<<" Type of node 2: "<<NodeTypeTmp2<<std::endl;
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
				BcMethods->ApplyWall(normal,Node.Get_connect(),f_in,Rho,U,V);
			}
			else
			{
			//Apply symmetry on the right side
				normal=3;
				BcMethods->ApplySymmetry(normal,Node.Get_connect(),Rho_def,UDef,f_in,Rho,U,V);
			//Apply Other boundary on the top side
				normal=4;
				BcMethods->ApplyWall(normal,Node.Get_connect(),f_in,Rho,U,V);
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
				BcMethods->ApplyWall(normal,Node.Get_connect(),f_in,Rho,U,V);
			}
			else
			{
			//Apply symmetry on the right side
				normal=3;
				BcMethods->ApplyPeriodic(normal,Node.Get_connect(),Rho_def,UDef,f_in,Rho,U,V);
			//Apply Other boundary on the top side
				normal=4;
				BcMethods->ApplyWall(normal,Node.Get_connect(),f_in,Rho,U,V);
			}

		}
		else
			if(NodeTypeTmp1==Wall)
				switch(NodeTypeTmp2)
				{
				case Pressure:
					BcMethods->ApplyPreVelSpecialWall(Node,f_in,Rho_def,BcMethods->Get_U(Node.Get_connect()[4]),BcMethods->Get_V(Node.Get_connect()[4]),Rho,U,V);
					break;
				case Velocity:
					BcMethods->ApplyPreVelSpecialWall(Node,f_in,BcMethods->Get_Rho(Node.Get_connect()[3]),UDef[0],UDef[1],Rho,U,V);
					break;
				case Wall:
					BcMethods->ApplyCornerSpecialWall(Node,f_in,Rho,U,V);
					break;
				default:
					std::cerr<<"Special Wall condition with one wall is not found. Direction of the special wall: "<<Node.Get_BcNormal()<<" Type of node 1: "<<NodeTypeTmp1<<" Type of node 2: "<<NodeTypeTmp2<<std::endl;
				}
			else
				switch(NodeTypeTmp1)
					{
				case Pressure:
					BcMethods->ApplyPreVelSpecialWall(Node,f_in,Rho_def,BcMethods->Get_U(Node.Get_connect()[3]),BcMethods->Get_V(Node.Get_connect()[3]),Rho,U,V);
					break;
				case Velocity:
					BcMethods->ApplyPreVelSpecialWall(Node,f_in,BcMethods->Get_Rho(Node.Get_connect()[4]),UDef[0],UDef[1],Rho,U,V);
					break;
				default:
					std::cerr<<"Special Wall condition with one wall is not found. Direction of the special wall: "<<Node.Get_BcNormal()<<" Type of node 1: "<<NodeTypeTmp1<<" Type of node 2: "<<NodeTypeTmp2<<std::endl;
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
				BcMethods->ApplyWall(normal,Node.Get_connect(),f_in,Rho,U,V);
			}
			else
			{
			//Apply symmetry on the left side
				normal=1;
				BcMethods->ApplySymmetry(normal,Node.Get_connect(),Rho_def,UDef,f_in,Rho,U,V);
			//Apply Other boundary on the top side
				normal=4;
				BcMethods->ApplyWall(normal,Node.Get_connect(),f_in,Rho,U,V);
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
				BcMethods->ApplyWall(normal,Node.Get_connect(),f_in,Rho,U,V);
			}
			else
			{
			//Apply symmetry on the left side
				normal=1;
				BcMethods->ApplyPeriodic(normal,Node.Get_connect(),Rho_def,UDef,f_in,Rho,U,V);
			//Apply Other boundary on the top side
				normal=4;
				BcMethods->ApplyWall(normal,Node.Get_connect(),f_in,Rho,U,V);
			}
		}
		else
			if(NodeTypeTmp1==Wall)
				switch(NodeTypeTmp2)
				{
				case Pressure:
					BcMethods->ApplyPreVelSpecialWall(Node,f_in,Rho_def,BcMethods->Get_U(Node.Get_connect()[4]),BcMethods->Get_V(Node.Get_connect()[4]),Rho,U,V);
					break;
				case Velocity:
					BcMethods->ApplyPreVelSpecialWall(Node,f_in,BcMethods->Get_Rho(Node.Get_connect()[1]),UDef[0],UDef[1],Rho,U,V);
					break;
				case Wall:
					BcMethods->ApplyCornerSpecialWall(Node,f_in,Rho,U,V);
					break;
				default:
					std::cerr<<"Special Wall condition with one wall is not found. Direction of the special wall: "<<Node.Get_BcNormal()<<" Type of node 1: "<<NodeTypeTmp1<<" Type of node 2: "<<NodeTypeTmp2<<std::endl;
				}
			else
				switch(NodeTypeTmp1)
					{
				case Pressure:
					BcMethods->ApplyPreVelSpecialWall(Node,f_in,Rho_def,BcMethods->Get_U(Node.Get_connect()[1]),BcMethods->Get_V(Node.Get_connect()[1]),Rho,U,V);
					break;
				case Velocity:
					BcMethods->ApplyPreVelSpecialWall(Node,f_in,BcMethods->Get_Rho(Node.Get_connect()[4]),UDef[0],UDef[1],Rho,U,V);
					break;
				default:
					std::cerr<<"Special Wall condition with one wall is not found. Direction of the special wall: "<<Node.Get_BcNormal()<<" Type of node 1: "<<NodeTypeTmp1<<" Type of node 2: "<<NodeTypeTmp2<<std::endl;
				}
		break;
	default:
		std::cerr<<"Direction: "<< Node.Get_BcNormal()<<" not found for Special Wall"<<std::endl;
	}


}


