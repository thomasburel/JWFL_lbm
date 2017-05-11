/*
 * NodeArray.cpp
 *
 *  Created on: 9 Dec 2015
 *      Author: thomas
 */

#include "NodeArrays.h"

NodeArrays::NodeArrays() {
	// TODO Auto-generated constructor stub
	nullNodeArray=0;
}

NodeArrays::~NodeArrays() {
//	TypeOfNode.clear();
//	if(!NodeIndexByType.empty())
//		NodeIndexByType.clear();
	// TODO Auto-generated destructor stub
}
NodeArrays2D::NodeArrays2D() {
	// TODO Auto-generated constructor stub

}

NodeArrays2D::~NodeArrays2D() {
	NodeInterior.clear();
	NodeSolid.clear();
	NodeGhost.clear();
	NodeCorner.clear();
	NodeGlobalCorner.clear();
	NodeWall.clear();
	NodeSpecialWall.clear();
	NodePeriodic.clear();
	NodeVelocity.clear();
	NodePressure.clear();
	NodeSymmetry.clear();
}


void NodeArrays2D::Get_coordinate(int index,double & x, double & y){
	if(index<TypeOfNode.size())
	{
		int localidx=NodeIndexByType[index];
		switch(TypeOfNode[index])
		{
		case Interior:
			x=NodeInterior[localidx].get_x();
			y=NodeInterior[localidx].get_y();
			break;
		case Solid:
			x=NodeSolid[localidx].get_x();
			y=NodeSolid[localidx].get_y();
			break;
		case Ghost:
			x=NodeGhost[localidx].get_x();
			y=NodeGhost[localidx].get_y();
			break;
		case Corner:
			x=NodeCorner[localidx].get_x();
			y=NodeCorner[localidx].get_y();
			break;
		case Wall:
			x=NodeWall[localidx].get_x();
			y=NodeWall[localidx].get_y();
			break;
		case Periodic:
			x=NodePeriodic[localidx].get_x();
			y=NodePeriodic[localidx].get_y();
			break;
		case Velocity:
			x=NodeVelocity[localidx].get_x();
			y=NodeVelocity[localidx].get_y();
			break;
		case Symmetry:
			x=NodeSymmetry[localidx].get_x();
			y=NodeSymmetry[localidx].get_y();
			break;
		case Pressure:
			x=NodePressure[localidx].get_x();
			y=NodePressure[localidx].get_y();
			break;
		case ConcaveCorner:
			x=NodeCorner[localidx].get_x();
			y=NodeCorner[localidx].get_y();
			break;
		case ConvexCorner:
			x=NodeCorner[localidx].get_x();
			y=NodeCorner[localidx].get_y();
			break;
		case SolidGhost:
			x=NodeGhost[localidx].get_x();
			y=NodeGhost[localidx].get_y();
			break;
		case GlobalCorner:
			x=NodeGlobalCorner[localidx].get_x();
			y=NodeGlobalCorner[localidx].get_y();
			break;
		case SpecialWall:
			x=NodeSpecialWall[localidx].get_x();
			y=NodeSpecialWall[localidx].get_y();
			break;
		}
	}
	else
	{
		std::cerr<<"Node number is higher than the number of nodes."<<std::endl;
		x=-10;
		y=-10;
	}
}
/*
void NodeArrays2D::ChangeNodeType(int NodeNumber,NodeType OldNodeType_, NodeType NewNodeType_){
		double x_tmp=0;
		double y_tmp=0;
		int index=NodeIndexByType[NodeNumber];
		Get_coordinate(index,x_tmp,y_tmp);
		unsigned int Connect_N_tmp=0;
		unsigned int Connect_S_tmp=0;
		unsigned int Connect_W_tmp=0;
		unsigned int Connect_E_tmp=0;
		unsigned int NbVelocity_tmp=0;
		double *U_tmp=0;
		double Rho_tmp=0;
		int* connect=0;

		switch(OldNodeType_)
		{
		case Interior:
			Connect_N_tmp=NodeInterior[index].Get_connect(2);
			Connect_S_tmp=NodeInterior[index].Get_connect(0);
			Connect_W_tmp=NodeInterior[index].Get_connect(3);
			Connect_E_tmp=NodeInterior[index].Get_connect(1);
			NbVelocity_tmp=NodeInterior[index].Get_NbVelocity();
			connect=NodeInterior[index].Get_connect();
			break;
		case Solid:
			Connect_N_tmp=NodeSolid[index].Get_connect(2);
			Connect_S_tmp=NodeSolid[index].Get_connect(0);
			Connect_W_tmp=NodeSolid[index].Get_connect(3);
			Connect_E_tmp=NodeSolid[index].Get_connect(1);
			NbVelocity_tmp=NodeSolid[index].Get_NbVelocity();
			connect=NodeSolid[index].Get_connect();
			break;
		case SolidGhost:
//			Node[NodeNumber]->Set_NodeType(SolidGhost);
			Connect_N_tmp=NodeGhost[index].Get_connect(2);
			Connect_S_tmp=NodeGhost[index].Get_connect(0);
			Connect_W_tmp=NodeGhost[index].Get_connect(3);
			Connect_E_tmp=NodeGhost[index].Get_connect(1);
			NbVelocity_tmp=NodeGhost[index].Get_NbVelocity();
			connect=NodeGhost[index].Get_connect();
			break;
		case Ghost:
			Connect_N_tmp=NodeGhost[index].Get_connect(2);
			Connect_S_tmp=NodeGhost[index].Get_connect(0);
			Connect_W_tmp=NodeGhost[index].Get_connect(3);
			Connect_E_tmp=NodeGhost[index].Get_connect(1);
			NbVelocity_tmp=NodeGhost[index].Get_NbVelocity();
			connect=NodeGhost[index].Get_connect();
			break;
		case Velocity:
			Connect_N_tmp=NodeVelocity[index].Get_connect(2);
			Connect_S_tmp=NodeVelocity[index].Get_connect(0);
			Connect_W_tmp=NodeVelocity[index].Get_connect(3);
			Connect_E_tmp=NodeVelocity[index].Get_connect(1);
			NbVelocity_tmp=NodeVelocity[index].Get_NbVelocity();
			connect=NodeVelocity[index].Get_connect();
			break;
		case Pressure:
			Connect_N_tmp=NodePressure[index].Get_connect(2);
			Connect_S_tmp=NodePressure[index].Get_connect(0);
			Connect_W_tmp=NodePressure[index].Get_connect(3);
			Connect_E_tmp=NodePressure[index].Get_connect(1);
			NbVelocity_tmp=NodePressure[index].Get_NbVelocity();
			connect=NodePressure[index].Get_connect();
			break;
		case Periodic:
			Connect_N_tmp=NodePeriodic[index].Get_connect(2);
			Connect_S_tmp=NodePeriodic[index].Get_connect(0);
			Connect_W_tmp=NodePeriodic[index].Get_connect(3);
			Connect_E_tmp=NodePeriodic[index].Get_connect(1);
			NbVelocity_tmp=NodePeriodic[index].Get_NbVelocity();
			connect=NodePeriodic[index].Get_connect();
			break;
		case Symmetry:
			Connect_N_tmp=NodeSymmetry[index].Get_connect(2);
			Connect_S_tmp=NodeSymmetry[index].Get_connect(0);
			Connect_W_tmp=NodeSymmetry[index].Get_connect(3);
			Connect_E_tmp=NodeSymmetry[index].Get_connect(1);
			NbVelocity_tmp=NodeSymmetry[index].Get_NbVelocity();
			connect=NodeSymmetry[index].Get_connect();
			break;
		case Corner:
			Connect_N_tmp=NodeCorner[index].Get_connect(2);
			Connect_S_tmp=NodeCorner[index].Get_connect(0);
			Connect_W_tmp=NodeCorner[index].Get_connect(3);
			Connect_E_tmp=NodeCorner[index].Get_connect(1);
			NbVelocity_tmp=NodeCorner[index].Get_NbVelocity();
			connect=NodeCorner[index].Get_connect();
			break;
		case ConcaveCorner:
			Connect_N_tmp=NodeCorner[index].Get_connect(2);
			Connect_S_tmp=NodeCorner[index].Get_connect(0);
			Connect_W_tmp=NodeCorner[index].Get_connect(3);
			Connect_E_tmp=NodeCorner[index].Get_connect(1);
			NbVelocity_tmp=NodeCorner[index].Get_NbVelocity();
			connect=NodeCorner[index].Get_connect();
			break;
		case ConvexCorner:
			Connect_N_tmp=NodeCorner[index].Get_connect(2);
			Connect_S_tmp=NodeCorner[index].Get_connect(0);
			Connect_W_tmp=NodeCorner[index].Get_connect(3);
			Connect_E_tmp=NodeCorner[index].Get_connect(1);
			NbVelocity_tmp=NodeCorner[index].Get_NbVelocity();
			connect=NodeCorner[index].Get_connect();
			break;
		case Wall:
			Connect_N_tmp=NodeWall[index].Get_connect(2);
			Connect_S_tmp=NodeWall[index].Get_connect(0);
			Connect_W_tmp=NodeWall[index].Get_connect(3);
			Connect_E_tmp=NodeWall[index].Get_connect(1);
			NbVelocity_tmp=NodeWall[index].Get_NbVelocity();
			connect=NodeWall[index].Get_connect();
			break;
		default:
			std::cerr<< "Node Type not found" << std::endl;
	        break;
		}


		switch(NewNodeType_)
		{
		case Interior:
			NodeInterior.push_back(NodeInterior2D(x_tmp, y_tmp));
			NodeInterior.back().Set_Connect(Connect_N_tmp,Connect_S_tmp,Connect_W_tmp,Connect_E_tmp);
			NodeInterior.back().Set_NbVelocity(NbVelocity_tmp);
			NodeInterior.back().Set_Index(index);
			NodeIndexByType[NodeNumber]=NodeInterior.size()-1;
			break;
		case Solid:
			Node[NodeNumber]=new NodeSolid2D(x_tmp, y_tmp);
			break;
		case SolidGhost:
			Node[NodeNumber]=new NodeGhost2D(x_tmp, y_tmp);
			Node[NodeNumber]->Set_NodeType(SolidGhost);
			break;
		case Ghost:
			Node[NodeNumber]=new NodeGhost2D(x_tmp, y_tmp);
			break;
		case Velocity:
			Node[NodeNumber]=new NodeVelocity2D(x_tmp, y_tmp);
			break;
		case Pressure:
			Node[NodeNumber]=new NodePressure2D(x_tmp, y_tmp);
			break;
		case Periodic:
			Node[NodeNumber]=new NodePeriodic2D(x_tmp, y_tmp);
			break;
		case Symmetry:
			Node[NodeNumber]=new NodeSymmetry2D(x_tmp, y_tmp);
			break;
		case Corner:
			Node[NodeNumber]=new NodeCorner2D(x_tmp, y_tmp);
			break;
		case ConcaveCorner:
			Node[NodeNumber]=new NodeCorner2D(x_tmp, y_tmp);
			break;
		case ConvexCorner:
			Node[NodeNumber]=new NodeCorner2D(x_tmp, y_tmp);
			break;
		case Wall:
			Node[NodeNumber]=new NodeWall2D(x_tmp, y_tmp);
			break;
		default:
			cerr<< "Node Type not found" << endl;
	        break;
		}
		Node[NodeNumber]->Set_Connect(Connect_N_tmp,Connect_S_tmp,Connect_W_tmp,Connect_E_tmp);
		Node[NodeNumber]->Set_NbVelocity(NbVelocity_tmp);
		Node[NodeNumber]->Set_Index(index);
}*/
void NodeArrays2D::Get_CoordinateNextNodeAtNormal(int index,double & x, double & y){
	if(index<TypeOfNode.size())
	{
		int localidx=NodeIndexByType[index];
		switch(TypeOfNode[index])
		{
		case Interior:
			Get_coordinate(NodeInterior[localidx].Get_connect()[0],x,y);
			break;
		case Solid:
			Get_coordinate(NodeSolid[localidx].Get_connect()[0],x,y);
			break;
		case Ghost:
			Get_coordinate(NodeGhost[localidx].Get_connect()[0],x,y);
			break;
		case Corner:
			Get_coordinate(NodeCorner[localidx].Get_connect()[NodeCorner[localidx].Get_BcNormal()],x,y);
			break;
		case Wall:
			Get_coordinate(NodeWall[localidx].Get_connect()[NodeWall[localidx].Get_BcNormal()],x,y);
			break;
		case Periodic:
			Get_coordinate(NodePeriodic[localidx].Get_connect()[NodePeriodic[localidx].Get_BcNormal()],x,y);
			break;
		case Velocity:
			Get_coordinate(NodeVelocity[localidx].Get_connect()[NodeVelocity[localidx].Get_BcNormal()],x,y);
			break;
		case Symmetry:
			Get_coordinate(NodeSymmetry[localidx].Get_connect()[NodeSymmetry[localidx].Get_BcNormal()],x,y);
			break;
		case Pressure:
			Get_coordinate(NodePressure[localidx].Get_connect()[NodePressure[localidx].Get_BcNormal()],x,y);
			break;
		case ConcaveCorner:
			Get_coordinate(NodeCorner[localidx].Get_connect()[NodeCorner[localidx].Get_BcNormal()],x,y);
			break;
		case ConvexCorner:
			Get_coordinate(NodeCorner[localidx].Get_connect()[NodeCorner[localidx].Get_BcNormal()],x,y);
			break;
		case SolidGhost:
			Get_coordinate(NodeGhost[localidx].Get_connect()[0],x,y);
			break;
		case GlobalCorner:
			Get_coordinate(NodeGlobalCorner[localidx].Get_connect()[NodeGlobalCorner[localidx].Get_BcNormal()],x,y);
			break;
		case SpecialWall:
			Get_coordinate(NodeSpecialWall[localidx].Get_connect()[NodeSpecialWall[localidx].Get_BcNormal()],x,y);
			break;
		}
	}
	else
	{
		std::cerr<<"Node number is higher than the number of nodes."<<std::endl;
		x=-10;
		y=-10;
	}
}
