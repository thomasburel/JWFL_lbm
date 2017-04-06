/*
 * NodeArray.cpp
 *
 *  Created on: 9 Dec 2015
 *      Author: thomas
 */

#include "NodeArrays.h"

NodeArrays::NodeArrays() {
	// TODO Auto-generated constructor stub

}

NodeArrays::~NodeArrays() {
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
