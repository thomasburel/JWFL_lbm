/*
 * Cell2D.cpp
 *
 *  Created on: 21 Apr 2015
 *      Author: thomas
 */

#include "Cell2D.h"

Cell2D::Cell2D()
{
	for (int i=0;i<4;i++)
	{
		NodeNumber[i]=i;
		Connect[i][0]=i;
		Connect[i][1]=0;
	}
	NodeNumber[3]=0;
	Connect_=new unsigned int[2];
}

Cell2D::~Cell2D() {
	delete [] Connect;
}

Node2D* Cell2D::NewNode(NodeType NodeType_, unsigned int x, unsigned int y)
{
	switch(NodeType_)
	{
	case Interior:
		return new NodeInterior2D(x, y);
		break;
	case Solid:
		return new NodeSolid2D(x, y);
		break;
	case Ghost:
		return new NodeGhost2D(x, y);
		break;
	case Velocity:
		return new NodeVelocity2D(x, y);
		break;
	case Pressure:
		return new NodePressure2D(x, y);
		break;
	case Periodic:
		return new NodePeriodic2D(x, y);
		break;
	case Corner:
		return new NodeCorner2D(x, y);
		break;
	case GlobalCorner:
		return new NodeCorner2D(x, y);
		break;
	case Wall:
		return new NodeWall2D(x, y);
		break;
	case Symmetry:
		return new NodeSymmetry2D(x, y);
		break;
	default:
		std::cerr<< "Node Type not found" << std::endl;
		return new NodeInterior2D(x, y);
        break;
	}

}

void Cell2D::Set_Face(int FaceNumber, int& node1, int& node2)
{
	switch(FaceNumber)
	{
	case 0:
		NodeNumber[0]=node1;
		NodeNumber[1]=node2;
		break;
	case 1:
		NodeNumber[1]=node1;
		NodeNumber[2]=node2;
		break;
	case 2:
		NodeNumber[2]=node1;
		NodeNumber[3]=node2;
		break;
	case 3:
		NodeNumber[3]=node1;
		NodeNumber[0]=node2;
		break;
	default:
		std::cerr<< "Wrong Face Number" << std::endl;
	        break;
	}
}
unsigned int* Cell2D::Get_Face(int FaceNumber)const
{
	unsigned int* Face_=new unsigned int[2];
	switch(FaceNumber)
	{
	case 0:
		Face_[0]=NodeNumber[0];
		Face_[1]=NodeNumber[1];
		break;
	case 1:
		Face_[0]=NodeNumber[1];
		Face_[1]=NodeNumber[2];
		break;
	case 2:
		Face_[0]=NodeNumber[2];
		Face_[1]=NodeNumber[3];
		break;
	case 3:
		Face_[0]=NodeNumber[3];
		Face_[1]=NodeNumber[0];
		break;
	default:
		std::cerr<< "Wrong Face Number" << std::endl;
		break;
	}
	return Face_;
}
void Cell2D::Set_Connect(int FaceNumber, int face_, int cell_)
{
	Connect[FaceNumber][0]=face_;
	Connect[FaceNumber][1]=cell_;
}
unsigned int* Cell2D::Get_Connect(int FaceNumber)const
{
	Connect_[0]=Connect[FaceNumber][0];
	Connect_[1]=Connect[FaceNumber][1];
	return Connect_;
}
unsigned int Cell2D::Get_NodeNumber(int NodeNumber_) const
{
	return NodeNumber[NodeNumber_];
}

void Cell2D::Set_NodeNumber(int NodeNumber_[4])
{
	NodeNumber[0]=NodeNumber_[0];
	NodeNumber[1]=NodeNumber_[1];
	NodeNumber[2]=NodeNumber_[2];
	NodeNumber[3]=NodeNumber_[3];
}
void Cell2D::Set_NodeNumber(int NodeNumber_, int IdNode)
{
	NodeNumber[NodeNumber_]=IdNode;
}
