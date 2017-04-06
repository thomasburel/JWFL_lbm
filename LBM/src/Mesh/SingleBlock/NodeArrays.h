/*
 * NodeArray.h
 *
 *  Created on: 9 Dec 2015
 *      Author: thomas
 */

#ifndef MESH_SINGLEBLOCK_NODEARRAYS_H_
#define MESH_SINGLEBLOCK_NODEARRAYS_H_
#include "Node2D.h"
#include <boost/serialization/vector.hpp>
#include <boost/serialization/map.hpp>
class NodeArrays {
public:
	NodeArrays();
	virtual ~NodeArrays();
	//Map variables
    std::map<int,NodeType> TypeOfNode;
    std::map<int,int> NodeIndexByType;
    //Save Index for models
	std::vector<int> Solid1stLayer;
	std::vector<int> Solid1stLayerInCornerArray;
	std::vector<int> CornerConcave;
	std::vector<int> CornerConvex;

    int Get_NodeIndex(int const IndexNodeInDomain){return NodeIndexByType[IndexNodeInDomain];};

};

class NodeArrays2D: public NodeArrays{
public:
	NodeArrays2D();
	virtual ~NodeArrays2D();
	//std::vector<NodeInterior2D>::iterator itBc;
	std::vector<NodeInterior2D> NodeInterior;
	std::vector<NodeSolid2D> NodeSolid;
	std::vector<NodeGhost2D> NodeGhost;
	std::vector<NodeCorner2D> NodeCorner;
	std::vector<NodeCorner2D> NodeGlobalCorner; //Special treatment must be applied for BC of the domain as symmetry-pressure
	std::vector<NodeWall2D> NodeWall;
	std::vector<NodeWall2D> NodeSpecialWall; //Special treatment must be applied for wall on the domain as Wall-Symmetry
	std::vector<NodePeriodic2D> NodePeriodic;
	std::vector<NodeVelocity2D> NodeVelocity;
	std::vector<NodePressure2D> NodePressure;
	std::vector<NodeSymmetry2D> NodeSymmetry;

	void Get_coordinate(int index,double & x, double & y);

};
#endif /* MESH_SINGLEBLOCK_NODEARRAYS_H_ */
