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
 //   virtual void ChangeNodeType(int NodeNumberbyType,NodeType OldNodeType_, NodeType NewNodeType_)=0;

    virtual int Get_SizeNodeIdInterior()=0;
	virtual int& Get_NodeIdInterior(int idx)=0;
    virtual int Get_SizeNodeIdSolid()=0;
	virtual int& Get_NodeIdSolid(int idx)=0;
    virtual int Get_SizeNodeIdGhost()=0;
	virtual int& Get_NodeIdGhost(int idx)=0;
    virtual int Get_SizeNodeIdCorner()=0;
	virtual int& Get_NodeIdCorner(int idx)=0;
    virtual int Get_SizeNodeIdCornerConcave()=0;
	virtual int& Get_NodeIdCornerConcave(int idx)=0;
    virtual int Get_SizeNodeIdCornerConvex()=0;
	virtual int& Get_NodeIdCornerConvex(int idx)=0;

    virtual int Get_SizeNodeIdGlobalCorner()=0;
	virtual int& Get_NodeIdGlobalCorner(int idx)=0;
    virtual int Get_SizeNodeIdWall()=0;
	virtual int& Get_NodeIdWall(int idx)=0;
    virtual int Get_SizeNodeIdSpecialWall()=0;
	virtual int& Get_NodeIdSpecialWall(int idx)=0;
    virtual int Get_SizeNodeIdPeriodic()=0;
	virtual int& Get_NodeIdPeriodic(int idx)=0;
    virtual int Get_SizeNodeIdVelocity()=0;
	virtual int& Get_NodeIdVelocity(int idx)=0;
    virtual int Get_SizeNodeIdPressure()=0;
	virtual int& Get_NodeIdPressure(int idx)=0;
    virtual int Get_SizeNodeIdSymmetry()=0;
	virtual int& Get_NodeIdSymmetry(int idx)=0;
// get normal
	virtual int& Get_NodeNormalInterior(int idx)=0;
	virtual int& Get_NodeNormalSolid(int idx)=0;
	virtual int& Get_NodeNormalGhost(int idx)=0;
	virtual int& Get_NodeNormalCorner(int idx)=0;
	virtual int& Get_NodeNormalCornerConcave(int idx)=0;
	virtual int& Get_NodeNormalCornerConvex(int idx)=0;
	virtual int& Get_NodeNormalGlobalCorner(int idx)=0;
	virtual int& Get_NodeNormalWall(int idx)=0;
	virtual int& Get_NodeNormalSpecialWall(int idx)=0;
	virtual int& Get_NodeNormalPeriodic(int idx)=0;
	virtual int& Get_NodeNormalVelocity(int idx)=0;
	virtual int& Get_NodeNormalPressure(int idx)=0;
	virtual int& Get_NodeNormalSymmetry(int idx)=0;
//get connection
	virtual int* Get_NodeConnectInterior(int idx)=0;
	virtual int* Get_NodeConnectSolid(int idx)=0;
	virtual int* Get_NodeConnectGhost(int idx)=0;
	virtual int* Get_NodeConnectCorner(int idx)=0;
	virtual int* Get_NodeConnectCornerConcave(int idx)=0;
	virtual int* Get_NodeConnectCornerConvex(int idx)=0;
	virtual int* Get_NodeConnectGlobalCorner(int idx)=0;
	virtual int* Get_NodeConnectWall(int idx)=0;
	virtual int* Get_NodeConnectSpecialWall(int idx)=0;
	virtual int* Get_NodeConnectPeriodic(int idx)=0;
	virtual int* Get_NodeConnectVelocity(int idx)=0;
	virtual int* Get_NodeConnectPressure(int idx)=0;
	virtual int* Get_NodeConnectSymmetry(int idx)=0;

	virtual void Get_coordinate(int index,double & x, double & y)=0;
	virtual void Get_CoordinateNextNodeAtNormal(int index,double & x, double & y)=0;
protected:
	int nullNodeArray;
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

	//virtual void ChangeNodeType(int NodeNumberbyType,NodeType OldNodeType_, NodeType NewNodeType_);

	virtual void Get_coordinate(int index,double & x, double & y);
	virtual void Get_CoordinateNextNodeAtNormal(int index,double & x, double & y);

	virtual int Get_SizeNodeIdInterior(){return NodeInterior.size();};
	virtual int& Get_NodeIdInterior(int idx){return NodeInterior[idx].Get_index();};
    virtual int Get_SizeNodeIdSolid(){return NodeSolid.size();};
	virtual int& Get_NodeIdSolid(int idx){return NodeSolid[idx].Get_index();};
    virtual int Get_SizeNodeIdGhost(){return NodeGhost.size();};
	virtual int& Get_NodeIdGhost(int idx){return NodeGhost[idx].Get_index();};
    virtual int Get_SizeNodeIdCorner(){return NodeCorner.size();};
	virtual int& Get_NodeIdCorner(int idx){return NodeCorner[idx].Get_index();};
    virtual int Get_SizeNodeIdCornerConcave(){return CornerConcave.size();};
	virtual int& Get_NodeIdCornerConcave(int idx){return NodeCorner[CornerConcave[idx]].Get_index();};
    virtual int Get_SizeNodeIdCornerConvex(){return CornerConvex.size();};
	virtual int& Get_NodeIdCornerConvex(int idx){return NodeCorner[CornerConvex[idx]].Get_index();};
    virtual int Get_SizeNodeIdGlobalCorner(){return NodeGlobalCorner.size();};
	virtual int& Get_NodeIdGlobalCorner(int idx){return NodeGlobalCorner[idx].Get_index();};
    virtual int Get_SizeNodeIdWall(){return NodeWall.size();};
	virtual int& Get_NodeIdWall(int idx){return NodeWall[idx].Get_index();};
    virtual int Get_SizeNodeIdSpecialWall(){return NodeSpecialWall.size();};
	virtual int& Get_NodeIdSpecialWall(int idx){return NodeSpecialWall[idx].Get_index();};
    virtual int Get_SizeNodeIdPeriodic(){return NodePeriodic.size();};
	virtual int& Get_NodeIdPeriodic(int idx){return NodePeriodic[idx].Get_index();};
    virtual int Get_SizeNodeIdVelocity(){return NodeVelocity.size();};
	virtual int& Get_NodeIdVelocity(int idx){return NodeVelocity[idx].Get_index();};
    virtual int Get_SizeNodeIdPressure(){return NodePressure.size();};
	virtual int& Get_NodeIdPressure(int idx){return NodePressure[idx].Get_index();};
    virtual int Get_SizeNodeIdSymmetry(){return NodeSymmetry.size();};
	virtual int& Get_NodeIdSymmetry(int idx){return NodeSymmetry[idx].Get_index();};

	// get normal
		virtual int& Get_NodeNormalInterior(int idx){return nullNodeArray;};
		virtual int& Get_NodeNormalSolid(int idx){return nullNodeArray;};
		virtual int& Get_NodeNormalGhost(int idx){return nullNodeArray;};
		virtual int& Get_NodeNormalCorner(int idx){return NodeCorner[idx].Get_BcNormal();};
		virtual int& Get_NodeNormalCornerConcave(int idx){return NodeCorner[CornerConcave[idx]].Get_BcNormal();};
		virtual int& Get_NodeNormalCornerConvex(int idx){return NodeCorner[CornerConvex[idx]].Get_BcNormal();};
		virtual int& Get_NodeNormalGlobalCorner(int idx){return NodeGlobalCorner[idx].Get_BcNormal();};
		virtual int& Get_NodeNormalWall(int idx){return NodeWall[idx].Get_BcNormal();};
		virtual int& Get_NodeNormalSpecialWall(int idx){return NodeSpecialWall[idx].Get_BcNormal();};
		virtual int& Get_NodeNormalPeriodic(int idx){return NodePeriodic[idx].Get_BcNormal();};
		virtual int& Get_NodeNormalVelocity(int idx){return NodeVelocity[idx].Get_BcNormal();};
		virtual int& Get_NodeNormalPressure(int idx){return NodePressure[idx].Get_BcNormal();};
		virtual int& Get_NodeNormalSymmetry(int idx){return NodeSymmetry[idx].Get_BcNormal();};
	//get connection
		virtual int* Get_NodeConnectInterior(int idx){return NodeInterior[idx].Get_connect();};
		virtual int* Get_NodeConnectSolid(int idx){return NodeSolid[idx].Get_connect();};
		virtual int* Get_NodeConnectGhost(int idx){return NodeGhost[idx].Get_connect();};
		virtual int* Get_NodeConnectCorner(int idx){return NodeCorner[idx].Get_connect();};
		virtual int* Get_NodeConnectCornerConcave(int idx){return NodeCorner[CornerConcave[idx]].Get_connect();};
		virtual int* Get_NodeConnectCornerConvex(int idx){return NodeCorner[CornerConvex[idx]].Get_connect();};
		virtual int* Get_NodeConnectGlobalCorner(int idx){return NodeGlobalCorner[idx].Get_connect();};
		virtual int* Get_NodeConnectWall(int idx){return NodeWall[idx].Get_connect();};
		virtual int* Get_NodeConnectSpecialWall(int idx){return NodeSpecialWall[idx].Get_connect();};
		virtual int* Get_NodeConnectPeriodic(int idx){return NodePeriodic[idx].Get_connect();};
		virtual int* Get_NodeConnectVelocity(int idx){return NodeVelocity[idx].Get_connect();};
		virtual int* Get_NodeConnectPressure(int idx){return NodePressure[idx].Get_connect();};
		virtual int* Get_NodeConnectSymmetry(int idx){return NodeSymmetry[idx].Get_connect();};
};
#endif /* MESH_SINGLEBLOCK_NODEARRAYS_H_ */
