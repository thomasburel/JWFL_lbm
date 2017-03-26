/*
 * Block2D.h
 *
 *  Created on: 17 Apr 2015
 *      Author: thomas
 */

#ifndef MESH_SINGLEBLOCK_BLOCK2D_H_
#define MESH_SINGLEBLOCK_BLOCK2D_H_

// include this header to serialize vectors
#include <boost/serialization/vector.hpp>
// include this header to serialize maps
#include <boost/serialization/map.hpp>
#include "Cell2D.h"
#include "Block.h"
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>// std::sort
#include "mpi.h"
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
using namespace std;

class Block2D : public Block {
public:
	Block2D(int dx_=1, int dy_=1);
	virtual ~Block2D();

	void AddBlock(int TotalNumberCells_x,int TotalNumberCells_y,int dims[2],int periodic_[2],int coord[2],int NyCell_G,int Nx_begin,int Ny_begin,int Nx_last,int Ny_last,NodeType bC[4],int x_last=0,int y_last=0,int start_GlobalNodes=1,int end_GlobalNodes=1, int start_GlobalElems=1, int end_GlobalElems=1);
//	void AddGhostCell();
	void WriteCoord();
	double* Get_Coord(int NodeNumber) const;
	const double* Get_X0() const;
	const double* Get_Y0() const;
	int Get_nnodes()const ;
	int* Get_Elems0();
	//virtual std::vector<Node2D*>* Get_PtrNode();
	virtual void Get_NbNodes(int & NbRealNodes_, int & NbTotalNodes_);
	virtual std::vector<int>& Get_PtrIdBc();
	virtual Node2D* Get_Node(int NodeNumber);
	virtual void Get_Connect_Node(std::vector<int> & IdNodeN_,std::vector<int> & IdNodeE_,std::vector<int> & IdNodeS_,std::vector<int> & IdNodeW_,
								  std::vector<int> & IdNodeSW_,std::vector<int> & IdNodeSE_,std::vector<int> & IdNodeNW_,std::vector<int> & IdNodeNE_);

	//virtual void ModifyMeshByUser(Parameters &Param);
	virtual void InitPatchBc(Parameters *PtrParam);
	virtual void Set_NodeIndexByTypeForPatchBc();
	virtual void reorganizeNodeByType();
	virtual void ConvertToPhysicalUnit(Parameters &Param);
	virtual NodeArrays2D* Get_NodeArrays2D();

	virtual void Set_Connect(Parameters& Param);
	virtual void Mark1stLayerSolid();

	void Correct_Solid_Ghost();
	void Get_GhostType(std::vector<int> & NodeTypeN,std::vector<int> & NodeTypeE,std::vector<int> & NodeTypeS,std::vector<int> & NodeTypeW,
			  std::vector<int> & NodeTypeSW,std::vector<int> & NodeTypeSE,std::vector<int> & NodeTypeNW,std::vector<int> & NodeTypeNE);
	void Set_GhostType(std::vector<int> & NodeTypeN,std::vector<int> & NodeTypeE,std::vector<int> & NodeTypeS,std::vector<int> & NodeTypeW,
			  std::vector<int> & NodeTypeSW,std::vector<int> & NodeTypeSE,std::vector<int> & NodeTypeNW,std::vector<int> & NodeTypeNE);
	void Remove_SolidTypeInCommunicators(std::vector<int> & NodeTypeN,std::vector<int> & NodeTypeE,std::vector<int> & NodeTypeS,std::vector<int> & NodeTypeW,
			  std::vector<int> & NodeTypeSW,std::vector<int> & NodeTypeSE,std::vector<int> & NodeTypeNW,std::vector<int> & NodeTypeNE);
	void Remove_SolidTypeInCommunicator(std::vector<int> & RealNodeType,std::vector<int> & GhostNodeType,
			  std::vector<int> & RealNodeId,std::vector<int> & GhostNodeId,std::vector<int> & RealIdToBeSaved,std::vector<int> & GhostIdToBeSaved);
	void Get_CommNodes(std::vector<int> & IdRNodeN,std::vector<int> & IdRNodeE,std::vector<int> & IdRNodeS,std::vector<int> & IdRNodeW,
			std::vector<int> & IdGNodeN,std::vector<int> & IdGNodeE,std::vector<int> & IdGNodeS,std::vector<int> & IdGNodeW,
			std::vector<int> & IdRNodeSW,std::vector<int> & IdRNodeSE,std::vector<int> & IdRNodeNW,std::vector<int> & IdRNodeNE,
			std::vector<int> & IdGNodeSW,std::vector<int> & IdGNodeSE,std::vector<int> & IdGNodeNW,std::vector<int> & IdGNodeNE);

public:
	void GenerateSolid(Parameters &Param);
	void SetSolidBoundaries();
	void RemoveUnphysicalSolid(int &nbTotalSolidRemoved,int &nbTotalSolidadded);
	void RemoveSolidInCommunicator();
	void RemoveSolid();
private:
	void NewCell(NodeType Nodes[4], bool OldNodes[4],bool GhostNodes[4], int Node2D_SubDomain[4],int Node2D_GlobalDomain[4],int x_[4],int y_[4], int FaceConnect[4], int CellConnect[4]);
	void NewGhostCell(NodeType Nodes[4], bool OldNodes[4], int Node2D_SubDomain[4],int x_[4],int y_[4], int FaceConnect[4], int CellConnect[4]);
	void ChangeNodeType(int NodeNumber, NodeType NewNodeType_);
	void ChangeCoord(int NodeNumber,unsigned int const x=0, unsigned int const y=0);
	NodeType Get_NodeType(int NodeNumber) const;
//	Node2D* Get_Node(int NodeNumber) ;
	void Correct_OrderingGhostNode();
	void Set_Connect();
	void Check_ID();
	void Correct_MarkNode();
	void Clear_MarkNode();
	void DefinedCornerType(int nodenumber);
	int Connect_lowOrder(int &NodeNumber,unsigned int& direction);
	void Connect_highOrder();



	bool DetectSolidBoundaries(int & nodeID);
	void DetectUnphysicalSolid(int & nodeID);
	void DetectDirectInterior(int & nodeID, int & nbinterior);
	void DetectSpecialWall(int & nodeID,int x,int y, bool & specialwall, unsigned int & directionwall);
	void CreateCorner(int &nodeID);
	void CreateCornerConcave(int &nodeID);
	void CreateCornerConvex(int &nodeID);
	void CreateWall(int &nodeID);
	void CreateSpecialWall(int &nodeID);
	void CreateWallandCorners(Parameters &Param, int &nodeID, int & nbinterior, bool SpecialNode);
	void AddGhostCells(int TotalNumberCells_x,int TotalNumberCells_y);
	void WriteCells(); //For debugging
	void WriteNodes(); //For debugging
	void Set_BcNormal();
	int  Get_BcNormal(NodeWall2D& Node);
	int  Get_BcNormal_SpecialWall(NodeWall2D& Node);
	int  Get_BcNormal(NodeCorner2D& Node);
	int  Get_BcNormal(NodeVelocity2D& Node);
	int  Get_BcNormal(NodePressure2D& Node);
	int  Get_BcNormal(NodeSymmetry2D& Node);
	int  Get_BcNormal(NodePeriodic2D& Node);
	void Remove_OneCell(int cellnumber);
	void Remove_ExtraGhostCells();
	void Remove_SolidCells();
	void Correct_GhostType(int  idNode, NodeType RealNodeType);
	void Set_CommNodes();

/// Private Variables
private:
	std::vector<Cell2D*> CellArray,GhostCellArraytmp; //Store cell in the Block
	std::vector<Node2D*> Node, Node_Ghosttmp; // Store Nodes in Vector with ghost node at the end but before Solid nodes
	NodeArrays2D NodeArrays;
	std::vector<int> IdBoundaries;
	std::vector<int>::iterator position; //find position to remove IDBoundaries
	std::vector<double> x,y,z,x_Ghosttmp,y_Ghosttmp;// Store Node coordinates in Vector with ghost node at the end
	std::vector<int> IdCellGhostNode; //Mark cell with ghost node to correct the node numbering before connectivity.
	std::vector<int> IdGhostCell; //Mark Ghost cell for parallel consideration and also for periodic boundary conditions
	std::vector<int> NdGhostNodeInCell,IdGhostNodeInCell; //
	int NbGhostNode;
	int NbRealNodes, NbTotalNodes;
	std::vector<int> Elems; //4 nodes to define each cell
    int dx,dy,nx,ny;
    std::map<int,int> LocalToGlobalNode;
	std::vector<int> IdNodeN,IdNodeE,IdNodeS,IdNodeW;
	std::vector<int> IdNodeSW,IdNodeSE,IdNodeNW,IdNodeNE;
//Mark real and ghost nodes for communication
	std::vector<int> IdRNodeN,IdRNodeE,IdRNodeS,IdRNodeW,IdGNodeN,IdGNodeE,IdGNodeS,IdGNodeW;
	std::vector<int> IdRNodeSW,IdRNodeSE,IdRNodeNW,IdRNodeNE,IdGNodeSW,IdGNodeSE,IdGNodeNW,IdGNodeNE;

	std::vector<int> IdSolidNode,Id_SolidGhost,IdSolidBc;//Mark solid node to sort the array
	std::vector<int> IdWalltmp,IdSpecialWalltmp,IdCornerConcavetmp,IdCornerConvextmp,IDRemoveSolidtmp,IDAddSolidtmp,IDRemoveSolidToBctmp;
	std::vector<int> IDRemoveSolidGhosttmp,IDAddSolidGhosttmp;
	std::vector<Node2D*> Node_Solidtmp,Node_tmp;
	int Coord[2],Dims[2];
	bool periodic[2];
	bool verbous;
	double* Coord_physical;
	int intTmpReturn;
	NodeType bC[4];

;

private:
	friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
    	ar & boost::serialization::base_object<Block>(*this);
    	ar & IdBoundaries & x & y & x_Ghosttmp & y_Ghosttmp & IdCellGhostNode & NdGhostNodeInCell & IdGhostNodeInCell;
    	ar & NbGhostNode & NbRealNodes & NbTotalNodes & dx & dy;
    	ar & LocalToGlobalNode;
    	ar & Elems & IdNodeN & IdNodeE & IdNodeS & IdNodeW & IdNodeSW & IdNodeSE & IdNodeNW & IdNodeNE;
    	ar & Coord & periodic & Dims & verbous;
    	ar & CellArray;
    	//ar & Node;
		ar & Node_Ghosttmp;
    }
};

#endif /* MESH_SINGLEBLOCK_BLOCK2D_H_ */
