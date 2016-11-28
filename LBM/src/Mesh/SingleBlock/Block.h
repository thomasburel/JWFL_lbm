/*
 * Block.h
 *
 *  Created on: 5 May 2015
 *      Author: thomas
 */

#ifndef MESH_SINGLEBLOCK_BLOCK_H_
#define MESH_SINGLEBLOCK_BLOCK_H_
/*
 *
 */
#include "NodeArrays.h"
#include "Cell2D.h"
#include "../../User/UserMesh.h"
#include "../../Core/Parameters.h"
// include this header to serialize vectors
#include <boost/serialization/vector.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
class Block : public UserMesh   {
	private:
	friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
    }
public:
	Block();
	virtual ~Block();
	//virtual std::vector<Node2D*>* Get_PtrNode()=0;
	virtual void Get_NbNodes(int & NbRealNodes_, int & NbTotalNodes_)=0;
	virtual std::vector<int>& Get_PtrIdBc()=0;
	virtual Node2D* Get_Node(int NodeNumber)=0;
	virtual void Get_Connect_Node(std::vector<int> & IdNodeN_,std::vector<int> & IdNodeE_,std::vector<int> & IdNodeS_,std::vector<int> & IdNodeW_,
								  std::vector<int> & IdNodeSW_,std::vector<int> & IdNodeSE_,std::vector<int> & IdNodeNW_,std::vector<int> & IdNodeNE_)=0;
	//virtual void ModifyMeshByUser(Parameters &Param)=0;
	virtual void reorganizeNodeByType()=0;
	virtual void ConvertToPhysicalUnit(Parameters &Param)=0;
	virtual NodeArrays2D* Get_NodeArrays2D()=0;
	virtual void Set_Connect(Parameters& Param)=0;
	virtual void Mark1stLayerSolid()=0;
};

#endif /* MESH_SINGLEBLOCK_BLOCK_H_ */
