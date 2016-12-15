/*
 * Cell2D.h
 *
 *  Created on: 21 Apr 2015
 *      Author: thomas
 */

#ifndef MESH_SINGLEBLOCK_CELL2D_H_
#define MESH_SINGLEBLOCK_CELL2D_H_
#include "Node2D.h"

#include <iostream>
// include this header to serialize vectors
#include <boost/serialization/vector.hpp>
//using namespace std;
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
/*
 *
 */
class Cell2D {
public:
	Cell2D();
//	Cell2D(unsigned int* NodeNumber_0);
	virtual ~Cell2D();
	Node2D* NewNode(NodeType const NodeType_, unsigned int const x=0, unsigned int const y=0);
	void Set_Face(int FaceNumber, int& node1, int& node2);
	short int* Get_Face(int FaceNumber)const;
	void Set_Connect(int FaceNumber, int face_, int cell_);
	short int* Get_Connect(int FaceNumber)const;
	void Set_NodeNumber(int NodeNumber_[4]);
	void Set_NodeNumber(int NodeNumber_, int IdNode);
	short int Get_NodeNumber(int NodeNumber_) const;
private:
	short int NodeNumber[4]; // Map the node number from the block to the cell
	short int Connect[4][2]; // Map the Face connect to the Face a in the cell b [Facenumber][Face/Cell] Face=0 and Cell=1
	short int* Connect_;
private:
	friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
    	ar & NodeNumber & Connect;
    }
};

#endif /* MESH_SINGLEBLOCK_CELL2D_H_ */
