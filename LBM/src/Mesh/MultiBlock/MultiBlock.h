/*
 * MultiBlock.h
 *
 *  Created on: 5 May 2015
 *      Author: thomas
 */

#ifndef MESH_MULTIBLOCK_MULTIBLOCK_H_
#define MESH_MULTIBLOCK_MULTIBLOCK_H_

/*
 *
 */
#include "../SingleBlock/Block.h"
#include "../../Core/Parameters.h"
#include "../../Parallelism/MpiManager.h"
class MultiBlock {
public:
	MultiBlock();
	virtual ~MultiBlock();
	virtual void Partitioning()=0;
	virtual void Modify_Block()=0;
	virtual Block* Get_Block()=0;
	virtual int Get_Nx()=0;
	virtual int Get_Ny()=0;
	virtual const double* Get_X0()=0;
	virtual const double* Get_Y0()=0;
	virtual int* Get_Elems0()=0;
	virtual int Get_nnodes()=0;
	virtual int Get_Start_Nodes()=0;
	virtual int Get_End_Nodes()=0;
	virtual int Get_Start_Elems()=0;
	virtual int Get_End_Elems()=0;
	virtual void Get_Connect_Node(std::vector<int> & IdNodeN_,std::vector<int> & IdNodeE_,std::vector<int> & IdNodeS_,std::vector<int> & IdNodeW_,
								  std::vector<int> & IdNodeSW_,std::vector<int> & IdNodeSE_,std::vector<int> & IdNodeNW_,std::vector<int> & IdNodeNE_)=0;
	virtual void Get_Connect_Node(std::vector<int> & IdRNodeN_,std::vector<int> & IdRNodeE_,std::vector<int> & IdRNodeS_,std::vector<int> & IdRNodeW_,
			std::vector<int> & IdGNodeN_,std::vector<int> & IdGNodeE_,std::vector<int> & IdGNodeS_,std::vector<int> & IdGNodeW_,
			std::vector<int> & IdRNodeSW_,std::vector<int> & IdRNodeSE_,std::vector<int> & IdRNodeNW_,std::vector<int> & IdRNodeNE_,
			std::vector<int> & IdGNodeSW_,std::vector<int> & IdGNodeSE_,std::vector<int> & IdGNodeNW_,std::vector<int> & IdGNodeNE_)=0;
	virtual void Communication(double **buf_send,double **buf_recv, int *size_buf)=0;
	virtual void CommunicationToGhost(double **buf_send,double **buf_recv, int *size_buf)=0;
	virtual void CommunicationFromGhost(double **buf_send,double **buf_recv, int *size_buf)=0;
	virtual void CommunicationFromGhost(double **buf_send,double **buf_recv, int *size_buf,MPI_Status * status,MPI_Request * request)=0;
	virtual void Send(double *buf_send, int size_buf,int blocDirection,int tag)=0;
	virtual void Recv(double *buf_recv, int size_buf,int blocDirection,int tag,MPI_Status status)=0;
	virtual void Get_NbNodes(int & NbRealNodes_, int & NbTotalNodes_)=0;
	virtual void ConvertToPhysicalUnit()=0;
	virtual void reorganizeNodeByType()=0;
	virtual NodeArrays2D* Get_NodeArrays2D()=0;
	virtual int* get_Block_Connect()=0;
protected:
	int ndims;
	Parameters * PtrParameters;
	ParallelManager* parallel;

private:
	friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
       ar & ndims ;
    }
};

#endif /* MESH_MULTIBLOCK_MULTIBLOCK_H_ */
