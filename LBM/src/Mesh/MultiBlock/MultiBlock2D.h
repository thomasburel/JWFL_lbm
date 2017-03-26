/*
 * MultiBlock2D.h
 *
 *  Created on: 20 Apr 2015
 *      Author: thomas
 */

#ifndef MESH_MULTIBLOCK_MULTIBLOCK2D_H_
#define MESH_MULTIBLOCK_MULTIBLOCK2D_H_
#include"MultiBlock.h"
#include "../SingleBlock/Block2D.h"

#include <numeric>      // std::accumulate
/*
 *
 */
class MultiBlock2D: public MultiBlock {
public:
	MultiBlock2D();
	MultiBlock2D(ParallelManager* parrallel_,Parameters * PtrParameters_);
	virtual ~MultiBlock2D();
	virtual void Partitioning();
	void Create_Block2D();
	virtual void Modify_Block();
	virtual void GeneratePatchBc();
	virtual Block2D* Get_Block();
	virtual int Get_Nx();
	virtual int Get_Ny();
	virtual const double* Get_X0();
	virtual const double* Get_Y0();
	virtual int Get_nnodes();
	virtual int* Get_Elems0();
	virtual int Get_Start_Nodes();
	virtual int Get_End_Nodes();
	virtual int Get_Start_Elems();
	virtual int Get_End_Elems();
	virtual void Get_Connect_Node(std::vector<int> & IdNodeN_,std::vector<int> & IdNodeE_,std::vector<int> & IdNodeS_,std::vector<int> & IdNodeW_,
								  std::vector<int> & IdNodeSW_,std::vector<int> & IdNodeSE_,std::vector<int> & IdNodeNW_,std::vector<int> & IdNodeNE_);
	virtual void Get_Connect_Node(std::vector<int> & IdRNodeN_,std::vector<int> & IdRNodeE_,std::vector<int> & IdRNodeS_,std::vector<int> & IdRNodeW_,
			std::vector<int> & IdGNodeN_,std::vector<int> & IdGNodeE_,std::vector<int> & IdGNodeS_,std::vector<int> & IdGNodeW_,
			std::vector<int> & IdRNodeSW_,std::vector<int> & IdRNodeSE_,std::vector<int> & IdRNodeNW_,std::vector<int> & IdRNodeNE_,
			std::vector<int> & IdGNodeSW_,std::vector<int> & IdGNodeSE_,std::vector<int> & IdGNodeNW_,std::vector<int> & IdGNodeNE_);
	virtual void Communication(double **buf_send,double **buf_recv, int *size_buf);
	virtual void CommunicationToGhost(double **buf_send,double **buf_recv, int *size_buf);
	virtual void CommunicationFromGhost(double **buf_send,double **buf_recv, int *size_buf);
	virtual void CommunicationFromGhost(double **buf_send,double **buf_recv, int *size_buf,MPI_Status * status,MPI_Request * request);
	virtual void Send(double *buf_send, int size_buf,int blocDirection,int tag);
	virtual void Recv(double *buf_recv, int size_buf,int blocDirection,int tag,MPI_Status status);
	virtual void Get_NbNodes(int & NbRealNodes_, int & NbTotalNodes_);
	virtual void ConvertToPhysicalUnit();
	virtual void reorganizeNodeByType();
	virtual NodeArrays2D* Get_NodeArrays2D();
	void Correct_SolidGhost();
	void Remove_SolidInComunicators();
	virtual int* get_Block_Connect(){return &BlockNeighbour[0];};
	virtual double SumBC(double *value);
	virtual double SumAllProcessors(double *value);
	virtual int NumberOfProcessors();

private:
	int Nx_G,Ny_G;
	int  Nx_size, Ny_size;//For serialization of Nx and Ny
	int *Nx;
	int *Ny;
	int start_nodes, end_nodes,start_elems, end_elems;
	Block2D Block2D_;
	MPI_Comm COMM_CART;
	MPI_Comm VertComm,HorizComm;
	int periods[2];
	int reorder;
	int coord[2];
	int dims[2];
	int BlockNeighbour[8];
	int    N,E,S,W, NE, SE, NW, SW;
	int rank_in_topo;
	bool verbous;
	bool BcW,BcE,BcN, BcS;//mark the block located to a border to calculate properties

private:
	friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
    	ar & boost::serialization::base_object<MultiBlock>(*this);
    	ar & Nx_G & Ny_G & Nx_size & Ny_size;
    	for (int i=0;i<Nx_size;i++)
    		ar & Nx[i];
    	for (int i=0;i<Ny_size;i++)
    		ar & Ny[i];
    	ar & start_nodes & end_nodes & start_elems & end_elems;
    	ar & Block2D_;
    	//ar & COMM_CART;
    	ar & periods & reorder & coord & dims & BlockNeighbour;
    	ar & N & E & S & W & NE & SE & NW & SW;
    	ar & rank_in_topo;
   //    ar & x & y & z;


/*       ar & outfile;
       ar & NbVariableOutput & NbVariableBreakpoint;
       ar & VariableOutput & VariableBreakpoint;
       ar & Dimension;
       ar & outputfilename & Breakpointfilename;*/
    }
};

#endif /* MESH_MULTIBLOCK_MULTIBLOCK2D_H_ */
