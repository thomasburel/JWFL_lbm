/*
 * MultiBlock2D.cpp
 *
 *  Created on: 20 Apr 2015
 *      Author: thomas
 */

#include "MultiBlock2D.h"

MultiBlock2D::MultiBlock2D() :
	Nx_G(1),Ny_G(1),
	start_nodes(1),end_nodes(0),start_elems(1),end_elems(0),
	reorder(false),COMM_CART(0),rank_in_topo(0),VertComm(0),HorizComm(0)
{
	for (int i=0;i<2;i++)
		{
			dims[i]=0;
			coord[i]=0;
			periods[i]=false;
		}

	ndims=2;
	dims[0]=0;
	dims[1]=0;
	Nx=new int [1];
	Ny=new int [1];
	Nx[0]=1;
	Ny[0]=1;
	for (int i=0;i<4;i++)
		BlockNeighbour[i]=MPI_PROC_NULL;
	N=0,E=1,S=2,W=3,NE=4, SE=5, NW=6, SW=7;
	Nx_size=0;
	Ny_size=0;
	verbous=false;
	BcW=false;
	BcS=false;
	BcE=false;
	BcN=false;
}

//MultiBlock2D::MultiBlock2D(ParrallelManager* parallel_, int Nx_, int Ny_) :
MultiBlock2D::MultiBlock2D(ParallelManager* parallel_, Parameters * PtrParameters_) :
		start_nodes(1),end_nodes(0),start_elems(1),end_elems(0),
		reorder(false),COMM_CART(0),VertComm(0),HorizComm(0)
{
	PtrParameters=PtrParameters_;
	for (int i=0;i<2;i++)
	{
		dims[i]=0;
		coord[i]=0;
		periods[i]=false;
	}
	if(PtrParameters->Get_GlobalBcType(0)==Periodic||PtrParameters->Get_GlobalBcType(2)==Periodic)
		periods[1]=true;
	else
		periods[1]=false;
	if(PtrParameters->Get_GlobalBcType(1)==Periodic||PtrParameters->Get_GlobalBcType(3)==Periodic)
		periods[0]=true;
	else
		periods[0]=false;

	parallel=parallel_;
	Nx_G=PtrParameters->Get_Nx();
	Ny_G=PtrParameters->Get_Ny();
	Nx_size=parallel->getSize();
	Ny_size=parallel->getSize();
	Nx=new int[Nx_size];
	Ny=new int[Ny_size];
	for (int i=0;i<parallel->getSize();i++)
	{
		Nx[i]=1;
		Ny[i]=1;
	}
	ndims=2;
	dims[0]=0;
	dims[1]=0;
	rank_in_topo=parallel->getRank();
	for (int i=0;i<4;i++)
		BlockNeighbour[i]=MPI_PROC_NULL;
	N=0,E=1,S=2,W=3,NE=4, SE=5, NW=6, SW=7;
	verbous=PtrParameters->Get_Verbous();
	BcW=false;
	BcS=false;
	BcE=false;
	BcN=false;
}

MultiBlock2D::~MultiBlock2D() {
	delete parallel;
	delete [] Nx,Ny;
}

void MultiBlock2D::Partitioning() {
	if(verbous)
		std::cout<<"My Rank is: "<<parallel->getRank()<< std::endl;
	int coordtmp[2];

	MPI_Dims_create(parallel->getSize(),ndims,dims);
	if(verbous)
		std::cout<<"dims is: "<< dims[0] << " " << dims[1] << std::endl;
	MPI_Cart_create(parallel->getGlobalCommunicator(),ndims,dims, periods, reorder, &COMM_CART);
	/* Create communicator Vertical and Horizontal blocks */
	int remain[2];
	remain[0]=0;remain[1]=1;
	MPI_Cart_sub(COMM_CART,remain,&VertComm);
	remain[0]=1;remain[1]=0;
	MPI_Cart_sub(COMM_CART,remain,&HorizComm);

	/* Get Block Neighbour at West and East */
	MPI_Cart_shift(COMM_CART,0,1,&BlockNeighbour[W],&BlockNeighbour[E]);

	  /* Get Block Neighbour at South and North */
	MPI_Cart_shift(COMM_CART,1,1,&BlockNeighbour[S],&BlockNeighbour[N]) ;

	/* Get coordinates of each process inside the grid */
	MPI_Comm_rank(COMM_CART,&rank_in_topo);

	MPI_Cart_coords(COMM_CART,parallel->getRank(),2,coord);
	/* Mark block on the limit of the domain   */

	if(coord[0]==0)
		BcW=true;
	if(coord[1]==0)
		BcS=true;
	if(coord[0]==dims[0]-1)
		BcE=true;
	if(coord[1]==dims[1]-1)
		BcN=true;
	/* Get Block Neighbour at corners*/

	coordtmp[0]=coord[0]-1;
	coordtmp[1]=coord[1]-1;
	if (coordtmp[0]<0 || coordtmp[0]>dims[0]-1||coordtmp[1]<0 || coordtmp[1]>dims[1]-1)
		if(((coord[0]>0 && periods[1]) || (coord[1]>0 && periods[0])) || (periods[0] && periods[1]))
		{
			if(periods[0] && periods[1] && coord[0]==0 && coord[1]==0)
			{
				coordtmp[0]=dims[0]-1;
				coordtmp[1]=dims[1]-1;
			}
		/*	else
			{
				coordtmp[0]=dims[0]-1;
			}*/
				MPI_Cart_rank (COMM_CART,coordtmp,&BlockNeighbour[SW]);
		}
		else
		{
			BlockNeighbour[SW]=-2;
		}
	else
		MPI_Cart_rank (COMM_CART,coordtmp,&BlockNeighbour[SW]);

	coordtmp[0]=coord[0]-1;
	coordtmp[1]=coord[1]+1;
	if (coordtmp[0]<0 || coordtmp[0]>dims[0]-1||coordtmp[1]<0 || coordtmp[1]>dims[1]-1)
		if(((coord[0]>0 && periods[1]) || (coord[1]<dims[1]-1 && periods[0])) || (periods[0] && periods[1]))
		{
			if(periods[0] && periods[1]&& coord[0]==0 && coord[1]==dims[1]-1)
			{
				coordtmp[0]=dims[0]-1;
				coordtmp[1]=0;
			}
	/*		else
			{
				coordtmp[0]=dims[0]-1;
			}*/
			MPI_Cart_rank (COMM_CART,coordtmp,&BlockNeighbour[NW]);
		}
		else
		{
			BlockNeighbour[NW]=-2;
		}
	else
		MPI_Cart_rank (COMM_CART,coordtmp,&BlockNeighbour[NW]);

	coordtmp[0]=coord[0]+1;
	coordtmp[1]=coord[1]-1;
	if (coordtmp[0]<0 || coordtmp[0]>dims[0]-1||coordtmp[1]<0 || coordtmp[1]>dims[1]-1)
		if(((coord[0]<dims[0]-1 && periods[1]) || (coord[1]>0 && periods[0])) || (periods[0] && periods[1]))
		{
			if(periods[0] && periods[1]&& coord[0]==dims[0]-1 && coord[1]==0)
			{
				coordtmp[0]=0;
				coordtmp[1]=dims[1]-1;
			}
	/*		else
			{
				coordtmp[0]=0;
			}*/
			MPI_Cart_rank (COMM_CART,coordtmp,&BlockNeighbour[SE]);
		}
		else
		{
			BlockNeighbour[SE]=-2;
		}
	else
		MPI_Cart_rank (COMM_CART,coordtmp,&BlockNeighbour[SE]);

	coordtmp[0]=coord[0]+1;
	coordtmp[1]=coord[1]+1;
	if (coordtmp[0]<0 || coordtmp[0]>dims[0]-1||coordtmp[1]<0 || coordtmp[1]>dims[1]-1)
		if(((coord[0]<dims[0]-1 && periods[1]) || (coord[1]<dims[1]-1 && periods[0])) || (periods[0] && periods[1]))
		{
			if(periods[0] && periods[1]&& coord[0]==dims[0]-1 && coord[1]==dims[1]-1)
			{
				coordtmp[0]=0;
				coordtmp[1]=0;
			}
	/*		else if(coord[0]<dims[0]-1 && periods[1])
			{
				coordtmp[0]=0;
			}
			else if(coord[1]>0 && periods[0])
			{

			}*/
			MPI_Cart_rank (COMM_CART,coordtmp,&BlockNeighbour[NE]);
		}
		else
		{
			BlockNeighbour[NE]=-2;
		}
	else
		MPI_Cart_rank (COMM_CART,coordtmp,&BlockNeighbour[NE]);
/*
 	int rank;
 	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
 	char buffer[50]; // make sure it's big enough
 	snprintf(buffer, sizeof(buffer), "BlockConnect_%d.txt", rank);
 	std::ofstream myFlux;
 	myFlux.open(buffer);
 	myFlux<< "Connection to block N is: "<<BlockNeighbour[N]<<std::endl;
 	myFlux<< "Connection to block S is: "<<BlockNeighbour[S]<<std::endl;
 	myFlux<< "Connection to block E is: "<<BlockNeighbour[E]<<std::endl;
 	myFlux<< "Connection to block W is: "<<BlockNeighbour[W]<<std::endl;
 	myFlux<< "Connection to block NE is: "<<BlockNeighbour[NE]<<std::endl;
 	myFlux<< "Connection to block NW is: "<<BlockNeighbour[NW]<<std::endl;
 	myFlux<< "Connection to block SE is: "<<BlockNeighbour[SE]<<std::endl;
 	myFlux<< "Connection to block SW is: "<<BlockNeighbour[SW]<<std::endl;
 	myFlux.close();
*/
	if(verbous)
	{
		printf ("processor : %d my processors neighbour are ", rank_in_topo);
		for (int i=0 ; i < 8; i++ )
			printf (" %d",BlockNeighbour[i]);
		printf (" coordinates are ( %d , %d )\n",coord[0],coord[1]);
	}
	MultiBlock2D::Create_Block2D();

}
void  MultiBlock2D::Create_Block2D() {
	if(verbous)
		std::cout<<"Check rank size "<<parallel->getRank()<<" "<<parallel->getSize()<<std::endl;
/// Calculate the dimensions in the block
		if (dims[0]==1)// one processor
			Nx[0]=Nx_G;
		else
			if(coord[0]==dims[0]-1)//last block in x direction
				if(Nx_G%dims[0]==0) //
					Nx[0]=Nx_G/dims[0];
				else
					Nx[0]=Nx_G-int(Nx_G/dims[0])*(dims[0]-1);
			else
					Nx[0]=int(Nx_G/dims[0]);

		if(dims[1]==1)
			Ny[0]=Ny_G;
		else
			if(coord[1]==dims[1]-1)
				if(Ny_G%dims[1]==0)
					Ny[0]=Ny_G/dims[1];
				else
					Ny[0]=Ny_G-int(Ny_G/dims[1])*(dims[1]-1);
			else
				Ny[0]=Ny_G/dims[1];

/// Share the dimensions in the block to all processors
	parallel->barrier();
	int *Nx_tmp=new int[parallel->getSize()];
	int *Ny_tmp=new int[parallel->getSize()];
	MPI_Allgather(&Nx[0],1,MPI_INT ,&Nx_tmp[0],1,MPI_INT , MPI_COMM_WORLD);
	MPI_Allgather(&Ny[0],1,MPI_INT ,&Ny_tmp[0],1,MPI_INT , MPI_COMM_WORLD);
	parallel->barrier();
	for(int i=0;i<parallel->getSize();i++)
	{
		Nx[i]=Nx_tmp[i];
		Ny[i]=Ny_tmp[i];
	}
	delete [] Nx_tmp;
	delete [] Ny_tmp;
	int *Nx_dim=new int[dims[0]];
	int *Ny_dim=new int[dims[1]];

/// Dimensions of the domains by dimensions of the MPI grid
	for (int i=0;i<dims[0];i++)
		if(i==dims[0]-1)
			Nx_dim[i]=Nx[parallel->getSize()-1];
		else
			Nx_dim[i]=Nx[0];
	for (int i=0;i<dims[1];i++)
		if(i==dims[1]-1)
			Ny_dim[i]=Ny[parallel->getSize()-1];
		else
			Ny_dim[i]=Ny[0];


/// Calculate start and end nodes of the global domain for each processor
/// Calculate start and end elements of the global domain for each processor
	if(dims[0]==1 && dims[1]==1) // Only one processor
	{
		start_nodes=1;
		start_elems=1;
		end_nodes+=(Nx[0]+1)*(Ny[0]+1);
		end_elems+=(Nx[0])*(Ny[0]);
	}

		if(dims[0]>1 && dims[1]==1)// Only processors in x direction
	{
		start_nodes=1;
		start_elems=1;
		if(parallel->getRank()>0)
		{
			start_nodes+=(Nx[0])*(Ny[0]+1);
			start_elems+=(Nx[0])*(Ny[0]);
			if(parallel->getRank()>1)
			{
				for (int i=2;i<=parallel->getRank();i++)
				{
					start_nodes+=(Nx[i-1])*(Ny[i-1]+1);
					start_elems+=(Nx[i-1])*(Ny[i-1]);
				}
			}
		}
		if(parallel->getRank()<parallel->getSize()-1)
		{
			for (int i=0;i<=parallel->getRank();i++)
			{
				end_nodes+=(Nx[i])*(Ny[i]+1);
				end_elems+=(Nx[i])*(Ny[i]);
			}
		}
		else
		{
			for (int i=0;i<=parallel->getRank()-1;i++)
			{
				end_nodes+=(Nx[i])*(Ny[i]+1);
				end_elems+=(Nx[i])*(Ny[i]);
			}
			end_nodes+=(Nx[parallel->getSize()-1]+1)*(Ny[parallel->getSize()-1]+1);
			end_elems+=(Nx[parallel->getSize()-1])*(Ny[parallel->getSize()-1]);
		}
	}
// 2D Grid (Don't touch)
	if(dims[1]>1)
	{
		if(coord[0]==0)
			if(coord[1]==0)
			{
				start_nodes=1;
				start_elems=1;
				end_nodes=(coord[1]+1) * Nx_dim[0] * Ny_dim[0];
				end_elems=(coord[1]+1) * Nx_dim[0] * Ny_dim[0];
			}
			else
			{
				start_nodes= coord[1] * Nx_dim[0] * Ny_dim[0]+1;
				start_elems= coord[1] * Nx_dim[0] * Ny_dim[0]+1;
				if(coord[1]==dims[1]-1)
				{
					end_nodes= Nx_dim[0]*((dims[1]-1)*Ny_dim[0]+Ny_dim[dims[1]-1]+1);
					end_elems= Nx_dim[0]*((dims[1]-1)*Ny_dim[0]+Ny_dim[dims[1]-1]);
				}
				else
				{
					end_nodes= (coord[1]+1) * Nx_dim[0] * Ny_dim[0];
					end_elems= (coord[1]+1) * Nx_dim[0] * Ny_dim[0];
				}
			}
		else
		{
			if(coord[1]==0)
			{
				start_nodes=coord[0]*Nx_dim[0]*((dims[1]-1)*Ny_dim[0]+Ny_dim[dims[1]-1]+1)+1;
				start_elems=coord[0]*Nx_dim[0]*((dims[1]-1)*Ny_dim[0]+Ny_dim[dims[1]-1])+1;
				end_elems=start_elems -1 + (coord[1]+1) * Nx_dim[coord[0]] * Ny_dim[0];
				if(coord[0]<dims[0]-1)
				{
					end_nodes=start_nodes -1 + (coord[1]+1) * Nx_dim[coord[0]] * Ny_dim[0];
				}
				else
				{
					end_nodes=start_nodes -1 + (coord[1]+1) * (Nx_dim[coord[0]]+1) * Ny_dim[0];
				}
			}
			else
			{
				start_nodes=coord[0]*Nx_dim[0]*((dims[1]-1)*Ny_dim[0]+Ny_dim[dims[1]-1]+1)+1;
				start_elems=coord[0]*Nx_dim[0]*((dims[1]-1)*Ny_dim[0]+Ny_dim[dims[1]-1])+1 + coord[1] * Nx_dim[coord[0]] * Ny_dim[0];
				end_elems=start_elems -1 + Nx_dim[coord[0]] * Ny_dim[coord[1]];
				if(coord[0]<dims[0]-1)
				{
					start_nodes+=coord[1] * Nx_dim[coord[0]] * Ny_dim[0];
					if(coord[1]==dims[1]-1)
					{
						end_nodes= start_nodes -1 + Nx_dim[0]*(Ny_dim[dims[1]-1]+1);
					}
					else
					{
						end_nodes= start_nodes -1 + Nx_dim[0] * Ny_dim[0];
					}
				}
				else
				{
					start_nodes+=coord[1] * (Nx_dim[coord[0]]+1) * Ny_dim[0];
					if(coord[1]==dims[1]-1)
					{
						end_nodes= start_nodes -1 + (Nx_dim[dims[0]-1]+1)*(Ny_dim[dims[1]-1]+1);
					}
					else
					{
						end_nodes= start_nodes -1 + (Nx_dim[dims[0]-1]+1) * Ny_dim[0];
					}
				}

			}
		}

	}

	//if(verbous)
		if(parallel->isMainProcessor())
		{
			for (int i=0;i<parallel->getSize();i++)
				std::cout<< " Nx["<<i<<"] = "<<Nx[i]<< " Ny["<<i<<"] = "<<Ny[i];

			std::cout<< std::endl;

			for (int i=0;i<dims[0];i++)
			{
				std::cout<< " Nx_dim["<<i<<"] = "<<Nx_dim[i];
			}
			for (int j=0;j<dims[1];j++)
			{
				std::cout<< " Ny_dim["<<j<<"] = "<<Ny_dim[j];
			}
			std::cout<< std::endl;
		}



		parallel->barrier();
		if(verbous)
			printf("For Rank: %d  the start nodes is %d and the end nodes is %d\nFor Rank: %d   the start elems is %d and end elems is %d\n",parallel->getRank(),start_nodes, end_nodes,parallel->getRank(),start_elems, end_elems);
//		printf("For Rank: %d  the start nodes is %d and the end nodes is %d\n",parallel->getRank(),start_nodes, end_nodes);
		parallel->barrier();
	int x_start=coord[0]*Nx_dim[0];
	int y_start=coord[1]*Ny_dim[0];
	NodeType bC[4];

	// South side 0
	// East side 1
	// North side 2
	// West side 3
	if(dims[1]==1) // One 1D Cart MPI map
		{
			bC[0]=PtrParameters->Get_GlobalBcType(0);
			bC[2]=PtrParameters->Get_GlobalBcType(2);
			if (dims[0]==1) // one processor
			{
				bC[1]=PtrParameters->Get_GlobalBcType(1);
				bC[3]=PtrParameters->Get_GlobalBcType(3);
			}
			else if(coord[0]==0)
			{
				bC[1]=Interior;
				bC[3]=PtrParameters->Get_GlobalBcType(3);
			}
			else if(coord[0]==dims[0]-1)
			{
				bC[1]=PtrParameters->Get_GlobalBcType(1);
				bC[3]=Interior;
			}
			else
			{
				bC[1]=Interior;
				bC[3]=Interior;
			}
		}
	else //2D Cart MPI Map
	if(coord[0]==0)
	{
		if(coord[1]==0)
		{
			bC[0]=PtrParameters->Get_GlobalBcType(0);
			bC[1]=Interior;
			bC[2]=Interior;
			bC[3]=PtrParameters->Get_GlobalBcType(3);
		}
		else if(coord[1]==dims[1]-1)
		{
			bC[0]=Interior;
			bC[1]=Interior;
			bC[2]=PtrParameters->Get_GlobalBcType(2);
			bC[3]=PtrParameters->Get_GlobalBcType(3);
		}
		else
		{
			bC[0]=Interior;
			bC[1]=Interior;
			bC[2]=Interior;
			bC[3]=PtrParameters->Get_GlobalBcType(3);

		}
	}
	else if(coord[0]==dims[0]-1)
	{
		if(coord[1]==0)
		{
			bC[0]=PtrParameters->Get_GlobalBcType(0);
			bC[1]=PtrParameters->Get_GlobalBcType(1);
			bC[2]=Interior;
			bC[3]=Interior;
		}
		else if(coord[1]==dims[1]-1)
		{
			bC[0]=Interior;
			bC[1]=PtrParameters->Get_GlobalBcType(1);
			bC[2]=PtrParameters->Get_GlobalBcType(2);
			bC[3]=Interior;
		}
		else
		{
			bC[0]=Interior;
			bC[1]=PtrParameters->Get_GlobalBcType(1);
			bC[2]=Interior;
			bC[3]=Interior;
		}
	}
	else if(coord[1]==0)
	{
		bC[0]=PtrParameters->Get_GlobalBcType(0);
		bC[1]=Interior;
		bC[2]=Interior;
		bC[3]=Interior;

	}
	else if(coord[1]==dims[1]-1)
	{
		bC[0]=Interior;
		bC[1]=Interior;
		bC[2]=PtrParameters->Get_GlobalBcType(2);
		bC[3]=Interior;
	}
	else
	{
		bC[0]=Interior;
		bC[1]=Interior;
		bC[2]=Interior;
		bC[3]=Interior;
	}
	 //MPI_Barrier(parallel->getGlobalCommunicator());
	if(verbous)
		if(parallel->isMainProcessor())
			std::cout<<"processor ID: "<<parallel->getRank() <<" Bc[0]: "<< bC[0]<<" Bc[1]: "<<bC[1]<<" Bc[2]: "<<bC[2]<<" Bc[3]: "<< bC[3]<<std::endl;

	Block2D_.AddBlock(Nx[parallel->getRank()],Ny[parallel->getRank()],dims,periods,coord,Ny_G,Nx_dim[0],Ny_dim[0],Nx_dim[dims[0]-1],Ny_dim[dims[1]-1],bC,x_start,y_start,start_nodes,end_nodes, start_elems);

}

Block2D* MultiBlock2D::Get_Block()
{
	return &Block2D_;

}
int MultiBlock2D::Get_Nx(){
	return Nx[parallel->getRank()];
}
int MultiBlock2D::Get_Ny(){
	return Ny[parallel->getRank()];
}
const double* MultiBlock2D::Get_X0(){
	return Block2D_.Get_X0();
}
const double* MultiBlock2D::Get_Y0(){
	return Block2D_.Get_Y0();
}
int MultiBlock2D::Get_nnodes(){
	return Block2D_.Get_nnodes();
}
int* MultiBlock2D::Get_Elems0() {
	return Block2D_.Get_Elems0();
}
int MultiBlock2D::Get_Start_Nodes(){
	return start_nodes;
}
int MultiBlock2D::Get_End_Nodes()	{
	return end_nodes;
}
int MultiBlock2D::Get_Start_Elems(){
	return start_elems;
}
int MultiBlock2D::Get_End_Elems(){
	return end_elems;
}
void MultiBlock2D::Get_Connect_Node(std::vector<int> & IdNodeN_,std::vector<int> & IdNodeE_,std::vector<int> & IdNodeS_,std::vector<int> & IdNodeW_,
							  std::vector<int> & IdNodeSW_,std::vector<int> & IdNodeSE_,std::vector<int> & IdNodeNW_,std::vector<int> & IdNodeNE_)
{
	Block2D_.Get_Connect_Node(IdNodeN_,IdNodeE_,IdNodeS_,IdNodeW_,IdNodeSW_,IdNodeSE_,IdNodeNW_,IdNodeNE_);

}
void MultiBlock2D::Get_Connect_Node(std::vector<int> & IdRNodeN_,std::vector<int> & IdRNodeE_,std::vector<int> & IdRNodeS_,std::vector<int> & IdRNodeW_,
		std::vector<int> & IdGNodeN_,std::vector<int> & IdGNodeE_,std::vector<int> & IdGNodeS_,std::vector<int> & IdGNodeW_,
		std::vector<int> & IdRNodeSW_,std::vector<int> & IdRNodeSE_,std::vector<int> & IdRNodeNW_,std::vector<int> & IdRNodeNE_,
		std::vector<int> & IdGNodeSW_,std::vector<int> & IdGNodeSE_,std::vector<int> & IdGNodeNW_,std::vector<int> & IdGNodeNE_)
{
	Block2D_.Get_CommNodes(IdRNodeN_,IdRNodeE_,IdRNodeS_,IdRNodeW_,
			IdGNodeN_,IdGNodeE_,IdGNodeS_,IdGNodeW_,
			IdRNodeSW_,IdRNodeSE_,IdRNodeNW_,IdRNodeNE_,
			IdGNodeSW_,IdGNodeSE_,IdGNodeNW_,IdGNodeNE_);
}
 void MultiBlock2D::Communication(double **buf_send,double **buf_recv, int *size_buf)
 {
	  /* Send to Neighbour E and receive from Neighbour W */
	 int tag_x_r=1;
	 int tag_x_l=2;
	 int tag_y_t=3;
	 int tag_y_b=4;
	 int tag_d_tr=5;
	 int tag_d_tl=6;
	 int tag_d_br=7;
	 int tag_d_bl=8;
	 MPI_Status status;

	 //MPI_Allgather(&Nx[0],1,MPI_INT ,&Nx_tmp[0],1,MPI_INT , MPI_COMM_WORLD); // static_cast<void*>(sendBuf)
	//  parallel->sendRecv(buf_send[0],buf_recv[0],size_buf[0],BlockNeighbour[E],size_buf[1],BlockNeighbour[W],1);
	 //if(BlockNeighbour[W]>=0)
	// MPI_Sendrecv(&buf_send[0][0],size_buf[0],MPI_DOUBLE,BlockNeighbour[E],tag_x_r,&buf_recv[0][0],size_buf[1],MPI_DOUBLE,BlockNeighbour[W],tag_x_r,MPI_COMM_WORLD,&status);
	// MPI_Sendrecv(&buf_send[1][0],size_buf[0],MPI_DOUBLE,BlockNeighbour[W],tag_x_r,&buf_recv[1][0],size_buf[0],MPI_DOUBLE,BlockNeighbour[E],tag_x_r,MPI_COMM_WORLD,&status);
	 if(BlockNeighbour[E]>=0)
	 {
		 MPI_Send(&buf_send[0][0],size_buf[0],MPI_DOUBLE,BlockNeighbour[E],tag_x_r,parallel->getGlobalCommunicator());
	 }
	 if(BlockNeighbour[W]>=0)
	 {
		 MPI_Recv(&buf_recv[0][0],size_buf[4],MPI_DOUBLE,BlockNeighbour[W],tag_x_r,parallel->getGlobalCommunicator(),&status);
		MPI_Send(&buf_send[1][0],size_buf[1],MPI_DOUBLE,BlockNeighbour[W],tag_x_l,parallel->getGlobalCommunicator());
	 }
	 if(BlockNeighbour[E]>=0)
	 {
		MPI_Recv(&buf_recv[1][0],size_buf[5],MPI_DOUBLE,BlockNeighbour[E],tag_x_l,parallel->getGlobalCommunicator(),&status);
	 }
	 if(BlockNeighbour[N]>=0)
	 {
		 MPI_Send(&buf_send[3][0],size_buf[3],MPI_DOUBLE,BlockNeighbour[N],tag_y_t,parallel->getGlobalCommunicator());
	 }
	 if(BlockNeighbour[S]>=0)
	 {
		 MPI_Recv(&buf_recv[3][0],size_buf[7],MPI_DOUBLE,BlockNeighbour[S],tag_y_t,parallel->getGlobalCommunicator(),&status);
		MPI_Send(&buf_send[2][0],size_buf[2],MPI_DOUBLE,BlockNeighbour[S],tag_y_b,parallel->getGlobalCommunicator());
	 }
	 if(BlockNeighbour[N]>=0)
	 {
		MPI_Recv(&buf_recv[2][0],size_buf[6],MPI_DOUBLE,BlockNeighbour[N],tag_y_b,parallel->getGlobalCommunicator(),&status);
	 }

	  /* Send to Neighbour W and receive from Neighbour E */
	  //parallel->sendRecv(buf_send[1],buf_recv[1],size_buf[1],BlockNeighbour[W],size_buf[0],BlockNeighbour[E],2);

	  /* Send to Neighbour S and receive from Neighbour N */
	  //parallel->sendRecv(buf_send[2],buf_recv[2],size_buf[3],BlockNeighbour[S],size_buf[2],BlockNeighbour[N],3);

	  /* Send to Neighbour N and receive from Neighbour S */
	  //parallel->sendRecv(buf_send[3],buf_recv[3],size_buf[2],BlockNeighbour[N],size_buf[3],BlockNeighbour[S],4);
 }
 void MultiBlock2D::CommunicationToGhost(double **buf_send,double **buf_recv, int *size_buf)
 {
	 int tag_x_r=1;
	 int tag_x_l=2;
	 int tag_y_t=3;
	 int tag_y_b=4;
	 int tag_d_tr=5;
	 int tag_d_tl=6;
	 int tag_d_br=7;
	 int tag_d_bl=8;
	 MPI_Status status;

	 if(BlockNeighbour[W]>=0)
	 {
		MPI_Send(&buf_send[1][0],size_buf[1],MPI_DOUBLE,BlockNeighbour[W],tag_x_l,parallel->getGlobalCommunicator());
	 }
	 if(BlockNeighbour[E]>=0)
	 {
		MPI_Recv(&buf_recv[1][0],size_buf[0],MPI_DOUBLE,BlockNeighbour[E],tag_x_l,parallel->getGlobalCommunicator(),&status);
	 }
	 if(BlockNeighbour[S]>=0)
	 {
		MPI_Send(&buf_send[2][0],size_buf[2],MPI_DOUBLE,BlockNeighbour[S],tag_y_b,parallel->getGlobalCommunicator());
	 }
	 if(BlockNeighbour[N]>=0)
	 {
		MPI_Recv(&buf_recv[2][0],size_buf[3],MPI_DOUBLE,BlockNeighbour[N],tag_y_b,parallel->getGlobalCommunicator(),&status);
	 }
}
 void MultiBlock2D::CommunicationFromGhost(double **buf_send,double **buf_recv, int *size_buf)
  {
 	 int tag_x_r=1;
 	 int tag_x_l=2;
 	 int tag_y_t=3;
 	 int tag_y_b=4;

 	 MPI_Status status;
 	 if(BlockNeighbour[E]>=0)
 	 {
 		 MPI_Send(&buf_send[0][0],size_buf[0],MPI_DOUBLE,BlockNeighbour[E],tag_x_r,parallel->getGlobalCommunicator());
 	 }
 	 if(BlockNeighbour[W]>=0)
 	 {
 		 MPI_Recv(&buf_recv[0][0],size_buf[1],MPI_DOUBLE,BlockNeighbour[W],tag_x_r,parallel->getGlobalCommunicator(),&status);
 	 }
 	 if(BlockNeighbour[N]>=0)
 	 {
 		 MPI_Send(&buf_send[3][0],size_buf[3],MPI_DOUBLE,BlockNeighbour[N],tag_y_t,parallel->getGlobalCommunicator());
 	 }
 	 if(BlockNeighbour[S]>=0)
 	 {
 		 MPI_Recv(&buf_recv[3][0],size_buf[2],MPI_DOUBLE,BlockNeighbour[S],tag_y_t,parallel->getGlobalCommunicator(),&status);
 	 }
  }
 void MultiBlock2D::CommunicationFromGhost(double **buf_send,double **buf_recv, int *size_buf,MPI_Status * status,MPI_Request * request)
 {
	 int tag_x_r=1;
	 int tag_x_l=2;
	 int tag_y_t=3;
	 int tag_y_b=4;

	 if(BlockNeighbour[E]>=0)
	 {
		 MPI_Isend(&buf_send[0][0],size_buf[0],MPI_DOUBLE,BlockNeighbour[E],tag_x_r,parallel->getGlobalCommunicator(),&request[0]);

	 }
	 if(BlockNeighbour[W]>=0)
	 {
		 MPI_Irecv(&buf_recv[0][0],size_buf[1],MPI_DOUBLE,BlockNeighbour[W],tag_x_r,parallel->getGlobalCommunicator(),&request[0]);
	 }
	 if(BlockNeighbour[N]>=0)
	 {
		 MPI_Isend(&buf_send[3][0],size_buf[3],MPI_DOUBLE,BlockNeighbour[N],tag_y_t,parallel->getGlobalCommunicator(),&request[1]);

	 }
	 if(BlockNeighbour[S]>=0)
	 {
		 MPI_Irecv(&buf_recv[3][0],size_buf[2],MPI_DOUBLE,BlockNeighbour[S],tag_y_t,parallel->getGlobalCommunicator(),&request[1]);
	 }

 }
 void MultiBlock2D::Get_NbNodes(int & NbRealNodes_, int & NbTotalNodes_){
	 Block2D_.Get_NbNodes(NbRealNodes_,NbTotalNodes_);
 }
 void MultiBlock2D::Send(double *buf_send, int size_buf,int blocDirection,int tag)
 {
	 if(BlockNeighbour[blocDirection]>=0)
		 MPI_Send(buf_send,size_buf,MPI_DOUBLE,BlockNeighbour[blocDirection],tag,parallel->getGlobalCommunicator());
 }
 void MultiBlock2D::Recv(double *buf_recv, int size_buf,int blocDirection,int tag,MPI_Status status)
 {
	 if(BlockNeighbour[blocDirection]>=0)
		 MPI_Recv(buf_recv,size_buf,MPI_DOUBLE,BlockNeighbour[blocDirection],tag,parallel->getGlobalCommunicator(),&status);
 }
 void MultiBlock2D::Modify_Block(){
		if(parallel->isMainProcessor())
			std::cout<<"******BEGGINING MODIFY THE MESH BY THE USERS******"<<std::endl;
		int NbRealNodes,NbTotalNodes;
		Block2D_.Get_NbNodes(NbRealNodes,NbTotalNodes);
		if(PtrParameters->Get_Verbous())
			std::cout<<"Processor ID: "<<parallel->getRank()<<" Number of node:"<<NbRealNodes<<std::endl;
		if(parallel->isMainProcessor())
			std::cout<<"******START GENERATING SOLID IN THE DOMAIN******"<<std::endl;
		Block2D_.GenerateSolid(*PtrParameters);
		if(parallel->isMainProcessor())
			std::cout<<"******END GENERATING SOLID IN THE DOMAIN******"<<std::endl
			<<"******START CHECKING THE MESH AND REMOVING/ADDING SOLID******"<<std::endl;
		int nbTotalSolidRemoved=1;
		int nbTotalSolidadded=1;
		int maxclean=10;int count=0;
		bool check_cleaning=true;
		while( count<maxclean && check_cleaning)
		{
			Block2D_.RemoveUnphysicalSolid(nbTotalSolidRemoved,nbTotalSolidadded);
			//std::cout<<"processor: "<<parallel->getRank()<<" number of solid removed: "<<nbTotalSolidRemoved<<" number of solid added: "<<nbTotalSolidadded<<std::endl;
			MultiBlock2D::Correct_SolidGhost();
			check_cleaning=Check_CleaningMesh(nbTotalSolidRemoved,nbTotalSolidadded);
			//std::cout<<"processor: "<<parallel->getRank()<<" check cleaning: "<<check_cleaning<<std::endl;
			count++;
		}
		//MPI_Barrier(parallel->getGlobalCommunicator());
		if(parallel->isMainProcessor())
			std::cout<<"Processor: "<<parallel->getRank() <<" Number of cleaning: "<<count<<std::endl
			<<"******END CHECKING THE MESH AND REMOVING/ADDING SOLID******"<<std::endl
			<<"******START GENERATING WALL AND CORNERS******"<<std::endl;
		Block2D_.SetSolidBoundaries();
//	 Block2D_.ModifyMeshByUser(*PtrParameters);
	 MultiBlock2D::Correct_SolidGhost();
//	 Block2D_.RemoveSolidInCommunicator();
	 if(parallel->isMainProcessor())
		 std::cout<<"******END GENERATING WALL AND CORNERS******"<<std::endl;
 }
 void MultiBlock2D::ConvertToPhysicalUnit(){
	 Block2D_.ConvertToPhysicalUnit(*PtrParameters);
 }
NodeArrays2D* MultiBlock2D::Get_NodeArrays2D(){
	return Block2D_.Get_NodeArrays2D();
}
void MultiBlock2D::reorganizeNodeByType(){
	Block2D_.reorganizeNodeByType();
	Block2D_.Set_Connect(*PtrParameters);
	Block2D_.Mark1stLayerSolid();

}

void MultiBlock2D::Correct_SolidGhost()
{
	 int tag_x_r=1;
	 int tag_x_l=2;
	 int tag_y_t=3;
	 int tag_y_b=4;
	 std::vector<int> IdRNodeN,IdRNodeE,IdRNodeS,IdRNodeW,IdRNodeSW,IdRNodeSE,IdRNodeNW,IdRNodeNE;
	 std::vector<int> IdGNodeN,IdGNodeE,IdGNodeS,IdGNodeW,IdGNodeSW,IdGNodeSE,IdGNodeNW,IdGNodeNE;
	 Block2D_.Get_GhostType(IdRNodeN,IdRNodeE,IdRNodeS,IdRNodeW,IdRNodeSW,IdRNodeSE,IdRNodeNW,IdRNodeNE);
	 IdGNodeN=IdRNodeN;
	 IdGNodeE= IdRNodeE;
	 IdGNodeS=IdRNodeS;
	 IdGNodeW=IdRNodeW;
	 IdGNodeSW=IdRNodeSW;
	 IdGNodeSE=IdRNodeSE;
	 IdGNodeNW=IdRNodeNW;
	 IdGNodeNE=IdRNodeNE;
	 MPI_Status status;
/*	 	int rank;
	 	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	 	char buffer[50]; // make sure it's big enough
	 	snprintf(buffer, sizeof(buffer), "ConnectionMarks_%d.txt", rank);
	 	std::ofstream myFlux;
	 	myFlux.open(buffer);
	 	myFlux<<std::endl<<" ****** after Removing nodes *******"<<std::endl;
	 	myFlux<<"West real nodes: "<<std::endl;
	 	for(int i=0;i<IdRNodeW.size();i++)
	 		myFlux<<IdRNodeW[i]<<" ";
	 	myFlux<<std::endl<<"West Ghost nodes: "<<std::endl;
	 	for(int i=0;i<IdGNodeW.size();i++)
	 		myFlux<<IdGNodeW[i]<<" ";
	 	myFlux<<std::endl<<"North real nodes: "<<std::endl;
	 	for(int i=0;i<IdRNodeN.size();i++)
	 		myFlux<<IdRNodeN[i]<<" ";
	 	myFlux<<std::endl<<"North Ghost nodes: "<<std::endl;
	 	for(int i=0;i<IdGNodeN.size();i++)
	 		myFlux<<IdGNodeN[i]<<" ";
	 	myFlux<<std::endl<<"South real nodes: "<<std::endl;
	 	for(int i=0;i<IdRNodeS.size();i++)
	 		myFlux<<IdRNodeS[i]<<" ";
	 	myFlux<<std::endl<<"South Ghost nodes: "<<std::endl;
	 	for(int i=0;i<IdGNodeS.size();i++)
	 		myFlux<<IdGNodeS[i]<<" ";
	 	myFlux<<std::endl<<"East real nodes: "<<std::endl;
	 	for(int i=0;i<IdRNodeE.size();i++)
	 		myFlux<<IdRNodeE[i]<<" ";
	 	myFlux<<std::endl<<"East Ghost nodes: "<<std::endl;
	 	for(int i=0;i<IdGNodeE.size();i++)
	 		myFlux<<IdGNodeE[i]<<" ";

	 	myFlux<<std::endl<<std::endl;
	 	myFlux<<"South West real nodes: "<<std::endl;

	 	for(int i=0;i<IdRNodeSW.size();i++)
	 		myFlux<<IdRNodeSW[i]<<" ";
	 	myFlux<<std::endl<<"South West Ghost nodes: "<<std::endl;
	 	for(int i=0;i<IdGNodeSW.size();i++)
	 		myFlux<<IdGNodeSW[i]<<" ";
	 	myFlux<<std::endl<<"North West real nodes: "<<std::endl;
	 	for(int i=0;i<IdRNodeNW.size();i++)
	 		myFlux<<IdRNodeNW[i]<<" ";
	 	myFlux<<std::endl<<"North West Ghost nodes: "<<std::endl;
	 	for(int i=0;i<IdGNodeNW.size();i++)
	 		myFlux<<IdGNodeNW[i]<<" ";
	 	myFlux<<std::endl<<"South East real nodes: "<<std::endl;
	 	for(int i=0;i<IdRNodeSE.size();i++)
	 		myFlux<<IdRNodeSE[i]<<" ";
	 	myFlux<<std::endl<<"South East Ghost nodes: "<<std::endl;
	 	for(int i=0;i<IdGNodeSE.size();i++)
	 		myFlux<<IdGNodeSE[i]<<" ";
	 	myFlux<<std::endl<<"North East real nodes: "<<std::endl;
	 	for(int i=0;i<IdRNodeNE.size();i++)
	 		myFlux<<IdRNodeNE[i]<<" ";
	 	myFlux<<std::endl<<"North East Ghost nodes: "<<std::endl;
	 	for(int i=0;i<IdGNodeNE.size();i++)
	 		myFlux<<IdGNodeNE[i]<<" ";
	 	myFlux<<std::endl<<" Size of NE Ghost: "<< IdGNodeNE.size()<<std::endl;*/
/*	 	int rank;
	 	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	 	char buffer[50]; // make sure it's big enough
	 	snprintf(buffer, sizeof(buffer), "ConnectionBlocks_%d.txt", rank);
	 	std::ofstream myFlux;
	 	myFlux.open(buffer);

	 	myFlux<<"Connection Bloc W: "<<BlockNeighbour[W]<<std::endl;
	 	myFlux<<"Connection Bloc E: "<<BlockNeighbour[E]<<std::endl;
	 	myFlux<<"Connection Bloc N: "<<BlockNeighbour[N]<<std::endl;
	 	myFlux<<"Connection Bloc S: "<<BlockNeighbour[S]<<std::endl;

	 	myFlux<<"Connection Bloc NW: "<<BlockNeighbour[NW]<<std::endl;
	 	myFlux<<"Connection Bloc NE: "<<BlockNeighbour[NE]<<std::endl;
	 	myFlux<<"Connection Bloc SW: "<<BlockNeighbour[SW]<<std::endl;
	 	myFlux<<"Connection Bloc SE: "<<BlockNeighbour[SE]<<std::endl;*/

	 if(BlockNeighbour[W]>=0)
	 {
		 MPI_Send(&IdRNodeW[0],IdRNodeW.size(),MPI_INT,BlockNeighbour[W],tag_x_l,parallel->getGlobalCommunicator());

	 }
	 if(BlockNeighbour[E]>=0)
	 {
		 MPI_Recv(&IdGNodeE[0],IdGNodeE.size(),MPI_INT,BlockNeighbour[E],tag_x_l,parallel->getGlobalCommunicator(),&status);
		 MPI_Send(&IdRNodeE[0],IdRNodeE.size(),MPI_INT,BlockNeighbour[E],tag_x_l,parallel->getGlobalCommunicator());
	 }
	 if(BlockNeighbour[W]>=0)
	 {
		 MPI_Recv(&IdGNodeW[0],IdGNodeW.size(),MPI_INT,BlockNeighbour[W],tag_x_l,parallel->getGlobalCommunicator(),&status);
	 }
	  if(BlockNeighbour[S]>=0)
	 {
		 MPI_Send(&IdRNodeS[0],IdRNodeS.size(),MPI_INT,BlockNeighbour[S],tag_y_b,parallel->getGlobalCommunicator());
	 }
	 if(BlockNeighbour[N]>=0)
	 {
		 MPI_Recv(&IdGNodeN[0],IdGNodeN.size(),MPI_INT,BlockNeighbour[N],tag_y_b,parallel->getGlobalCommunicator(),&status);
		 MPI_Send(&IdRNodeN[0],IdRNodeN.size(),MPI_INT,BlockNeighbour[N],tag_y_b,parallel->getGlobalCommunicator());
	 }
	 if(BlockNeighbour[S]>=0)
	 {
		 MPI_Recv(&IdGNodeS[0],IdGNodeS.size(),MPI_INT,BlockNeighbour[S],tag_y_b,parallel->getGlobalCommunicator(),&status);
	 }
	 int tag_d_bl=5;
	//N=0,E=1,S=2,W=3,NE=4, SE=5, NW=6, SW=7;


	if(BlockNeighbour[SW]>=0)
	{
		MPI_Send(&IdRNodeSW[0],IdRNodeSW.size(),MPI_INT,BlockNeighbour[SW],tag_d_bl,parallel->getGlobalCommunicator());
	}
	if(BlockNeighbour[NE]>=0)
	{
		MPI_Recv(&IdGNodeNE[0],IdGNodeNE.size(),MPI_INT,BlockNeighbour[NE],tag_d_bl,parallel->getGlobalCommunicator(),&status);
		MPI_Send(&IdRNodeNE[0],IdRNodeNE.size(),MPI_INT,BlockNeighbour[NE],tag_d_bl,parallel->getGlobalCommunicator());
	}
	if(BlockNeighbour[SW]>=0)
	{
		MPI_Recv(&IdGNodeSW[0],IdGNodeSW.size(),MPI_INT,BlockNeighbour[SW],tag_d_bl,parallel->getGlobalCommunicator(),&status);
	}

	if(BlockNeighbour[SE]>=0)
	{
		MPI_Send(&IdRNodeSE[0],IdRNodeSE.size(),MPI_INT,BlockNeighbour[SE],tag_d_bl,parallel->getGlobalCommunicator());
	}
	if(BlockNeighbour[NW]>=0)
	{
		MPI_Recv(&IdGNodeNW[0],IdGNodeNW.size(),MPI_INT,BlockNeighbour[NW],tag_d_bl,parallel->getGlobalCommunicator(),&status);
		MPI_Send(&IdRNodeNW[0],IdRNodeNW.size(),MPI_INT,BlockNeighbour[NW],tag_d_bl,parallel->getGlobalCommunicator());
	}
	if(BlockNeighbour[SE]>=0)
	{
		MPI_Recv(&IdGNodeSE[0],IdGNodeSE.size(),MPI_INT,BlockNeighbour[SE],tag_d_bl,parallel->getGlobalCommunicator(),&status);
	}
	 	Block2D_.Set_GhostType(IdGNodeN,IdGNodeE,IdGNodeS,IdGNodeW,IdGNodeSW,IdGNodeSE,IdGNodeNW,IdGNodeNE);

}
double MultiBlock2D::SumBC(double *value){
	double sum=0;
	if(BcW||BcE)
		MPI_Reduce(value,&sum,1, MPI_DOUBLE , MPI_SUM ,0, VertComm);
	if(BcN||BcS)
		MPI_Reduce(value,&sum,1, MPI_DOUBLE , MPI_SUM ,0, HorizComm);
	return sum;
}
double MultiBlock2D::SumAllProcessors(double *value){
	double sum=0;
	MPI_Allreduce(&value[0],&sum,1, MPI_DOUBLE , MPI_SUM ,parallel->getGlobalCommunicator());
	return sum;
}
int MultiBlock2D::NumberOfProcessors(){
	return parallel->getSize();
}
