/*
 * Block2D.cpp
 *
 *  This Object creates the nodes, cells and connections in local and global
 *  starting cell bottom left
 *  directions: [0] down, [1] right, [2] up, [3] left
 *  Created on: 17 Apr 2015
 *      Author: thomas
 */

#include "Block2D.h"

Block2D::Block2D(int dx_, int dy_) : dx(dx_),dy(dy_)
{
	verbous=false;
	NbGhostNode=0;
	NbRealNodes=0;
	NbTotalNodes=0;
	Coord_physical=new double[2];
	periodic[0]=false;
	periodic[1]=false;
	nx=0;ny=0;
}

Block2D::~Block2D() {
	delete [] Coord_physical;
}
///Correct Ordering Ghost node for real cells
void Block2D::Correct_OrderingGhostNode()
{
	bool verbous=false;
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

	char buffer[50]; // make sure it's big enough
	snprintf(buffer, sizeof(buffer), "connections_%d.txt", rank);
	std::ofstream myFlux;
	if(verbous)
		myFlux.open(buffer);

	snprintf(buffer, sizeof(buffer), "Ghost_%d.txt", rank);
	std::ofstream myFlux2;
	if(verbous)
		myFlux2.open(buffer);
	std::vector<int>::iterator it;
	int count1=0; //Use to look through the variable NdGhostNodeInCell (number of Ghost node)
	int count2=Node.size()-NbGhostNode+1; // position of the ghost node in the node array
	int count3=0;//Use to look through the variable IdGhostNodeInCell (ID Ghost node)
	if(verbous)
		myFlux2<< "Node Size "<< Node.size() << " Ghost Size "<< NbGhostNode<<std::endl;

	for(it = IdCellGhostNode.begin() ; it != IdCellGhostNode.end() ; ++it)
	{

		for (int i=0;i<NdGhostNodeInCell[count1];i++)
		{
			CellArray[*it]->Set_NodeNumber(IdGhostNodeInCell[count3],count2);
			// To set the node number number 1 of the next cell shared with the node number 2 of the actual cell
			if (IdGhostNodeInCell[count3]==2 && *it!=CellArray[*it]->Get_Connect(2)[1])
				CellArray[*it+1]->Set_NodeNumber(1,count2);
			// To set the node number number 2 of the next cell shared with the node number 3 of the actual cell
			if (IdGhostNodeInCell[count3]==3 && *it!=CellArray[*it]->Get_Connect(3)[1])
				CellArray[*it+1]->Set_NodeNumber(2,count2);
			if(verbous)
			{
				myFlux2<<"Processor number is: "<<rank<<" Cell with ghost node, Id is: "<<*it<<" Number in cell is: "<<IdGhostNodeInCell[count3]<<" Node number is: "<<count2-1<<" Number of Ghost node in the cell: "<<NdGhostNodeInCell[count1]<<std::endl;
			}
			count3++;
			count2++;
		}
		count1++;


	}

}

///Add Ghost nodes in the Node array and connect or update connections
void Block2D::Set_Connect() {

	bool verbous=false;
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

	char buffer[50]; // make sure it's big enough
	snprintf(buffer, sizeof(buffer), "connections_%d.txt", rank);
	std::ofstream myFlux;
	if(verbous)
		myFlux.open(buffer);

	snprintf(buffer, sizeof(buffer), "Ghost_%d.txt", rank);
	std::ofstream myFlux2;
	if(verbous)
		myFlux2.open(buffer);





	for (int i=0;i<CellArray.size();i++)
	{
		//CellArray[CellArray[i]->Get_Connect(0)[1]]->Get_NodeNumber(0);
		if(verbous)
		{
			myFlux<<"Processor number is: "<<rank<<" Cell number is: "<<i<<" South Connect: "<<CellArray[i]->Get_Connect(0)[1]<<" East Connect: "<<CellArray[i]->Get_Connect(1)[1]<<" North Connect: "<<CellArray[i]->Get_Connect(2)[1]<<" West Connect: "<<CellArray[i]->Get_Connect(3)[1] <<std::endl;
			myFlux<<"Processor number is: "<<rank<<" Node number 0 is: "<<CellArray[i]->Get_NodeNumber(0)-1;
		}
		Node[CellArray[i]->Get_NodeNumber(0)-1]->Set_Connect(
				CellArray[i]->Get_NodeNumber(3),
				CellArray[CellArray[i]->Get_Connect(0)[1]]->Get_NodeNumber(0),
				CellArray[CellArray[i]->Get_Connect(3)[1]]->Get_NodeNumber(0),
				CellArray[i]->Get_NodeNumber(1));

		if(verbous)
		{
			myFlux<<" Node number 1 is: "<<CellArray[i]->Get_NodeNumber(1)-1;
		}
		Node[CellArray[i]->Get_NodeNumber(1)-1]->Set_Connect(
				CellArray[i]->Get_NodeNumber(2),
				CellArray[CellArray[i]->Get_Connect(0)[1]]->Get_NodeNumber(1),
				CellArray[i]->Get_NodeNumber(0),
				CellArray[CellArray[i]->Get_Connect(1)[1]]->Get_NodeNumber(1));

		if(verbous)
		{
			myFlux<<" Node number 2 is: "<<CellArray[i]->Get_NodeNumber(2)-1;
		}
		Node[CellArray[i]->Get_NodeNumber(2)-1]->Set_Connect(
				CellArray[CellArray[i]->Get_Connect(2)[1]]->Get_NodeNumber(2),
				CellArray[i]->Get_NodeNumber(1),
				CellArray[i]->Get_NodeNumber(3),
				CellArray[CellArray[i]->Get_Connect(1)[1]]->Get_NodeNumber(2));

		if(verbous)
		{
			myFlux<<" Node number 3 is: "<<CellArray[i]->Get_NodeNumber(3)-1<<std::endl;
		}
		Node[CellArray[i]->Get_NodeNumber(3)-1]->Set_Connect(
				CellArray[CellArray[i]->Get_Connect(2)[1]]->Get_NodeNumber(3),
				CellArray[i]->Get_NodeNumber(0),
				CellArray[CellArray[i]->Get_Connect(3)[1]]->Get_NodeNumber(3),
				CellArray[i]->Get_NodeNumber(2));

	}


}
/// Method to create a new Cell of Quad type
void Block2D::NewCell(NodeType Nodes[4], bool NewNodes[4],bool GhostNodes[4], int Node2D_SubDomain[4],int Node2D_GlobalDomain[4], int x_[4],int y_[4], int FaceConnect[4], int CellConnect[4])
{
/// Create a new element by 4 nodes (global numbering)
	Elems.push_back(Node2D_GlobalDomain[0]);
	Elems.push_back(Node2D_GlobalDomain[1]);
	Elems.push_back(Node2D_GlobalDomain[2]);
	Elems.push_back(Node2D_GlobalDomain[3]);

/// Create a cell with old/new nodes and separate ghost nodes
	Cell2D* cell=new Cell2D();

	int NbGhostNode=0;
	for(int i=0;i<4;i++)
	{

		if (NewNodes[i])
		{
			if(GhostNodes[i])
			{
				Node_Ghosttmp.push_back(cell->NewNode(Ghost,x_[i],y_[i]));
				x_Ghosttmp.push_back(x_[i]);
				y_Ghosttmp.push_back(y_[i]);
				IdGhostNodeInCell.push_back(i);
				NbGhostNode++;
			}
			else
			{
				Node.push_back(cell->NewNode(Nodes[i],x_[i],y_[i]));
				Node.back()->Set_NodeType(Nodes[i]);
				x.push_back(x_[i]);
				y.push_back(y_[i]);
				if(verbous)
				std::cout<<"Node number is: "<<Node2D_GlobalDomain[i] <<" Node type: "<<Nodes[i]<<std::endl;
				if(Nodes[i]!=Interior)
					IdBoundaries.push_back(Node.size()-1);
			}
			Node2D_SubDomain[i]=Node.size();
		}
	}

	if (NbGhostNode>0)
		NdGhostNodeInCell.push_back(NbGhostNode);

// Set local node numbering in the cell
	cell->Set_NodeNumber(Node2D_SubDomain);
// Set connections (Faces and Cells)
	cell->Set_Connect(0,FaceConnect[0],CellConnect[0]);
	cell->Set_Connect(1,FaceConnect[1],CellConnect[1]);
	cell->Set_Connect(2,FaceConnect[2],CellConnect[2]);
	cell->Set_Connect(3,FaceConnect[3],CellConnect[3]);
// Store the new cell
	CellArray.push_back(cell);


}

void Block2D::NewGhostCell(NodeType Nodes[4], bool NewNodes[4], int Node2D_SubDomain[4], int x_[4],int y_[4], int FaceConnect[4], int CellConnect[4])
{
/// Create a cell with old/new nodes and separate ghost nodes
	Cell2D* cell=new Cell2D();

	int NbGhostNode=0;
	for(int i=0;i<4;i++)
	{
		if (NewNodes[i])
		{
			Node.push_back(cell->NewNode(Ghost,x_[i],y_[i]));
			x.push_back(x_[i]);
			y.push_back(y_[i]);
			IdGhostNodeInCell.push_back(i);
			NbGhostNode++;
			Node2D_SubDomain[i]=Node.size();
		}
	}

	if (NbGhostNode>0)
		NdGhostNodeInCell.push_back(NbGhostNode);

// Set local node numbering in the cell
	cell->Set_NodeNumber(Node2D_SubDomain);
// Set connections (Faces and Cells)
	cell->Set_Connect(0,FaceConnect[0],CellConnect[0]);
	cell->Set_Connect(1,FaceConnect[1],CellConnect[1]);
	cell->Set_Connect(2,FaceConnect[2],CellConnect[2]);
	cell->Set_Connect(3,FaceConnect[3],CellConnect[3]);
// Store the new cell
	CellArray.push_back(cell);
}

/// Method to generate a mesh block inside a MPI Block.
void Block2D::AddBlock(int TotalNumberCells_x,int TotalNumberCells_y,int dims[2],int periodic_[2],int coord[2],int NyCell_G,int Nx_begin,int Ny_begin,int Nx_last,int Ny_last,NodeType bC[4],int x_last,int y_last,int start_GlobalNodes,int end_GlobalNodes, int start_GlobalElems, int end_GlobalElems)
{
	verbous=false;
	NodeType Nodes[4];
	int x_[4];
	int y_[4];
	int x_tmp=x_last,y_tmp=y_last;
	int FaceConnect[4];
	int CellConnect[4];
	bool NewNodes[4];
	bool GhostNodes[4];
	Coord[0]=coord[0];
	Coord[1]=coord[1];
	Dims[0]=dims[0];
	Dims[1]=dims[1];
	periodic[0]=periodic_[0];
	periodic[1]=periodic_[1];
	bool dimsM1=dims[0]-2==coord[0];
	int Node2D_SubDomain[4];
	int Node2D_GlobalDomain[4];
	int nNode2D_SubDomain=1;
	int nComputationalNode=1;
	int nNodes_Y=(NyCell_G+1)*TotalNumberCells_x;
	GhostNodes[0]=false;
	GhostNodes[1]=false;
	GhostNodes[2]=false;
	GhostNodes[3]=false;
/// Calculation the global numbering
	int start_GlobalNodes_Domain_N;
	int start_GlobalNodes_Domain_E;

	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

	if(coord[1]==dims[1]-1)// Treat the MPI blocks at the top of the physical domain
	{
		start_GlobalNodes_Domain_N=start_GlobalNodes+TotalNumberCells_x*(TotalNumberCells_y+1);
		if(dimsM1)// only previous last block in x directions
			start_GlobalNodes_Domain_E=start_GlobalNodes_Domain_N+(dims[1]-1)*(Nx_last+1)*Ny_begin;
		else if(coord[0]==dims[0]-1) //treat the last block in the x direction
		{
			start_GlobalNodes_Domain_N=start_GlobalNodes+(Nx_last+1)*(Ny_last+1);
			start_GlobalNodes_Domain_E=start_GlobalNodes+Nx_last*(Ny_last);
		}
		else // Other top blocks
			start_GlobalNodes_Domain_E=start_GlobalNodes+(dims[1]-1)*Nx_begin*Ny_begin+Nx_begin*(Ny_last+1);

	}
	else
	{
		if(coord[0]==dims[0]-1)//treat the last blocks in the x direction
		{
			start_GlobalNodes_Domain_N=start_GlobalNodes+(TotalNumberCells_x+1)*TotalNumberCells_y;
			start_GlobalNodes_Domain_E=start_GlobalNodes+Nx_last*Ny_begin;
		}
		else// Other blocks in the domain
		{
			start_GlobalNodes_Domain_N=start_GlobalNodes+TotalNumberCells_x*TotalNumberCells_y;
			if(dimsM1)// only block on x directions
				start_GlobalNodes_Domain_E=start_GlobalNodes+(dims[1]-1-coord[1])*Nx_begin*Ny_begin+Nx_begin*(Ny_last+1)+coord[1]*(Nx_last+1)*(Ny_begin);
			else
				start_GlobalNodes_Domain_E=start_GlobalNodes+(dims[1]-1)*Nx_begin*Ny_begin+Nx_begin*(Ny_last+1);
		}
	}
	if(verbous)
		std::cout<<"start Global Nodes Next Domain is: "<<start_GlobalNodes_Domain_N<<" start Global Nodes East Domain is: "<<start_GlobalNodes_Domain_E<<std::endl;

/// Create Cells, nodes and connections starting from the West bottom cell (0,0)
	if (CellArray.empty())
	{

		int nbCells=0;
		// Cell (0,0)
		if (bC[3]==Interior)
			Nodes[0]=bC[0];
		else
			if(bC[0]==Interior)
				Nodes[0]=bC[3];
			else
				Nodes[0]=GlobalCorner;
		Nodes[1]=bC[0];
		Nodes[2]=Interior;
		Nodes[3]=bC[3];

		// Coordinate of the cell
		x_[0]=x_tmp;
		x_[1]=x_tmp+dx;
		x_[2]=x_tmp+dx;
		x_[3]=x_tmp;

		y_[0]=y_tmp;
		y_[1]=y_tmp;
		y_[2]=y_tmp+dy;
		y_[3]=y_tmp+dy;

		//Face connections (local)
		FaceConnect[0]=2;
		FaceConnect[1]=3;
		FaceConnect[2]=0;
		FaceConnect[3]=1;

		// Cell connections (local)
		CellConnect[0]=nbCells; //no connection
		CellConnect[1]=TotalNumberCells_y-1;
		CellConnect[2]=1;
		CellConnect[3]=nbCells; //no connection

		//Determine if we need a new node or already exist
		NewNodes[0]=true;
		NewNodes[1]=true;
		NewNodes[2]=true;
		NewNodes[3]=true;

		// Determine local numbering

		Node2D_SubDomain[0]=nNode2D_SubDomain++;
		Node2D_SubDomain[1]=nNode2D_SubDomain++;
		Node2D_SubDomain[2]=nNode2D_SubDomain++;
		Node2D_SubDomain[3]=nNode2D_SubDomain++;

		//Calculate global numbering
		Node2D_GlobalDomain[0]=start_GlobalNodes+Node2D_SubDomain[0]-1;
		Node2D_GlobalDomain[1]=start_GlobalNodes+Node2D_SubDomain[1]-1;
		Node2D_GlobalDomain[2]=start_GlobalNodes+Node2D_SubDomain[2]-1;
		Node2D_GlobalDomain[3]=start_GlobalNodes+Node2D_SubDomain[3]-1;

		//Map Local numbering to Global numbering
		for (int i=0;i<4;i++)
		LocalToGlobalNode[Node2D_SubDomain[i]]=Node2D_GlobalDomain[i];

		//Mark Nodes for parallel
		IdNodeSW.push_back(Node2D_SubDomain[0]-1);
		IdNodeS.push_back(Node2D_SubDomain[1]-1);
		IdNodeW.push_back(Node2D_SubDomain[3]-1);

		//Create the cell
		NewCell(Nodes,NewNodes,GhostNodes,Node2D_SubDomain,Node2D_GlobalDomain,x_,y_,FaceConnect,CellConnect);
		y_tmp=y_tmp+dy;
		nbCells++;


		// Cell (0,1:TotalNumberCells_y-1)
	for (int j=1;j<TotalNumberCells_y-1;j++)
		{
			Nodes[0]=bC[3];
			Nodes[1]=Interior;
			Nodes[2]=Interior;
			Nodes[3]=bC[3];

			x_[0]=x_tmp;
			x_[1]=x_tmp+dx;
			x_[2]=x_tmp+dx;
			x_[3]=x_tmp;

			y_[0]=y_tmp;
			y_[1]=y_tmp;
			y_[2]=y_tmp+dy;
			y_[3]=y_tmp+dy;

			FaceConnect[0]=2;
			FaceConnect[1]=3;
			FaceConnect[2]=0;
			FaceConnect[3]=1;

			CellConnect[0]=nbCells-1;
			CellConnect[1]=nbCells+TotalNumberCells_y-1;
			CellConnect[2]=nbCells+1;
			CellConnect[3]=nbCells;

			NewNodes[0]=false;
			NewNodes[1]=false;
			NewNodes[2]=true;
			NewNodes[3]=true;

			Node2D_SubDomain[0]=CellArray[CellConnect[0]]->Get_NodeNumber(3);
			Node2D_SubDomain[1]=CellArray[CellConnect[0]]->Get_NodeNumber(2);
			Node2D_SubDomain[2]=nNode2D_SubDomain++;
			Node2D_SubDomain[3]=nNode2D_SubDomain++;

			Node2D_GlobalDomain[0]=start_GlobalNodes+Node2D_SubDomain[0]-1;
			Node2D_GlobalDomain[1]=start_GlobalNodes+Node2D_SubDomain[1]-1;
			Node2D_GlobalDomain[2]=start_GlobalNodes+Node2D_SubDomain[2]-1;
			Node2D_GlobalDomain[3]=start_GlobalNodes+Node2D_SubDomain[3]-1;
			//Map Local numbering to Global numbering
			for (int i=0;i<4;i++)
			LocalToGlobalNode[Node2D_SubDomain[i]]=Node2D_GlobalDomain[i];

			//Mark Nodes for parallel
			IdNodeW.push_back(Node2D_SubDomain[3]-1);

			NewCell(Nodes,NewNodes,GhostNodes,Node2D_SubDomain,Node2D_GlobalDomain,x_,y_,FaceConnect,CellConnect);
			y_tmp=y_tmp+dy;
			nbCells++;
		}



		x_tmp=x_tmp+dx;

		// Cell (1:TotalNumberCells_x-1,1:TotalNumberCells_y-1)
		if(TotalNumberCells_x>1)
		{
			for (int i=1;i<TotalNumberCells_x-1;i++)
			{
				y_tmp=y_last;
				Nodes[0]=bC[0];
				Nodes[1]=bC[0];
				Nodes[2]=Interior;
				Nodes[3]=Interior;

				x_[0]=x_tmp;
				x_[1]=x_tmp+dx;
				x_[2]=x_tmp+dx;
				x_[3]=x_tmp;

				y_[0]=y_tmp;
				y_[1]=y_tmp;
				y_[2]=y_tmp+dy;
				y_[3]=y_tmp+dy;

				FaceConnect[0]=2;
				FaceConnect[1]=3;
				FaceConnect[2]=0;
				FaceConnect[3]=1;

				CellConnect[0]=nbCells;
				CellConnect[1]=nbCells+TotalNumberCells_y-1;
				CellConnect[2]=nbCells+1;
				CellConnect[3]=nbCells-TotalNumberCells_y+1;

				NewNodes[0]=false;
				NewNodes[1]=true;
				NewNodes[2]=true;
				NewNodes[3]=false;

				Node2D_SubDomain[0]=CellArray[CellConnect[3]]->Get_NodeNumber(1);
				Node2D_SubDomain[1]=nNode2D_SubDomain++;
				Node2D_SubDomain[2]=nNode2D_SubDomain++;
				Node2D_SubDomain[3]=CellArray[CellConnect[3]]->Get_NodeNumber(2);

				Node2D_GlobalDomain[0]=start_GlobalNodes+Node2D_SubDomain[0]-1;
				Node2D_GlobalDomain[1]=start_GlobalNodes+Node2D_SubDomain[1]-1;
				Node2D_GlobalDomain[2]=start_GlobalNodes+Node2D_SubDomain[2]-1;
				Node2D_GlobalDomain[3]=start_GlobalNodes+Node2D_SubDomain[3]-1;
				//Map Local numbering to Global numbering
				for (int k=0;k<4;k++)
				LocalToGlobalNode[Node2D_SubDomain[k]]=Node2D_GlobalDomain[k];

				//Mark Nodes for parallel
				/*if(coord[0]!=dims[0]-1 && i==TotalNumberCells_x-2)
				{
					IdNodeSE.push_back(Node2D_SubDomain[1]-1);
					IdNodeE.push_back(Node2D_SubDomain[2]-1);
				}*/
				if(i==TotalNumberCells_x-2)
				{
					IdNodeS.push_back(Node2D_SubDomain[1]-1);
					IdNodeSE.push_back(Node2D_SubDomain[1]-1);
					IdNodeE.push_back(Node2D_SubDomain[2]-1);
				}
				else
				{
					IdNodeS.push_back(Node2D_SubDomain[1]-1);
				}

				NewCell(Nodes,NewNodes,GhostNodes,Node2D_SubDomain,Node2D_GlobalDomain,x_,y_,FaceConnect,CellConnect);
				y_tmp=y_tmp+dy;
				nbCells++;

				for (int j=1;j<TotalNumberCells_y-1;j++)
				{
					Nodes[0]=Interior;
					Nodes[1]=Interior;
					Nodes[2]=Interior;
					Nodes[3]=Interior;

					x_[0]=x_tmp;
					x_[1]=x_tmp+dx;
					x_[2]=x_tmp+dx;
					x_[3]=x_tmp;

					y_[0]=y_tmp;
					y_[1]=y_tmp;
					y_[2]=y_tmp+dy;
					y_[3]=y_tmp+dy;

					FaceConnect[0]=2;
					FaceConnect[1]=3;
					FaceConnect[2]=0;
					FaceConnect[3]=1;

					CellConnect[0]=nbCells-1;
					CellConnect[1]=nbCells+TotalNumberCells_y-1;
					CellConnect[2]=nbCells+1;
					CellConnect[3]=nbCells-TotalNumberCells_y+1;

					NewNodes[0]=false;
					NewNodes[1]=false;
					NewNodes[2]=true;
					NewNodes[3]=false;

					Node2D_SubDomain[0]=CellArray[CellConnect[0]]->Get_NodeNumber(3);
					Node2D_SubDomain[1]=CellArray[CellConnect[0]]->Get_NodeNumber(2);
					Node2D_SubDomain[2]=nNode2D_SubDomain++;
					Node2D_SubDomain[3]=CellArray[CellConnect[3]]->Get_NodeNumber(2);

					Node2D_GlobalDomain[0]=start_GlobalNodes+Node2D_SubDomain[0]-1;
					Node2D_GlobalDomain[1]=start_GlobalNodes+Node2D_SubDomain[1]-1;
					Node2D_GlobalDomain[2]=start_GlobalNodes+Node2D_SubDomain[2]-1;
					Node2D_GlobalDomain[3]=start_GlobalNodes+Node2D_SubDomain[3]-1;
					//Map Local numbering to Global numbering
					for (int k=0;k<4;k++)
					LocalToGlobalNode[Node2D_SubDomain[k]]=Node2D_GlobalDomain[k];

					//Mark Nodes for parallel
					//if(coord[0]!=dims[0]-1 && i==TotalNumberCells_x-2)
					if(i==TotalNumberCells_x-2)
					{
						IdNodeE.push_back(Node2D_SubDomain[2]-1);
					}

					NewCell(Nodes,NewNodes,GhostNodes,Node2D_SubDomain,Node2D_GlobalDomain,x_,y_,FaceConnect,CellConnect);
					nbCells++;
					y_tmp=y_tmp+dy;
				}
				x_tmp=x_tmp+dx;
			}

			// Cell (TotalNumberCells_x,0)
			nComputationalNode=nNode2D_SubDomain-1;
			y_tmp=y_last;
			Nodes[0]=bC[0];
			if (bC[1]==Interior)
				Nodes[1]=bC[0];
			else
				if(bC[0]==Interior)
					Nodes[1]=bC[1];
				else
					Nodes[1]=GlobalCorner;
			Nodes[2]=bC[1];
			Nodes[3]=Interior;
			if(verbous)
				std::cout<< "Corner check. Type of node 0: "<< Nodes[0]<< " type of node 1: "<< Nodes[1]<< " type of node 2: "<< Nodes[2]<< " type of node 3: "<< Nodes[2]<<std::endl;

			x_[0]=x_tmp;
			x_[1]=x_tmp+dx;
			x_[2]=x_tmp+dx;
			x_[3]=x_tmp;

			y_[0]=y_tmp;
			y_[1]=y_tmp;
			y_[2]=y_tmp+dy;
			y_[3]=y_tmp+dy;

			FaceConnect[0]=2;
			FaceConnect[1]=3;
			FaceConnect[2]=0;
			FaceConnect[3]=1;

			CellConnect[0]=nbCells;
			CellConnect[1]=nbCells;
			CellConnect[2]=nbCells+1;
			CellConnect[3]=nbCells-TotalNumberCells_y+1;

			NewNodes[0]=false;
			NewNodes[1]=true;
			NewNodes[2]=true;
			NewNodes[3]=false;

			Node2D_SubDomain[0]=CellArray[CellConnect[3]]->Get_NodeNumber(1);
			Node2D_SubDomain[1]=nNode2D_SubDomain++;
			Node2D_SubDomain[2]=nNode2D_SubDomain++;
			Node2D_SubDomain[3]=CellArray[CellConnect[3]]->Get_NodeNumber(2);

			if(bC[1]!=Interior)
			{

				GhostNodes[0]=false;
				GhostNodes[1]=false;
				GhostNodes[2]=false;
				GhostNodes[3]=false;
				Node2D_GlobalDomain[0]=start_GlobalNodes+Node2D_SubDomain[0]-1;
				Node2D_GlobalDomain[1]=start_GlobalNodes_Domain_E;
				Node2D_GlobalDomain[2]=start_GlobalNodes_Domain_E+1;
				Node2D_GlobalDomain[3]=start_GlobalNodes+Node2D_SubDomain[3]-1;

				//Map Local numbering to Global numbering
				for (int i=0;i<4;i++)
				LocalToGlobalNode[Node2D_SubDomain[i]]=Node2D_GlobalDomain[i];

			}
			else
			{

				GhostNodes[0]=false;
				GhostNodes[1]=true;
				GhostNodes[2]=true;
				GhostNodes[3]=false;

				//Mark cell with ghost node to correct the node numbering before connectivity
				IdCellGhostNode.push_back(nbCells);

				Node2D_GlobalDomain[0]=start_GlobalNodes+Node2D_SubDomain[0]-1;
				Node2D_GlobalDomain[1]=start_GlobalNodes_Domain_E;
				Node2D_GlobalDomain[2]=start_GlobalNodes_Domain_E+3;
				Node2D_GlobalDomain[3]=start_GlobalNodes+Node2D_SubDomain[3]-1;

				//Map Local numbering to Global numbering
				LocalToGlobalNode[Node2D_SubDomain[0]]=Node2D_GlobalDomain[0];
				LocalToGlobalNode[Node2D_SubDomain[1]]=-1;
				LocalToGlobalNode[Node2D_SubDomain[2]]=-1;
				LocalToGlobalNode[Node2D_SubDomain[3]]=Node2D_GlobalDomain[3];
			}

			//Mark Nodes for parallel
			/*if(coord[0]==dims[0]-1)
			{
				IdNodeSE.push_back(Node2D_SubDomain[1]-1);
				IdNodeE.push_back(Node2D_SubDomain[2]-1);
			}*/
			NewCell(Nodes,NewNodes,GhostNodes,Node2D_SubDomain,Node2D_GlobalDomain,x_,y_,FaceConnect,CellConnect);
			y_tmp=y_tmp+dy;
			nbCells++;
			// Cell (TotalNumberCells_x,1:TotalNumberCells_y-1)
			for (int j=1;j<TotalNumberCells_y-1;j++)
			{
				Nodes[0]=Interior;
				Nodes[1]=bC[1];//Wall;
				Nodes[2]=bC[1];//Wall;
				Nodes[3]=Interior;

				x_[0]=x_tmp;
				x_[1]=x_tmp+dx;
				x_[2]=x_tmp+dx;
				x_[3]=x_tmp;

				y_[0]=y_tmp;
				y_[1]=y_tmp;
				y_[2]=y_tmp+dy;
				y_[3]=y_last+dy;

				FaceConnect[0]=2;
				FaceConnect[1]=3;
				FaceConnect[2]=0;
				FaceConnect[3]=1;

				CellConnect[0]=nbCells-1;
				CellConnect[1]=nbCells;
				CellConnect[2]=nbCells+1;
				CellConnect[3]=nbCells-TotalNumberCells_y+1;

				NewNodes[0]=false;
				NewNodes[1]=false;
				NewNodes[2]=true;
				NewNodes[3]=false;

				Node2D_SubDomain[0]=CellArray[CellConnect[0]]->Get_NodeNumber(3);
				Node2D_SubDomain[1]=CellArray[CellConnect[0]]->Get_NodeNumber(2);
				Node2D_SubDomain[2]=nNode2D_SubDomain++;
				Node2D_SubDomain[3]=CellArray[CellConnect[3]]->Get_NodeNumber(2);

//				std::cout<<"bc outlet "<<bC[1]<<std::endl;

				if(bC[1]!=Interior)
				{

					GhostNodes[0]=false;
					GhostNodes[1]=false;
					GhostNodes[2]=false;
					GhostNodes[3]=false;

					Node2D_GlobalDomain[0]=start_GlobalNodes+Node2D_SubDomain[0]-1;
					Node2D_GlobalDomain[1]=start_GlobalNodes_Domain_E+j;
					Node2D_GlobalDomain[2]=start_GlobalNodes_Domain_E+1+j;
					Node2D_GlobalDomain[3]=start_GlobalNodes+Node2D_SubDomain[3]-1;

					//Map Local numbering to Global numbering
					for (int i=0;i<4;i++)
					LocalToGlobalNode[Node2D_SubDomain[i]]=Node2D_GlobalDomain[i];
				}
				else
				{

					GhostNodes[0]=false;
					GhostNodes[1]=true;
					GhostNodes[2]=true;
					GhostNodes[3]=false;

					//Mark cell with ghost node to correct the node numbering before connectivity
					IdCellGhostNode.push_back(nbCells);


					Node2D_GlobalDomain[0]=start_GlobalNodes+Node2D_SubDomain[0]-1;
					Node2D_GlobalDomain[1]=start_GlobalNodes_Domain_E+2*j+1;
					Node2D_GlobalDomain[2]=start_GlobalNodes_Domain_E+2*j+3;
					Node2D_GlobalDomain[3]=start_GlobalNodes+Node2D_SubDomain[3]-1;


					//Map Local numbering to Global numbering
					LocalToGlobalNode[Node2D_SubDomain[0]]=Node2D_GlobalDomain[0];
					LocalToGlobalNode[Node2D_SubDomain[1]]=-1;
					LocalToGlobalNode[Node2D_SubDomain[2]]=-1;
					LocalToGlobalNode[Node2D_SubDomain[3]]=Node2D_GlobalDomain[3];
				}

				//Mark Nodes for parallel
				/*if(coord[0]==dims[0]-1)
				IdNodeE.push_back(Node2D_SubDomain[2]-1);*/

				NewCell(Nodes,NewNodes,GhostNodes,Node2D_SubDomain,Node2D_GlobalDomain,x_,y_,FaceConnect,CellConnect);
				y_tmp=y_tmp+dy;
				nbCells++;
			}


		}



		// Cell (TotalNumberCells_x,TotalNumberCells_y)
		Nodes[0]=Interior;
		Nodes[1]=bC[1];
		if (bC[1]==Interior)
			Nodes[2]=bC[2];
		else
			Nodes[2]=GlobalCorner;
		Nodes[3]=bC[2];

		x_[0]=x_tmp;
		x_[1]=x_tmp+dx;
		x_[2]=x_tmp+dx;
		x_[3]=x_tmp;

		y_[0]=y_tmp;
		y_[1]=y_tmp;
		y_[2]=y_tmp+dy;
		y_[3]=y_tmp+dy;

		FaceConnect[0]=2;
		FaceConnect[1]=3;
		FaceConnect[2]=0;
		FaceConnect[3]=1;

		CellConnect[0]=nbCells-1;
		CellConnect[1]=nbCells;
		CellConnect[2]=nbCells;
		CellConnect[3]=nbCells+1;

		NewNodes[0]=false;
		NewNodes[1]=false;
		NewNodes[2]=true;
		NewNodes[3]=true;

		Node2D_SubDomain[0]=CellArray[CellConnect[0]]->Get_NodeNumber(3);
		Node2D_SubDomain[1]=CellArray[CellConnect[0]]->Get_NodeNumber(2);
		Node2D_SubDomain[2]=nNode2D_SubDomain++;
		Node2D_SubDomain[3]=nNode2D_SubDomain++;

		if(bC[1]!=Interior) // Block at East side
		{
			if(bC[2]!=Interior) // Block at Top East side
			{

				GhostNodes[0]=false;
				GhostNodes[1]=false;
				GhostNodes[2]=false;
				GhostNodes[3]=false;

				Node2D_GlobalDomain[0]=start_GlobalNodes+Node2D_SubDomain[0]-1;
				Node2D_GlobalDomain[1]=start_GlobalNodes_Domain_E+Ny_last-1;
				Node2D_GlobalDomain[2]=start_GlobalNodes_Domain_E+Ny_last;
				Node2D_GlobalDomain[3]=start_GlobalNodes_Domain_E+Ny_last+1;

				//Map Local numbering to Global numbering
				for (int i=0;i<4;i++)
				LocalToGlobalNode[Node2D_SubDomain[i]]=Node2D_GlobalDomain[i];
			}
			else // Block at East side but not Top
			{
				GhostNodes[0]=false;
				GhostNodes[1]=false;
				GhostNodes[2]=true;
				GhostNodes[3]=true;

				//Mark cell with ghost node to correct the node numbering before connectivity
				IdCellGhostNode.push_back(nbCells);

				Node2D_GlobalDomain[0]=start_GlobalNodes+Node2D_SubDomain[0]-1;
				Node2D_GlobalDomain[1]=start_GlobalNodes_Domain_E+TotalNumberCells_y-1;
				if (dims[1]-2==coord[1])
				{
					Node2D_GlobalDomain[2]=start_GlobalNodes_Domain_N+Nx_last*Ny_last;
					Node2D_GlobalDomain[3]=start_GlobalNodes_Domain_N+(Nx_last-1)*Ny_last;
				}
				else
				{
					Node2D_GlobalDomain[2]=start_GlobalNodes_Domain_N+TotalNumberCells_x*TotalNumberCells_y;
					Node2D_GlobalDomain[3]=start_GlobalNodes_Domain_N+(TotalNumberCells_x-1)*TotalNumberCells_y;
				}
				//Map Local numbering to Global numbering
				LocalToGlobalNode[Node2D_SubDomain[0]]=Node2D_GlobalDomain[0];
				LocalToGlobalNode[Node2D_SubDomain[1]]=Node2D_GlobalDomain[1];
				LocalToGlobalNode[Node2D_SubDomain[2]]=-1;
				LocalToGlobalNode[Node2D_SubDomain[3]]=-1;
			}
		}
		else// Block NOT at East side
		{
			if(bC[2]!=Interior) // Block NOT at East side but in top
			{
				GhostNodes[0]=false;
				GhostNodes[1]=true;
				GhostNodes[2]=true;
				GhostNodes[3]=false;

				//Mark cell with ghost node to correct the node numbering before connectivity
				IdCellGhostNode.push_back(nbCells);

				if(dimsM1)
				{
					Node2D_GlobalDomain[0]=start_GlobalNodes+Node2D_SubDomain[0]-1;
					Node2D_GlobalDomain[1]=start_GlobalNodes_Domain_E+(Ny_last-1)*2+1;
					Node2D_GlobalDomain[2]=start_GlobalNodes_Domain_E+(Nx_last+1)*(Ny_last+1)-1;
					Node2D_GlobalDomain[3]=start_GlobalNodes+Node2D_SubDomain[0];
				}
				else
				{
					Node2D_GlobalDomain[0]=start_GlobalNodes+Node2D_SubDomain[0]-1;
					Node2D_GlobalDomain[1]=start_GlobalNodes_Domain_E+(Ny_last-1)*2+1;
					Node2D_GlobalDomain[2]=start_GlobalNodes_Domain_E+TotalNumberCells_x*(TotalNumberCells_y+1)-1;
					Node2D_GlobalDomain[3]=start_GlobalNodes+Node2D_SubDomain[0];
				}
				//Map Local numbering to Global numbering
				LocalToGlobalNode[Node2D_SubDomain[0]]=Node2D_GlobalDomain[0];
				LocalToGlobalNode[Node2D_SubDomain[1]]=-1;
				LocalToGlobalNode[Node2D_SubDomain[2]]=-1;
				LocalToGlobalNode[Node2D_SubDomain[3]]=Node2D_GlobalDomain[3];
			}
			else
			{

				GhostNodes[0]=false;
				GhostNodes[1]=true;
				GhostNodes[2]=true;
				GhostNodes[3]=true;

				//Mark cell with ghost node to correct the node numbering before connectivity
				IdCellGhostNode.push_back(nbCells);

				if(dimsM1)
				{
					Node2D_GlobalDomain[0]=start_GlobalNodes+Node2D_SubDomain[0]-1;
					Node2D_GlobalDomain[1]=start_GlobalNodes_Domain_E+(TotalNumberCells_y-1)*2+1;
					Node2D_GlobalDomain[2]=start_GlobalNodes_Domain_E+(Nx_last+1)*TotalNumberCells_y;
					if (dims[1]-2==coord[1])
						Node2D_GlobalDomain[3]=start_GlobalNodes_Domain_N+Ny_last*(2+TotalNumberCells_x-3);
					else
						Node2D_GlobalDomain[3]=start_GlobalNodes_Domain_N+TotalNumberCells_y*(2+TotalNumberCells_x-3);
				}
				else if (dims[1]-2==coord[1])
				{
					Node2D_GlobalDomain[0]=start_GlobalNodes+Node2D_SubDomain[0]-1;
					Node2D_GlobalDomain[1]=start_GlobalNodes_Domain_E+(TotalNumberCells_y-1)*2+1;
					Node2D_GlobalDomain[2]=start_GlobalNodes_Domain_E+TotalNumberCells_x*TotalNumberCells_y;
					Node2D_GlobalDomain[3]=start_GlobalNodes_Domain_N+Ny_last*(2+TotalNumberCells_x-3);
				}
				else
				{
					Node2D_GlobalDomain[0]=start_GlobalNodes+Node2D_SubDomain[0]-1;
					Node2D_GlobalDomain[1]=start_GlobalNodes_Domain_E+(TotalNumberCells_y-1)*2+1;
					Node2D_GlobalDomain[2]=start_GlobalNodes_Domain_E+TotalNumberCells_x*TotalNumberCells_y;
					Node2D_GlobalDomain[3]=start_GlobalNodes_Domain_N+TotalNumberCells_y*(2+TotalNumberCells_x-3);
				}
				//Map Local numbering to Global numbering
				LocalToGlobalNode[Node2D_SubDomain[0]]=Node2D_GlobalDomain[0];
				LocalToGlobalNode[Node2D_SubDomain[1]]=-1;
				LocalToGlobalNode[Node2D_SubDomain[2]]=-1;
				LocalToGlobalNode[Node2D_SubDomain[3]]=-1;
			}
		}

		//Mark Nodes for parallel
		//if(coord[0]==dims[0]-1)
		//{
			IdNodeNE.push_back(Node2D_SubDomain[0]-1);
			IdNodeN.push_back(Node2D_SubDomain[0]-1);
		//}


		NewCell(Nodes,NewNodes,GhostNodes,Node2D_SubDomain,Node2D_GlobalDomain,x_,y_,FaceConnect,CellConnect);
		nbCells++;



		// Cell (TotalNumberCells_x-1:1,TotalNumberCells_y)
		for (int i=TotalNumberCells_x-1;i>1;i--)
		{
			x_tmp=x_tmp-dx;
			Nodes[0]=Interior;
			Nodes[1]=Interior;
			Nodes[2]=bC[2];//Wall;
			Nodes[3]=bC[2];//Wall;

			x_[0]=x_tmp;
			x_[1]=x_tmp+dx;
			x_[2]=x_tmp+dx;
			x_[3]=x_tmp;

			y_[0]=y_tmp;
			y_[1]=y_tmp;
			y_[2]=y_tmp+dy;
			y_[3]=y_tmp+dy;

			FaceConnect[0]=2;
			FaceConnect[1]=3;
			FaceConnect[2]=0;
			FaceConnect[3]=1;

			CellConnect[0]=CellConnect[0]-TotalNumberCells_y+1;
			CellArray[CellConnect[0]]->Set_Connect(2,0,nbCells); //Correct the connectivity for the cell in south side
			CellConnect[1]=nbCells-1;
			CellConnect[2]=nbCells;
			CellConnect[3]=nbCells+1;

			NewNodes[0]=false;
			NewNodes[1]=false;
			NewNodes[2]=false;
			NewNodes[3]=true;

			Node2D_SubDomain[0]=CellArray[CellConnect[0]]->Get_NodeNumber(3);
			Node2D_SubDomain[1]=CellArray[CellConnect[0]]->Get_NodeNumber(2);
			Node2D_SubDomain[2]=CellArray[CellConnect[1]]->Get_NodeNumber(3);
			Node2D_SubDomain[3]=nNode2D_SubDomain++;

			if(bC[1]!=Interior)
			{
				if(bC[2]!=Interior)
				{
					GhostNodes[0]=false;
					GhostNodes[1]=false;
					GhostNodes[2]=false;
					GhostNodes[3]=false;

					Node2D_GlobalDomain[0]=start_GlobalNodes+Node2D_SubDomain[0]-1;
					Node2D_GlobalDomain[1]=start_GlobalNodes+Node2D_SubDomain[1]-1;
					Node2D_GlobalDomain[2]=start_GlobalNodes_Domain_E+TotalNumberCells_y+(TotalNumberCells_x-i);
					Node2D_GlobalDomain[3]=start_GlobalNodes_Domain_E+TotalNumberCells_y+(TotalNumberCells_x-i)+1;

					//Map Local numbering to Global numbering
					for (int i=0;i<4;i++)
					LocalToGlobalNode[Node2D_SubDomain[i]]=Node2D_GlobalDomain[i];
				}
				else
				{
					GhostNodes[0]=false;
					GhostNodes[1]=false;
					GhostNodes[2]=true;
					GhostNodes[3]=true;

					//Mark cell with ghost node to correct the node numbering before connectivity
					IdCellGhostNode.push_back(nbCells);

					Node2D_GlobalDomain[0]=start_GlobalNodes+Node2D_SubDomain[0]-1;
					Node2D_GlobalDomain[1]=start_GlobalNodes+Node2D_SubDomain[1]-1;
					if (dims[1]-2==coord[1])
					{
						if(i>2)
						{
							Node2D_GlobalDomain[2]=start_GlobalNodes_Domain_N+Ny_last*(2+i-2);
							Node2D_GlobalDomain[3]=start_GlobalNodes_Domain_N+Ny_last*(2+i-3);
						}
						else
						{
							Node2D_GlobalDomain[2]=start_GlobalNodes_Domain_N+(Ny_last-1)*2+2;
							Node2D_GlobalDomain[3]=start_GlobalNodes_Domain_N+1;
						}
					}
					else
					{
						if(i>2)
						{
							Node2D_GlobalDomain[2]=start_GlobalNodes_Domain_N+TotalNumberCells_y*(2+i-2);
							Node2D_GlobalDomain[3]=start_GlobalNodes_Domain_N+TotalNumberCells_y*(2+i-3);
						}
						else
						{
							Node2D_GlobalDomain[2]=start_GlobalNodes_Domain_N+(Ny_begin-1)*2+2;
							Node2D_GlobalDomain[3]=start_GlobalNodes_Domain_N+1;
						}
					}
					//Map Local numbering to Global numbering
					LocalToGlobalNode[Node2D_SubDomain[0]]=Node2D_GlobalDomain[0];
					LocalToGlobalNode[Node2D_SubDomain[1]]=Node2D_GlobalDomain[1];
					LocalToGlobalNode[Node2D_SubDomain[2]]=-1;
					LocalToGlobalNode[Node2D_SubDomain[3]]=-1;
				}
			}
			else
			{
				Node2D_GlobalDomain[0]=start_GlobalNodes+Node2D_SubDomain[0]-1;
				Node2D_GlobalDomain[1]=start_GlobalNodes+Node2D_SubDomain[1]-1;

				if(bC[2]!=Interior)
				{
					GhostNodes[0]=false;
					GhostNodes[1]=false;
					GhostNodes[2]=false;
					GhostNodes[3]=false;

					Node2D_GlobalDomain[2]=start_GlobalNodes+TotalNumberCells_x*TotalNumberCells_y+(TotalNumberCells_x-1-i);
					Node2D_GlobalDomain[3]=start_GlobalNodes+TotalNumberCells_x*TotalNumberCells_y+(TotalNumberCells_x-1-i)+1;

					//Map Local numbering to Global numbering
					for (int i=0;i<4;i++)
					LocalToGlobalNode[Node2D_SubDomain[i]]=Node2D_GlobalDomain[i];
				}
				else
				{
					GhostNodes[0]=false;
					GhostNodes[1]=false;
					GhostNodes[2]=true;
					GhostNodes[3]=true;

					//Mark cell with ghost node to correct the node numbering before connectivity
					IdCellGhostNode.push_back(nbCells);

					if (dims[1]-2==coord[1])
					{
						if(i>2)
						{
							Node2D_GlobalDomain[2]=start_GlobalNodes_Domain_N+Ny_last*(2+i-2);
							Node2D_GlobalDomain[3]=start_GlobalNodes_Domain_N+Ny_last*(2+i-3);
						}
						else
						{
							Node2D_GlobalDomain[2]=start_GlobalNodes_Domain_N+(Ny_last-1)*2+2;
							Node2D_GlobalDomain[3]=start_GlobalNodes_Domain_N+1;
						}
					}
					else
					{
						if(i>2)
						{
							Node2D_GlobalDomain[2]=start_GlobalNodes_Domain_N+Ny_begin*(2+i-2);
							Node2D_GlobalDomain[3]=start_GlobalNodes_Domain_N+Ny_begin*(2+i-3);
						}
						else
						{
							Node2D_GlobalDomain[2]=start_GlobalNodes_Domain_N+(Ny_begin-1)*2+2;
							Node2D_GlobalDomain[3]=start_GlobalNodes_Domain_N+1;
						}
					}
					//Map Local numbering to Global numbering
					LocalToGlobalNode[Node2D_SubDomain[0]]=Node2D_GlobalDomain[0];
					LocalToGlobalNode[Node2D_SubDomain[1]]=Node2D_GlobalDomain[1];
					LocalToGlobalNode[Node2D_SubDomain[2]]=-1;
					LocalToGlobalNode[Node2D_SubDomain[3]]=-1;
				}
			}

			//Mark Nodes for parallel
			/*if(coord[0]!=dims[0]-1 && i==TotalNumberCells_x-1	)
			{
				IdNodeNE.push_back(Node2D_GlobalDomain[2]-start_GlobalNodes);
				IdNodeN.push_back(Node2D_GlobalDomain[3]-start_GlobalNodes);

			}
			else
				IdNodeN.push_back(Node2D_GlobalDomain[3]-start_GlobalNodes);*/
			IdNodeN.push_back(Node2D_SubDomain[0]-1);

			NewCell(Nodes,NewNodes,GhostNodes,Node2D_SubDomain,Node2D_GlobalDomain,x_,y_,FaceConnect,CellConnect);
			nbCells++;

		}


// Cell (0,TotalNumberCells_y)
		x_tmp=x_tmp-dx;
		Nodes[0]=bC[1];
		Nodes[1]=Interior;
		Nodes[2]=bC[2];
		if (bC[3]==Interior)
			Nodes[3]=bC[2];
		else
			Nodes[3]=GlobalCorner;

		x_[0]=x_tmp;
		x_[1]=x_tmp+dx;
		x_[2]=x_tmp+dx;
		x_[3]=x_tmp;

		y_[0]=y_tmp;
		y_[1]=y_tmp;
		y_[2]=y_tmp+dy;
		y_[3]=y_tmp+dy;

		FaceConnect[0]=2;
		FaceConnect[1]=3;
		FaceConnect[2]=0;
		FaceConnect[3]=1;

		CellConnect[0]=CellConnect[0]-TotalNumberCells_y+1;
		CellArray[CellConnect[0]]->Set_Connect(2,0,nbCells); //Correct the connectivity for the cell in south side
		CellConnect[1]=nbCells-1;
		CellConnect[2]=nbCells;
		CellConnect[3]=nbCells;

		NewNodes[0]=false;
		NewNodes[1]=false;
		NewNodes[2]=false;
		NewNodes[3]=true;

		Node2D_SubDomain[0]=CellArray[CellConnect[0]]->Get_NodeNumber(3);
		Node2D_SubDomain[1]=CellArray[CellConnect[0]]->Get_NodeNumber(2);
		Node2D_SubDomain[2]=CellArray[CellConnect[1]]->Get_NodeNumber(3);
		Node2D_SubDomain[3]=nNode2D_SubDomain++;

		if(bC[1]!=Interior)
		{
			if(bC[2]!=Interior)
			{

				GhostNodes[0]=false;
				GhostNodes[1]=false;
				GhostNodes[2]=false;
				GhostNodes[3]=false;

				Node2D_GlobalDomain[0]=start_GlobalNodes+Node2D_SubDomain[0]-1;
				Node2D_GlobalDomain[1]=start_GlobalNodes+Node2D_SubDomain[1]-1;
				Node2D_GlobalDomain[2]=start_GlobalNodes+nComputationalNode+1+TotalNumberCells_y+(TotalNumberCells_x-2);
				Node2D_GlobalDomain[3]=start_GlobalNodes+nComputationalNode+1+TotalNumberCells_y+(TotalNumberCells_x-1);

				//Map Local numbering to Global numbering
				for (int i=0;i<4;i++)
				LocalToGlobalNode[Node2D_SubDomain[i]]=Node2D_GlobalDomain[i];
			}
			else
			{

				GhostNodes[0]=false;
				GhostNodes[1]=false;
				GhostNodes[2]=true;
				GhostNodes[3]=true;

				//Mark cell with ghost node to correct the node numbering before connectivity
				IdCellGhostNode.push_back(nbCells);

				Node2D_GlobalDomain[0]=start_GlobalNodes+Node2D_SubDomain[0]-1;
				Node2D_GlobalDomain[1]=start_GlobalNodes+Node2D_SubDomain[1]-1;
				Node2D_GlobalDomain[2]=start_GlobalNodes+nComputationalNode+TotalNumberCells_y+1;
				Node2D_GlobalDomain[3]=start_GlobalNodes+nComputationalNode+TotalNumberCells_y;

				//Map Local numbering to Global numbering
				LocalToGlobalNode[Node2D_SubDomain[0]]=Node2D_GlobalDomain[0];
				LocalToGlobalNode[Node2D_SubDomain[1]]=Node2D_GlobalDomain[1];
				LocalToGlobalNode[Node2D_SubDomain[2]]=-1;
				LocalToGlobalNode[Node2D_SubDomain[3]]=-1;
			}
		}
		else
		{
			Node2D_GlobalDomain[0]=start_GlobalNodes+Node2D_SubDomain[0]-1;
			Node2D_GlobalDomain[1]=start_GlobalNodes+Node2D_SubDomain[1]-1;
			if(bC[2]!=Interior)
			{

				GhostNodes[0]=false;
				GhostNodes[1]=false;
				GhostNodes[2]=false;
				GhostNodes[3]=false;

				Node2D_GlobalDomain[2]=start_GlobalNodes+nComputationalNode+(TotalNumberCells_x-2);
				Node2D_GlobalDomain[3]=start_GlobalNodes+nComputationalNode+(TotalNumberCells_x-1);

				//Map Local numbering to Global numbering
				for (int i=0;i<4;i++)
				LocalToGlobalNode[Node2D_SubDomain[i]]=Node2D_GlobalDomain[i];
			}
			else
			{

				GhostNodes[0]=false;
				GhostNodes[1]=false;
				GhostNodes[2]=true;
				GhostNodes[3]=true;

				//Mark cell with ghost node to correct the node numbering before connectivity
				IdCellGhostNode.push_back(nbCells);

				Node2D_GlobalDomain[2]=start_GlobalNodes+nComputationalNode+1;
				Node2D_GlobalDomain[3]=start_GlobalNodes+nComputationalNode;

				//Map Local numbering to Global numbering
				LocalToGlobalNode[Node2D_SubDomain[0]]=Node2D_GlobalDomain[0];
				LocalToGlobalNode[Node2D_SubDomain[1]]=Node2D_GlobalDomain[1];
				LocalToGlobalNode[Node2D_SubDomain[2]]=-1;
				LocalToGlobalNode[Node2D_SubDomain[3]]=-1;
			}
		}

		//Mark Nodes for parallel
		//IdNodeNW.push_back(Node2D_GlobalDomain[3]-start_GlobalNodes);
		IdNodeNW.push_back(Node2D_SubDomain[0]-1);

		NewCell(Nodes,NewNodes,GhostNodes,Node2D_SubDomain,Node2D_GlobalDomain,x_,y_,FaceConnect,CellConnect);
		nbCells++;
//		std::cout<<"Compare size "<<nNode2D_SubDomain<<" "<<Node.size()+Node_Ghosttmp.size()<<std::endl;





/// Add Ghost nodes at the end
		NbRealNodes=Node.size();
		NbGhostNode=x_Ghosttmp.size();
		if(verbous)
			std::cout<<"Number of computation cells is: "<<NbRealNodes<<std::endl;
		Node.insert(Node.end(),Node_Ghosttmp.begin(),Node_Ghosttmp.end());
		x.insert(x.end(),x_Ghosttmp.begin(),x_Ghosttmp.end());
		y.insert(y.end(),y_Ghosttmp.begin(),y_Ghosttmp.end());
		Node_Ghosttmp.clear();
		x_Ghosttmp.clear();
		y_Ghosttmp.clear();



/// Set Connections at each node
		Block2D::Correct_OrderingGhostNode();
///Add Ghost Cell for Parallel and periodic boundary conditions
		AddGhostCells(TotalNumberCells_x,TotalNumberCells_y);

		Block2D::Set_Connect();
		Correct_MarkNode();
		Block2D::Set_CommNodes();
		Clear_MarkNode();
		if(verbous)
			Block2D::Check_ID();
		NbTotalNodes=Node.size();
		//WriteCells(); WriteNodes();


	}
	else
	{
		std::cout<<"Error: a block has already been created"<<endl;

	}

}
void Block2D::AddGhostCells(int TotalNumberCells_x,int TotalNumberCells_y)
{
	int x_[4],y_[4], FaceConnect[4], CellConnect[4],Node2D_SubDomain[4];
	int nbCells=CellArray.size();
	int x_tmp=Node[CellArray[nbCells-1]->Get_NodeNumber(0)-1]->get_x()-dx;
	int y_tmp=Node[CellArray[nbCells-1]->Get_NodeNumber(0)-1]->get_y();
	NodeType Nodes[4];
	int GhostCellConnectToReal;

	bool NewNodes[4];
	int nNode2D_SubDomain=Node.size()+Node_Ghosttmp.size()+1;


// Cell (-1,TotalNumberCells_y)
	GhostCellConnectToReal=nbCells-1;
	Nodes[0]=Ghost;
	Nodes[1]=Node[CellArray[GhostCellConnectToReal]->Get_NodeNumber(0)-1]->get_NodeType();
	Nodes[2]=Node[CellArray[GhostCellConnectToReal]->Get_NodeNumber(3)-1]->get_NodeType();
	Nodes[3]=Ghost;

	x_[0]=x_tmp;
	x_[1]=x_tmp+dx;
	x_[2]=x_tmp+dx;
	x_[3]=x_tmp;

	y_[0]=y_tmp;
	y_[1]=y_tmp;
	y_[2]=y_tmp+dy;
	y_[3]=y_tmp+dy;

	FaceConnect[0]=2;
	FaceConnect[1]=3;
	FaceConnect[2]=0;
	FaceConnect[3]=1;

	CellConnect[0]=nbCells+1;
	CellConnect[1]=GhostCellConnectToReal;
	CellConnect[2]=nbCells;
	CellConnect[3]=nbCells;

	NewNodes[0]=true;
	NewNodes[1]=false;
	NewNodes[2]=false;
	NewNodes[3]=true;

	Node2D_SubDomain[0]=nNode2D_SubDomain++;
	Node2D_SubDomain[1]=CellArray[CellConnect[1]]->Get_NodeNumber(0);
	Node2D_SubDomain[2]=CellArray[CellConnect[1]]->Get_NodeNumber(3);
	Node2D_SubDomain[3]=nNode2D_SubDomain++;

	//Mark cell with ghost node to correct the node numbering before connectivity
	IdGhostCell.push_back(nbCells);


	NewGhostCell(Nodes,NewNodes,Node2D_SubDomain,x_,y_,FaceConnect,CellConnect);
	//Correct connectivity of neighbour cells
	CellArray[CellConnect[1]]->Set_Connect(3,1,nbCells);

	nbCells++;
	y_tmp=y_tmp-dy;

// Cell (-1,TotalNumberCells_y-2:0)
	for (int i=TotalNumberCells_y-1;i>0;i--)
	{
		GhostCellConnectToReal=CellArray[CellArray[nbCells-1]->Get_Connect(1)[1]]->Get_Connect(0)[1];
		Nodes[0]=Ghost;
		Nodes[1]=Node[CellArray[GhostCellConnectToReal]->Get_NodeNumber(0)-1]->get_NodeType();
		Nodes[2]=Node[CellArray[GhostCellConnectToReal]->Get_NodeNumber(3)-1]->get_NodeType();
		Nodes[3]=Ghost;

		x_[0]=x_tmp;
		x_[1]=x_tmp+dx;
		x_[2]=x_tmp+dx;
		x_[3]=x_tmp;

		y_[0]=y_tmp;
		y_[1]=y_tmp;
		y_[2]=y_tmp+dy;
		y_[3]=y_tmp+dy;

		FaceConnect[0]=2;
		FaceConnect[1]=3;
		FaceConnect[2]=0;
		FaceConnect[3]=1;

		CellConnect[0]=nbCells+1;
		CellConnect[1]=GhostCellConnectToReal;
		CellConnect[2]=nbCells-1;
		CellConnect[3]=nbCells;

		NewNodes[0]=true;
		NewNodes[1]=false;
		NewNodes[2]=false;
		NewNodes[3]=false;

		Node2D_SubDomain[0]=nNode2D_SubDomain++;
		Node2D_SubDomain[1]=CellArray[CellConnect[1]]->Get_NodeNumber(0);
		Node2D_SubDomain[2]=CellArray[CellConnect[1]]->Get_NodeNumber(3);
		Node2D_SubDomain[3]=CellArray[CellConnect[2]]->Get_NodeNumber(0);

		//Mark cell with ghost node to correct the node numbering before connectivity
		IdGhostCell.push_back(nbCells);


		NewGhostCell(Nodes,NewNodes,Node2D_SubDomain,x_,y_,FaceConnect,CellConnect);
		//Correct connectivity of neighbour cells
		CellArray[CellConnect[1]]->Set_Connect(3,1,nbCells);

		nbCells++;
		y_tmp=y_tmp-dy;
	}
// Cell (-1,-1)
	GhostCellConnectToReal=nbCells-1;
	Nodes[0]=Ghost;
	Nodes[1]=Ghost;
	Nodes[2]=Node[0]->get_NodeType();
	Nodes[3]=Ghost;

	x_[0]=x_tmp;
	x_[1]=x_tmp+dx;
	x_[2]=x_tmp+dx;
	x_[3]=x_tmp;

	y_[0]=y_tmp;
	y_[1]=y_tmp;
	y_[2]=y_tmp+dy;
	y_[3]=y_tmp+dy;

	FaceConnect[0]=2;
	FaceConnect[1]=3;
	FaceConnect[2]=0;
	FaceConnect[3]=1;

	CellConnect[0]=nbCells;
	CellConnect[1]=nbCells+1;
	CellConnect[2]=nbCells-1;
	CellConnect[3]=nbCells;

	NewNodes[0]=true;
	NewNodes[1]=true;
	NewNodes[2]=false;
	NewNodes[3]=false;

	Node2D_SubDomain[0]=nNode2D_SubDomain++;
	Node2D_SubDomain[1]=nNode2D_SubDomain++;
	Node2D_SubDomain[2]=CellArray[CellConnect[2]]->Get_NodeNumber(1);
	Node2D_SubDomain[3]=CellArray[CellConnect[2]]->Get_NodeNumber(0);

	//Mark cell with ghost node to correct the node numbering before connectivity
	IdGhostCell.push_back(nbCells);

	NewGhostCell(Nodes,NewNodes,Node2D_SubDomain,x_,y_,FaceConnect,CellConnect);

	nbCells++;
	x_tmp=x_tmp+dx;
// Cell (0:TotalNumberCells_x,-1)
	for (int i=0;i<TotalNumberCells_x;i++)
	{
		GhostCellConnectToReal=CellArray[CellArray[nbCells-1]->Get_Connect(2)[1]]->Get_Connect(1)[1];
		Nodes[0]=Ghost;
		Nodes[1]=Ghost;
		Nodes[2]=Node[CellArray[GhostCellConnectToReal]->Get_NodeNumber(1)-1]->get_NodeType();
		Nodes[3]=Node[CellArray[GhostCellConnectToReal]->Get_NodeNumber(0)-1]->get_NodeType();

		x_[0]=x_tmp;
		x_[1]=x_tmp+dx;
		x_[2]=x_tmp+dx;
		x_[3]=x_tmp;

		y_[0]=y_tmp;
		y_[1]=y_tmp;
		y_[2]=y_tmp+dy;
		y_[3]=y_tmp+dy;

		FaceConnect[0]=2;
		FaceConnect[1]=3;
		FaceConnect[2]=0;
		FaceConnect[3]=1;

		CellConnect[0]=nbCells;
		CellConnect[1]=nbCells+1;
		CellConnect[2]=GhostCellConnectToReal;
		CellConnect[3]=nbCells-1;

		NewNodes[0]=false;
		NewNodes[1]=true;
		NewNodes[2]=false;
		NewNodes[3]=false;

		Node2D_SubDomain[0]=CellArray[CellConnect[3]]->Get_NodeNumber(1);
		Node2D_SubDomain[1]=nNode2D_SubDomain++;
		Node2D_SubDomain[2]=CellArray[CellConnect[2]]->Get_NodeNumber(1);
		Node2D_SubDomain[3]=CellArray[CellConnect[2]]->Get_NodeNumber(0);

		//Mark cell with ghost node to correct the node numbering before connectivity
		IdGhostCell.push_back(nbCells);


		NewGhostCell(Nodes,NewNodes,Node2D_SubDomain,x_,y_,FaceConnect,CellConnect);
		//Correct connectivity of neighbour cells
		CellArray[CellConnect[2]]->Set_Connect(0,2,nbCells);

		nbCells++;
		x_tmp=x_tmp+dx;
	}
// Cell (TotalNumberCells_x+1,-1)
	GhostCellConnectToReal=CellArray[nbCells-1]->Get_Connect(2)[1];
	Nodes[0]=Ghost;
	Nodes[1]=Ghost;
	Nodes[2]=Ghost;
	Nodes[3]=Node[CellArray[GhostCellConnectToReal]->Get_NodeNumber(1)-1]->get_NodeType();

	x_[0]=x_tmp;
	x_[1]=x_tmp+dx;
	x_[2]=x_tmp+dx;
	x_[3]=x_tmp;

	y_[0]=y_tmp;
	y_[1]=y_tmp;
	y_[2]=y_tmp+dy;
	y_[3]=y_tmp+dy;

	FaceConnect[0]=2;
	FaceConnect[1]=3;
	FaceConnect[2]=0;
	FaceConnect[3]=1;

	CellConnect[0]=nbCells;
	CellConnect[1]=nbCells;
	CellConnect[2]=nbCells+1;
	CellConnect[3]=nbCells-1;

	NewNodes[0]=false;
	NewNodes[1]=true;
	NewNodes[2]=true;
	NewNodes[3]=false;

	Node2D_SubDomain[0]=CellArray[CellConnect[3]]->Get_NodeNumber(1);
	Node2D_SubDomain[1]=nNode2D_SubDomain++;
	Node2D_SubDomain[2]=nNode2D_SubDomain++;
	Node2D_SubDomain[3]=CellArray[CellConnect[3]]->Get_NodeNumber(2);

	//Mark cell with ghost node to correct the node numbering before connectivity
	IdGhostCell.push_back(nbCells);


	NewGhostCell(Nodes,NewNodes,Node2D_SubDomain,x_,y_,FaceConnect,CellConnect);
	//Correct connectivity of neighbour cells
	//CellArray[CellConnect[1]]->Set_Connect(3,1,nbCells);

	nbCells++;
	y_tmp=y_tmp+dy;
// Cell (TotalNumberCells_x+1,0:TotalNumberCells_y)
	for (int i=0;i<TotalNumberCells_y;i++)
	{
		GhostCellConnectToReal=CellArray[CellArray[nbCells-1]->Get_Connect(3)[1]]->Get_Connect(2)[1];
		Nodes[0]=Node[CellArray[GhostCellConnectToReal]->Get_NodeNumber(1)-1]->get_NodeType();
		Nodes[1]=Ghost;
		Nodes[2]=Ghost;
		Nodes[3]=Node[CellArray[GhostCellConnectToReal]->Get_NodeNumber(2)-1]->get_NodeType();

		x_[0]=x_tmp;
		x_[1]=x_tmp+dx;
		x_[2]=x_tmp+dx;
		x_[3]=x_tmp;

		y_[0]=y_tmp;
		y_[1]=y_tmp;
		y_[2]=y_tmp+dy;
		y_[3]=y_tmp+dy;

		FaceConnect[0]=2;
		FaceConnect[1]=3;
		FaceConnect[2]=0;
		FaceConnect[3]=1;

		CellConnect[0]=nbCells-1;
		CellConnect[1]=nbCells;
		CellConnect[2]=nbCells+1;
		CellConnect[3]=GhostCellConnectToReal;

		NewNodes[0]=false;
		NewNodes[1]=false;
		NewNodes[2]=true;
		NewNodes[3]=false;

		Node2D_SubDomain[0]=CellArray[CellConnect[3]]->Get_NodeNumber(1);
		Node2D_SubDomain[1]=CellArray[CellConnect[0]]->Get_NodeNumber(2);
		Node2D_SubDomain[2]=nNode2D_SubDomain++;
		Node2D_SubDomain[3]=CellArray[CellConnect[3]]->Get_NodeNumber(2);

		//Mark cell with ghost node to correct the node numbering before connectivity
		IdGhostCell.push_back(nbCells);

		NewGhostCell(Nodes,NewNodes,Node2D_SubDomain,x_,y_,FaceConnect,CellConnect);
		//Correct connectivity of neighbour cells
		CellArray[CellConnect[3]]->Set_Connect(1,3,nbCells);

		nbCells++;
		y_tmp=y_tmp+dy;
	}
// Cell (TotalNumberCells_x+1,TotalNumberCells_y+1)
	GhostCellConnectToReal=CellArray[CellArray[nbCells-1]->Get_Connect(3)[1]]->Get_Connect(2)[1];
	Nodes[0]=Node[CellArray[nbCells-1]->Get_NodeNumber(3)-1]->get_NodeType();
	Nodes[1]=Ghost;
	Nodes[2]=Ghost;
	Nodes[3]=Ghost;

	x_[0]=x_tmp;
	x_[1]=x_tmp+dx;
	x_[2]=x_tmp+dx;
	x_[3]=x_tmp;

	y_[0]=y_tmp;
	y_[1]=y_tmp;
	y_[2]=y_tmp+dy;
	y_[3]=y_tmp+dy;

	FaceConnect[0]=2;
	FaceConnect[1]=3;
	FaceConnect[2]=0;
	FaceConnect[3]=1;

	CellConnect[0]=nbCells-1;
	CellConnect[1]=nbCells;
	CellConnect[2]=nbCells;
	CellConnect[3]=nbCells+1;

	NewNodes[0]=false;
	NewNodes[1]=false;
	NewNodes[2]=true;
	NewNodes[3]=true;

	Node2D_SubDomain[0]=CellArray[CellConnect[0]]->Get_NodeNumber(3);
	Node2D_SubDomain[1]=CellArray[CellConnect[0]]->Get_NodeNumber(2);
	Node2D_SubDomain[2]=nNode2D_SubDomain++;
	Node2D_SubDomain[3]=nNode2D_SubDomain++;

	//Mark cell with ghost node to correct the node numbering before connectivity
	IdGhostCell.push_back(nbCells);

	NewGhostCell(Nodes,NewNodes,Node2D_SubDomain,x_,y_,FaceConnect,CellConnect);
	//Correct connectivity of neighbour cells
	//CellArray[CellConnect[3]]->Set_Connect(1,3,nbCells);

	nbCells++;
	x_tmp=x_tmp-dx;

// Cell (0:TotalNumberCells_x+1,TotalNumberCells_y+1)
	for (int i=TotalNumberCells_x;i>0;i--)
	{
		GhostCellConnectToReal=CellArray[CellArray[nbCells-1]->Get_Connect(0)[1]]->Get_Connect(3)[1];
		Nodes[0]=Node[CellArray[GhostCellConnectToReal]->Get_NodeNumber(3)-1]->get_NodeType();
		Nodes[1]=Node[CellArray[GhostCellConnectToReal]->Get_NodeNumber(2)-1]->get_NodeType();
		Nodes[2]=Ghost;
		Nodes[3]=Ghost;

		x_[0]=x_tmp;
		x_[1]=x_tmp+dx;
		x_[2]=x_tmp+dx;
		x_[3]=x_tmp;

		y_[0]=y_tmp;
		y_[1]=y_tmp;
		y_[2]=y_tmp+dy;
		y_[3]=y_tmp+dy;

		FaceConnect[0]=2;
		FaceConnect[1]=3;
		FaceConnect[2]=0;
		FaceConnect[3]=1;

		CellConnect[0]=GhostCellConnectToReal;
		CellConnect[1]=nbCells-1;
		CellConnect[2]=nbCells;
		CellConnect[3]=nbCells+1;

		NewNodes[0]=false;
		NewNodes[1]=false;
		NewNodes[2]=false;
		NewNodes[3]=true;

		Node2D_SubDomain[0]=CellArray[CellConnect[0]]->Get_NodeNumber(3);
		Node2D_SubDomain[1]=CellArray[CellConnect[0]]->Get_NodeNumber(2);
		Node2D_SubDomain[2]=CellArray[CellConnect[1]]->Get_NodeNumber(3);
		Node2D_SubDomain[3]=nNode2D_SubDomain++;

		//Mark cell with ghost node to correct the node numbering before connectivity
		IdGhostCell.push_back(nbCells);


		NewGhostCell(Nodes,NewNodes,Node2D_SubDomain,x_,y_,FaceConnect,CellConnect);
		//Correct connectivity of neighbour cells
		CellArray[CellConnect[0]]->Set_Connect(2,0,nbCells);

		nbCells++;
		x_tmp=x_tmp-dx;
	}
// Cell (-1,TotalNumberCells_y+1)
	GhostCellConnectToReal=CellArray[CellArray[nbCells-1]->Get_Connect(0)[1]]->Get_Connect(3)[1];
	Nodes[0]=Ghost;
	Nodes[1]=Node[CellArray[nbCells-1]->Get_NodeNumber(0)-1]->get_NodeType();
	Nodes[2]=Ghost;
	Nodes[3]=Ghost;

	x_[0]=x_tmp;
	x_[1]=x_tmp+dx;
	x_[2]=x_tmp+dx;
	x_[3]=x_tmp;

	y_[0]=y_tmp;
	y_[1]=y_tmp;
	y_[2]=y_tmp+dy;
	y_[3]=y_tmp+dy;

	FaceConnect[0]=2;
	FaceConnect[1]=3;
	FaceConnect[2]=0;
	FaceConnect[3]=1;

	CellConnect[0]=GhostCellConnectToReal;
	CellConnect[1]=nbCells-1;
	CellConnect[2]=nbCells;
	CellConnect[3]=nbCells;

	NewNodes[0]=false;
	NewNodes[1]=false;
	NewNodes[2]=false;
	NewNodes[3]=true;

	Node2D_SubDomain[0]=CellArray[CellConnect[0]]->Get_NodeNumber(3);
	Node2D_SubDomain[1]=CellArray[CellConnect[0]]->Get_NodeNumber(2);
	Node2D_SubDomain[2]=CellArray[CellConnect[1]]->Get_NodeNumber(3);
	Node2D_SubDomain[3]=nNode2D_SubDomain++;

	//Mark cell with ghost node to correct the node numbering before connectivity
	IdGhostCell.push_back(nbCells);

	NewGhostCell(Nodes,NewNodes,Node2D_SubDomain,x_,y_,FaceConnect,CellConnect);
	//Correct connectivity of neighbour cells
	CellArray[CellConnect[0]]->Set_Connect(2,0,nbCells);



	nbCells++;
}
void Block2D::WriteCells()
{
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

	char buffer[50]; // make sure it's big enough
	snprintf(buffer, sizeof(buffer), "Cells_%d.txt", rank);
	std::ofstream myFlux;
	myFlux.open(buffer);
	int it;
	for (it=0; it<CellArray.size(); it++)
	{
	myFlux<<"Cell: "<<it<<std::endl;
	myFlux<<"Nodes: "<<CellArray[it]->Get_NodeNumber(0)<<" "<<CellArray[it]->Get_NodeNumber(1)<<" "<<CellArray[it]->Get_NodeNumber(2)<<" "<<CellArray[it]->Get_NodeNumber(3)<<std::endl;
	myFlux<<"Connections: "<<CellArray[it]->Get_Connect(0)[1]<<" "<<CellArray[it]->Get_Connect(1)[1]<<" "<<CellArray[it]->Get_Connect(2)[1]<<" "<<CellArray[it]->Get_Connect(3)[1]<<std::endl;
	myFlux<<std::endl;
	}
}
void Block2D::WriteNodes()
{
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

	char buffer[50]; // make sure it's big enough
	snprintf(buffer, sizeof(buffer), "Nodes_%d.txt", rank);
	std::ofstream myFlux;
	myFlux.open(buffer);
	int it;
	for (it=0; it<Node.size(); it++)
	{
	myFlux<<"Node: "<<it<<std::endl;
	myFlux<<"X Y: "<<Node[it]->get_x()<<" "<<Node[it]->get_y()<<std::endl;
	myFlux<<"Connection: "<<Node[it]->Get_connect(0)<<" "<<Node[it]->Get_connect(1)<<" "<<Node[it]->Get_connect(2)<<" "<<Node[it]->Get_connect(3)<<std::endl;
	myFlux<<std::endl;
	}
}
void Block2D::WriteCoord()
{
/*    string const FileName("Coordinates.vtk");
    std::ofstream myFlux(FileName.c_str());

    if(myFlux)
    {
    	myFlux << "# vtk DataFile Version 3.0" << endl << "vtk output" << endl << "ASCII" << endl << "DATASET UNSTRUCTURED_GRID" << endl;
    	myFlux << "POINTS "<<  Node.size() << " float"<< endl;

        for (int i=0;i<Node.size();i++)
        {
        	myFlux << Get_Coord(i)[0] << " " << Get_Coord(i)[1] << " " << 0 << endl;
        }

        myFlux << "CELLS "<<CellArray.size()<< " "<< 5*CellArray.size()<< endl;
		for (int i=0;i<CellArray.size();i++)
		{
			myFlux << 4 << " " <<  CellArray[i]->Get_NodeNumber(0)-1 << " " << CellArray[i]->Get_NodeNumber(1)-1 << " " << CellArray[i]->Get_NodeNumber(2)-1 << " " << CellArray[i]->Get_NodeNumber(3)-1 << endl;
		}

        myFlux << "CELL_TYPES "<<CellArray.size()<< endl;
        if(CellArray.size()>=1)
        {
			for (int i=0;i<CellArray.size()-1;i++)
			{
				myFlux << 9 << " " ;
			}
        }
		myFlux << 9 << endl ;

        myFlux << "POINT_DATA " << Node.size() << endl;
		myFlux << "SCALARS " << "XCoordinate" << " FLOAT" << endl;
		myFlux << "LOOKUP_TABLE default" << endl;
        for (int i=0;i<Node.size();i++)
        {
        	myFlux << Get_Coord(i)[0]<< endl;

        }
		myFlux << "SCALARS " << "YCoordinate" << " FLOAT" << endl;
		myFlux << "LOOKUP_TABLE default" << endl;
        for (int i=0;i<Node.size();i++)
        {
        	myFlux << Get_Coord(i)[1] << endl;
        }

    }
    else
    {
        cout << "ERREUR: Impossible d'ouvrir le fichier." << endl;
    }*/
}


void Block2D::ChangeNodeType(int NodeNumber, NodeType NewNodeType_)
{

	unsigned int x_tmp=(unsigned int)x[NodeNumber];
	unsigned int y_tmp=(unsigned int)y[NodeNumber];
	unsigned int Connect_N_tmp=Node[NodeNumber]->Get_connect(2)+1;
	unsigned int Connect_S_tmp=Node[NodeNumber]->Get_connect(0)+1;
	unsigned int Connect_W_tmp=Node[NodeNumber]->Get_connect(3)+1;
	unsigned int Connect_E_tmp=Node[NodeNumber]->Get_connect(1)+1;
	unsigned int NbVelocity_tmp=Node[NodeNumber]->Get_NbVelocity();
	NodeType OldNodeType_=Node[NodeNumber]->get_NodeType();

	delete Node[NodeNumber];

	switch(NewNodeType_)
	{
	case Interior:
		Node[NodeNumber]=new NodeInterior2D(x_tmp, y_tmp);
		break;
	case Solid:
		if(OldNodeType_==Ghost)
			Id_SolidGhost.push_back(NodeNumber);
			//Node[NodeNumber]=new NodeGhost2D(x_tmp, y_tmp,SolidGhost);
		//else
		Node[NodeNumber]=new NodeSolid2D(x_tmp, y_tmp);
		break;
	case Ghost:
		Node[NodeNumber]=new NodeGhost2D(x_tmp, y_tmp);
		break;
	case Velocity:
		Node[NodeNumber]=new NodeVelocity2D(x_tmp, y_tmp);
		break;
	case Pressure:
		Node[NodeNumber]=new NodePressure2D(x_tmp, y_tmp);
		break;
	case Periodic:
		Node[NodeNumber]=new NodePeriodic2D(x_tmp, y_tmp);
		break;
	case Corner:
		Node[NodeNumber]=new NodeCorner2D(x_tmp, y_tmp);
		break;
	case Wall:
		Node[NodeNumber]=new NodeWall2D(x_tmp, y_tmp);
		break;
	default:
		cerr<< "Node Type not found" << endl;
        break;
	}
	Node[NodeNumber]->Set_Connect(Connect_N_tmp,Connect_S_tmp,Connect_W_tmp,Connect_E_tmp);
	Node[NodeNumber]->Set_NbVelocity(NbVelocity_tmp);

}


void Block2D::ChangeCoord(int NodeNumber,unsigned int const x_, unsigned int const y_)
{
	x[NodeNumber]=x_ ;
	y[NodeNumber]=y_ ;
}


NodeType Block2D::Get_NodeType(int NodeNumber) const
{
	return Node[NodeNumber]->get_NodeType();
}

double* Block2D::Get_Coord(int NodeNumber) const
{
	Coord_physical[0]=x[NodeNumber];
	Coord_physical[1]=y[NodeNumber];
	return Coord_physical;
}

const double* Block2D::Get_X0() const
{
	const double* ptr=0;
	ptr=&x[0];
	return ptr;
}

const double* Block2D::Get_Y0() const
{
	const double* ptr=0;
	ptr=&y[0];
	return ptr;
}

int Block2D::Get_nnodes()const {
	return NbTotalNodes;
}

int* Block2D::Get_Elems0()
{
	return &Elems[0];
}

/*std::vector<Node2D*>* Block2D::Get_PtrNode() {
	return &Node;
}*/

Node2D* Block2D::Get_Node(int NodeNumber)
{
	return Node[NodeNumber];
}
std::vector<int>& Block2D::Get_PtrIdBc() {return IdBoundaries;}

void Block2D::Check_ID(){
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	char buffer[50]; // make sure it's big enough
	snprintf(buffer, sizeof(buffer), "check_connnect_%d.txt", rank);
	std::ofstream myFlux;
	myFlux.open(buffer);
	std::vector<int>::iterator it;
	myFlux<<"Processor number is: "<<rank<<" South Nodes are: ";
	for(it = IdNodeS.begin() ; it != IdNodeS.end() ; ++it)
		myFlux<<*it<<" ";
	myFlux<<std::endl;
	myFlux<<"Processor number is: "<<rank<<" East Nodes are: ";
	for(it = IdNodeE.begin() ; it != IdNodeE.end() ; ++it)
		myFlux<<*it<<" ";
	myFlux<<std::endl;
	myFlux<<"Processor number is: "<<rank<<" North Nodes are: ";
	for(it = IdNodeN.begin() ; it != IdNodeN.end() ; ++it)
		myFlux<<*it<<" ";
	myFlux<<std::endl;
	myFlux<<"Processor number is: "<<rank<<" West Nodes are: ";
	for(it = IdNodeW.begin() ; it != IdNodeW.end() ; ++it)
		myFlux<<*it<<" ";
	myFlux<<std::endl;
	myFlux<<"Processor number is: "<<rank<<" Corner South West Nodes are: ";
	for(it = IdNodeSW.begin() ; it != IdNodeSW.end() ; ++it)
		myFlux<<*it<<" ";
	myFlux<<std::endl;
	myFlux<<"Processor number is: "<<rank<<" Corner South East Nodes are: ";
	for(it = IdNodeSE.begin() ; it != IdNodeSE.end() ; ++it)
		myFlux<<*it<<" ";
	myFlux<<std::endl;
	myFlux<<"Processor number is: "<<rank<<" Corner North East Nodes are: ";
	for(it = IdNodeNE.begin() ; it != IdNodeNE.end() ; ++it)
		myFlux<<*it<<" ";
	myFlux<<std::endl;
	myFlux<<"Processor number is: "<<rank<<" Corner North West Nodes are: ";
	for(it = IdNodeNW.begin() ; it != IdNodeNW.end() ; ++it)
		myFlux<<*it<<" ";
	myFlux<<std::endl;

//	myFlux<<x[0]<<std::endl;

	myFlux<<"Processor number is: "<<rank<<" South Nodes coordinates are: ";
	for(it = IdNodeS.begin() ; it != IdNodeS.end() ; ++it)
		myFlux<<x[*it]<<" "<<y[*it]<<" ";
	myFlux<<std::endl;
	myFlux<<"Processor number is: "<<rank<<" East Nodes coordinates are: ";
	for(it = IdNodeE.begin() ; it != IdNodeE.end() ; ++it)
		myFlux<<x[*it]<<" "<<y[*it]<<" ";
	myFlux<<std::endl;
	myFlux<<"Processor number is: "<<rank<<" North Nodes coordinates are: ";
	for(it = IdNodeN.begin() ; it != IdNodeN.end() ; ++it)
		myFlux<<x[*it]<<" "<<y[*it]<<" ";
	myFlux<<std::endl;
	myFlux<<"Processor number is: "<<rank<<" West Nodes coordinates are: ";
	for(it = IdNodeW.begin() ; it != IdNodeW.end() ; ++it)
		myFlux<<x[*it]<<" "<<y[*it]<<" ";
	myFlux<<std::endl;
	myFlux<<"Processor number is: "<<rank<<" Corner South West Nodes coordinates are: ";
	for(it = IdNodeSW.begin() ; it != IdNodeSW.end() ; ++it)
		myFlux<<x[*it]<<" "<<y[*it]<<" ";
	myFlux<<std::endl;
	myFlux<<"Processor number is: "<<rank<<" Corner South East Nodes coordinates are: ";
	for(it = IdNodeSE.begin() ; it != IdNodeSE.end() ; ++it)
		myFlux<<x[*it]<<" "<<y[*it]<<" ";
	myFlux<<std::endl;
	myFlux<<"Processor number is: "<<rank<<" Corner North East Nodes coordinates are: ";
	for(it = IdNodeNE.begin() ; it != IdNodeNE.end() ; ++it)
		myFlux<<x[*it]<<" "<<y[*it]<<" ";
	myFlux<<std::endl;
	myFlux<<"Processor number is: "<<rank<<" Corner North West Nodes coordinates are: ";
	for(it = IdNodeNW.begin() ; it != IdNodeNW.end() ; ++it)
		myFlux<<x[*it]<<" "<<y[*it]<<" ";
	myFlux<<std::endl;
}

/// Function to put mark node on the borders
void Block2D::Correct_MarkNode(){

	std::vector<int>::iterator it;


	for(it = IdNodeE.begin() ; it != IdNodeE.end() ; ++it)
		*it=Node[*it]->Get_connect(1);

	for(it = IdNodeN.begin() ; it != IdNodeN.end() ; ++it)
		*it=Node[*it]->Get_connect(2);

	for(it = IdNodeSE.begin() ; it != IdNodeSE.end() ; ++it)
		*it=Node[*it]->Get_connect(1);

	for(it = IdNodeNE.begin() ; it != IdNodeNE.end() ; ++it)
		*it=Node[Node[*it]->Get_connect(1)]->Get_connect(2);

	for(it = IdNodeNW.begin() ; it != IdNodeNW.end() ; ++it)
		*it=Node[*it]->Get_connect(2);

}

void Block2D::Get_Connect_Node(std::vector<int> & IdNodeN_,std::vector<int> & IdNodeE_,std::vector<int> & IdNodeS_,std::vector<int> & IdNodeW_,
							  std::vector<int> & IdNodeSW_,std::vector<int> & IdNodeSE_,std::vector<int> & IdNodeNW_,std::vector<int> & IdNodeNE_){
	IdNodeN_=IdNodeN;
	IdNodeE_=IdNodeE;
	IdNodeS_=IdNodeS;
	IdNodeW_=IdNodeW;
	IdNodeSW_=IdNodeSW;
	IdNodeSE_=IdNodeSE;
	IdNodeNW_=IdNodeNW;
	IdNodeNE_=IdNodeNE;
}
void Block2D::Get_NbNodes(int & NbRealNodes_, int & NbTotalNodes_){
	NbRealNodes_=NbRealNodes;
	NbTotalNodes_=NbTotalNodes;
}
void Block2D::GenerateSolid(Parameters &Param, int &ndSolidNode,int &firstSolid){
	bool testSolid;
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

	//for (int i=0;i<NbRealNodes;i++)
	for (int i=0;i<Node.size();i++)
	{
		testSolid=false;//Fluid by default
		Node[i]->Set_Index(i);
// ******** Generation of Solid in the domain*********
		ChangeNode(*Node[i],testSolid);
		if(testSolid==true)
		{
			if(Param.Get_Verbous())
				std::cout<<"Node set as Solid: "<<i<< " in processor: "<<rank<<std::endl;
			//Save the first Solid node in the processor to sort array later on
			if(ndSolidNode==0)
				firstSolid=i;
			//Store the number of Solid node
			ndSolidNode++;
			//store id of solid nodes
			if(i<NbRealNodes)
				IdSolidNode.push_back(i);
			//create the solid node}
			Block2D::ChangeNodeType(i,Solid);//convert node to solid
			//Extend the solid to the ghost of the domain for the automatic wall detection functions. Otherwise, the wall can be consider as a corner
			if(Node[i]->get_y()==0)//Convert ghost node to solid at the border of the domain (bottom)
				Block2D::ChangeNodeType(Node[i]->Get_connect(0),Solid);//convert node to solid
			if(Node[i]->get_y()==Param.Get_Ny())//Convert ghost node to solid at the border of the domain (North)
				Block2D::ChangeNodeType(Node[i]->Get_connect(2),Solid);//convert node to solid
			if(Node[i]->get_x()==0)//Convert ghost node to solid at the border of the domain (West)
				Block2D::ChangeNodeType(Node[i]->Get_connect(3),Solid);//convert node to solid
			if(Node[i]->get_x()==Param.Get_Nx())//Convert ghost node to solid at the border of the domain (East)
				Block2D::ChangeNodeType(Node[i]->Get_connect(1),Solid);//convert node to solid
		}
		//Save nodes which are not solid. To reduce memory consumption, only nodes after the first solid are stored
		/*if(ndSolidNode>0 && Solid==false)
			Node_tmp.push_back(Node[i]);
		ndSolidNodetoremove.push_back(ndSolidNode);*/
		//}
	}
	//Correct_Solid_Ghost();
	std::cout<<"Processor: "<<rank<<" Number of Solid: "<<IdSolidNode.size()<<" Number of Computation nodes: "<<NbRealNodes-IdSolidNode.size()<<std::endl;
}
/*void Block2D::Correct_Solid_Ghost()
{
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

	char buffer[50]; // make sure it's big enough
	snprintf(buffer, sizeof(buffer), "solid_%d.txt", rank);
	std::ofstream myFlux;
	myFlux.open(buffer);

	int ndNotSolid=0;
	for(int i=0;i<Id_SolidGhost.size();i++)
	{
		ndNotSolid=0;
		for (int j=0;j<4;j++)
		{

			if(Node[Node[Id_SolidGhost[i]]->Get_connect(j)]->get_NodeType()!=Solid)
				ndNotSolid++;
		}
		if(ndNotSolid>0)
		{
			myFlux<<" Solid converted to Ghost is: "<<Id_SolidGhost[i]<<" x: "<<Node[Id_SolidGhost[i]]->get_x()<<" y: "<<Node[Id_SolidGhost[i]]->get_y()<<std::endl ;
			for (int j=0;j<4;j++)
				myFlux<<Node[Node[Id_SolidGhost[i]]->Get_connect(j)]->get_NodeType()<<" ";
			myFlux<<std::endl ;
		}
	}

}*/
void Block2D::DetectDirectInterior(int & nodeID, int & nbinterior)
{
	nbinterior=0;// detect corners
	for (unsigned int j=0;j<4;j++)
	{
		if(Node[Node[nodeID]->Get_connect(j)]->get_NodeType()==Interior
				||Node[Node[nodeID]->Get_connect(j)]->get_NodeType()==Ghost)
		{
			nbinterior++;
		}
	}
}
void Block2D::DetectSpecialWall(int & nodeID,int x,int y, bool & specialwall, unsigned int & directionwall)
{
	directionwall=-1;//Negative value to generate an error if it is not a special wall
	int nbspecialwall=0;// check if only one special wall is detected
	specialwall=false;
//	int nx,ny,x,y;
	if(x==0 && (Node[Node[nodeID]->Get_connect(3)]->get_NodeType()==Solid||Node[Node[nodeID]->Get_connect(3)]->get_NodeType()==Corner||Node[Node[nodeID]->Get_connect(3)]->get_NodeType()==ConcaveCorner||Node[Node[nodeID]->Get_connect(3)]->get_NodeType()==ConvexCorner||Node[Node[nodeID]->Get_connect(3)]->get_NodeType()==Wall))
	{
		if(Node[Node[nodeID]->Get_connect(0)]->get_NodeType()==Wall||Node[Node[nodeID]->Get_connect(2)]->get_NodeType()==Wall
			||Node[Node[nodeID]->Get_connect(0)]->get_NodeType()==SpecialWall||Node[Node[nodeID]->Get_connect(2)]->get_NodeType()==SpecialWall
			||(Node[Node[nodeID]->Get_connect(0)]->get_NodeType()==Solid&&Node[Node[nodeID]->Get_connect(2)]->get_NodeType()==Solid))//Corner
			specialwall=false;
		else
		{
			specialwall=true;
			if(Node[Node[nodeID]->Get_connect(0)]->get_NodeType()!=Solid)
				directionwall=4;
			else
				directionwall=2;
		}
	}
	else if(x==nx && (Node[Node[nodeID]->Get_connect(1)]->get_NodeType()==Solid||Node[Node[nodeID]->Get_connect(1)]->get_NodeType()==Corner||Node[Node[nodeID]->Get_connect(1)]->get_NodeType()==ConcaveCorner||Node[Node[nodeID]->Get_connect(1)]->get_NodeType()==ConvexCorner||Node[Node[nodeID]->Get_connect(1)]->get_NodeType()==Wall))
	{
		if(Node[Node[nodeID]->Get_connect(0)]->get_NodeType()==Wall||Node[Node[nodeID]->Get_connect(2)]->get_NodeType()==Wall
			||Node[Node[nodeID]->Get_connect(0)]->get_NodeType()==SpecialWall||Node[Node[nodeID]->Get_connect(2)]->get_NodeType()==SpecialWall
			||(Node[Node[nodeID]->Get_connect(0)]->get_NodeType()==Solid&&Node[Node[nodeID]->Get_connect(2)]->get_NodeType()==Solid))//Corner
			specialwall=false;
		else
		{
			specialwall=true;
			if(Node[Node[nodeID]->Get_connect(0)]->get_NodeType()!=Solid)
				directionwall=4;
			else
				directionwall=2;
		}
	}
	else
		if(y==0 && (Node[Node[nodeID]->Get_connect(2)]->get_NodeType()==Solid||Node[Node[nodeID]->Get_connect(2)]->get_NodeType()==Corner||Node[Node[nodeID]->Get_connect(2)]->get_NodeType()==ConcaveCorner||Node[Node[nodeID]->Get_connect(2)]->get_NodeType()==ConvexCorner||Node[Node[nodeID]->Get_connect(2)]->get_NodeType()==Wall))
	{
		if(Node[Node[nodeID]->Get_connect(1)]->get_NodeType()==Wall||Node[Node[nodeID]->Get_connect(3)]->get_NodeType()==Wall
			||Node[Node[nodeID]->Get_connect(1)]->get_NodeType()==SpecialWall||Node[Node[nodeID]->Get_connect(3)]->get_NodeType()==SpecialWall
			||(Node[Node[nodeID]->Get_connect(1)]->get_NodeType()==Solid&&Node[Node[nodeID]->Get_connect(3)]->get_NodeType()==Solid))//Corner or inside solid
			specialwall=false;
		else
		{
			specialwall=true;
		}
	}
	else if(y==ny && (Node[Node[nodeID]->Get_connect(0)]->get_NodeType()==Solid||Node[Node[nodeID]->Get_connect(0)]->get_NodeType()==Corner||Node[Node[nodeID]->Get_connect(0)]->get_NodeType()==ConcaveCorner||Node[Node[nodeID]->Get_connect(0)]->get_NodeType()==ConvexCorner||Node[Node[nodeID]->Get_connect(0)]->get_NodeType()==Wall))
	{
		if(Node[Node[nodeID]->Get_connect(1)]->get_NodeType()==Wall||Node[Node[nodeID]->Get_connect(3)]->get_NodeType()==Wall
			||Node[Node[nodeID]->Get_connect(1)]->get_NodeType()==SpecialWall||Node[Node[nodeID]->Get_connect(3)]->get_NodeType()==SpecialWall
			||(Node[Node[nodeID]->Get_connect(1)]->get_NodeType()==Solid&&Node[Node[nodeID]->Get_connect(3)]->get_NodeType()==Solid))//Corner
			specialwall=false;
		else
		{
			specialwall=true;
			if(Node[Node[nodeID]->Get_connect(1)]->get_NodeType()!=Solid)
				directionwall=1;
			else
				directionwall=3;
		}
	}
	else
		specialwall=false;

/*		for (unsigned int j=0;j<4;j++)
	{
		if(Node[Node[nodeID]->Get_connect(j)]->get_NodeType()!=Interior && Node[Node[nodeID]->Get_connect(j)]->get_NodeType()!=Corner && Node[Node[nodeID]->Get_connect(j)]->get_NodeType()!=ConcaveCorner && Node[Node[nodeID]->Get_connect(j)]->get_NodeType()!=ConvexCorner && Node[Node[nodeID]->Get_connect(j)]->get_NodeType()!=Wall && Node[Node[nodeID]->Get_connect(j)]->get_NodeType()!=Solid)
		{
			nbspecialwall++;
			directionwall=j;
		}
	}
	if(nbspecialwall==0)
		specialwall=false;
	if(nbspecialwall==1)
	{
		specialwall=true;
		//std::cout<<"**** Detect Special Wall *****    nodeID: "<< nodeID<<" type connect 0: "<<Node[Node[nodeID]->Get_connect(0)]->get_NodeType()<<" type connect 1: "<<Node[Node[nodeID]->Get_connect(1)]->get_NodeType()<<" type connect 2: "<<Node[Node[nodeID]->Get_connect(2)]->get_NodeType()<<" type connect 3: "<<Node[Node[nodeID]->Get_connect(3)]->get_NodeType()<<std::endl;
	}
	if(nbspecialwall>1)
		std::cerr<<" Error detection special wall    nodeID: "<< nodeID<<std::endl;*/
}
void Block2D::CreateCorner(int &nodeID){
	Block2D::ChangeNodeType(nodeID,Corner);
	DefinedCornerType(nodeID);
	IdBoundaries.push_back(nodeID);
}
void Block2D::CreateWall(int &nodeID){
	Block2D::ChangeNodeType(nodeID,Wall);
	IdBoundaries.push_back(nodeID);
}
void Block2D::CreateSpecialWall(int &nodeID){
	Block2D::ChangeNodeType(nodeID,Wall);
	Node[nodeID]->Set_NodeType(SpecialWall);
	IdBoundaries.push_back(nodeID);
}
void Block2D::CreateWallandCorners(Parameters &Param, int &nodeID, int & nbinterior, bool SpecialNode)
{
	if (nbinterior>0)//else solid or concave corner
	{
			if(nbinterior>1)//else convex corner
			{
				if(Param.Get_Verbous())
					std::cout<<"Corner number: "<<nodeID<<std::endl;
				CreateCorner(nodeID);
			}
			else
			{
				CreateWall(nodeID);
			}
	}
	else
	{
			nbinterior=0;// detect corners
			//North East direction
			if(Node[Node[Node[nodeID]->Get_connect(1)]->Get_connect(2)]->get_NodeType()==Interior||Node[Node[Node[nodeID]->Get_connect(1)]->Get_connect(2)]->get_NodeType()==Ghost)
			{
				if(Param.Get_Verbous())
					std::cout<<"Corner number: "<<nodeID<<std::endl;
				CreateCorner(nodeID);
				nbinterior++;
			}
			//North West direction
			if(Node[Node[Node[nodeID]->Get_connect(3)]->Get_connect(2)]->get_NodeType()==Interior||Node[Node[Node[nodeID]->Get_connect(3)]->Get_connect(2)]->get_NodeType()==Ghost)
			{
				if(Param.Get_Verbous())
					std::cout<<"Corner number: "<<nodeID<<std::endl;
				CreateCorner(nodeID);
				nbinterior++;
			}
			//South West direction
			if(Node[Node[Node[nodeID]->Get_connect(3)]->Get_connect(0)]->get_NodeType()==Interior||Node[Node[Node[nodeID]->Get_connect(3)]->Get_connect(0)]->get_NodeType()==Ghost)
			{
				if(Param.Get_Verbous())
					std::cout<<"Corner number: "<<nodeID<<std::endl;
				CreateCorner(nodeID);
				nbinterior++;
			}
			//South East direction
			if(Node[Node[Node[nodeID]->Get_connect(1)]->Get_connect(0)]->get_NodeType()==Interior||Node[Node[Node[nodeID]->Get_connect(1)]->Get_connect(0)]->get_NodeType()==Ghost)
			{
				if(Param.Get_Verbous())
					std::cout<<"Corner number: "<<nodeID<<std::endl;
				CreateCorner(nodeID);
				nbinterior++;
			}
			if(nbinterior++>1 &&SpecialNode==false)
				std::cerr<<"Error detection concave corners; Node ID: "<<nodeID<<" x coordinate: "<<Node[nodeID]->get_x()<<" y coordinate: "<<Node[nodeID]->get_y()<<std::endl;
	}

}
void Block2D::ModifyMeshByUser(Parameters &Param){
	int ndSolidNode=0;
	int firstSolid=0;
	std::vector<int> ndSolidNodetoremove;

	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

	if(Param.Get_Verbous())
		std::cout<<"Processor ID: "<<rank<<" Number of node:"<<NbRealNodes<<std::endl;
	nx=Param.Get_Nx();ny=Param.Get_Ny();
	GenerateSolid(Param, ndSolidNode,firstSolid);
	int computenode=Node.size()-NbGhostNode-ndSolidNode;
	int sizeIdBcBeforeSolid=IdBoundaries.size();
	int nbinterior=0;// detect corners
	bool specialwall=false;
	std::vector<int> IdSpecialNode;
	unsigned int directionspecialwall;

	for (int i=0;i<IdSolidNode.size();i++)
	{

		DetectSpecialWall(IdSolidNode[i], Node[IdSolidNode[i]]->get_x(),Node[IdSolidNode[i]]->get_y(),specialwall,directionspecialwall);
		if(specialwall)
			IdSpecialNode.push_back(IdSolidNode[i]);

		DetectDirectInterior(IdSolidNode[i], nbinterior);
		CreateWallandCorners(Param,IdSolidNode[i] , nbinterior,specialwall);
	}
	//Remove Solid Ghost node convert to wall or corner
	for(int i=0;i<Id_SolidGhost.size();i++)
	{
		Block2D::ChangeNodeType(Id_SolidGhost[i],Solid);
	}
	for(int i=0;i<IdSpecialNode.size();i++)
	{
		DetectSpecialWall(IdSpecialNode[i], Node[IdSpecialNode[i]]->get_x(),Node[IdSpecialNode[i]]->get_y(),specialwall,directionspecialwall);
		if(specialwall)
			CreateSpecialWall(IdSpecialNode[i]);
	}
	if(Param.Get_Verbous())
		for(int i=sizeIdBcBeforeSolid;i<IdBoundaries.size();i++)
			std::cout<<"Processor: "<<rank<<" Node id: "<<IdBoundaries[i]<<" x coordinate: " <<Node[IdBoundaries[i]]->get_x() <<" y coordinate: " << Node[IdBoundaries[i]]->get_y()<<" type of node: "<<Node[IdBoundaries[i]]->get_NodeType()<<std::endl;
	Node_tmp.clear();
	Node_Solidtmp.clear();

}
void Block2D::DefinedCornerType(int nodenumber){
	int nbinterior=0;
	for (unsigned int i=5;i<9;i++)
	{
		if (Node[Connect_lowOrder(nodenumber,i)]->get_NodeType()==Interior ||Node[Connect_lowOrder(nodenumber,i)]->get_NodeType()==Ghost)
		{
			nbinterior++;
		}
	}
///normal convex corners
	if(nbinterior>1)
	{
		Node[nodenumber]->Set_CornerType(ConvexCorner);
		if(verbous)
			std::cout<<"Node number "<<nodenumber <<" is convex and type: "<<Node[nodenumber]->get_NodeType()<<std::endl;
	}
///Normal concave corners
	else
	{
		if(verbous)
			std::cout<<"Node number "<<nodenumber <<" is concave"<<std::endl;
		Node[nodenumber]->Set_CornerType(ConcaveCorner);
	}
}
void Block2D::ConvertToPhysicalUnit(Parameters &Param){
	std::vector<double>::iterator it;
	double convertx=Param.Get_deltax();//local variable is quicker
	double converty=Param.Get_deltax();//local variable is quicker
	for (int i=0;i<x.size() ; ++i)
	{
		x[i]=convertx*x[i];
		y[i]=converty*y[i];
	}

}
void Block2D::reorganizeNodeByType(){
/*	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	char buffer[50]; // make sure it's big enough
	snprintf(buffer, sizeof(buffer), "CheckNodes_%d.txt", rank);
	std::ofstream myFlux;
	myFlux.open(buffer);*/
	unsigned int Connect_N_tmp=0;
	unsigned int Connect_S_tmp=0;
	unsigned int Connect_W_tmp=0;
	unsigned int Connect_E_tmp=0;
	unsigned int NbVelocity_tmp=0;
	for (int i=0;i<Node.size();i++)
	{
		NodeArrays.TypeOfNode[i]=Node[i]->get_NodeType();
		unsigned int x_tmp=(unsigned int)Node[i]->get_x();
		unsigned int y_tmp=(unsigned int)Node[i]->get_y();
		Connect_N_tmp=Node[i]->Get_connect(2)+1;
		Connect_S_tmp=Node[i]->Get_connect(0)+1;
		Connect_W_tmp=Node[i]->Get_connect(3)+1;
		Connect_E_tmp=Node[i]->Get_connect(1)+1;
		NbVelocity_tmp=Node[i]->Get_NbVelocity();
		//if(x_tmp==3 && y_tmp==0)
//			myFlux<<"Processor: "<<rank<<" x: "<<x_tmp<<" y: "<<y_tmp<<" Node type is: "<<NodeArrays.TypeOfNode[i]<<std::endl;
		/*		if(x_tmp==20 && y_tmp==20)
			myFlux<<"Processor: "<<rank<<" x: "<<x_tmp<<" y: "<<y_tmp<<" Node type is: "<<NodeArrays.TypeOfNode[i]<<std::endl;
		if(x_tmp==20 && y_tmp==21)
			myFlux<<"Processor: "<<rank<<" x: "<<x_tmp<<" y: "<<y_tmp<<" Node type is: "<<NodeArrays.TypeOfNode[i]<<std::endl;
		if(x_tmp==19 && y_tmp==19)
			myFlux<<"Processor: "<<rank<<" x: "<<x_tmp<<" y: "<<y_tmp<<" Node type is: "<<NodeArrays.TypeOfNode[i]<<std::endl;
		if(x_tmp==19 && y_tmp==20)
			myFlux<<"Processor: "<<rank<<" x: "<<x_tmp<<" y: "<<y_tmp<<" Node type is: "<<NodeArrays.TypeOfNode[i]<<std::endl;
		if(x_tmp==19 && y_tmp==21)
			myFlux<<"Processor: "<<rank<<" x: "<<x_tmp<<" y: "<<y_tmp<<" Node type is: "<<NodeArrays.TypeOfNode[i]<<std::endl;
*/
		switch(NodeArrays.TypeOfNode[i])
		{
		case Interior:
			NodeArrays.NodeIndexByType[i]=NodeArrays.NodeInterior.size();
			NodeArrays.NodeInterior.push_back(NodeInterior2D(x_tmp,y_tmp));
			NodeArrays.NodeInterior[NodeArrays.NodeInterior.size()-1].Set_Connect(Connect_N_tmp,Connect_S_tmp,Connect_W_tmp,Connect_E_tmp);
			NodeArrays.NodeInterior[NodeArrays.NodeInterior.size()-1].Set_NbVelocity(NbVelocity_tmp);
			NodeArrays.NodeInterior[NodeArrays.NodeInterior.size()-1].Set_Index(i);
			break;
		case Solid:
			NodeArrays.NodeIndexByType[i]=NodeArrays.NodeSolid.size();
			NodeArrays.NodeSolid.push_back(NodeSolid2D(x_tmp,y_tmp));
			NodeArrays.NodeSolid[NodeArrays.NodeSolid.size()-1].Set_Connect(Connect_N_tmp,Connect_S_tmp,Connect_W_tmp,Connect_E_tmp);
			NodeArrays.NodeSolid[NodeArrays.NodeSolid.size()-1].Set_NbVelocity(NbVelocity_tmp);
			NodeArrays.NodeSolid[NodeArrays.NodeSolid.size()-1].Set_Index(i);
			break;
		case Ghost:
			NodeArrays.NodeIndexByType[i]=NodeArrays.NodeGhost.size();
			NodeArrays.NodeGhost.push_back(NodeGhost2D(x_tmp,y_tmp));
			NodeArrays.NodeGhost[NodeArrays.NodeGhost.size()-1].Set_Connect(Connect_N_tmp,Connect_S_tmp,Connect_W_tmp,Connect_E_tmp);
			NodeArrays.NodeGhost[NodeArrays.NodeGhost.size()-1].Set_NbVelocity(NbVelocity_tmp);
			NodeArrays.NodeGhost[NodeArrays.NodeGhost.size()-1].Set_Index(i);
			break;
		case GlobalCorner:
			NodeArrays.NodeIndexByType[i]=NodeArrays.NodeGlobalCorner.size();
			NodeArrays.NodeGlobalCorner.push_back(NodeCorner2D(x_tmp,y_tmp,Node[i]->Get_RhoDef(),Node[i]->Get_UDef()));
			NodeArrays.NodeGlobalCorner[NodeArrays.NodeGlobalCorner.size()-1].Set_Connect(Connect_N_tmp,Connect_S_tmp,Connect_W_tmp,Connect_E_tmp);
			NodeArrays.NodeGlobalCorner[NodeArrays.NodeGlobalCorner.size()-1].Set_NbVelocity(NbVelocity_tmp);
			NodeArrays.NodeGlobalCorner[NodeArrays.NodeGlobalCorner.size()-1].Set_Index(i);
			NodeArrays.NodeGlobalCorner[NodeArrays.NodeGlobalCorner.size()-1].Set_CornerType(Concave);
			NodeArrays.NodeGlobalCorner[NodeArrays.NodeGlobalCorner.size()-1].Set_NodeType(GlobalCorner);
			break;
		case Corner:
			NodeArrays.NodeIndexByType[i]=NodeArrays.NodeCorner.size();
			NodeArrays.NodeCorner.push_back(NodeCorner2D(x_tmp,y_tmp,Node[i]->Get_RhoDef(),Node[i]->Get_UDef()));
			NodeArrays.NodeCorner[NodeArrays.NodeCorner.size()-1].Set_Connect(Connect_N_tmp,Connect_S_tmp,Connect_W_tmp,Connect_E_tmp);
			NodeArrays.NodeCorner[NodeArrays.NodeCorner.size()-1].Set_NbVelocity(NbVelocity_tmp);
			NodeArrays.NodeCorner[NodeArrays.NodeCorner.size()-1].Set_Index(i);
			NodeArrays.NodeCorner[NodeArrays.NodeCorner.size()-1].Set_CornerType(Concave);
			break;
		case SpecialWall:
			NodeArrays.NodeIndexByType[i]=NodeArrays.NodeSpecialWall.size();
			NodeArrays.NodeSpecialWall.push_back(NodeWall2D(x_tmp,y_tmp));
			NodeArrays.NodeSpecialWall[NodeArrays.NodeSpecialWall.size()-1].Set_Connect(Connect_N_tmp,Connect_S_tmp,Connect_W_tmp,Connect_E_tmp);
			NodeArrays.NodeSpecialWall[NodeArrays.NodeSpecialWall.size()-1].Set_NbVelocity(NbVelocity_tmp);
			NodeArrays.NodeSpecialWall[NodeArrays.NodeSpecialWall.size()-1].Set_Index(i);
			NodeArrays.NodeSpecialWall[NodeArrays.NodeSpecialWall.size()-1].Set_NodeType(SpecialWall);
			break;
		case Wall:
			NodeArrays.NodeIndexByType[i]=NodeArrays.NodeWall.size();
			NodeArrays.NodeWall.push_back(NodeWall2D(x_tmp,y_tmp));
			NodeArrays.NodeWall[NodeArrays.NodeWall.size()-1].Set_Connect(Connect_N_tmp,Connect_S_tmp,Connect_W_tmp,Connect_E_tmp);
			NodeArrays.NodeWall[NodeArrays.NodeWall.size()-1].Set_NbVelocity(NbVelocity_tmp);
			NodeArrays.NodeWall[NodeArrays.NodeWall.size()-1].Set_Index(i);
			break;
		case Velocity:
			NodeArrays.NodeIndexByType[i]=NodeArrays.NodeVelocity.size();
			NodeArrays.NodeVelocity.push_back(NodeVelocity2D(x_tmp,y_tmp,Node[i]->Get_UDef()));
			NodeArrays.NodeVelocity[NodeArrays.NodeVelocity.size()-1].Set_Connect(Connect_N_tmp,Connect_S_tmp,Connect_W_tmp,Connect_E_tmp);
			NodeArrays.NodeVelocity[NodeArrays.NodeVelocity.size()-1].Set_NbVelocity(NbVelocity_tmp);
			NodeArrays.NodeVelocity[NodeArrays.NodeVelocity.size()-1].Set_Index(i);
			break;
		case Pressure:
			NodeArrays.NodeIndexByType[i]=NodeArrays.NodePressure.size();
			NodeArrays.NodePressure.push_back(NodePressure2D(x_tmp,y_tmp,Node[i]->Get_RhoDef()));
			NodeArrays.NodePressure[NodeArrays.NodePressure.size()-1].Set_Connect(Connect_N_tmp,Connect_S_tmp,Connect_W_tmp,Connect_E_tmp);
			NodeArrays.NodePressure[NodeArrays.NodePressure.size()-1].Set_NbVelocity(NbVelocity_tmp);
			NodeArrays.NodePressure[NodeArrays.NodePressure.size()-1].Set_Index(i);
			break;
		case ConcaveCorner:
			NodeArrays.NodeIndexByType[i]=NodeArrays.NodeCorner.size();
			NodeArrays.NodeCorner.push_back(NodeCorner2D(x_tmp,y_tmp,Node[i]->Get_RhoDef(),Node[i]->Get_UDef()));
			NodeArrays.NodeCorner[NodeArrays.NodeCorner.size()-1].Set_Connect(Connect_N_tmp,Connect_S_tmp,Connect_W_tmp,Connect_E_tmp);
			NodeArrays.NodeCorner[NodeArrays.NodeCorner.size()-1].Set_NbVelocity(NbVelocity_tmp);
			NodeArrays.NodeCorner[NodeArrays.NodeCorner.size()-1].Set_CornerType(Concave);
			NodeArrays.NodeCorner[NodeArrays.NodeCorner.size()-1].Set_Index(i);
			break;
		case ConvexCorner:
			NodeArrays.NodeIndexByType[i]=NodeArrays.NodeCorner.size();
			NodeArrays.NodeCorner.push_back(NodeCorner2D(x_tmp,y_tmp,Node[i]->Get_RhoDef(),Node[i]->Get_UDef()));
			NodeArrays.NodeCorner[NodeArrays.NodeCorner.size()-1].Set_Connect(Connect_N_tmp,Connect_S_tmp,Connect_W_tmp,Connect_E_tmp);
			NodeArrays.NodeCorner[NodeArrays.NodeCorner.size()-1].Set_NbVelocity(NbVelocity_tmp);
			NodeArrays.NodeCorner[NodeArrays.NodeCorner.size()-1].Set_CornerType(Convex);
			NodeArrays.NodeCorner[NodeArrays.NodeCorner.size()-1].Set_Index(i);
			break;
		case SolidGhost:
			NodeArrays.NodeIndexByType[i]=NodeArrays.NodeGhost.size();
			NodeArrays.NodeGhost.push_back(NodeGhost2D(x_tmp,y_tmp));
			NodeArrays.NodeGhost[NodeArrays.NodeGhost.size()-1].Set_Connect(Connect_N_tmp,Connect_S_tmp,Connect_W_tmp,Connect_E_tmp);
			NodeArrays.NodeGhost[NodeArrays.NodeGhost.size()-1].Set_NbVelocity(NbVelocity_tmp);
			NodeArrays.NodeGhost[NodeArrays.NodeGhost.size()-1].Set_CornerType(SolidGhost);
			NodeArrays.NodeGhost[NodeArrays.NodeGhost.size()-1].Set_Index(i);
			break;
		case Symmetry:
			NodeArrays.NodeIndexByType[i]=NodeArrays.NodeSymmetry.size();
			NodeArrays.NodeSymmetry.push_back(NodeSymmetry2D(x_tmp,y_tmp));
			NodeArrays.NodeSymmetry[NodeArrays.NodeSymmetry.size()-1].Set_Connect(Connect_N_tmp,Connect_S_tmp,Connect_W_tmp,Connect_E_tmp);
			NodeArrays.NodeSymmetry[NodeArrays.NodeSymmetry.size()-1].Set_NbVelocity(NbVelocity_tmp);
			NodeArrays.NodeSymmetry[NodeArrays.NodeSymmetry.size()-1].Set_Index(i);
			// ******** Set Symmetry Type ********
			SetSymmetryType(NodeArrays.NodeSymmetry[NodeArrays.NodeSymmetry.size()-1].Get_SymmetryType(),NodeArrays.NodeSymmetry[NodeArrays.NodeSymmetry.size()-1].get_x() , NodeArrays.NodeSymmetry[NodeArrays.NodeSymmetry.size()-1].get_y());
			break;
		case Periodic:
			NodeArrays.NodeIndexByType[i]=NodeArrays.NodePeriodic.size();
			NodeArrays.NodePeriodic.push_back(NodePeriodic2D(x_tmp,y_tmp));
			NodeArrays.NodePeriodic[NodeArrays.NodePeriodic.size()-1].Set_Connect(Connect_N_tmp,Connect_S_tmp,Connect_W_tmp,Connect_E_tmp);
			NodeArrays.NodePeriodic[NodeArrays.NodePeriodic.size()-1].Set_NbVelocity(NbVelocity_tmp);
			NodeArrays.NodePeriodic[NodeArrays.NodePeriodic.size()-1].Set_Index(i);
			break;
		default:
			std::cerr<<"Type "<< NodeArrays.TypeOfNode[i]<<" of node not find in the sort of Node array and type "<<std::endl;
			break;
		}
	}
}
NodeArrays2D* Block2D::Get_NodeArrays2D(){
	return &NodeArrays;
}

///
void Block2D::Set_Connect(Parameters& Param){
	int* tmpConnect=0;
	tmpConnect=new int[Param.Get_NbVelocities()];
	for(int i=0;i<NbTotalNodes;i++)
	{
		tmpConnect[0]=i;
		for (unsigned int j=1;j<Param.Get_NbVelocities();j++)
			tmpConnect[j]=Connect_lowOrder(i,j);
		switch(NodeArrays.TypeOfNode[i])
		{
		case Interior:
			NodeArrays.NodeInterior[NodeArrays.NodeIndexByType[i]].Set_Connect(tmpConnect,Param.Get_NbVelocities());
			break;
		case Solid:
			NodeArrays.NodeSolid[NodeArrays.NodeIndexByType[i]].Set_Connect(tmpConnect,Param.Get_NbVelocities());
			break;
		case Ghost:
			NodeArrays.NodeGhost[NodeArrays.NodeIndexByType[i]].Set_Connect(tmpConnect,Param.Get_NbVelocities());
			break;
		case Corner:
			NodeArrays.NodeCorner[NodeArrays.NodeIndexByType[i]].Set_Connect(tmpConnect,Param.Get_NbVelocities());
			break;
		case GlobalCorner:
			NodeArrays.NodeGlobalCorner[NodeArrays.NodeIndexByType[i]].Set_Connect(tmpConnect,Param.Get_NbVelocities());
			break;
		case Wall:
			NodeArrays.NodeWall[NodeArrays.NodeIndexByType[i]].Set_Connect(tmpConnect,Param.Get_NbVelocities());
			break;
		case SpecialWall:
			NodeArrays.NodeSpecialWall[NodeArrays.NodeIndexByType[i]].Set_Connect(tmpConnect,Param.Get_NbVelocities());
			break;
		case Velocity:
			NodeArrays.NodeVelocity[NodeArrays.NodeIndexByType[i]].Set_Connect(tmpConnect,Param.Get_NbVelocities());
			break;
		case Pressure:
			NodeArrays.NodePressure[NodeArrays.NodeIndexByType[i]].Set_Connect(tmpConnect,Param.Get_NbVelocities());
			break;
		case ConcaveCorner:
			NodeArrays.NodeCorner[NodeArrays.NodeIndexByType[i]].Set_Connect(tmpConnect,Param.Get_NbVelocities());
			break;
		case ConvexCorner:
			NodeArrays.NodeCorner[NodeArrays.NodeIndexByType[i]].Set_Connect(tmpConnect,Param.Get_NbVelocities());
			break;
		case SolidGhost:
			NodeArrays.NodeGhost[NodeArrays.NodeIndexByType[i]].Set_Connect(tmpConnect,Param.Get_NbVelocities());
			break;
		case Symmetry:
			NodeArrays.NodeSymmetry[NodeArrays.NodeIndexByType[i]].Set_Connect(tmpConnect,Param.Get_NbVelocities());
			break;
		case Periodic:
			NodeArrays.NodePeriodic[NodeArrays.NodeIndexByType[i]].Set_Connect(tmpConnect,Param.Get_NbVelocities());
			break;
		default:
			std::cerr<<"Type of node not find in the sort of Node array"<<std::endl;
			break;
		}
	}
	delete [] tmpConnect;
	for(int i=0;i<NbTotalNodes;i++)
		delete Node[i];
	Node.clear();
	Block2D::Set_BcNormal();
}
int Block2D::Connect_lowOrder(int &NodeNumber,unsigned int& direction){
	switch (direction)
	{
	case 1:
		return Node[NodeNumber]->Get_connect(1);
		break;
	case 2:
		return Node[NodeNumber]->Get_connect(2);
		break;
	case 3:
		return Node[NodeNumber]->Get_connect(3);
		break;
	case 4:
		return Node[NodeNumber]->Get_connect(0);
		break;
	case 5:
		if ((Node[NodeNumber]->Get_connect(1)==(unsigned int)NodeNumber)||(Node[NodeNumber]->Get_connect(2)==(unsigned int)NodeNumber))
			return (unsigned int)NodeNumber;
		else
			return Node[Node[NodeNumber]->Get_connect(1)]->Get_connect(2);
		break;
	case 6:
		if ((Node[NodeNumber]->Get_connect(3)==(unsigned int)NodeNumber)||(Node[NodeNumber]->Get_connect(2)==(unsigned int)NodeNumber))
			return (unsigned int)NodeNumber;
		else
			return Node[Node[NodeNumber]->Get_connect(3)]->Get_connect(2);
		break;
	case 7:
		if ((Node[NodeNumber]->Get_connect(3)==(unsigned int)NodeNumber)||(Node[NodeNumber]->Get_connect(0)==(unsigned int)NodeNumber))
			return (unsigned int)NodeNumber;
		else
			return Node[Node[NodeNumber]->Get_connect(3)]->Get_connect(0);
		break;
	case 8:
		if ((Node[NodeNumber]->Get_connect(1)==(unsigned int)NodeNumber)||(Node[NodeNumber]->Get_connect(0)==(unsigned int)NodeNumber))
			return (unsigned int)NodeNumber;
		else
			return Node[Node[NodeNumber]->Get_connect(1)]->Get_connect(0);
		break;
	default:
		std::cerr<<"Wrong Connection ask. Direction "<<direction<<" is not defined."<<std::endl;
		return direction;
		break;
	}
}
void Block2D::Set_BcNormal()
{
/*	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	char buffer[50]; // make sure it's big enough
	snprintf(buffer, sizeof(buffer), "normalWall_%d.txt", rank);
	std::ofstream myFlux;
	myFlux.open(buffer);*/
	for (int i=0;i<NodeArrays.NodeWall.size();i++)
	{
			NodeArrays.NodeWall[i].Set_BcNormal(Get_BcNormal(NodeArrays.NodeWall[i]));
	}
	for (int i=0;i<NodeArrays.NodeSpecialWall.size();i++)
	{
			NodeArrays.NodeSpecialWall[i].Set_BcNormal(Get_BcNormal_SpecialWall(NodeArrays.NodeSpecialWall[i]));
//			std::cout<<std::endl<<"x: "<<NodeArrays.NodeSpecialWall[i].get_x()<<"y: "<<NodeArrays.NodeSpecialWall[i].get_y()<<" Normal: "<<Get_BcNormal_SpecialWall(NodeArrays.NodeSpecialWall[i])<<std::endl;
	}
	for (int i=0;i<NodeArrays.NodeCorner.size();i++)
	{
			NodeArrays.NodeCorner[i].Set_BcNormal(Get_BcNormal(NodeArrays.NodeCorner[i]));
	}
	for (int i=0;i<NodeArrays.NodeGlobalCorner.size();i++)
	{
			NodeArrays.NodeGlobalCorner[i].Set_BcNormal(Get_BcNormal(NodeArrays.NodeGlobalCorner[i]));
	}
	for (int i=0;i<NodeArrays.NodeSymmetry.size();i++)
	{
			NodeArrays.NodeSymmetry[i].Set_BcNormal(Get_BcNormal(NodeArrays.NodeSymmetry[i]));
	}
	for (int i=0;i<NodeArrays.NodePeriodic.size();i++)
	{
			NodeArrays.NodePeriodic[i].Set_BcNormal(Get_BcNormal(NodeArrays.NodePeriodic[i]));
	}
	for (int i=0;i<NodeArrays.NodePressure.size();i++)
	{
			NodeArrays.NodePressure[i].Set_BcNormal(Get_BcNormal(NodeArrays.NodePressure[i]));
	}
	for (int i=0;i<NodeArrays.NodeVelocity.size();i++)
	{
			NodeArrays.NodeVelocity[i].Set_BcNormal(Get_BcNormal(NodeArrays.NodeVelocity[i]));
	}
}

int Block2D::Get_BcNormal(NodeCorner2D & nodeIn)
{
	int nbinterior=0, count=0;
	int tmpdirection,intTmpReturn=0;
	if(nodeIn.get_NodeType()==GlobalCorner)
	{
		///Normal concave corners
			for (unsigned int i=5;i<9;i++)
			{
				if (NodeArrays.TypeOfNode[nodeIn.Get_connect()[i]]==Interior )
				{
					intTmpReturn=i;
					nbinterior++;
				}
				if(nbinterior>1)
				{
					std::cerr<<"Wrong corner normal detection for Global corner. "<<std::endl;
				}
			}
	}
	else
	{
///Normal concave corners
	for (unsigned int i=5;i<9;i++)
	{
//		if (!(nodeIn.Get_connect()[i]==nodeIn.Get_index()||nodeIn.get_NodeType()==Solid||NodeArrays.TypeOfNode[nodeIn.Get_connect()[i]]==Solid) && (NodeArrays.TypeOfNode[nodeIn.Get_connect()[i]]==Interior ||NodeArrays.TypeOfNode[nodeIn.Get_connect()[i]]==Ghost))
		if (NodeArrays.TypeOfNode[nodeIn.Get_connect()[i]]==Interior ||NodeArrays.TypeOfNode[nodeIn.Get_connect()[i]]==Ghost)
		{
			intTmpReturn=i;
			nbinterior++;
		}
	}
///normal convex corners
		if(nbinterior>1)
		{
/*			for(int i=0;i<4;i++)
				tmpdirection[i]=0;
			int count=0;
			for (unsigned int i=1;i<5;i++)
			{
				if (NodeArrays.TypeOfNode[nodeIn.Get_connect()[i]]!=Interior && NodeArrays.TypeOfNode[nodeIn.Get_connect()[i]]!=Ghost )
				{
					tmpdirection[count]=i;
					count++;
				}
			}
			if(tmpdirection[0]==3 && tmpdirection[1]==4)
				intTmpReturn=5;
			if(tmpdirection[0]==1 && tmpdirection[1]==4)
				intTmpReturn=6;
			if(tmpdirection[0]==1 && tmpdirection[1]==2)
				intTmpReturn=7;
			if(tmpdirection[0]==2 && tmpdirection[1]==3)
				intTmpReturn=8;
*/
			for (unsigned int i=5;i<9;i++)
			{
				if (NodeArrays.TypeOfNode[nodeIn.Get_connect()[i]]==Solid )
				{
					tmpdirection=i;
					count++;
				}			}
			switch(tmpdirection)
			{
				case 5:	
					intTmpReturn=7;
					break;	
				case 6:	
					intTmpReturn=8;
					break;	
				case 7:	
					intTmpReturn=5;
					break;	
				case 8:	
					intTmpReturn=6;
					break;	
				default:
					std::cerr<<"Wrong convex corner normal detection. Detection of direction: "<<tmpdirection<<" x: "<<nodeIn.get_x()<<" y: "<<nodeIn.get_y()<<" node index: "<<nodeIn.Get_index()<<std::endl;
					break;
			}

		}
	}

	return intTmpReturn;
}
int Block2D::Get_BcNormal(NodeWall2D & nodeIn)
{
	intTmpReturn=0;
	for (unsigned int i=1;i<5;i++)
		if (NodeArrays.TypeOfNode[nodeIn.Get_connect()[i]]==Interior)
			intTmpReturn=i;
	if(intTmpReturn==0)
		for (unsigned int i=1;i<5;i++)
			if (NodeArrays.TypeOfNode[nodeIn.Get_connect()[i]]==Ghost)
			{
				if(i==1)
				{
					if (NodeArrays.TypeOfNode[nodeIn.Get_connect()[2]]==Wall||NodeArrays.TypeOfNode[nodeIn.Get_connect()[4]]==Wall
						||NodeArrays.TypeOfNode[nodeIn.Get_connect()[2]]==SpecialWall||NodeArrays.TypeOfNode[nodeIn.Get_connect()[4]]==SpecialWall
						||NodeArrays.TypeOfNode[nodeIn.Get_connect()[2]]==Corner||NodeArrays.TypeOfNode[nodeIn.Get_connect()[4]]==Corner
						||NodeArrays.TypeOfNode[nodeIn.Get_connect()[2]]==ConcaveCorner||NodeArrays.TypeOfNode[nodeIn.Get_connect()[4]]==ConcaveCorner
						||NodeArrays.TypeOfNode[nodeIn.Get_connect()[2]]==ConvexCorner||NodeArrays.TypeOfNode[nodeIn.Get_connect()[4]]==ConvexCorner)
						intTmpReturn=i;
				}
				else if(i==4)
				{
					if (NodeArrays.TypeOfNode[nodeIn.Get_connect()[1]]==Wall||NodeArrays.TypeOfNode[nodeIn.Get_connect()[3]]==Wall
						||NodeArrays.TypeOfNode[nodeIn.Get_connect()[1]]==SpecialWall||NodeArrays.TypeOfNode[nodeIn.Get_connect()[3]]==SpecialWall
						||NodeArrays.TypeOfNode[nodeIn.Get_connect()[1]]==Corner||NodeArrays.TypeOfNode[nodeIn.Get_connect()[3]]==Corner
						||NodeArrays.TypeOfNode[nodeIn.Get_connect()[1]]==ConcaveCorner||NodeArrays.TypeOfNode[nodeIn.Get_connect()[3]]==ConcaveCorner
						||NodeArrays.TypeOfNode[nodeIn.Get_connect()[1]]==ConvexCorner||NodeArrays.TypeOfNode[nodeIn.Get_connect()[3]]==ConvexCorner)
					intTmpReturn=i;
				}
				else
				{
					if (NodeArrays.TypeOfNode[nodeIn.Get_connect()[i+1]]==Wall||NodeArrays.TypeOfNode[nodeIn.Get_connect()[i-1]]==Wall
						||NodeArrays.TypeOfNode[nodeIn.Get_connect()[i+1]]==SpecialWall||NodeArrays.TypeOfNode[nodeIn.Get_connect()[i-1]]==SpecialWall
						||NodeArrays.TypeOfNode[nodeIn.Get_connect()[i+1]]==Corner||NodeArrays.TypeOfNode[nodeIn.Get_connect()[i-1]]==Corner
						||NodeArrays.TypeOfNode[nodeIn.Get_connect()[i+1]]==ConcaveCorner||NodeArrays.TypeOfNode[nodeIn.Get_connect()[i-1]]==ConcaveCorner
						||NodeArrays.TypeOfNode[nodeIn.Get_connect()[i+1]]==ConvexCorner||NodeArrays.TypeOfNode[nodeIn.Get_connect()[i-1]]==ConvexCorner)
					intTmpReturn=i;
				}
			}
	return intTmpReturn;
}
int Block2D::Get_BcNormal_SpecialWall(NodeWall2D & nodeIn)
{
	if(nodeIn.get_x()==0 )
	{
		if(NodeArrays.TypeOfNode[nodeIn.Get_connect()[2]]!=Solid)
			intTmpReturn=5;
		else
			intTmpReturn=8;
	}
	if(nodeIn.get_x()==nx )
	{
		if(NodeArrays.TypeOfNode[nodeIn.Get_connect()[2]]!=Solid)
			intTmpReturn=6;
		else
			intTmpReturn=7;
	}
	if(nodeIn.get_y()==0 )
	{
		if(NodeArrays.TypeOfNode[nodeIn.Get_connect()[1]]!=Solid)
			intTmpReturn=5;
		else
			intTmpReturn=6;
	}
	if(nodeIn.get_y()==ny )
	{
		if(NodeArrays.TypeOfNode[nodeIn.Get_connect()[1]]!=Solid)
			intTmpReturn=8;
		else
			intTmpReturn=7;
	}
	return intTmpReturn;
}

int Block2D::Get_BcNormal(NodeSymmetry2D & nodeIn)
{
	int nbinterior=0;
	for (unsigned int i=1;i<5;i++)
//		if (nodeIn.stream()[i] && NodeArrays.TypeOfNode[nodeIn.Get_connect()[i]]==Interior)
		if (NodeArrays.TypeOfNode[nodeIn.Get_connect()[i]]==Interior)
			intTmpReturn=i;
	return intTmpReturn;
}
int Block2D::Get_BcNormal(NodePeriodic2D & nodeIn)
{
	int nbinterior=0;
	for (unsigned int i=1;i<5;i++)
//		if (nodeIn.stream()[i] && NodeArrays.TypeOfNode[nodeIn.Get_connect()[i]]==Interior)
		if (NodeArrays.TypeOfNode[nodeIn.Get_connect()[i]]==Interior)
			intTmpReturn=i;
	return intTmpReturn;
}
int Block2D::Get_BcNormal(NodePressure2D & nodeIn)
{
	int nbinterior=0;
	for (unsigned int i=1;i<5;i++)
		if (NodeArrays.TypeOfNode[nodeIn.Get_connect()[i]]==Interior)
			intTmpReturn=i;
	return intTmpReturn;
}
int Block2D::Get_BcNormal(NodeVelocity2D & nodeIn)
{
	int nbinterior=0;
	for (unsigned int i=1;i<5;i++)
		if (NodeArrays.TypeOfNode[nodeIn.Get_connect()[i]]==Interior)
			intTmpReturn=i;
	return intTmpReturn;
}
void Block2D::Get_GhostType(std::vector<int> & NodeTypeN,std::vector<int> & NodeTypeE,std::vector<int> & NodeTypeS,std::vector<int> & NodeTypeW,
		  std::vector<int> & NodeTypeSW,std::vector<int> & NodeTypeSE,std::vector<int> & NodeTypeNW,std::vector<int> & NodeTypeNE)
{

// Get Node type of the sides of the 2D block for the real nodes
	for(int i=0;i<IdRNodeN.size();i++)
		NodeTypeN.push_back((int)Node[IdRNodeN[i]]->get_NodeType());
	for(int i=0;i<IdRNodeE.size();i++)
		NodeTypeE.push_back((int)Node[IdRNodeE[i]]->get_NodeType());
	for(int i=0;i<IdRNodeS.size();i++)
		NodeTypeS.push_back((int)Node[IdRNodeS[i]]->get_NodeType());
	for(int i=0;i<IdRNodeW.size();i++)
		NodeTypeW.push_back((int)Node[IdRNodeW[i]]->get_NodeType());
// Get Node type of the corners of the 2D block
	for(int i=0;i<IdRNodeSW.size();i++)
		NodeTypeSW.push_back((int)Node[IdRNodeSW[i]]->get_NodeType());
	for(int i=0;i<IdRNodeSE.size();i++)
		NodeTypeSE.push_back((int)Node[IdRNodeSE[i]]->get_NodeType());
	for(int i=0;i<IdRNodeNW.size();i++)
		NodeTypeNW.push_back((int)Node[IdRNodeNW[i]]->get_NodeType());
	for(int i=0;i<IdRNodeNE.size();i++)
		NodeTypeNE.push_back((int)Node[IdRNodeNE[i]]->get_NodeType());
/*
 	int rank;
 	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
 	char buffer[50]; // make sure it's big enough
 	snprintf(buffer, sizeof(buffer), "NodeTypeReal_%d.txt", rank);
 	std::ofstream myFlux;
 	myFlux.open(buffer);
 	myFlux<<std::endl<<" ****** after Removing nodes *******"<<std::endl;
 	myFlux<<"West real nodes: "<<std::endl;
 	for(int i=0;i<IdRNodeW.size();i++)
 		myFlux<<IdRNodeW[i]<<" NodeType: "<<NodeTypeW[i]<<" ";
 	myFlux<<std::endl<<"West Ghost nodes: "<<std::endl;
 	for(int i=0;i<IdGNodeW.size();i++)
 		myFlux<<IdGNodeW[i]<<" ";
 	myFlux<<std::endl<<"North real nodes: "<<std::endl;
 	for(int i=0;i<IdRNodeN.size();i++)
 		myFlux<<IdRNodeN[i]<<" NodeType: "<<NodeTypeN[i]<<" ";
 	myFlux<<std::endl<<"North Ghost nodes: "<<std::endl;
 	for(int i=0;i<IdGNodeN.size();i++)
 		myFlux<<IdGNodeN[i]<<" ";
 	myFlux<<std::endl<<"South real nodes: "<<std::endl;
 	for(int i=0;i<IdRNodeS.size();i++)
 		myFlux<<IdRNodeS[i]<<" NodeType: "<<NodeTypeS[i]<<" ";
 	myFlux<<std::endl<<"South Ghost nodes: "<<std::endl;
 	for(int i=0;i<IdGNodeS.size();i++)
 		myFlux<<IdGNodeS[i]<<" ";
 	myFlux<<std::endl<<"East real nodes: "<<std::endl;
 	for(int i=0;i<IdRNodeE.size();i++)
 		myFlux<<IdRNodeE[i]<<" NodeType: "<<NodeTypeE[i]<<" ";
 	myFlux<<std::endl<<"East Ghost nodes: "<<std::endl;
 	for(int i=0;i<IdGNodeE.size();i++)
 		myFlux<<IdGNodeE[i]<<" ";

 	myFlux<<std::endl<<std::endl;
 	myFlux<<"South West real nodes: "<<std::endl;

 	for(int i=0;i<IdRNodeSW.size();i++)
 		myFlux<<IdRNodeSW[i]<<" NodeType: "<<NodeTypeSW[i]<<" ";
 	myFlux<<std::endl<<"South West Ghost nodes: "<<std::endl;
 	for(int i=0;i<IdGNodeSW.size();i++)
 		myFlux<<IdGNodeSW[i]<<" ";
 	myFlux<<std::endl<<"North West real nodes: "<<std::endl;
 	for(int i=0;i<IdRNodeNW.size();i++)
 		myFlux<<IdRNodeNW[i]<<" NodeType: "<<NodeTypeNW[i]<<" ";
 	myFlux<<std::endl<<"North West Ghost nodes: "<<std::endl;
 	for(int i=0;i<IdGNodeNW.size();i++)
 		myFlux<<IdGNodeNW[i]<<" ";
 	myFlux<<std::endl<<"South East real nodes: "<<std::endl;
 	for(int i=0;i<IdRNodeSE.size();i++)
 		myFlux<<IdRNodeSE[i]<<" NodeType: "<<NodeTypeSE[i]<<" ";
 	myFlux<<std::endl<<"South East Ghost nodes: "<<std::endl;
 	for(int i=0;i<IdGNodeSE.size();i++)
 		myFlux<<IdGNodeSE[i]<<" ";
 	myFlux<<std::endl<<"North East real nodes: "<<std::endl;
 	for(int i=0;i<IdRNodeNE.size();i++)
 		myFlux<<IdRNodeNE[i]<<" NodeType: "<<NodeTypeNE[i]<<" ";
 	myFlux<<std::endl<<"North East Ghost nodes: "<<std::endl;
 	for(int i=0;i<IdGNodeNE.size();i++)
 		myFlux<<IdGNodeNE[i]<<" ";
 	myFlux<<std::endl<<" Size of NE Ghost: "<< IdGNodeNE.size()<<std::endl;*/

}
void Block2D::Set_GhostType(std::vector<int> & NodeTypeN,std::vector<int> & NodeTypeE,std::vector<int> & NodeTypeS,std::vector<int> & NodeTypeW,
		  std::vector<int> & NodeTypeSW,std::vector<int> & NodeTypeSE,std::vector<int> & NodeTypeNW,std::vector<int> & NodeTypeNE)
{

	// Get Node type of the sides of the 2D block
	if(IdGNodeN.size()>0)
	for(int i=0;i<IdGNodeN.size();i++)
		Correct_GhostType(IdGNodeN[i], (NodeType)NodeTypeN[i]);
	if(IdGNodeE.size()>0)
	for(int i=0;i<IdGNodeE.size();i++)
		Correct_GhostType(IdGNodeE[i], (NodeType)NodeTypeE[i]);
	if(IdGNodeS.size()>0)
	for(int i=0;i<IdGNodeS.size();i++)
		Correct_GhostType(IdGNodeS[i], (NodeType)NodeTypeS[i]);
	if(IdGNodeW.size()>0)
	for(int i=0;i<IdGNodeW.size();i++)
		Correct_GhostType(IdGNodeW[i], (NodeType)NodeTypeW[i]);
// Get Node type of the corners of the 2D block
	if(IdGNodeSW.size()>0)
	for(int i=0;i<IdGNodeSW.size();i++)
		Correct_GhostType(IdGNodeSW[i], (NodeType)NodeTypeSW[i]);
	if(IdGNodeSE.size()>0)
	for(int i=0;i<IdGNodeSE.size();i++)
		Correct_GhostType(IdGNodeSE[i], (NodeType)NodeTypeSE[i]);
	if(IdGNodeNW.size()>0)
	for(int i=0;i<IdGNodeNW.size();i++)
		Correct_GhostType(IdGNodeNW[i], (NodeType)NodeTypeNW[i]);
	if(IdGNodeNE.size()>0)
	for(int i=0;i<IdGNodeNE.size();i++)
		Correct_GhostType(IdGNodeNE[i], (NodeType)NodeTypeNE[i]);
/*
	int rank;
 	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
 	char buffer[50]; // make sure it's big enough
 	snprintf(buffer, sizeof(buffer), "NodeTypeGhost_%d.txt", rank);
 	std::ofstream myFlux;
 	myFlux.open(buffer);
 	myFlux<<std::endl<<" ****** after Removing nodes *******"<<std::endl;
 	myFlux<<"West real nodes: "<<std::endl;
 	for(int i=0;i<IdRNodeW.size();i++)
 		myFlux<<IdRNodeW[i]<<" ";
 	myFlux<<std::endl<<"West Ghost nodes: "<<std::endl;
 	for(int i=0;i<IdGNodeW.size();i++)
 		myFlux<<IdGNodeW[i]<<" NodeType: "<<NodeTypeW[i]<<" ";
 	myFlux<<std::endl<<"North real nodes: "<<std::endl;
 	for(int i=0;i<IdRNodeN.size();i++)
 		myFlux<<IdRNodeN[i]<<" ";
 	myFlux<<std::endl<<"North Ghost nodes: "<<std::endl;
 	for(int i=0;i<IdGNodeN.size();i++)
 		myFlux<<IdGNodeN[i]<<" NodeType: "<<NodeTypeN[i]<<" ";
 	myFlux<<std::endl<<"South real nodes: "<<std::endl;
 	for(int i=0;i<IdRNodeS.size();i++)
 		myFlux<<IdRNodeS[i]<<" ";
 	myFlux<<std::endl<<"South Ghost nodes: "<<std::endl;
 	for(int i=0;i<IdGNodeS.size();i++)
 		myFlux<<IdGNodeS[i]<<" NodeType: "<<NodeTypeS[i]<<" ";
 	myFlux<<std::endl<<"East real nodes: "<<std::endl;
 	for(int i=0;i<IdRNodeE.size();i++)
 		myFlux<<IdRNodeE[i]<<" ";
 	myFlux<<std::endl<<"East Ghost nodes: "<<std::endl;
 	for(int i=0;i<IdGNodeE.size();i++)
 		myFlux<<IdGNodeE[i]<<" NodeType: "<<NodeTypeE[i]<<" ";

 	myFlux<<std::endl<<std::endl;
 	myFlux<<"South West real nodes: "<<std::endl;

 	for(int i=0;i<IdRNodeSW.size();i++)
 		myFlux<<IdRNodeSW[i]<<" ";
 	myFlux<<std::endl<<"South West Ghost nodes: "<<std::endl;
 	for(int i=0;i<IdGNodeSW.size();i++)
 		myFlux<<IdGNodeSW[i]<<" NodeType: "<<NodeTypeSW[i]<<" ";
 	myFlux<<std::endl<<"North West real nodes: "<<std::endl;
 	for(int i=0;i<IdRNodeNW.size();i++)
 		myFlux<<IdRNodeNW[i]<<" ";
 	myFlux<<std::endl<<"North West Ghost nodes: "<<std::endl;
 	for(int i=0;i<IdGNodeNW.size();i++)
 		myFlux<<IdGNodeNW[i]<<" NodeType: "<<NodeTypeNW[i]<<" ";
 	myFlux<<std::endl<<"South East real nodes: "<<std::endl;
 	for(int i=0;i<IdRNodeSE.size();i++)
 		myFlux<<IdRNodeSE[i]<<" ";
 	myFlux<<std::endl<<"South East Ghost nodes: "<<std::endl;
 	for(int i=0;i<IdGNodeSE.size();i++)
 		myFlux<<IdGNodeSE[i]<<" NodeType: "<<NodeTypeSE[i]<<" ";
 	myFlux<<std::endl<<"North East real nodes: "<<std::endl;
 	for(int i=0;i<IdRNodeNE.size();i++)
 		myFlux<<IdRNodeNE[i]<<" ";
 	myFlux<<std::endl<<"North East Ghost nodes: "<<std::endl;
 	for(int i=0;i<IdGNodeNE.size();i++)
 		myFlux<<IdGNodeNE[i]<<" NodeType: "<<NodeTypeNE[i]<<" ";
 	myFlux<<std::endl<<" Size of NE Ghost: "<< IdGNodeNE.size()<<std::endl;*/
}
void Block2D::Correct_GhostType(int  idNode, NodeType RealNodeType)
{
	if(RealNodeType!=Solid)
		Block2D::ChangeNodeType(idNode, Ghost);
	/*else
		Block2D::ChangeNodeType(idNode, Solid);*/
}
//Remove unused mark nodes
void Block2D::Clear_MarkNode(){
	IdNodeN.clear();	IdNodeS.clear();	IdNodeE.clear();	IdNodeW.clear();
	IdNodeNW.clear();	IdNodeNE.clear();	IdNodeSW.clear();	IdNodeSE.clear();
}
void Block2D::Get_CommNodes(std::vector<int> & IdRNodeN_,std::vector<int> & IdRNodeE_,std::vector<int> & IdRNodeS_,std::vector<int> & IdRNodeW_,
		std::vector<int> & IdGNodeN_,std::vector<int> & IdGNodeE_,std::vector<int> & IdGNodeS_,std::vector<int> & IdGNodeW_,
		std::vector<int> & IdRNodeSW_,std::vector<int> & IdRNodeSE_,std::vector<int> & IdRNodeNW_,std::vector<int> & IdRNodeNE_,
		std::vector<int> & IdGNodeSW_,std::vector<int> & IdGNodeSE_,std::vector<int> & IdGNodeNW_,std::vector<int> & IdGNodeNE_)
{
	IdGNodeN_=IdGNodeN;IdRNodeN_=IdRNodeN;	IdGNodeE_=IdGNodeE;IdRNodeE_=IdRNodeE;	IdGNodeW_=IdGNodeW;IdRNodeW_=IdRNodeW;	IdGNodeS_=IdGNodeS;IdRNodeS_=IdRNodeS;
	IdGNodeSW_=IdGNodeSW;IdRNodeSW_=IdRNodeSW;	IdGNodeSE_=IdGNodeSE;IdRNodeSE_=IdRNodeSE;	IdGNodeNE_=IdGNodeNE;IdRNodeNE_=IdRNodeNE;	IdGNodeNW_=IdGNodeNW;IdRNodeNW_=IdRNodeNW;
}
void Block2D::Set_CommNodes()
{
	//Create the Real and Ghost mark nodes
	IdRNodeW.push_back(IdNodeSW[0]);
	IdRNodeW.insert(IdRNodeW.end(),IdNodeW.begin(),IdNodeW.end());
	IdGNodeW.push_back(Node[IdNodeSW[0]]->Get_connect(3));
	for(int i=0;i<IdNodeW.size();i++)
		IdGNodeW.push_back(Node[IdNodeW[i]]->Get_connect(3));
	if(Node[IdNodeNW[0]]->get_NodeType()!=Ghost)
	{
		IdGNodeNW.push_back(Node[Node[IdNodeNW[0]]->Get_connect(3)]->Get_connect(2));
		IdRNodeNW.push_back(IdNodeNW[0]);
		IdRNodeW.push_back(IdNodeNW[0]);
		IdGNodeW.push_back(Node[IdNodeNW[0]]->Get_connect(3));
	}
	else
	{
		IdGNodeNW.push_back(Node[IdNodeNW[0]]->Get_connect(3));
		IdRNodeNW.push_back(Node[IdNodeNW[0]]->Get_connect(0));
	}

	IdRNodeS.push_back(IdNodeSW[0]);
	IdRNodeS.insert(IdRNodeS.end(),IdNodeS.begin(),IdNodeS.end());
	IdGNodeS.push_back(Node[IdNodeSW[0]]->Get_connect(0));
	for(int i=0;i<IdNodeS.size();i++)
		IdGNodeS.push_back(Node[IdNodeS[i]]->Get_connect(0));
	if(Node[IdNodeSE[0]]->get_NodeType()!=Ghost)
	{
		IdGNodeSE.push_back(Node[Node[IdNodeSE[0]]->Get_connect(0)]->Get_connect(1));
		IdRNodeSE.push_back(IdNodeSE[0]);
		IdRNodeS.push_back(IdNodeSE[0]);
		IdGNodeS.push_back(Node[IdNodeSE[0]]->Get_connect(0));
	}
	else
	{
		IdGNodeSE.push_back(Node[IdNodeSE[0]]->Get_connect(0));
		IdRNodeSE.push_back(Node[IdNodeSE[0]]->Get_connect(3));
	}
	IdRNodeSW.push_back(IdNodeSW[0]);
	IdGNodeSW.push_back(Node[Node[IdNodeSW[0]]->Get_connect(0)]->Get_connect(3));

	if(Node[IdNodeSE[0]]->get_NodeType()!=Ghost)//right side of the domain
	{
		if(Node[IdNodeNE[0]]->get_NodeType()!=Ghost)//top side of the domain
		{
			IdGNodeNE.push_back(Node[Node[IdNodeNE[0]]->Get_connect(1)]->Get_connect(2));
			IdRNodeNE.push_back(IdNodeNE[0]);

			IdRNodeN.push_back(IdNodeNW[0]);
			IdRNodeN.insert(IdRNodeN.end(),IdNodeN.rbegin(),IdNodeN.rend());
			IdRNodeN.push_back(IdNodeNE[0]);

			for(int i=0;i<IdRNodeN.size();i++)
				IdGNodeN.push_back(Node[IdRNodeN[i]]->Get_connect(2));


			IdRNodeE.push_back(IdNodeSE[0]);
			IdRNodeE.insert(IdRNodeE.end(),IdNodeE.begin(),IdNodeE.end());
			IdRNodeE.push_back(IdNodeNE[0]);
			for(int i=0;i<IdRNodeE.size();i++)
				IdGNodeE.push_back(Node[IdRNodeE[i]]->Get_connect(1));
		}
		else
		{
			IdGNodeNE.push_back(Node[IdNodeNE[0]]->Get_connect(1));
			IdRNodeNE.push_back(Node[IdNodeNE[0]]->Get_connect(0));

			IdGNodeN.push_back(IdNodeNW[0]);
			IdGNodeN.insert(IdGNodeN.end(),IdNodeN.rbegin(),IdNodeN.rend());
			IdGNodeN.push_back(IdNodeNE[0]);
			for(int i=0;i<IdGNodeN.size();i++)
				IdRNodeN.push_back(Node[IdGNodeN[i]]->Get_connect(0));


			IdRNodeE.push_back(IdNodeSE[0]);
			IdRNodeE.insert(IdRNodeE.end(),IdNodeE.begin(),IdNodeE.end());
			for(int i=0;i<IdRNodeE.size();i++)
				IdGNodeE.push_back(Node[IdRNodeE[i]]->Get_connect(1));
		}
	}
	else
	{
		if(Node[IdNodeNW[0]]->get_NodeType()!=Ghost)//top side of the domain
		{
			IdGNodeNE.push_back(Node[IdNodeNE[0]]->Get_connect(2));
			IdRNodeNE.push_back(Node[IdNodeNE[0]]->Get_connect(3));

			IdRNodeN.push_back(IdNodeNW[0]);
			IdRNodeN.insert(IdRNodeN.end(),IdNodeN.rbegin(),IdNodeN.rend());
			for(int i=0;i<IdRNodeN.size();i++)
				IdGNodeN.push_back(Node[IdRNodeN[i]]->Get_connect(2));

			IdGNodeE.push_back(IdNodeSE[0]);
			IdGNodeE.insert(IdGNodeE.end(),IdNodeE.begin(),IdNodeE.end());
			IdGNodeE.push_back(IdNodeNE[0]);
			for(int i=0;i<IdGNodeE.size();i++)
				IdRNodeE.push_back(Node[IdGNodeE[i]]->Get_connect(3));
		}
		else
		{
			IdGNodeNE.push_back(IdNodeNE[0]);
			IdRNodeNE.push_back(Node[Node[IdNodeNE[0]]->Get_connect(0)]->Get_connect(3));

			IdGNodeN.push_back(IdNodeNW[0]);
			IdGNodeN.insert(IdGNodeN.end(),IdNodeN.rbegin(),IdNodeN.rend());
			for(int i=0;i<IdGNodeN.size();i++)
				IdRNodeN.push_back(Node[IdGNodeN[i]]->Get_connect(0));

			IdGNodeE.push_back(IdNodeSE[0]);
			IdGNodeE.insert(IdGNodeE.end(),IdNodeE.begin(),IdNodeE.end());
			for(int i=0;i<IdGNodeE.size();i++)
				IdRNodeE.push_back(Node[IdGNodeE[i]]->Get_connect(3));
		}
	}

// if periodic boundary conditions, we need to shift the real by one node to keep the code working in the same way.
// On the periodic boundary conditions, by this way, the distribution will be collide stream in each side of the periodic conditions.
// If no periodic boundary, we remove the unused marking nodes on the boundary of the domain

// Boundaries treatment
	if(periodic[0]==true)
	{
		if(Coord[0]==0)
		{
			for(int i=0;i<IdRNodeW.size();i++)
				IdRNodeW[i]=Node[IdRNodeW[i]]->Get_connect(1);
		}
		if(Coord[0]==Dims[0]-1)
		{
			for(int i=0;i<IdRNodeE.size();i++)
				IdRNodeE[i]=Node[IdRNodeE[i]]->Get_connect(3);
		}
	}
	else
	{
		if(Coord[0]==0)
		{IdRNodeW.clear();IdGNodeW.clear();}
		if(Coord[0]==Dims[0]-1)
		{IdRNodeE.clear();IdGNodeE.clear();}
	}
	if(periodic[1]==true)
	{
		if(Coord[1]==0)
		{
			for(int i=0;i<IdRNodeS.size();i++)
				IdRNodeS[i]=Node[IdRNodeS[i]]->Get_connect(2);
		}
		if(Coord[1]==Dims[1]-1)
		{
			for(int i=0;i<IdRNodeN.size();i++)
				IdRNodeN[i]=Node[IdRNodeN[i]]->Get_connect(0);
		}
	}
	else
	{
		if(Coord[1]==0)
		{IdRNodeS.clear();IdGNodeS.clear();}
		if(Coord[1]==Dims[1]-1)
		{IdRNodeN.clear();IdGNodeN.clear();}
	}
//Corners treament
	if(periodic[0]||periodic[1])
	{
		//Treat the 4 corners of the domain
		if(Coord[0]==0 && Coord[1]==0)
			if(periodic[0]&&periodic[1])
			{
				IdRNodeSW[0]=Node[Node[IdRNodeSW[0]]->Get_connect(1)]->Get_connect(2);
				IdRNodeNW[0]=Node[IdRNodeNW[0]]->Get_connect(1);
				IdRNodeSE[0]=Node[IdRNodeSE[0]]->Get_connect(2);
			}
			else
			{
				if(periodic[0])
				{
					IdRNodeNW[0]=Node[IdRNodeNW[0]]->Get_connect(1);
					IdRNodeSW.clear();IdGNodeSW.clear();
					IdRNodeSE.clear();IdGNodeSE.clear();
				}
				else
				{
					IdRNodeSE[0]=Node[IdRNodeSE[0]]->Get_connect(2);
					IdRNodeSW.clear();IdGNodeSW.clear();
					IdRNodeNW.clear();IdGNodeNW.clear();
				}
			}
		else
			if(Coord[0]==0 && Coord[1]==Dims[1]-1)
						if(periodic[0]&&periodic[1])
						{
							IdRNodeNW[0]=Node[Node[IdRNodeNW[0]]->Get_connect(1)]->Get_connect(0);
							IdRNodeNE[0]=Node[IdRNodeNE[0]]->Get_connect(0);
							IdRNodeSW[0]=Node[IdRNodeSW[0]]->Get_connect(1);
						}
						else
						{
							if(periodic[0])
							{
								IdRNodeSW[0]=Node[IdRNodeSW[0]]->Get_connect(1);
								IdRNodeNW.clear();IdGNodeNW.clear();
								IdRNodeNE.clear();IdGNodeNE.clear();
							}
							else
							{
								IdRNodeNE[0]=Node[IdRNodeNE[0]]->Get_connect(0);
								IdRNodeNW.clear();IdGNodeNW.clear();
								IdRNodeSW.clear();IdGNodeSW.clear();

							}
						}
			else
				if(Coord[0]==Dims[0]-1 && Coord[1]==0)
							if(periodic[0]&&periodic[1])
							{
								IdRNodeSE[0]=Node[Node[IdRNodeSE[0]]->Get_connect(3)]->Get_connect(2);
								IdRNodeNE[0]=Node[IdRNodeNE[0]]->Get_connect(3);
								IdRNodeSW[0]=Node[IdRNodeSW[0]]->Get_connect(2);
							}
							else
							{
								if(periodic[0])
								{
									IdRNodeNE[0]=Node[IdRNodeNE[0]]->Get_connect(3);
									IdRNodeSW.clear();IdGNodeSW.clear();
									IdRNodeSE.clear();IdGNodeSE.clear();
								}
								else
								{
									IdRNodeSW[0]=Node[IdRNodeSW[0]]->Get_connect(2);
									IdRNodeSE.clear();IdGNodeSE.clear();
									IdRNodeNE.clear();IdGNodeNE.clear();

								}
							}
				else
					if(Coord[0]==Dims[0]-1 && Coord[1]==Dims[1]-1)
								if(periodic[0]&&periodic[1])
								{
									IdRNodeNE[0]=Node[Node[IdRNodeNE[0]]->Get_connect(3)]->Get_connect(0);
									IdRNodeNW[0]=Node[IdRNodeNW[0]]->Get_connect(0);
									IdRNodeSE[0]=Node[IdRNodeSE[0]]->Get_connect(3);
								}
								else
								{
									if(periodic[0])
									{
										IdRNodeSE[0]=Node[IdRNodeSE[0]]->Get_connect(3);
										IdRNodeNW.clear();IdGNodeNW.clear();
										IdRNodeNE.clear();IdGNodeNE.clear();
									}
									else
									{
										IdRNodeNW[0]=Node[IdRNodeNW[0]]->Get_connect(0);
										IdRNodeNE.clear();IdGNodeNE.clear();
										IdRNodeSE.clear();IdGNodeSE.clear();

									}
								}
		//Treat the boundary of the domain excluding the corners
					else
					{
						if(Coord[0]==0)
						{
							if(periodic[0])
							{
								IdRNodeNW[0]=Node[IdRNodeNW[0]]->Get_connect(1);
								IdRNodeSW[0]=Node[IdRNodeSW[0]]->Get_connect(1);
							}
							else
							{
								IdRNodeNW.clear();IdGNodeNW.clear();
								IdRNodeSW.clear();IdGNodeSW.clear();
							}
						}
						if(Coord[0]==Dims[0]-1)
						{
							if(periodic[0])
							{
								IdRNodeNE[0]=Node[IdRNodeNE[0]]->Get_connect(3);
								IdRNodeSE[0]=Node[IdRNodeSE[0]]->Get_connect(3);
							}
							else
							{
								IdRNodeNE.clear();IdGNodeNE.clear();
								IdRNodeSE.clear();IdGNodeSE.clear();
							}
						}
						if(Coord[1]==0)
						{
							if(periodic[1])
							{
								IdRNodeSE[0]=Node[IdRNodeSE[0]]->Get_connect(2);
								IdRNodeSW[0]=Node[IdRNodeSW[0]]->Get_connect(2);
							}
							else
							{
								IdRNodeSE.clear();IdGNodeSE.clear();
								IdRNodeSW.clear();IdGNodeSW.clear();
							}
						}
						if(Coord[1]==Dims[1]-1)
						{
							if(periodic[1])
							{
								IdRNodeNE[0]=Node[IdRNodeNE[0]]->Get_connect(0);
								IdRNodeNW[0]=Node[IdRNodeNW[0]]->Get_connect(0);
							}
							else
							{
								IdRNodeNE.clear();IdGNodeNE.clear();
								IdRNodeNW.clear();IdGNodeNW.clear();
							}
						}
					}
	}
	else
	{
		//Remove unused mark nodes on the boundary of the domain
		if (Dims[1]==1)
		{
			IdRNodeNW.clear();IdGNodeNW.clear();
			IdRNodeNE.clear();IdGNodeNE.clear();
			IdRNodeSW.clear();IdGNodeSW.clear();
			IdRNodeSE.clear();IdGNodeSE.clear();
		}
		if (Coord[0]==0)
		{
			IdRNodeNW.clear();IdGNodeNW.clear();
			IdRNodeSW.clear();IdGNodeSW.clear();
		}
		if (Coord[0]==Dims[0]-1)
		{
			IdRNodeNE.clear();IdGNodeNE.clear();
			IdRNodeSE.clear();IdGNodeSE.clear();
		}

		if (Coord[1]==Dims[1]-1)
		{
			IdRNodeNW.clear();IdGNodeNW.clear();
			IdRNodeNE.clear();IdGNodeNE.clear();
		}
		if (Coord[1]==0)
		{
			IdRNodeSW.clear();IdGNodeSW.clear();
			IdRNodeSE.clear();IdGNodeSE.clear();
		}
	}
/*
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	char buffer[50]; // make sure it's big enough
	snprintf(buffer, sizeof(buffer), "ConnectionMarks_%d.txt", rank);
	std::ofstream myFlux;
	myFlux.open(buffer);
	myFlux<<std::endl<<" ****** Removing nodes *******"<<std::endl;
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
		myFlux<<IdRNodeSW[i]<<" x: "<< Node[IdRNodeSW[i]]->get_x()<<" y: "<< Node[IdRNodeSW[i]]->get_y()<<" ";
	myFlux<<std::endl<<"South West Ghost nodes: "<<std::endl;
	for(int i=0;i<IdGNodeSW.size();i++)
		myFlux<<IdGNodeSW[i]<<" x: "<< Node[IdGNodeSW[i]]->get_x()<<" y: "<< Node[IdGNodeSW[i]]->get_y()<<" ";
	myFlux<<std::endl<<"North West real nodes: "<<std::endl;
	for(int i=0;i<IdRNodeNW.size();i++)
		myFlux<<IdRNodeNW[i]<<" x: "<< Node[IdRNodeNW[i]]->get_x()<<" y: "<< Node[IdRNodeNW[i]]->get_y()<<" ";
	myFlux<<std::endl<<"North West Ghost nodes: "<<std::endl;
	for(int i=0;i<IdGNodeNW.size();i++)
		myFlux<<IdGNodeNW[i]<<" x: "<< Node[IdGNodeNW[i]]->get_x()<<" y: "<< Node[IdGNodeNW[i]]->get_y()<<" ";
	myFlux<<std::endl<<"South East real nodes: "<<std::endl;
	for(int i=0;i<IdRNodeSE.size();i++)
		myFlux<<IdRNodeSE[i]<<" x: "<< Node[IdRNodeSE[i]]->get_x()<<" y: "<< Node[IdRNodeSE[i]]->get_y()<<" ";
	myFlux<<std::endl<<"South East Ghost nodes: "<<std::endl;
	for(int i=0;i<IdGNodeSE.size();i++)
		myFlux<<IdGNodeSE[i]<<" x: "<< Node[IdGNodeSE[i]]->get_x()<<" y: "<< Node[IdGNodeSE[i]]->get_y()<<" ";
	myFlux<<std::endl<<"North East real nodes: "<<std::endl;
	for(int i=0;i<IdRNodeNE.size();i++)
		myFlux<<IdRNodeNE[i]<<" x: "<< Node[IdRNodeNE[i]]->get_x()<<" y: "<< Node[IdRNodeNE[i]]->get_y()<<" ";
	myFlux<<std::endl<<"North East Ghost nodes: "<<std::endl;
	for(int i=0;i<IdGNodeNE.size();i++)
		myFlux<<IdGNodeNE[i]<<" x: "<< Node[IdGNodeNE[i]]->get_x()<<" y: "<< Node[IdGNodeNE[i]]->get_y()<<" ";
		*/
}
void Block2D::Mark1stLayerSolid(){
	for(int i=0;i<NodeArrays.NodeWall.size();i++)
		switch(NodeArrays.NodeWall[i].Get_BcNormal())
		{
		case 1:
			NodeArrays.Solid1stLayer.push_back(NodeArrays.NodeWall[i].Get_connect()[3]);
			break;
		case 2:
			NodeArrays.Solid1stLayer.push_back(NodeArrays.NodeWall[i].Get_connect()[4]);
			break;
		case 3:
			NodeArrays.Solid1stLayer.push_back(NodeArrays.NodeWall[i].Get_connect()[1]);
			break;
		case 4:
			NodeArrays.Solid1stLayer.push_back(NodeArrays.NodeWall[i].Get_connect()[2]);
			break;
		}
	for(int i=0;i<NodeArrays.NodeCorner.size();i++)
		switch(NodeArrays.NodeCorner[i].Get_BcNormal())
		{
		case 5:
			if(NodeArrays.NodeCorner[i].Get_CornerType()==Convex)
			{
				NodeArrays.Solid1stLayer.push_back(NodeArrays.NodeCorner[i].Get_connect()[7]);
				NodeArrays.CornerConvex.push_back(i);
			}
			else
			{
				NodeArrays.Solid1stLayer.push_back(NodeArrays.NodeCorner[i].Get_connect()[7]);
				NodeArrays.Solid1stLayer.push_back(NodeArrays.NodeCorner[i].Get_connect()[3]);
				NodeArrays.Solid1stLayer.push_back(NodeArrays.NodeCorner[i].Get_connect()[4]);
				NodeArrays.CornerConcave.push_back(i);
			}
			break;
		case 6:
			if(NodeArrays.NodeCorner[i].Get_CornerType()==Convex)
			{
				NodeArrays.Solid1stLayer.push_back(NodeArrays.NodeCorner[i].Get_connect()[8]);
				NodeArrays.CornerConvex.push_back(i);
			}
			else
			{
				NodeArrays.Solid1stLayer.push_back(NodeArrays.NodeCorner[i].Get_connect()[8]);
				NodeArrays.Solid1stLayer.push_back(NodeArrays.NodeCorner[i].Get_connect()[1]);
				NodeArrays.Solid1stLayer.push_back(NodeArrays.NodeCorner[i].Get_connect()[4]);
				NodeArrays.CornerConcave.push_back(i);
			}			break;
		case 7:
			if(NodeArrays.NodeCorner[i].Get_CornerType()==Convex)
			{
				NodeArrays.Solid1stLayer.push_back(NodeArrays.NodeCorner[i].Get_connect()[5]);
				NodeArrays.CornerConvex.push_back(i);
			}
			else
			{
				NodeArrays.Solid1stLayer.push_back(NodeArrays.NodeCorner[i].Get_connect()[5]);
				NodeArrays.Solid1stLayer.push_back(NodeArrays.NodeCorner[i].Get_connect()[1]);
				NodeArrays.Solid1stLayer.push_back(NodeArrays.NodeCorner[i].Get_connect()[2]);
				NodeArrays.CornerConcave.push_back(i);
			}			break;
		case 8:
			if(NodeArrays.NodeCorner[i].Get_CornerType()==Convex)
			{
				NodeArrays.Solid1stLayer.push_back(NodeArrays.NodeCorner[i].Get_connect()[6]);
				NodeArrays.CornerConvex.push_back(i);
			}
			else
			{
				NodeArrays.Solid1stLayer.push_back(NodeArrays.NodeCorner[i].Get_connect()[6]);
				NodeArrays.Solid1stLayer.push_back(NodeArrays.NodeCorner[i].Get_connect()[3]);
				NodeArrays.	Solid1stLayer.push_back(NodeArrays.NodeCorner[i].Get_connect()[2]);
				NodeArrays.CornerConcave.push_back(i);
			}			break;
		}
	std::sort( NodeArrays.Solid1stLayer.begin(), NodeArrays.Solid1stLayer.end() );
	NodeArrays.Solid1stLayer.erase( std::unique( NodeArrays.Solid1stLayer.begin(), NodeArrays.Solid1stLayer.end() ), NodeArrays.Solid1stLayer.end() );

}
