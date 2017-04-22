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
	for(int i=0;i<CellArray.size();i++) delete CellArray[i];
	for(int i=0;i<GhostCellArraytmp.size();i++) delete GhostCellArraytmp[i];
	for(int i=0;i<Node.size();i++) delete Node[i];
	for(int i=0;i<Node_Ghosttmp.size();i++) delete Node_Ghosttmp[i];
	for(int i=0;i<Node_Solidtmp.size();i++) delete Node_Solidtmp[i];
	for(int i=0;i<Node_tmp.size();i++) delete Node_tmp[i];

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
void Block2D::AddBlock(int TotalNumberCells_x,int TotalNumberCells_y,int dims[2],int periodic_[2],int coord[2],int NyCell_G,int Nx_begin,int Ny_begin,int Nx_last,int Ny_last,NodeType bCIn[4],int x_last,int y_last,int start_GlobalNodes,int end_GlobalNodes, int start_GlobalElems, int end_GlobalElems)
{
	verbous=false;
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	for(int i=0;i<4;i++)
		{bC[i]=bCIn[i];
		//std::cout<<"Processor: "<<rank<<" Type: "<<bC[i]<<" in: "<<i<<std::endl;;
		}

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
//Remove Second layer of ghost generated
	//	RemoveWrongGhostCells();

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
void Block2D::ConvertSolidInGhostToGhostnode(){
	for(int i=NbRealNodes;i<Node.size();i++)
		if(Node[i]->get_NodeType()==Solid)
			ChangeNodeType(i,SolidGhost);
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
///Deos not work. it need to renumbering node and cells...
void Block2D::RemoveWrongGhostCells(){
	std::vector<int> CellToBeRemoved,NodeToBeRemoved;
	//Get Cells and node to be removed
	for(int i=0;i<IdGhostCell.size();i++)
	{
		for(int j=0;j<4;j++)
			if(IsGhostConnectToGhost(CellArray[IdGhostCell[i]]->Get_NodeNumber(j)-1))
			{
				CellToBeRemoved.push_back(IdGhostCell[i]);
				NodeToBeRemoved.push_back(CellArray[IdGhostCell[i]]->Get_NodeNumber(j)-1);
			}
	}
	//Sort to remove in order and remove duplicates for cells and nodes
	std::sort( CellToBeRemoved.begin(), CellToBeRemoved.end() );
	CellToBeRemoved.erase( std::unique( CellToBeRemoved.begin(), CellToBeRemoved.end() ), CellToBeRemoved.end() );

	std::sort( NodeToBeRemoved.begin(), NodeToBeRemoved.end() );
	NodeToBeRemoved.erase( std::unique( NodeToBeRemoved.begin(), NodeToBeRemoved.end() ), NodeToBeRemoved.end() );


	if(!CellArray.empty())
	{
		for(int i=CellToBeRemoved.size()-1;i>=0;i--)
		{
			Remove_CellConnections(CellToBeRemoved[i]);
		}
		for(int i=CellToBeRemoved.size()-1;i>=0;i--)
			CellArray.erase(CellArray.begin()+CellToBeRemoved[i]);
	}

	if(!Node.empty())
	{
		for(int i=NodeToBeRemoved.size()-1;i>=0;i--)
		{
			Remove_NodeConnections(NodeToBeRemoved[i]);
		}

		for(int i=NodeToBeRemoved.size()-1;i>=0;i--)
			Node.erase(Node.begin()+NodeToBeRemoved[i]);
	}

}
void Block2D::Remove_CellConnections(int cellid){
	int idCellConnect,idFaceConnect;
	//Remove connection to neighbours
	for(int i=0;i<4;i++)
	{
		idFaceConnect=CellArray[cellid]->Get_Connect(i)[0];
		idCellConnect=CellArray[cellid]->Get_Connect(i)[1];
		CellArray[idCellConnect]->Remove_Connect(idCellConnect,idFaceConnect);
	}
	//Remove connection of the cell
	for(int i=0;i<4;i++)
	{
		CellArray[cellid]->Remove_Connect(cellid,i);
	}
}
void Block2D::Remove_NodeConnections(int nodeid){
	int idnodeConnect;
	//Remove connection to neighbours
	for(int i=0;i<4;i++)
	{
		idnodeConnect=Node[nodeid]->Get_connect(i);
		Node[idnodeConnect]->Remove_Connect(OppositeDirection(i));
	}
	for(int i=0;i<4;i++)
	{
		Node[nodeid]->Remove_Connect(i);
	}
}
int Block2D::OppositeDirection(int direction){
	switch (direction)
	{
	case 0:
		return 2;
		break;
	case 1:
		return 3;
		break;
	case 2:
		return 0;
		break;
	case 3:
		return 1;
		break;
	default:
		std::cerr<<"Wrong Connection ask. Default opposite Connection is north node "<<std::endl;
		return 2;
		break;
	}
}
bool Block2D::IsGhostConnectToGhost(int idx){
	int ndDiffGhost=0;
	for(unsigned int j=1;j<9;j++)
	{

		if(Node[Connect_lowOrder(idx,j)]->get_NodeType()!=Ghost)
			ndDiffGhost++;
	}
	if (ndDiffGhost>0)
		return false;
	else
		return true;
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

	signed short int x_tmp=(signed short int)x[NodeNumber];
	signed short int y_tmp=(signed short int)y[NodeNumber];
	unsigned int Connect_N_tmp=Node[NodeNumber]->Get_connect(2)+1;
	unsigned int Connect_S_tmp=Node[NodeNumber]->Get_connect(0)+1;
	unsigned int Connect_W_tmp=Node[NodeNumber]->Get_connect(3)+1;
	unsigned int Connect_E_tmp=Node[NodeNumber]->Get_connect(1)+1;
	unsigned int NbVelocity_tmp=Node[NodeNumber]->Get_NbVelocity();
	NodeType OldNodeType_=Node[NodeNumber]->get_NodeType();
	int index=Node[NodeNumber]->Get_index();

	delete Node[NodeNumber];

	switch(NewNodeType_)
	{
	case Interior:
		Node[NodeNumber]=new NodeInterior2D(x_tmp, y_tmp);
		break;
	case Solid:
		Node[NodeNumber]=new NodeSolid2D(x_tmp, y_tmp);
		break;
	case SolidGhost:
		Node[NodeNumber]=new NodeGhost2D(x_tmp, y_tmp);
		Node[NodeNumber]->Set_NodeType(SolidGhost);
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
	case Symmetry:
		Node[NodeNumber]=new NodeSymmetry2D(x_tmp, y_tmp);
		break;
	case Corner:
		Node[NodeNumber]=new NodeCorner2D(x_tmp, y_tmp);
		break;
	case ConcaveCorner:
		Node[NodeNumber]=new NodeCorner2D(x_tmp, y_tmp);
		break;
	case ConvexCorner:
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
	Node[NodeNumber]->Set_Index(index);
}


void Block2D::ChangeCoord(int NodeNumber,signed short int const x_, signed short int const y_)
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
//void Block2D::GenerateSolid(Parameters &Param, int &ndSolidNode,int &firstSolid)
void Block2D::GenerateSolid(Parameters &Param)
{
	bool testSolid;
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	nx=Param.Get_Nx();ny=Param.Get_Ny();
	PtrParametersUserMesh=&Param;
	UserMesh::SetUserMeshVariables();

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
/*			if(ndSolidNode==0)
				firstSolid=i;*/
			//Store the number of Solid node
//			ndSolidNode++;
			//store id of solid nodes
			if(i<NbRealNodes)
	//		if(Node[i]->get_NodeType()!=Ghost)
				IdSolidNode.push_back(i);
			//create the solid node}
			Block2D::ChangeNodeType(i,Solid);//convert node to solid
		}
		//Extend the solid to the ghost of the domain for the automatic wall detection functions. Otherwise, the wall can be consider as a corner
		if(Node[i]->get_y()==0 && (testSolid||Node[i]->get_NodeType()==Wall))//Convert ghost node to solid at the border of the domain (bottom)
		{
			Block2D::ChangeNodeType(Node[i]->Get_connect(0),Solid);//convert node to solid
			if(Node[i]->get_x()==0)//Correct Global corner
				Block2D::ChangeNodeType(Node[Node[i]->Get_connect(0)]->Get_connect(3),Solid);//convert node to solid
			if(Node[i]->get_x()==Param.Get_Nx())//Correct Global corner
				Block2D::ChangeNodeType(Node[Node[i]->Get_connect(0)]->Get_connect(1),Solid);//convert node to solid
		}
		if(Node[i]->get_y()==Param.Get_Ny() && (testSolid||Node[i]->get_NodeType()==Wall))//Convert ghost node to solid at the border of the domain (North)
		{
			Block2D::ChangeNodeType(Node[i]->Get_connect(2),Solid);//convert node to solid
			if(Node[i]->get_x()==0)//Correct Global corner
				Block2D::ChangeNodeType(Node[Node[i]->Get_connect(2)]->Get_connect(3),Solid);//convert node to solid
			if(Node[i]->get_x()==Param.Get_Nx())//Correct Global corner
				Block2D::ChangeNodeType(Node[Node[i]->Get_connect(2)]->Get_connect(1),Solid);//convert node to solid
		}
		if(Node[i]->get_x()==0 && (testSolid||Node[i]->get_NodeType()==Wall))//Convert ghost node to solid at the border of the domain (West)
			Block2D::ChangeNodeType(Node[i]->Get_connect(3),Solid);//convert node to solid
		if(Node[i]->get_x()==Param.Get_Nx() && (testSolid||Node[i]->get_NodeType()==Wall))//Convert ghost node to solid at the border of the domain (East)
			Block2D::ChangeNodeType(Node[i]->Get_connect(1),Solid);//convert node to solid

		//Save nodes which are not solid. To reduce memory consumption, only nodes after the first solid are stored
		/*if(ndSolidNode>0 && Solid==false)
			Node_tmp.push_back(Node[i]);
		ndSolidNodetoremove.push_back(ndSolidNode);*/
		//}
	}
	//Correct_Solid_Ghost();
	std::cout<<"Processor: "<<rank<<" Number of Solid: "<<IdSolidNode.size()<<" Number of Computation nodes: "<<NbRealNodes-IdSolidNode.size()<<std::endl;
//	if(rank==0)
//		for(int i=0;i<IdCornerConvextmp.size();i++)
//			std::cout<<"x: "<<Node[IdCornerConvextmp[i]]->get_x()<<" y: "<<Node[IdCornerConvextmp[i]]->get_y()<<std::endl;
}
void Block2D::SetSolidBoundaries(){
//Detect boundaries
	for(int i=0;i<IdSolidNode.size();i++)
	{
		if(DetectSolidBoundaries(IdSolidNode[i]))
			IdSolidBc.push_back(IdSolidNode[i]);
	}
//Create wall and corners
	for (int i=0;i<IdWalltmp.size();i++)
		CreateWall(IdWalltmp[i]);
	for (int i=0;i<IdSpecialWalltmp.size();i++)
		CreateSpecialWall(IdSpecialWalltmp[i]);
	for (int i=0;i<IdCornerConcavetmp.size();i++)
		CreateCornerConcave(IdCornerConcavetmp[i]);
	for (int i=0;i<IdCornerConvextmp.size();i++)
		CreateCornerConvex(IdCornerConvextmp[i]);
//Remove wrong solid
/*	for (int i=0;i<IDRemoveSolidtmp.size();i++)
	{
		Block2D::ChangeNodeType(IDRemoveSolidtmp[i],Interior);
		std::cout<<"ID : "<<IDRemoveSolidtmp[i]<<std::endl;
	}*/

	Node_tmp.clear();
	Node_Solidtmp.clear();
}
void Block2D::RemoveUnphysicalSolid(int &nbTotalSolidRemoved,int &nbTotalSolidadded){
	std::vector<int>::iterator position;
//	int nbTotalSolidRemoved=1;
//	int nbTotalSolidadded=1;
//	int maxclean=10;int count=0;
//	int rank;
//	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
//	while( count<maxclean && (nbTotalSolidRemoved>0 || nbTotalSolidadded>0))//(count<maxclean ||(nbTotalSolidRemoved>0 && nbTotalSolidadded>0))// &&
//	{
		for(int i=0;i<IdSolidNode.size();i++)
		{
			DetectUnphysicalSolid(IdSolidNode[i]);
		}

	//**********Remove wrong solid************
		//Remove duplicate and sort the array
		std::sort( IDRemoveSolidtmp.begin(), IDRemoveSolidtmp.end() );
		IDRemoveSolidtmp.erase( std::unique( IDRemoveSolidtmp.begin(), IDRemoveSolidtmp.end() ), IDRemoveSolidtmp.end() );

		for (int i=0;i<IDRemoveSolidtmp.size();i++)
		{
			position = std::find(IdSolidNode.begin(), IdSolidNode.end(), IDRemoveSolidtmp[i]);
			if (position != IdSolidNode.end()) // == myVector.end() means the element was not found
				IdSolidNode.erase(position);
		// Convert node to interior
			Block2D::ChangeNodeType(IDRemoveSolidtmp[i],Interior);
		}

		std::sort( IDRemoveSolidToBctmp.begin(), IDRemoveSolidToBctmp.end() );
		IDRemoveSolidToBctmp.erase( std::unique( IDRemoveSolidToBctmp.begin(), IDRemoveSolidToBctmp.end() ), IDRemoveSolidToBctmp.end() );

		for (int i=0;i<IDRemoveSolidToBctmp.size();i++)
		{
			position = std::find(IdSolidNode.begin(), IdSolidNode.end(), IDRemoveSolidToBctmp[i]);
			if (position != IdSolidNode.end()) // == myVector.end() means the element was not found
				IdSolidNode.erase(position);
		// Convert node to boundary condition

			if(Node[IDRemoveSolidToBctmp[i]]->get_x()==0)
				{
					Block2D::ChangeNodeType(IDRemoveSolidToBctmp[i],bC[3]);
				//std::cout<<"Processor: "<<rank<<" Type: "<<bC[3]<<" in: "<<3<<std::endl;
				}
			if(Node[IDRemoveSolidToBctmp[i]]->get_x()==nx)
				{
					Block2D::ChangeNodeType(IDRemoveSolidToBctmp[i],bC[1]);
				//	std::cout<<"Processor: "<<rank<<" Type: "<<bC[1]<<" in: "<<1<<std::endl;
				}
			if(Node[IDRemoveSolidToBctmp[i]]->get_y()==0)
				{
					Block2D::ChangeNodeType(IDRemoveSolidToBctmp[i],bC[0]);
				//	std::cout<<"Processor: "<<rank<<" Type: "<<bC[0]<<" in: "<<0<<std::endl;
				}
			if(Node[IDRemoveSolidToBctmp[i]]->get_y()==ny)
				{
					Block2D::ChangeNodeType(IDRemoveSolidToBctmp[i],bC[2]);
				//	std::cout<<"Processor: "<<rank<<" Type: "<<bC[2]<<" in: "<<2<<std::endl;
				}
		}
		nbTotalSolidRemoved=IDRemoveSolidtmp.size()+IDRemoveSolidToBctmp.size();
//		std::cout<<"Processor: "<<rank<<"Number of solid removed: "<<nbTotalSolidRemoved<<std::endl;
		//clear the array
		IDRemoveSolidtmp.clear();
		IDRemoveSolidToBctmp.clear();
		//**********Add solids************
		//Remove duplicate and sort the array
		std::sort( IDAddSolidtmp.begin(), IDAddSolidtmp.end() );
		IDAddSolidtmp.erase( std::unique( IDAddSolidtmp.begin(), IDAddSolidtmp.end() ), IDAddSolidtmp.end() );
		// Convert node to interior
		for (int i=0;i<IDAddSolidtmp.size();i++)
			{Block2D::ChangeNodeType(IDAddSolidtmp[i],Solid);}
		nbTotalSolidadded=IDAddSolidtmp.size();
//		std::cout<<"Processor: "<<rank<<"Number of solid added: "<<nbTotalSolidadded<<std::endl;
		IdSolidNode.insert(IdSolidNode.end(),IDAddSolidtmp.begin(), IDAddSolidtmp.end());
		//clear the array
		IDAddSolidtmp.clear();
//		count++;
//	}
	//std::cout<<"number of cleaning: "<<count<<std::endl;
}
void Block2D::DetectUnphysicalSolid(int & nodeID){
	int nbSolidDirect=0;
	int nbSolidDiagonal=0;
	int nbWall=0;
	for(unsigned int j=1;j<5;j++)
	{
		if(Node[Connect_lowOrder(nodeID,j)]->get_NodeType()==Solid)
		{	nbSolidDirect++;	}
	}
	for(unsigned int j=5;j<9;j++)
	{
		if(Node[Connect_lowOrder(nodeID,j)]->get_NodeType()==Solid)
		{	nbSolidDiagonal++;	}
	}
	if((nbSolidDirect+nbSolidDiagonal)<8)
	{
		if(nbSolidDirect==3 &&nbSolidDiagonal<2)
		{
			//Add solid on the diagonal for wall boundary in order to avoid only solid nodes
			unsigned int j=0;unsigned int k=0;unsigned int l=0;
			for(j=2;j<5;j++)
			{
				k=j-1;
				l=j+3;
					if(Node[Connect_lowOrder(nodeID,j)]->get_NodeType()==Solid &&
							Node[Connect_lowOrder(nodeID,k)]->get_NodeType()==Solid &&
							Node[Connect_lowOrder(nodeID,l)]->get_NodeType()!=Solid)
						IDAddSolidtmp.push_back(Connect_lowOrder(nodeID,l));
			}
			j=1;k=4;l=8;
			if(Node[Connect_lowOrder(nodeID,j)]->get_NodeType()==Solid &&
					Node[Connect_lowOrder(nodeID,k)]->get_NodeType()==Solid &&
					Node[Connect_lowOrder(nodeID,l)]->get_NodeType()!=Solid)
				IDAddSolidtmp.push_back(Connect_lowOrder(nodeID,l));
		}
		else if(nbSolidDirect==4 && nbSolidDiagonal<3)
		{
			unsigned int j=5;
			while(nbSolidDiagonal<3)
			{
				if(Node[Connect_lowOrder(nodeID,j)]->get_NodeType()!=Solid)
				{
					IDAddSolidtmp.push_back(Connect_lowOrder(nodeID,j));
					nbSolidDiagonal++;
				}
				j++;
				if(j>8) {nbSolidDiagonal=3; std::cerr<<"Error for adding solids"<<std::endl;}
			}

		}
		else if(nbSolidDirect==2)
		{
			if(Node[nodeID]->get_x()==0||Node[nodeID]->get_x()==nx
			 ||Node[nodeID]->get_y()==0||Node[nodeID]->get_y()==ny)
				IDRemoveSolidToBctmp.push_back(nodeID);


			else if(	Node[Node[nodeID]->Get_connect(0)]->get_NodeType()==Node[Node[nodeID]->Get_connect(1)]->get_NodeType()||
				Node[Node[nodeID]->Get_connect(0)]->get_NodeType()==Node[Node[nodeID]->Get_connect(3)]->get_NodeType()||
				Node[Node[nodeID]->Get_connect(2)]->get_NodeType()==Node[Node[nodeID]->Get_connect(1)]->get_NodeType()||
				Node[Node[nodeID]->Get_connect(2)]->get_NodeType()==Node[Node[nodeID]->Get_connect(3)]->get_NodeType()	)
			{
				if(nbSolidDiagonal==0)
					IDRemoveSolidtmp.push_back(nodeID);
			}
			else
				IDRemoveSolidtmp.push_back(nodeID);

		}
		else if(nbSolidDirect==0 ||nbSolidDirect==1)//nbSolidDirect==0 or ==1
		{
			IDRemoveSolidtmp.push_back(nodeID);

		}
	}
}

bool Block2D::DetectSolidBoundaries(int & nodeID){
	int nbSolidDirect=0;
	int nbSolidDiagonal=0;
	int nbWall=0;
	for(unsigned int j=1;j<5;j++)
	{
		if(Node[Connect_lowOrder(nodeID,j)]->get_NodeType()==Solid)
		{	nbSolidDirect++;	}
	}
	for(unsigned int j=5;j<9;j++)
	{
		if(Node[Connect_lowOrder(nodeID,j)]->get_NodeType()==Solid)
		{	nbSolidDiagonal++;	}
	}

	if((nbSolidDirect+nbSolidDiagonal)<8)
	{
		if(nbSolidDirect==3)
		{
			if(Node[nodeID]->get_x()==0||Node[nodeID]->get_x()==nx
			 ||Node[nodeID]->get_y()==0||Node[nodeID]->get_y()==ny)
			{
				if(Node[nodeID]->get_x()==0)
				{
					if(Node[Node[nodeID]->Get_connect(1)]->get_NodeType()!=Solid)
						IdWalltmp.push_back(nodeID);
					else if((Node[Node[nodeID]->Get_connect(2)]->get_NodeType()==Solid || Node[Node[nodeID]->Get_connect(2)]->get_NodeType()==Wall) &&
							(Node[Node[nodeID]->Get_connect(0)]->get_NodeType()==Solid || Node[Node[nodeID]->Get_connect(0)]->get_NodeType()==Wall))
						IdCornerConcavetmp.push_back(nodeID);
					else
						IdSpecialWalltmp.push_back(nodeID);
				}
				else if(Node[nodeID]->get_x()==nx)
				{
					if(Node[Node[nodeID]->Get_connect(3)]->get_NodeType()!=Solid)
						IdWalltmp.push_back(nodeID);
					else if((Node[Node[nodeID]->Get_connect(2)]->get_NodeType()==Solid || Node[Node[nodeID]->Get_connect(2)]->get_NodeType()==Wall) &&
							(Node[Node[nodeID]->Get_connect(0)]->get_NodeType()==Solid || Node[Node[nodeID]->Get_connect(0)]->get_NodeType()==Wall))
						IdCornerConcavetmp.push_back(nodeID);
					else
						IdSpecialWalltmp.push_back(nodeID);
				}
				else if(Node[nodeID]->get_y()==0)
				{
					if(Node[Node[nodeID]->Get_connect(2)]->get_NodeType()!=Solid)
						IdWalltmp.push_back(nodeID);
					else if((Node[Node[nodeID]->Get_connect(1)]->get_NodeType()==Solid || Node[Node[nodeID]->Get_connect(1)]->get_NodeType()==Wall) &&
							(Node[Node[nodeID]->Get_connect(3)]->get_NodeType()==Solid || Node[Node[nodeID]->Get_connect(3)]->get_NodeType()==Wall))
						IdCornerConcavetmp.push_back(nodeID);
					else
						IdSpecialWalltmp.push_back(nodeID);
				}
				else if(Node[nodeID]->get_y()==ny)
				{
					if(Node[Node[nodeID]->Get_connect(0)]->get_NodeType()!=Solid)
						IdWalltmp.push_back(nodeID);
					else if((Node[Node[nodeID]->Get_connect(1)]->get_NodeType()==Solid || Node[Node[nodeID]->Get_connect(1)]->get_NodeType()==Wall) &&
							(Node[Node[nodeID]->Get_connect(3)]->get_NodeType()==Solid || Node[Node[nodeID]->Get_connect(3)]->get_NodeType()==Wall))
						IdCornerConcavetmp.push_back(nodeID);
					else
						IdSpecialWalltmp.push_back(nodeID);
				}
/*
				for(unsigned int j=1;j<5;j++)
					if(Node[Connect_lowOrder(nodeID,j)]->get_NodeType()==Wall)
						nbWall++;
				if(nbWall>0)
					IdCornerConcavetmp.push_back(nodeID);
				else
					IdSpecialWalltmp.push_back(nodeID);
					*/
			}

			else
				IdWalltmp.push_back(nodeID);
		}
		else if(nbSolidDirect==4)
		{
			IdCornerConcavetmp.push_back(nodeID);
		}
		else if(nbSolidDirect==2)
		{
			IdCornerConvextmp.push_back(nodeID);
		}
		else
		{
			std::cout<<"Solid boundary not found. Number of Solid connected: "<<nbSolidDirect+nbSolidDiagonal<<std::endl;
			IDRemoveSolidtmp.push_back(nodeID);
		}
		return true;

	}
	else
		return false;

}

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


}
void Block2D::CreateCorner(int &nodeID){
	Block2D::ChangeNodeType(nodeID,Corner);
	DefinedCornerType(nodeID);
	IdBoundaries.push_back(nodeID);
}
void Block2D::CreateCornerConcave(int &nodeID){
	Block2D::ChangeNodeType(nodeID,Corner);
	Node[nodeID]->Set_CornerType(ConcaveCorner);
	IdBoundaries.push_back(nodeID);
}
void Block2D::CreateCornerConvex(int &nodeID){
	Block2D::ChangeNodeType(nodeID,Corner);
	Node[nodeID]->Set_CornerType(ConvexCorner);
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
bool Block2D::IsSpecialWallAtGlocalCorner(int idx){

	unsigned int Connect_S_tmp=Node[idx]->Get_connect(0);
	unsigned int Connect_E_tmp=Node[idx]->Get_connect(1);
	unsigned int Connect_N_tmp=Node[idx]->Get_connect(2);
	unsigned int Connect_W_tmp=Node[idx]->Get_connect(3);

	if(Node[Connect_N_tmp]->get_NodeType()==Wall || Node[Connect_S_tmp]->get_NodeType()==Wall || Node[Connect_W_tmp]->get_NodeType()==Wall || Node[Connect_E_tmp]->get_NodeType()==Wall)
		return true;
	else
		return false;
}
void Block2D::reorganizeNodeByType(){

	unsigned int Connect_N_tmp=0;
	unsigned int Connect_S_tmp=0;
	unsigned int Connect_W_tmp=0;
	unsigned int Connect_E_tmp=0;
	unsigned int NbVelocity_tmp=0;
	int count=0;

	for (int i=0;i<Node.size();i++)
	{
		NodeArrays.TypeOfNode[i]=Node[i]->get_NodeType();
		signed short int x_tmp=(signed short int)Node[i]->get_x();
		signed short int y_tmp=(signed short int)Node[i]->get_y();
		Connect_N_tmp=Node[i]->Get_connect(2)+1;
		Connect_S_tmp=Node[i]->Get_connect(0)+1;
		Connect_W_tmp=Node[i]->Get_connect(3)+1;
		Connect_E_tmp=Node[i]->Get_connect(1)+1;
		NbVelocity_tmp=Node[i]->Get_NbVelocity();

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
			NodeArrays.NodeGhost[NodeArrays.NodeGhost.size()-1].Set_GhostType(InteriorGhostType);
			NodeArrays.NodeGhost[NodeArrays.NodeGhost.size()-1].Set_Index(i);
			break;
		case GlobalCorner:
			if(IsSpecialWallAtGlocalCorner(i))
			{
				NodeArrays.TypeOfNode[i]=SpecialWall;
				NodeArrays.NodeIndexByType[i]=NodeArrays.NodeSpecialWall.size();
				NodeArrays.NodeSpecialWall.push_back(NodeWall2D(x_tmp,y_tmp));
				NodeArrays.NodeSpecialWall[NodeArrays.NodeSpecialWall.size()-1].Set_Connect(Connect_N_tmp,Connect_S_tmp,Connect_W_tmp,Connect_E_tmp);
				NodeArrays.NodeSpecialWall[NodeArrays.NodeSpecialWall.size()-1].Set_NbVelocity(NbVelocity_tmp);
				NodeArrays.NodeSpecialWall[NodeArrays.NodeSpecialWall.size()-1].Set_Index(i);
				NodeArrays.NodeSpecialWall[NodeArrays.NodeSpecialWall.size()-1].Set_NodeType(SpecialWall);
			}
			else
			{
				NodeArrays.NodeIndexByType[i]=NodeArrays.NodeGlobalCorner.size();
				NodeArrays.NodeGlobalCorner.push_back(NodeCorner2D(x_tmp,y_tmp,Node[i]->Get_RhoDef(),Node[i]->Get_UDef()));
				NodeArrays.NodeGlobalCorner[NodeArrays.NodeGlobalCorner.size()-1].Set_Connect(Connect_N_tmp,Connect_S_tmp,Connect_W_tmp,Connect_E_tmp);
				NodeArrays.NodeGlobalCorner[NodeArrays.NodeGlobalCorner.size()-1].Set_NbVelocity(NbVelocity_tmp);
				NodeArrays.NodeGlobalCorner[NodeArrays.NodeGlobalCorner.size()-1].Set_Index(i);
				NodeArrays.NodeGlobalCorner[NodeArrays.NodeGlobalCorner.size()-1].Set_CornerType(Concave);
				NodeArrays.NodeGlobalCorner[NodeArrays.NodeGlobalCorner.size()-1].Set_NodeType(GlobalCorner);
			}
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
			NodeArrays.NodeGhost[NodeArrays.NodeGhost.size()-1].Set_GhostType(SolidGhostType);
			NodeArrays.NodeGhost[NodeArrays.NodeGhost.size()-1].Set_Index(i);
			count++;
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
	if(!Node.empty())
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
	else
	{
		std::cerr<<"Node Array is empty."<<std::endl;
		return direction;
	}
}
void Block2D::Set_BcNormal()
{
	for (int i=0;i<NodeArrays.NodeWall.size();i++)
	{
			NodeArrays.NodeWall[i].Set_BcNormal(Get_BcNormal(NodeArrays.NodeWall[i]));
	}
	for (int i=0;i<NodeArrays.NodeSpecialWall.size();i++)
	{
			NodeArrays.NodeSpecialWall[i].Set_BcNormal(Get_BcNormal_SpecialWall(NodeArrays.NodeSpecialWall[i]));
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
	int tmpdirection[4];
	int intTmpReturn=0;

	if(nodeIn.get_NodeType()==GlobalCorner)
	{
		for (unsigned int i=5;i<9;i++)
		{
			if (NodeArrays.TypeOfNode[nodeIn.Get_connect()[i]]==Interior )
			{
				intTmpReturn=i;
				nbinterior++;
			}
		}
		if(nbinterior>1 || intTmpReturn==0)
		{
			std::cerr<<"Wrong corner normal detection for Global corner. "<<"x: "<< nodeIn.get_x()<<" y: "<< nodeIn.get_y()<<std::endl;
		}
	}
	else
	if(nodeIn.Get_CornerType()==Concave)
	{
		for (unsigned int i=5;i<9;i++)
		{
			if (NodeArrays.TypeOfNode[nodeIn.Get_connect()[i]]==Interior ||NodeArrays.TypeOfNode[nodeIn.Get_connect()[i]]==Ghost
					|| NodeArrays.TypeOfNode[nodeIn.Get_connect()[i]]==Velocity || NodeArrays.TypeOfNode[nodeIn.Get_connect()[i]]==Pressure
					|| NodeArrays.TypeOfNode[nodeIn.Get_connect()[i]]==Symmetry || NodeArrays.TypeOfNode[nodeIn.Get_connect()[i]]==Periodic)
			{
				intTmpReturn=i;
				nbinterior++;
			}
		}
		if(nbinterior>1 || intTmpReturn==0)
		{
			std::cerr<<"Wrong corner normal detection for Concave corner. "<<"x: "<< nodeIn.get_x()<<" y: "<< nodeIn.get_y()<<std::endl;
			for (unsigned int i=5;i<9;i++)
			{
				if (NodeArrays.TypeOfNode[nodeIn.Get_connect()[i]]==Interior ||NodeArrays.TypeOfNode[nodeIn.Get_connect()[i]]==Ghost)
				{
					std::cerr<<"Direction "<<i<<" is an interior"<<std::endl;
				}
			}

		}
	}
	else
	{
		for(int i=0;i<4;i++)
			tmpdirection[i]=0;
		int count=0;
		for (unsigned int i=1;i<5;i++)
		{
			if (NodeArrays.TypeOfNode[nodeIn.Get_connect()[i]]!=Interior && NodeArrays.TypeOfNode[nodeIn.Get_connect()[i]]!=Ghost
					&& NodeArrays.TypeOfNode[nodeIn.Get_connect()[i]]!=Velocity && NodeArrays.TypeOfNode[nodeIn.Get_connect()[i]]!=Pressure
					&& NodeArrays.TypeOfNode[nodeIn.Get_connect()[i]]!=Symmetry && NodeArrays.TypeOfNode[nodeIn.Get_connect()[i]]!=Periodic)
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
		if(intTmpReturn<5 )
			std::cerr<<"Wrong corner normal detection for Convex corner. "<<"x: "<< nodeIn.get_x()<<" y: "<< nodeIn.get_y()<<" Direction detected: "<<intTmpReturn<<std::endl;
	}

	return intTmpReturn;
}

int Block2D::Get_BcNormal(NodeWall2D & nodeIn)
{
	intTmpReturn=0;
	if(nodeIn.get_x()==0)
	{
		intTmpReturn=1;
	}
	else if(nodeIn.get_x()==nx)
	{
		intTmpReturn=3;
	}
	else if(nodeIn.get_y()==0)
	{
		intTmpReturn=2;
	}
	else if(nodeIn.get_y()==ny)
	{
		intTmpReturn=4;
	}
	else
	{
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
		if(intTmpReturn==0)
			for (unsigned int i=1;i<5;i++)
				if (NodeArrays.TypeOfNode[nodeIn.Get_connect()[i]]==Pressure ||NodeArrays.TypeOfNode[nodeIn.Get_connect()[i]]==Velocity||NodeArrays.TypeOfNode[nodeIn.Get_connect()[i]]==Symmetry ||NodeArrays.TypeOfNode[nodeIn.Get_connect()[i]]==Periodic )
				{
					intTmpReturn=i;
				}
	}
	if(intTmpReturn==0)
		std::cerr<<"Wrong normal detection for a Wall node. "<<"x: "<< nodeIn.get_x()<<" y: "<< nodeIn.get_y()<<" Direction detected: "<<intTmpReturn<<std::endl;

	return intTmpReturn;
}
int Block2D::Get_BcNormal_SpecialWall(NodeWall2D & nodeIn)

{NodeType check;
	intTmpReturn=0;
	//treat global corner convert to special node
	if((nodeIn.get_x()==0 && nodeIn.get_y()==0)  || (nodeIn.get_x()==0 && nodeIn.get_y()==ny)  || (nodeIn.get_x()==nx && nodeIn.get_y()==0)  || (nodeIn.get_x()==nx && nodeIn.get_y()==ny))
	{
		if((nodeIn.get_x()==0 && nodeIn.get_y()==0))
		{
			intTmpReturn=5;
			check=NodeArrays.TypeOfNode[nodeIn.Get_connect()[3]];
			check=NodeArrays.TypeOfNode[nodeIn.Get_connect()[4]];
			check=NodeArrays.TypeOfNode[nodeIn.Get_connect()[7]];
		}
		else if((nodeIn.get_x()==nx && nodeIn.get_y()==0))
		{
			intTmpReturn=6;
			check=NodeArrays.TypeOfNode[nodeIn.Get_connect()[1]];
			check=NodeArrays.TypeOfNode[nodeIn.Get_connect()[4]];
			check=NodeArrays.TypeOfNode[nodeIn.Get_connect()[8]];
		}
		else if((nodeIn.get_x()==nx && nodeIn.get_y()==ny))
		{
			intTmpReturn=7;
			check=NodeArrays.TypeOfNode[nodeIn.Get_connect()[1]];
			check=NodeArrays.TypeOfNode[nodeIn.Get_connect()[2]];
			check=NodeArrays.TypeOfNode[nodeIn.Get_connect()[5]];
		}
		else if((nodeIn.get_x()==0 && nodeIn.get_y()==ny))
		{
			intTmpReturn=8;
			check=NodeArrays.TypeOfNode[nodeIn.Get_connect()[2]];
			check=NodeArrays.TypeOfNode[nodeIn.Get_connect()[3]];
			check=NodeArrays.TypeOfNode[nodeIn.Get_connect()[6]];
		}
	}
	else
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
	}
	if(intTmpReturn==0)
		std::cerr<<"Wrong normal detection for a Special wall node. "<<"x: "<< nodeIn.get_x()<<" y: "<< nodeIn.get_y()<<" Direction detected: "<<intTmpReturn<<std::endl;

	return intTmpReturn;
}

int Block2D::Get_BcNormal(NodeSymmetry2D & nodeIn)
{
	intTmpReturn=0;
	for (unsigned int i=1;i<5;i++)
		if (NodeArrays.TypeOfNode[nodeIn.Get_connect()[i]]==Interior)
			intTmpReturn=i;
	if(intTmpReturn==0)
		for (unsigned int i=1;i<5;i++)
			if (NodeArrays.TypeOfNode[nodeIn.Get_connect()[i]]==Wall||NodeArrays.TypeOfNode[nodeIn.Get_connect()[i]]==Corner ||
					NodeArrays.TypeOfNode[nodeIn.Get_connect()[i]]==ConvexCorner||NodeArrays.TypeOfNode[nodeIn.Get_connect()[i]]==ConcaveCorner)
				intTmpReturn=i;
	if(intTmpReturn==0)
		std::cerr<<"Wrong normal detection for a Symmetry node. "<<"x: "<< nodeIn.get_x()<<" y: "<< nodeIn.get_y()<<" Direction detected: "<<intTmpReturn<<std::endl;

	return intTmpReturn;
}
int Block2D::Get_BcNormal(NodePeriodic2D & nodeIn)
{
	intTmpReturn=0;
	for (unsigned int i=1;i<5;i++)
		if (NodeArrays.TypeOfNode[nodeIn.Get_connect()[i]]==Interior)
			intTmpReturn=i;
	if(intTmpReturn==0)
		for (unsigned int i=1;i<5;i++)
			if (NodeArrays.TypeOfNode[nodeIn.Get_connect()[i]]==Wall||NodeArrays.TypeOfNode[nodeIn.Get_connect()[i]]==Corner ||
					NodeArrays.TypeOfNode[nodeIn.Get_connect()[i]]==ConvexCorner||NodeArrays.TypeOfNode[nodeIn.Get_connect()[i]]==ConcaveCorner)
				intTmpReturn=i;
	if(intTmpReturn==0)
		std::cerr<<"Wrong normal detection for a Periodic node. "<<"x: "<< nodeIn.get_x()<<" y: "<< nodeIn.get_y()<<" Direction detected: "<<intTmpReturn<<std::endl;

	return intTmpReturn;
}
int Block2D::Get_BcNormal(NodePressure2D & nodeIn)
{
	intTmpReturn=0;
	for (unsigned int i=1;i<5;i++)
		if (NodeArrays.TypeOfNode[nodeIn.Get_connect()[i]]==Interior)
			intTmpReturn=i;
	if(intTmpReturn==0)
		for (unsigned int i=1;i<5;i++)
			if (NodeArrays.TypeOfNode[nodeIn.Get_connect()[i]]==Wall||NodeArrays.TypeOfNode[nodeIn.Get_connect()[i]]==Corner ||
					NodeArrays.TypeOfNode[nodeIn.Get_connect()[i]]==ConvexCorner||NodeArrays.TypeOfNode[nodeIn.Get_connect()[i]]==ConcaveCorner)
				intTmpReturn=i;
	if(intTmpReturn==0)
		std::cerr<<"Wrong normal detection for a Pressure node. "<<"x: "<< nodeIn.get_x()<<" y: "<< nodeIn.get_y()<<" Direction detected: "<<intTmpReturn<<std::endl;
	return intTmpReturn;
}
int Block2D::Get_BcNormal(NodeVelocity2D & nodeIn)
{
	intTmpReturn=0;
	for (unsigned int i=1;i<5;i++)
		if (NodeArrays.TypeOfNode[nodeIn.Get_connect()[i]]==Interior)
			intTmpReturn=i;
	if(intTmpReturn==0)
		for (unsigned int i=1;i<5;i++)
			if (NodeArrays.TypeOfNode[nodeIn.Get_connect()[i]]==Wall||NodeArrays.TypeOfNode[nodeIn.Get_connect()[i]]==Corner||
					NodeArrays.TypeOfNode[nodeIn.Get_connect()[i]]==ConvexCorner||NodeArrays.TypeOfNode[nodeIn.Get_connect()[i]]==ConcaveCorner)
				intTmpReturn=i;
	if(intTmpReturn==0)
		std::cerr<<"Wrong normal detection for a Velocity node. "<<"x: "<< nodeIn.get_x()<<" y: "<< nodeIn.get_y()<<" Direction detected: "<<intTmpReturn<<std::endl;

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


}
void Block2D::Get_SolidGhost(std::vector<int> & GhostSolidIdN,std::vector<int> & GhostSolidIdE,std::vector<int> & GhostSolidIdS,std::vector<int> & GhostSolidIdW,
		  std::vector<int> & GhostSolidIdSW,std::vector<int> & GhostSolidIdSE,std::vector<int> & GhostSolidIdNW,std::vector<int> & GhostSolidIdNE)
{
	GhostSolidIdN.clear();GhostSolidIdE.clear();GhostSolidIdS.clear();GhostSolidIdW.clear();
	GhostSolidIdSW.clear(); GhostSolidIdSE.clear(); GhostSolidIdNW.clear(); GhostSolidIdNE.clear();

	if(!NodeArrays.TypeOfNode.empty()){
	// Get Node type of the sides of the 2D block for the real nodes
		for(int i=0;i<IdRNodeN.size();i++)
		{
			if(NodeArrays.TypeOfNode[IdRNodeN[i]]==Solid)
			{
				GhostSolidIdN.push_back(i);
			}
			else
				GhostSolidIdN.push_back(-1);
		}
		for(int i=0;i<IdRNodeE.size();i++)
		{
			if(NodeArrays.TypeOfNode[IdRNodeE[i]]==Solid)
			{
				GhostSolidIdE.push_back(i);
			}
			else
				GhostSolidIdE.push_back(-1);
		}
		for(int i=0;i<IdRNodeS.size();i++)
		{
			if(NodeArrays.TypeOfNode[IdRNodeS[i]]==Solid)
			{
				GhostSolidIdS.push_back(i);
			}
			else
				GhostSolidIdS.push_back(-1);
		}
		for(int i=0;i<IdRNodeW.size();i++)
		{
			if(NodeArrays.TypeOfNode[IdRNodeW[i]]==Solid)
			{
				GhostSolidIdW.push_back(i);
			}
			else
				GhostSolidIdW.push_back(-1);
		}
		// Get Node type of the corners of the 2D block
		for(int i=0;i<IdRNodeSW.size();i++)
		{
			if(NodeArrays.TypeOfNode[IdRNodeSW[i]]==Solid)
			{
				GhostSolidIdSW.push_back(i);
			}
			else
				GhostSolidIdSW.push_back(-1);
		}
		for(int i=0;i<IdRNodeSE.size();i++)
		{
			if(NodeArrays.TypeOfNode[IdRNodeSE[i]]==Solid)
			{
				GhostSolidIdSE.push_back(i);
			}
			else
				GhostSolidIdSE.push_back(-1);
		}
		for(int i=0;i<IdRNodeNW.size();i++)
		{
			if(NodeArrays.TypeOfNode[IdRNodeNW[i]]==Solid)
			{
				GhostSolidIdNW.push_back(i);
			}
			else
				GhostSolidIdNW.push_back(-1);
		}
		for(int i=0;i<IdRNodeNE.size();i++)
		{
			if(NodeArrays.TypeOfNode[IdRNodeNE[i]]==Solid)
			{
				GhostSolidIdNE.push_back(i);
			}
			else
				GhostSolidIdNE.push_back(-1);
		}
	}
	else
		std::cerr<<"NodeArrays is empty."<<std::endl;

}
void Block2D::Get_GhostFirstLayer(std::vector<int> & GhostSolidIdN,std::vector<int> & GhostSolidIdE,std::vector<int> & GhostSolidIdS,std::vector<int> & GhostSolidIdW,
		  std::vector<int> & GhostSolidIdSW,std::vector<int> & GhostSolidIdSE,std::vector<int> & GhostSolidIdNW,std::vector<int> & GhostSolidIdNE)
{
	GhostSolidIdN.clear();GhostSolidIdE.clear();GhostSolidIdS.clear();GhostSolidIdW.clear();
	GhostSolidIdSW.clear(); GhostSolidIdSE.clear(); GhostSolidIdNW.clear(); GhostSolidIdNE.clear();

	if(!NodeArrays.TypeOfNode.empty()){
	// Get Node type of the sides of the 2D block for the real nodes
		for(int i=0;i<IdRNodeN.size();i++)
		{
			if(NodeArrays.TypeOfNode[IdRNodeN[i]]==Solid)
			{
				if(NodeArrays.NodeSolid[NodeArrays.NodeIndexByType[IdRNodeN[i]]].IsFirstLayer())
					GhostSolidIdN.push_back(i);
				else
					GhostSolidIdN.push_back(-1);
			}
			else
				GhostSolidIdN.push_back(-1);
		}
		for(int i=0;i<IdRNodeE.size();i++)
		{
			if(NodeArrays.TypeOfNode[IdRNodeE[i]]==Solid)
			{
				if(NodeArrays.NodeSolid[NodeArrays.NodeIndexByType[IdRNodeE[i]]].IsFirstLayer())
					GhostSolidIdE.push_back(i);
				else
					GhostSolidIdE.push_back(-1);
			}
			else
				GhostSolidIdE.push_back(-1);
		}
		for(int i=0;i<IdRNodeS.size();i++)
		{
			if(NodeArrays.TypeOfNode[IdRNodeS[i]]==Solid)
			{
				if(NodeArrays.NodeSolid[NodeArrays.NodeIndexByType[IdRNodeS[i]]].IsFirstLayer())
					GhostSolidIdS.push_back(i);
				else
					GhostSolidIdS.push_back(-1);
			}
			else
				GhostSolidIdS.push_back(-1);
		}
		for(int i=0;i<IdRNodeW.size();i++)
		{
			if(NodeArrays.TypeOfNode[IdRNodeW[i]]==Solid)
			{
				if(NodeArrays.NodeSolid[NodeArrays.NodeIndexByType[IdRNodeW[i]]].IsFirstLayer())
					GhostSolidIdW.push_back(i);
				else
					GhostSolidIdW.push_back(-1);
			}
			else
				GhostSolidIdW.push_back(-1);
		}
		// Get Node type of the corners of the 2D block
		for(int i=0;i<IdRNodeSW.size();i++)
		{
			if(NodeArrays.TypeOfNode[IdRNodeSW[i]]==Solid)
			{
				if(NodeArrays.NodeSolid[NodeArrays.NodeIndexByType[IdRNodeSW[i]]].IsFirstLayer())
					GhostSolidIdSW.push_back(i);
				else
					GhostSolidIdSW.push_back(-1);
			}
			else
				GhostSolidIdSW.push_back(-1);
		}
		for(int i=0;i<IdRNodeSE.size();i++)
		{
			if(NodeArrays.TypeOfNode[IdRNodeSE[i]]==Solid)
			{
				if(NodeArrays.NodeSolid[NodeArrays.NodeIndexByType[IdRNodeSE[i]]].IsFirstLayer())
					GhostSolidIdSE.push_back(i);
				else
					GhostSolidIdSE.push_back(-1);
			}
			else
				GhostSolidIdSE.push_back(-1);
		}
		for(int i=0;i<IdRNodeNW.size();i++)
		{
			if(NodeArrays.TypeOfNode[IdRNodeNW[i]]==Solid)
			{
				if(NodeArrays.NodeSolid[NodeArrays.NodeIndexByType[IdRNodeNW[i]]].IsFirstLayer())
					GhostSolidIdNW.push_back(i);
				else
					GhostSolidIdNW.push_back(-1);
			}
			else
				GhostSolidIdNW.push_back(-1);
		}
		for(int i=0;i<IdRNodeNE.size();i++)
		{
			if(NodeArrays.TypeOfNode[IdRNodeNE[i]]==Solid)
			{
				if(NodeArrays.NodeSolid[NodeArrays.NodeIndexByType[IdRNodeNE[i]]].IsFirstLayer())
					GhostSolidIdNE.push_back(i);
				else
					GhostSolidIdNE.push_back(-1);
			}
			else
				GhostSolidIdNE.push_back(-1);
		}
	}
	else
		std::cerr<<"NodeArrays is empty."<<std::endl;

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

}
void Block2D::Set_GhostTypeAsGhost(std::vector<int> & NodeTypeN,std::vector<int> & NodeTypeE,std::vector<int> & NodeTypeS,std::vector<int> & NodeTypeW,
		  std::vector<int> & NodeTypeSW,std::vector<int> & NodeTypeSE,std::vector<int> & NodeTypeNW,std::vector<int> & NodeTypeNE)
{

	// Get Node type of the sides of the 2D block
	if(IdGNodeN.size()>0)
	for(int i=0;i<IdGNodeN.size();i++)
		Correct_SolidGhostType(IdGNodeN[i], (NodeType)NodeTypeN[i]);
	if(IdGNodeE.size()>0)
	for(int i=0;i<IdGNodeE.size();i++)
		Correct_SolidGhostType(IdGNodeE[i], (NodeType)NodeTypeE[i]);
	if(IdGNodeS.size()>0)
	for(int i=0;i<IdGNodeS.size();i++)
		Correct_SolidGhostType(IdGNodeS[i], (NodeType)NodeTypeS[i]);
	if(IdGNodeW.size()>0)
	for(int i=0;i<IdGNodeW.size();i++)
		Correct_SolidGhostType(IdGNodeW[i], (NodeType)NodeTypeW[i]);
// Get Node type of the corners of the 2D block
	if(IdGNodeSW.size()>0)
	for(int i=0;i<IdGNodeSW.size();i++)
		Correct_SolidGhostType(IdGNodeSW[i], (NodeType)NodeTypeSW[i]);
	if(IdGNodeSE.size()>0)
	for(int i=0;i<IdGNodeSE.size();i++)
		Correct_SolidGhostType(IdGNodeSE[i], (NodeType)NodeTypeSE[i]);
	if(IdGNodeNW.size()>0)
	for(int i=0;i<IdGNodeNW.size();i++)
		Correct_SolidGhostType(IdGNodeNW[i], (NodeType)NodeTypeNW[i]);
	if(IdGNodeNE.size()>0)
	for(int i=0;i<IdGNodeNE.size();i++)
		Correct_SolidGhostType(IdGNodeNE[i], (NodeType)NodeTypeNE[i]);

}
void Block2D::Remove_SolidTypeInCommunicatorsForRealNodes(){

	std::vector<int> RealNodeType,RealIdToSaved,GhostIdToSaved;
	// Get Node type of the sides of the 2D block
	for(int i=0;i<IdRNodeN.size();i++)
		RealNodeType.push_back((int)NodeArrays.TypeOfNode[IdRNodeN[i]]);
	Remove_SolidTypeInCommunicatorForRealNodes(RealNodeType,IdRNodeN,RealIdToSaved);
	SetCommunicatorinSolidForRealNodes(SolidIdRNodeN,RealIdToSaved);
	RealNodeType.clear();RealIdToSaved.clear();

	for(int i=0;i<IdRNodeE.size();i++)
		RealNodeType.push_back((int)NodeArrays.TypeOfNode[IdRNodeE[i]]);
	Remove_SolidTypeInCommunicatorForRealNodes(RealNodeType,IdRNodeE,RealIdToSaved);
	SetCommunicatorinSolidForRealNodes(SolidIdRNodeE,RealIdToSaved);
	RealNodeType.clear();RealIdToSaved.clear();

	for(int i=0;i<IdRNodeS.size();i++)
		RealNodeType.push_back((int)NodeArrays.TypeOfNode[IdRNodeS[i]]);
	Remove_SolidTypeInCommunicatorForRealNodes(RealNodeType,IdRNodeS,RealIdToSaved);
	SetCommunicatorinSolidForRealNodes(SolidIdRNodeS,RealIdToSaved);
	RealNodeType.clear();RealIdToSaved.clear();

	for(int i=0;i<IdRNodeW.size();i++)
		RealNodeType.push_back((int)NodeArrays.TypeOfNode[IdRNodeW[i]]);
	Remove_SolidTypeInCommunicatorForRealNodes(RealNodeType,IdRNodeW,RealIdToSaved);
	SetCommunicatorinSolidForRealNodes(SolidIdRNodeW,RealIdToSaved);
	RealNodeType.clear();RealIdToSaved.clear();
// Get Node type of the corners of the 2D block
	for(int i=0;i<IdRNodeSW.size();i++)
		RealNodeType.push_back((int)NodeArrays.TypeOfNode[IdRNodeSW[i]]);
	Remove_SolidTypeInCommunicatorForRealNodes(RealNodeType,IdRNodeSW,RealIdToSaved);
	SetCommunicatorinSolidForRealNodes(SolidIdRNodeSW,RealIdToSaved);
	RealNodeType.clear();RealIdToSaved.clear();

	for(int i=0;i<IdRNodeSE.size();i++)
		RealNodeType.push_back((int)NodeArrays.TypeOfNode[IdRNodeSE[i]]);
	Remove_SolidTypeInCommunicatorForRealNodes(RealNodeType,IdRNodeSE,RealIdToSaved);
	SetCommunicatorinSolidForRealNodes(SolidIdRNodeSE,RealIdToSaved);
	RealNodeType.clear();RealIdToSaved.clear();

	for(int i=0;i<IdRNodeNW.size();i++)
		RealNodeType.push_back((int)NodeArrays.TypeOfNode[IdRNodeNW[i]]);
	Remove_SolidTypeInCommunicatorForRealNodes(RealNodeType,IdRNodeNW,RealIdToSaved);
	SetCommunicatorinSolidForRealNodes(SolidIdRNodeNW,RealIdToSaved);
	RealNodeType.clear();RealIdToSaved.clear();

	for(int i=0;i<IdRNodeNE.size();i++)
		RealNodeType.push_back((int)NodeArrays.TypeOfNode[IdRNodeNE[i]]);
	Remove_SolidTypeInCommunicatorForRealNodes(RealNodeType,IdRNodeNE,RealIdToSaved);
	SetCommunicatorinSolidForRealNodes(SolidIdRNodeNE,RealIdToSaved);
	RealNodeType.clear();RealIdToSaved.clear();


}
void Block2D::Remove_SolidTypeInCommunicatorForRealNodes(std::vector<int> & RealNodeType,
		  std::vector<int> & RealNodeId,std::vector<int> & RealIdToBeSaved){
	std::vector<int> RealIdToBeRemoved;
	RealIdToBeSaved.clear();
	for(int i=0;i<RealNodeType.size();i++)
	{
		if ((NodeType)RealNodeType[i]==Solid)
			RealIdToBeRemoved.push_back(i);
	}

	if(!RealNodeId.empty())
	for(int i=RealIdToBeRemoved.size()-1;i>=0;i--)
	{
		RealIdToBeSaved.push_back(RealNodeId[RealIdToBeRemoved[i]]);
		RealNodeId.erase(RealNodeId.begin()+RealIdToBeRemoved[i]);
	}
}
void Block2D::SetCommunicatorinSolidForRealNodes(std::vector<int> & RealNodeId,std::vector<int> & RealIdToBeSaved){

//	RealNodeId=RealIdToBeSaved;
	//remove node who are not in the first layer
	int idx=0;
	RealNodeId.clear();
	for(int i=0;i<RealIdToBeSaved.size();i++)
	{
		idx=NodeArrays.NodeIndexByType[RealIdToBeSaved[i]];
		if(idx<NodeArrays.NodeSolid.size())
		{
			if(NodeArrays.NodeSolid[idx].IsFirstLayer())
							 RealNodeId.push_back(NodeArrays.NodeSolid[idx].Get_index());
		}
	}

}
void Block2D::Remove_SolidTypeInCommunicatorsForGhostNodes(std::vector<int> & GhostSolidIdN,std::vector<int> & GhostSolidIdE,std::vector<int> & GhostSolidIdS,std::vector<int> & GhostSolidIdW,
		  std::vector<int> & GhostSolidIdSW,std::vector<int> & GhostSolidIdSE,std::vector<int> & GhostSolidIdNW,std::vector<int> & GhostSolidIdNE,
		  std::vector<int> & GhostFirstLayerSolidIdN,std::vector<int> & GhostFirstLayerSolidIdE,std::vector<int> & GhostFirstLayerSolidIdS,std::vector<int> & GhostFirstLayerSolidIdW,
		  std::vector<int> & GhostFirstLayerSolidIdSW,std::vector<int> & GhostFirstLayerSolidIdSE,std::vector<int> & GhostFirstLayerSolidIdNW,std::vector<int> & GhostFirstLayerSolidIdNE){
	std::vector<int> GhostIdToSaved;

	// Get Node type of the sides of the 2D block
	Remove_SolidTypeInCommunicatorForGhostNodes(GhostSolidIdN,GhostFirstLayerSolidIdN,IdGNodeN,GhostIdToSaved);
	SetCommunicatorinSolidForGhostNodes(SolidIdGNodeN,GhostIdToSaved);
	GhostIdToSaved.clear();

	Remove_SolidTypeInCommunicatorForGhostNodes(GhostSolidIdE,GhostFirstLayerSolidIdE,IdGNodeE,GhostIdToSaved);
	SetCommunicatorinSolidForGhostNodes(SolidIdGNodeE,GhostIdToSaved);
	GhostIdToSaved.clear();

	Remove_SolidTypeInCommunicatorForGhostNodes(GhostSolidIdS,GhostFirstLayerSolidIdS,IdGNodeS,GhostIdToSaved);
	SetCommunicatorinSolidForGhostNodes(SolidIdGNodeS,GhostIdToSaved);
	GhostIdToSaved.clear();

	Remove_SolidTypeInCommunicatorForGhostNodes(GhostSolidIdW,GhostFirstLayerSolidIdW,IdGNodeW,GhostIdToSaved);
	SetCommunicatorinSolidForGhostNodes(SolidIdGNodeW,GhostIdToSaved);
	GhostIdToSaved.clear();
// Get Node type of the corners of the 2D block
	Remove_SolidTypeInCommunicatorForGhostNodes(GhostSolidIdSW,GhostFirstLayerSolidIdSW,IdGNodeSW,GhostIdToSaved);
	SetCommunicatorinSolidForGhostNodes(SolidIdGNodeSW,GhostIdToSaved);
	GhostIdToSaved.clear();

	Remove_SolidTypeInCommunicatorForGhostNodes(GhostSolidIdSE,GhostFirstLayerSolidIdSE,IdGNodeSE,GhostIdToSaved);
	SetCommunicatorinSolidForGhostNodes(SolidIdGNodeSE,GhostIdToSaved);
	GhostIdToSaved.clear();

	Remove_SolidTypeInCommunicatorForGhostNodes(GhostSolidIdNW,GhostFirstLayerSolidIdNW,IdGNodeNW,GhostIdToSaved);
	SetCommunicatorinSolidForGhostNodes(SolidIdGNodeNW,GhostIdToSaved);
	GhostIdToSaved.clear();

	Remove_SolidTypeInCommunicatorForGhostNodes(GhostSolidIdNE,GhostFirstLayerSolidIdNE,IdGNodeNE,GhostIdToSaved);
	SetCommunicatorinSolidForGhostNodes(SolidIdGNodeNE,GhostIdToSaved);
	GhostIdToSaved.clear();


}
void Block2D::Remove_SolidTypeInCommunicatorForGhostNodes(std::vector<int> & GhostSolidId,std::vector<int> & GhostSolidFirstLayerId,
		  std::vector<int> & GhostNodeId,std::vector<int> & GhostIdToBeSaved){
	std::vector<int> GhostIdToBeRemoved,GhostSolidIdToBeRemoved;
	GhostIdToBeSaved.clear();
	//mark solids
	for(int i=0;i<GhostSolidId.size();i++)
	{
		if (GhostSolidId[i]>-1)
			GhostIdToBeRemoved.push_back(GhostSolidId[i]);
	}
	//Filtering the first layer
	for(int i=0;i<GhostSolidFirstLayerId.size();i++)
	{
		if (GhostSolidFirstLayerId[i]>-1)
			GhostSolidIdToBeRemoved.push_back(GhostSolidFirstLayerId[i]);
	}

	if(!GhostNodeId.empty())
	{
		//Save the solid needed it
		for(int i=GhostSolidIdToBeRemoved.size()-1;i>=0;i--)
		{
			GhostIdToBeSaved.push_back(GhostNodeId[GhostSolidIdToBeRemoved[i]]);
		}
		//Remove solids
		for(int i=GhostIdToBeRemoved.size()-1;i>=0;i--)
		{
			GhostNodeId.erase(GhostNodeId.begin()+GhostIdToBeRemoved[i]);
		}
	}
}
void Block2D::SetCommunicatorinSolidForGhostNodes(std::vector<int> & GhostNodeId,std::vector<int> & GhostIdToBeSaved){
	GhostNodeId=GhostIdToBeSaved;
}

void Block2D::Correct_GhostType(int  idNode, NodeType RealNodeType)
{
	if(RealNodeType!=Solid)
	{
		if(RealNodeType==Interior)
		{
			Block2D::ChangeNodeType(idNode, Ghost);
		}
		else
		{
			Block2D::ChangeNodeType(idNode, SolidGhost);
		}
	}
	else
	{
		if(Node[idNode]->get_NodeType()!=Solid)
		{
			Block2D::ChangeNodeType(idNode, Solid);
		}
	}
}
void Block2D::Correct_SolidGhostType(int  idNode, NodeType RealNodeType)
{
	if(RealNodeType!=Solid)
	{
		if(RealNodeType==Interior)
		{
			Block2D::ChangeNodeType(idNode, Ghost);
		}
		else
		{
			Block2D::ChangeNodeType(idNode, SolidGhost);
		}
	}
	else
	{
		if(Node[idNode]->get_NodeType()!=Solid)
		{
			Block2D::ChangeNodeType(idNode, Solid);
		}
		else
		{
			Block2D::ChangeNodeType(idNode, SolidGhost);
		}
	}
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
void Block2D::Get_SolidCommNodes(std::vector<int> & IdRNodeN_,std::vector<int> & IdRNodeE_,std::vector<int> & IdRNodeS_,std::vector<int> & IdRNodeW_,
		std::vector<int> & IdGNodeN_,std::vector<int> & IdGNodeE_,std::vector<int> & IdGNodeS_,std::vector<int> & IdGNodeW_,
		std::vector<int> & IdRNodeSW_,std::vector<int> & IdRNodeSE_,std::vector<int> & IdRNodeNW_,std::vector<int> & IdRNodeNE_,
		std::vector<int> & IdGNodeSW_,std::vector<int> & IdGNodeSE_,std::vector<int> & IdGNodeNW_,std::vector<int> & IdGNodeNE_)
{
	IdGNodeN_=SolidIdGNodeN;IdRNodeN_=SolidIdRNodeN;	IdGNodeE_=SolidIdGNodeE;IdRNodeE_=SolidIdRNodeE;	IdGNodeW_=SolidIdGNodeW;IdRNodeW_=SolidIdRNodeW;	IdGNodeS_=SolidIdGNodeS;IdRNodeS_=SolidIdRNodeS;
	IdGNodeSW_=SolidIdGNodeSW;IdRNodeSW_=SolidIdRNodeSW;	IdGNodeSE_=SolidIdGNodeSE;IdRNodeSE_=SolidIdRNodeSE;	IdGNodeNE_=SolidIdGNodeNE;IdRNodeNE_=SolidIdRNodeNE;	IdGNodeNW_=SolidIdGNodeNW;IdRNodeNW_=SolidIdRNodeNW;
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
	for(int i=0;i<NodeArrays.Solid1stLayer.size();i++)
	{
		if(NodeArrays.TypeOfNode[NodeArrays.Solid1stLayer[i]]==Solid)
			NodeArrays.NodeSolid[NodeArrays.NodeIndexByType[NodeArrays.Solid1stLayer[i]]].Set_FirstLayer(true);
	}
}

void  Block2D::RemoveSolid(){
	std::vector<Cell2D*> CellArraytmp;
	std::vector<Node2D*> NodeArraytmp,NodeFirstWallLayerArraytmp;
	std::map<int,int> MapCellArraytmp,MapNodeArraytmp,MapFirstWallLayerArraytmp;
	std::vector<int> IdCellSolid,IdCellGhost,IdCellFirstWallLayer,IdCellFirstGhostLayer;
	int countSolid,countGhost;
	for(int i=0;i<CellArray.size();i++)
	{
		countSolid=0;countGhost=0;
		for(int j=0;j<4;j++)
		{
			if(Node[CellArray[i]->Get_NodeNumber(j)]->get_NodeType()==Solid)
				countSolid++;
			if(Node[CellArray[i]->Get_NodeNumber(j)]->get_NodeType()==Ghost)
				countGhost++;
		}
		if(countSolid>0)
			if(countSolid>=4)
			{
				IdCellSolid.push_back(i);
				//Remove connection between the Solid cell and others cell (can be also a solid cell)
				for(int j=0;j<4;j++)
				{
					//Remove the connection in the solid cell
					CellArray[i]->Set_Connect(j,j,i);
					//Remove the connection with the connected cells
					CellArray[CellArray[i]->Get_Connect(j)[1]]->Set_Connect(CellArray[i]->Get_Connect(j)[0],CellArray[i]->Get_Connect(j)[0],CellArray[i]->Get_Connect(j)[1]);

				}
			}
			else
				IdCellFirstWallLayer.push_back(i);
		if(countGhost>0)
			if(countGhost>=4)
			{
				IdCellGhost.push_back(i);
				for(int j=0;j<4;j++)
				{
					//Remove the connection in the solid cell
					CellArray[i]->Set_Connect(j,j,i);
					//Remove the connection with the connected cells
					CellArray[CellArray[i]->Get_Connect(j)[1]]->Set_Connect(CellArray[i]->Get_Connect(j)[0],CellArray[i]->Get_Connect(j)[0],CellArray[i]->Get_Connect(j)[1]);
				}
			}
			else
				IdCellFirstGhostLayer.push_back(i);
		//Save computational cells
		if(countSolid==0&&countGhost==0)
		{
			MapCellArraytmp[CellArraytmp.size()]=i;
			CellArraytmp.push_back(CellArray[i]);
		}
	}
	int nbSolidToRemove=0; int nbGhostToRemove=0; int notSolid=0; int nbGhostConnect=0;
	//Removing Solid nodes in the node array
	for(int i=0;i<Node.size();i++)
	{
		if(Node[i]->get_NodeType()==Ghost)
		{
			nbGhostConnect=0;
			for(unsigned int j=1;j<9;j++)
			{
				//Check if ghost node is used or not
				if(Node[Connect_lowOrder(i,j)]->get_NodeType()==Ghost || Node[Connect_lowOrder(i,j)]->get_NodeType()==Solid)
					nbGhostConnect++;
			}
			//Keep ghost
			if(nbGhostConnect!=8)
			{
				//map the new node array to old node array to keep node informations
				MapNodeArraytmp[i]=NodeArraytmp.size();
				//Update the index by removing solid
				Node[i]->Set_Index(Node[i]->Get_index()-nbSolidToRemove-nbGhostToRemove);
				//Store the nodes
				NodeArraytmp.push_back(Node[i]);
			}//Remove ghost
			else
			{
				nbGhostToRemove++;
			}
		}
		else if(Node[i]->get_NodeType()!=Solid)
		{
			//map the new node array to old node array to keep node informations
			MapNodeArraytmp[i]=NodeArraytmp.size();
			//Update the index by removing solid
			Node[i]->Set_Index(Node[i]->Get_index()-nbSolidToRemove);
			//Store the nodes
			NodeArraytmp.push_back(Node[i]);
		}
		else
		{
			//Count number of solid to remove it
			nbSolidToRemove++;
			//Keep the first layer
			notSolid=0;
			for(unsigned int j=1;j<9;j++)
			{
				//Get first layer to do not remove it
				if(Node[Connect_lowOrder(i,j)]->get_NodeType()==Wall && Node[Connect_lowOrder(i,j)]->get_NodeType()==Corner)
					notSolid++;
			}
			if(notSolid>0)
			{
				//map the new node array to old node array to keep node informations
				MapNodeArraytmp[i]=NodeArraytmp.size();
				//Update the index by removing solid
				Node[i]->Set_Index(Node[i]->Get_index()-nbSolidToRemove);
				//Store the nodes
				NodeArraytmp.push_back(Node[i]);
			}


		}

	}
	//Correct numbering for the connections between cells
	for(int i=0;i<CellArraytmp.size();i++)
	{
		for(int j=0;j<4;j++)
		{
			//Set the connection for the j face to cell
			CellArraytmp[i]->Set_Connect(j,CellArraytmp[i]->Get_Connect(j)[0],MapCellArraytmp[i]);
			//Adjust the numbering
			CellArraytmp[i]->Set_NodeNumber(j,MapNodeArraytmp[CellArraytmp[i]->Get_NodeNumber(j)]);
		}
	}
}
void Block2D::InitPatchBc(Parameters *PtrParam){
	PatchsBc.initPatchBc(Node, PtrParam);
	std::vector<int> PatchIdInType=PatchsBc.Get_PatchIdInType();
	std::vector<SolverEnum::PatchType> PatchTypeInType=PatchsBc.Get_PatchTypeInType();
	int NumberOfPatchBc=PatchsBc.Get_NumberOfPatchBc();
	std::vector<int> NodeIdx;
	for(int j=0;j<NumberOfPatchBc;j++)
	{
		switch(PatchTypeInType[j])
		{
			case SolverEnum::Pressure:
				NodeIdx=PatchsBc.Get_PressurePatch()[PatchIdInType[j]].Get_NodeIndex();
				for (int i=0;i<NodeIdx.size();i++)
					if(Node[NodeIdx[i]]->get_NodeType()!=GlobalCorner && Node[NodeIdx[i]]->get_NodeType()!=Ghost)
						ChangeNodeType(NodeIdx[i], Pressure);
			break;
			case SolverEnum::Velocity:
				NodeIdx=PatchsBc.Get_VelocityPatch()[PatchIdInType[j]].Get_NodeIndex();
				for (int i=0;i<NodeIdx.size();i++)
					if(Node[NodeIdx[i]]->get_NodeType()!=GlobalCorner && Node[NodeIdx[i]]->get_NodeType()!=Ghost)
						ChangeNodeType(NodeIdx[i], Velocity);
			break;
			case SolverEnum::Periodic:
				NodeIdx=PatchsBc.Get_PeriodicPatch()[PatchIdInType[j]].Get_NodeIndex();
				for (int i=0;i<NodeIdx.size();i++)
					if(Node[NodeIdx[i]]->get_NodeType()!=GlobalCorner && Node[NodeIdx[i]]->get_NodeType()!=Ghost)
						ChangeNodeType(NodeIdx[i], Periodic);
			break;
			case SolverEnum::Symmetry:
				NodeIdx=PatchsBc.Get_SymmetryPatch()[PatchIdInType[j]].Get_NodeIndex();
				for (int i=0;i<NodeIdx.size();i++)
					if(Node[NodeIdx[i]]->get_NodeType()!=GlobalCorner && Node[NodeIdx[i]]->get_NodeType()!=Ghost)
						ChangeNodeType(NodeIdx[i], Symmetry);
			break;
			case SolverEnum::Wall:
				NodeIdx=PatchsBc.Get_WallPatch()[PatchIdInType[j]].Get_NodeIndex();
				for (int i=0;i<NodeIdx.size();i++)
					if(Node[NodeIdx[i]]->get_NodeType()!=GlobalCorner && Node[NodeIdx[i]]->get_NodeType()!=Ghost)
						ChangeNodeType(NodeIdx[i], Wall);
			break;
		}
	}
}
///Set the node arrays in Patch by type and also reorganise (and remove solid) by local index (all node in the processor)
void Block2D::Set_NodeIndexByTypeForPatchBc(){
	std::vector<int> PatchIdInType=PatchsBc.Get_PatchIdInType();
	std::vector<SolverEnum::PatchType> PatchTypeInType=PatchsBc.Get_PatchTypeInType();
	int NumberOfPatchBc=PatchsBc.Get_NumberOfPatchBc();
	std::vector<int> NodeIdx,NodeIdxtmp,NodeIdxSpecialWallstmp,NodeIdxGlobalCornertmp,NodeIdxbyTypetmp,NodeIdxbyTypeSpecialWallstmp,NodeIdxbyTypeGlobalCornertmp;
	for(int j=0;j<NumberOfPatchBc;j++)
	{
		switch(PatchTypeInType[j])
		{
			case SolverEnum::Pressure:
				NodeIdx=PatchsBc.Get_PressurePatch()[PatchIdInType[j]].Get_NodeIndex();
				NodeIdxtmp.clear();
				NodeIdxSpecialWallstmp.clear();
				NodeIdxGlobalCornertmp.clear();
				NodeIdxbyTypetmp.clear();
				NodeIdxbyTypeSpecialWallstmp.clear();
				NodeIdxbyTypeGlobalCornertmp.clear();
				for (int i=0;i<NodeIdx.size();i++)
				{
					if(NodeArrays.TypeOfNode[NodeIdx[i]]==Pressure)//Removing solid and special wall
					{
						NodeIdxtmp.push_back(NodeIdx[i]);
						NodeIdxbyTypetmp.push_back(NodeArrays.NodeIndexByType[NodeIdx[i]]);
					}
					else if(NodeArrays.TypeOfNode[NodeIdx[i]]==SpecialWall)//Mark Special wall
					{
						NodeIdxSpecialWallstmp.push_back(NodeIdx[i]);
						NodeIdxbyTypeSpecialWallstmp.push_back(NodeArrays.NodeIndexByType[NodeIdx[i]]);
					}
					else if(NodeArrays.TypeOfNode[NodeIdx[i]]==GlobalCorner)//Mark Global Corner
					{
						NodeIdxGlobalCornertmp.push_back(NodeIdx[i]);
						NodeIdxbyTypeGlobalCornertmp.push_back(NodeArrays.NodeIndexByType[NodeIdx[i]]);
					}
				}
				PatchsBc.Set_NodeIndex(PatchTypeInType[j],PatchIdInType[j],NodeIdxtmp);
				PatchsBc.Set_NodeIndexSpecialWalls(PatchTypeInType[j],PatchIdInType[j],NodeIdxSpecialWallstmp);
				PatchsBc.Set_NodeIndexGlobalCorner(PatchTypeInType[j],PatchIdInType[j],NodeIdxGlobalCornertmp);
				PatchsBc.Set_NodeIndexByType(PatchTypeInType[j],PatchIdInType[j],NodeIdxbyTypetmp);
				PatchsBc.Set_NodeIndexByTypeSpecialWalls(PatchTypeInType[j],PatchIdInType[j],NodeIdxbyTypeSpecialWallstmp);
				PatchsBc.Set_NodeIndexByTypeGlobalCorner(PatchTypeInType[j],PatchIdInType[j],NodeIdxbyTypeGlobalCornertmp);

				break;
			case SolverEnum::Velocity:
				NodeIdx=PatchsBc.Get_VelocityPatch()[PatchIdInType[j]].Get_NodeIndex();
				NodeIdxtmp.clear();
				NodeIdxSpecialWallstmp.clear();
				NodeIdxGlobalCornertmp.clear();
				NodeIdxbyTypetmp.clear();
				NodeIdxbyTypeSpecialWallstmp.clear();
				NodeIdxbyTypeGlobalCornertmp.clear();
				for (int i=0;i<NodeIdx.size();i++)
				{
					if(NodeArrays.TypeOfNode[NodeIdx[i]]==Velocity)//Removing solid and special wall
					{
						NodeIdxtmp.push_back(NodeIdx[i]);
						NodeIdxbyTypetmp.push_back(NodeArrays.NodeIndexByType[NodeIdx[i]]);
					}
					if(NodeArrays.TypeOfNode[NodeIdx[i]]==SpecialWall)//Mark Special wall
					{
						NodeIdxSpecialWallstmp.push_back(NodeIdx[i]);
						NodeIdxbyTypeSpecialWallstmp.push_back(NodeArrays.NodeIndexByType[NodeIdx[i]]);
					}
					if(NodeArrays.TypeOfNode[NodeIdx[i]]==GlobalCorner)//Mark Global Corner
					{
						NodeIdxGlobalCornertmp.push_back(NodeIdx[i]);
						NodeIdxbyTypeGlobalCornertmp.push_back(NodeArrays.NodeIndexByType[NodeIdx[i]]);
					}
				}
				PatchsBc.Set_NodeIndex(PatchTypeInType[j],PatchIdInType[j],NodeIdxtmp);
				PatchsBc.Set_NodeIndexSpecialWalls(PatchTypeInType[j],PatchIdInType[j],NodeIdxSpecialWallstmp);
				PatchsBc.Set_NodeIndexGlobalCorner(PatchTypeInType[j],PatchIdInType[j],NodeIdxGlobalCornertmp);
				PatchsBc.Set_NodeIndexByType(PatchTypeInType[j],PatchIdInType[j],NodeIdxbyTypetmp);
				PatchsBc.Set_NodeIndexByTypeSpecialWalls(PatchTypeInType[j],PatchIdInType[j],NodeIdxbyTypeSpecialWallstmp);
				PatchsBc.Set_NodeIndexByTypeGlobalCorner(PatchTypeInType[j],PatchIdInType[j],NodeIdxbyTypeGlobalCornertmp);
			break;
			case SolverEnum::Periodic:
				NodeIdx=PatchsBc.Get_PeriodicPatch()[PatchIdInType[j]].Get_NodeIndex();
				NodeIdxtmp.clear();
				NodeIdxSpecialWallstmp.clear();
				NodeIdxGlobalCornertmp.clear();
				NodeIdxbyTypetmp.clear();
				NodeIdxbyTypeSpecialWallstmp.clear();
				NodeIdxbyTypeGlobalCornertmp.clear();
				for (int i=0;i<NodeIdx.size();i++)
				{
					if(NodeArrays.TypeOfNode[NodeIdx[i]]==Periodic)//Removing solid and special wall
					{
						NodeIdxtmp.push_back(NodeIdx[i]);
						NodeIdxbyTypetmp.push_back(NodeArrays.NodeIndexByType[NodeIdx[i]]);
					}
					if(NodeArrays.TypeOfNode[NodeIdx[i]]==SpecialWall)//Mark Special wall
					{
						NodeIdxSpecialWallstmp.push_back(NodeIdx[i]);
						NodeIdxbyTypeSpecialWallstmp.push_back(NodeArrays.NodeIndexByType[NodeIdx[i]]);
					}
					if(NodeArrays.TypeOfNode[NodeIdx[i]]==GlobalCorner)//Mark Global Corner
					{
						NodeIdxGlobalCornertmp.push_back(NodeIdx[i]);
						NodeIdxbyTypeGlobalCornertmp.push_back(NodeArrays.NodeIndexByType[NodeIdx[i]]);
					}
				}
				PatchsBc.Set_NodeIndex(PatchTypeInType[j],PatchIdInType[j],NodeIdxtmp);
				PatchsBc.Set_NodeIndexSpecialWalls(PatchTypeInType[j],PatchIdInType[j],NodeIdxSpecialWallstmp);
				PatchsBc.Set_NodeIndexGlobalCorner(PatchTypeInType[j],PatchIdInType[j],NodeIdxGlobalCornertmp);
				PatchsBc.Set_NodeIndexByType(PatchTypeInType[j],PatchIdInType[j],NodeIdxbyTypetmp);
				PatchsBc.Set_NodeIndexByTypeSpecialWalls(PatchTypeInType[j],PatchIdInType[j],NodeIdxbyTypeSpecialWallstmp);
				PatchsBc.Set_NodeIndexByTypeGlobalCorner(PatchTypeInType[j],PatchIdInType[j],NodeIdxbyTypeGlobalCornertmp);
			break;
			case SolverEnum::Symmetry:
				NodeIdx=PatchsBc.Get_SymmetryPatch()[PatchIdInType[j]].Get_NodeIndex();
				NodeIdxtmp.clear();
				NodeIdxSpecialWallstmp.clear();
				NodeIdxGlobalCornertmp.clear();
				NodeIdxbyTypetmp.clear();
				NodeIdxbyTypeSpecialWallstmp.clear();
				NodeIdxbyTypeGlobalCornertmp.clear();
				for (int i=0;i<NodeIdx.size();i++)
				{
					if(NodeArrays.TypeOfNode[NodeIdx[i]]==Symmetry)//Removing solid and special wall
					{
						NodeIdxtmp.push_back(NodeIdx[i]);
						NodeIdxbyTypetmp.push_back(NodeArrays.NodeIndexByType[NodeIdx[i]]);
					}
					if(NodeArrays.TypeOfNode[NodeIdx[i]]==SpecialWall)//Mark Special wall
					{
						NodeIdxSpecialWallstmp.push_back(NodeIdx[i]);
						NodeIdxbyTypeSpecialWallstmp.push_back(NodeArrays.NodeIndexByType[NodeIdx[i]]);
					}
					if(NodeArrays.TypeOfNode[NodeIdx[i]]==GlobalCorner)//Mark Global Corner
					{
						NodeIdxGlobalCornertmp.push_back(NodeIdx[i]);
						NodeIdxbyTypeGlobalCornertmp.push_back(NodeArrays.NodeIndexByType[NodeIdx[i]]);
					}
				}
				PatchsBc.Set_NodeIndex(PatchTypeInType[j],PatchIdInType[j],NodeIdxtmp);
				PatchsBc.Set_NodeIndexSpecialWalls(PatchTypeInType[j],PatchIdInType[j],NodeIdxSpecialWallstmp);
				PatchsBc.Set_NodeIndexGlobalCorner(PatchTypeInType[j],PatchIdInType[j],NodeIdxGlobalCornertmp);
				PatchsBc.Set_NodeIndexByType(PatchTypeInType[j],PatchIdInType[j],NodeIdxbyTypetmp);
				PatchsBc.Set_NodeIndexByTypeSpecialWalls(PatchTypeInType[j],PatchIdInType[j],NodeIdxbyTypeSpecialWallstmp);
				PatchsBc.Set_NodeIndexByTypeGlobalCorner(PatchTypeInType[j],PatchIdInType[j],NodeIdxbyTypeGlobalCornertmp);
			break;
			case SolverEnum::Wall:
				NodeIdx=PatchsBc.Get_WallPatch()[PatchIdInType[j]].Get_NodeIndex();
				NodeIdxtmp.clear();
				NodeIdxSpecialWallstmp.clear();
				NodeIdxGlobalCornertmp.clear();
				NodeIdxbyTypetmp.clear();
				NodeIdxbyTypeSpecialWallstmp.clear();
				NodeIdxbyTypeGlobalCornertmp.clear();
	/*			for (int i=0;i<NodeIdx.size();i++)
				{
					if(NodeArrays.TypeOfNode[NodeIdx[i]]==Wall)//Removing solid and special wall
						NodeIdxbyTypetmp.push_back(NodeArrays.NodeIndexByType[NodeIdx[i]]);
					if(NodeArrays.TypeOfNode[NodeIdx[i]]==SpecialWall||NodeArrays.TypeOfNode[NodeIdx[i]]==ConcaveCorner)//Mark Special wall
						NodeIdxbyTypeSpecialWallstmp.push_back(NodeArrays.NodeIndexByType[NodeIdx[i]]);
					if(NodeArrays.TypeOfNode[NodeIdx[i]]==GlobalCorner)//Mark Global Corner
						NodeIdxbyTypeGlobalCornertmp.push_back(NodeArrays.NodeIndexByType[NodeIdx[i]]);
				}
				*/
				PatchsBc.Set_NodeIndex(PatchTypeInType[j],PatchIdInType[j],NodeIdxtmp);
				PatchsBc.Set_NodeIndexSpecialWalls(PatchTypeInType[j],PatchIdInType[j],NodeIdxSpecialWallstmp);
				PatchsBc.Set_NodeIndexGlobalCorner(PatchTypeInType[j],PatchIdInType[j],NodeIdxGlobalCornertmp);
				PatchsBc.Set_NodeIndexByType(PatchTypeInType[j],PatchIdInType[j],NodeIdxbyTypetmp);
				PatchsBc.Set_NodeIndexByTypeSpecialWalls(PatchTypeInType[j],PatchIdInType[j],NodeIdxbyTypeSpecialWallstmp);
				PatchsBc.Set_NodeIndexByTypeGlobalCorner(PatchTypeInType[j],PatchIdInType[j],NodeIdxbyTypeGlobalCornertmp);			break;
		}
	}
}
//enum NodeType {Interior, Solid, Ghost, Corner, Wall, Periodic, Velocity, Symmetry, Pressure,ConcaveCorner,ConvexCorner, SolidGhost,GlobalCorner, SpecialWall};

void Block2D::Write_CommunicationNodes(){
	int rank;
	double x_tmp=0;double y_tmp=0;int index=0;
 	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

 	char buffer[50]; // make sure it's big enough
 	snprintf(buffer, sizeof(buffer), "IdCommm_InFluid_%d.txt", rank);
 	std::ofstream myFlux;
 	myFlux.open(buffer);
 	myFlux<<std::endl<<" ****** Communication in Fluid *******"<<std::endl;
 	myFlux<<"West real nodes: "<<std::endl;
 	for(int i=0;i<IdRNodeW.size();i++)
 	{
 		index=IdRNodeW[i];
 		NodeArrays.Get_coordinate(index,x_tmp, y_tmp);
 		myFlux<<"Node ID: "<<index<<" x: "<<x_tmp <<" y: "<<y_tmp<<std::endl;
 	}
 	myFlux<<std::endl<<"West Ghost nodes: "<<std::endl;
 	for(int i=0;i<IdGNodeW.size();i++)
 	{
 		index=IdGNodeW[i];
 		NodeArrays.Get_coordinate(index,x_tmp, y_tmp);
 		myFlux<<"Node ID: "<<index<<" x: "<<x_tmp <<" y: "<<y_tmp<<std::endl;
 	}
 	myFlux<<std::endl<<"North real nodes: "<<std::endl;
 	for(int i=0;i<IdRNodeN.size();i++)
 	{
 		index=IdRNodeN[i];
 		NodeArrays.Get_coordinate(index,x_tmp, y_tmp);
 		myFlux<<"Node ID: "<<index<<" x: "<<x_tmp <<" y: "<<y_tmp<<std::endl;
 	}
 	myFlux<<std::endl<<"North Ghost nodes: "<<std::endl;
 	for(int i=0;i<IdGNodeN.size();i++)
 	{
 		index=IdGNodeN[i];
 		NodeArrays.Get_coordinate(index,x_tmp, y_tmp);
 		myFlux<<"Node ID: "<<index<<" x: "<<x_tmp <<" y: "<<y_tmp<<std::endl;
 	}
 	myFlux<<std::endl<<"South real nodes: "<<std::endl;
 	for(int i=0;i<IdRNodeS.size();i++)
 	{
 		index=IdRNodeS[i];
 		NodeArrays.Get_coordinate(index,x_tmp, y_tmp);
 		myFlux<<"Node ID: "<<index<<" x: "<<x_tmp <<" y: "<<y_tmp<<std::endl;
 	}
 	myFlux<<std::endl<<"South Ghost nodes: "<<std::endl;
 	for(int i=0;i<IdGNodeS.size();i++)
 	{
 		index=IdGNodeS[i];
 		NodeArrays.Get_coordinate(index,x_tmp, y_tmp);
 		myFlux<<"Node ID: "<<index<<" x: "<<x_tmp <<" y: "<<y_tmp<<std::endl;
 	}
 	myFlux<<std::endl<<"East real nodes: "<<std::endl;
 	for(int i=0;i<IdRNodeE.size();i++)
 	{
 		index=IdRNodeE[i];
 		NodeArrays.Get_coordinate(index,x_tmp, y_tmp);
 		myFlux<<"Node ID: "<<index<<" x: "<<x_tmp <<" y: "<<y_tmp<<std::endl;
 	}
 	myFlux<<std::endl<<"East Ghost nodes: "<<std::endl;
 	for(int i=0;i<IdGNodeE.size();i++)
 	{
 		index=IdGNodeE[i];
 		NodeArrays.Get_coordinate(index,x_tmp, y_tmp);
 		myFlux<<"Node ID: "<<index<<" x: "<<x_tmp <<" y: "<<y_tmp<<std::endl;
 	}

 	myFlux<<std::endl<<std::endl;
 	myFlux<<"South West real nodes: "<<std::endl;
 	for(int i=0;i<IdRNodeSW.size();i++)
 	{
 		index=IdRNodeSW[i];
 		NodeArrays.Get_coordinate(index,x_tmp, y_tmp);
 		myFlux<<"Node ID: "<<index<<" x: "<<x_tmp <<" y: "<<y_tmp<<std::endl;
 	}
 	myFlux<<std::endl<<"South West Ghost nodes: "<<std::endl;
 	for(int i=0;i<IdGNodeSW.size();i++)
 	{
 		index=IdGNodeSW[i];
 		NodeArrays.Get_coordinate(index,x_tmp, y_tmp);
 		myFlux<<"Node ID: "<<index<<" x: "<<x_tmp <<" y: "<<y_tmp<<std::endl;
 	}
 	myFlux<<std::endl<<"North West real nodes: "<<std::endl;
 	for(int i=0;i<IdRNodeNW.size();i++)
 	{
 		index=IdRNodeNW[i];
 		NodeArrays.Get_coordinate(index,x_tmp, y_tmp);
 		myFlux<<"Node ID: "<<index<<" x: "<<x_tmp <<" y: "<<y_tmp<<std::endl;
 	}
 	myFlux<<std::endl<<"North West Ghost nodes: "<<std::endl;
 	for(int i=0;i<IdGNodeNW.size();i++)
 	{
 		index=IdGNodeNW[i];
 		NodeArrays.Get_coordinate(index,x_tmp, y_tmp);
 		myFlux<<"Node ID: "<<index<<" x: "<<x_tmp <<" y: "<<y_tmp<<std::endl;
 	}
 	myFlux<<std::endl<<"South East real nodes: "<<std::endl;
 	for(int i=0;i<IdRNodeSE.size();i++)
 	{
 		index=IdRNodeSE[i];
 		NodeArrays.Get_coordinate(index,x_tmp, y_tmp);
 		myFlux<<"Node ID: "<<index<<" x: "<<x_tmp <<" y: "<<y_tmp<<std::endl;
 	}
 	myFlux<<std::endl<<"South East Ghost nodes: "<<std::endl;
 	for(int i=0;i<IdGNodeSE.size();i++)
 	{
 		index=IdGNodeSE[i];
 		NodeArrays.Get_coordinate(index,x_tmp, y_tmp);
 		myFlux<<"Node ID: "<<index<<" x: "<<x_tmp <<" y: "<<y_tmp<<std::endl;
 	}
 	myFlux<<std::endl<<"North East real nodes: "<<std::endl;
 	for(int i=0;i<IdRNodeNE.size();i++)
 	{
 		index=IdRNodeNE[i];
 		NodeArrays.Get_coordinate(index,x_tmp, y_tmp);
 		myFlux<<"Node ID: "<<index<<" x: "<<x_tmp <<" y: "<<y_tmp<<std::endl;
 	}
 	myFlux<<std::endl<<"North East Ghost nodes: "<<std::endl;
 	for(int i=0;i<IdGNodeNE.size();i++)
 	{
 		index=IdGNodeNE[i];
 		NodeArrays.Get_coordinate(index,x_tmp, y_tmp);
 		myFlux<<"Node ID: "<<index<<" x: "<<x_tmp <<" y: "<<y_tmp<<std::endl;
 	}

}
void Block2D::Write_CommunicationSolidNodes(){
	int rank;
	double x_tmp=0;double y_tmp=0;int index=0;
 	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

 	char buffer[50]; // make sure it's big enough
 	snprintf(buffer, sizeof(buffer), "IdCommm_InSolid_%d.txt", rank);
 	std::ofstream myFlux;
 	myFlux.open(buffer);
 	myFlux<<std::endl<<" ****** Communication in Solid *******"<<std::endl;
 	myFlux<<"West real nodes: "<<std::endl;
 	for(int i=0;i<SolidIdRNodeW.size();i++)
 	{
 		index=SolidIdRNodeW[i];
 		NodeArrays.Get_coordinate(index,x_tmp, y_tmp);
 		myFlux<<"Node ID: "<<index<<" x: "<<x_tmp <<" y: "<<y_tmp<<std::endl;
 	}
 	myFlux<<std::endl<<"West Ghost nodes: "<<std::endl;
 	for(int i=0;i<SolidIdGNodeW.size();i++)
 	{
 		index=SolidIdGNodeW[i];
 		NodeArrays.Get_coordinate(index,x_tmp, y_tmp);
 		myFlux<<"Node ID: "<<index<<" x: "<<x_tmp <<" y: "<<y_tmp<<std::endl;
 	}
 	myFlux<<std::endl<<"North real nodes: "<<std::endl;
 	for(int i=0;i<SolidIdRNodeN.size();i++)
 	{
 		index=SolidIdRNodeN[i];
 		NodeArrays.Get_coordinate(index,x_tmp, y_tmp);
 		myFlux<<"Node ID: "<<index<<" x: "<<x_tmp <<" y: "<<y_tmp<<std::endl;
 	}
 	myFlux<<std::endl<<"North Ghost nodes: "<<std::endl;
 	for(int i=0;i<SolidIdGNodeN.size();i++)
 	{
 		index=SolidIdGNodeN[i];
 		NodeArrays.Get_coordinate(index,x_tmp, y_tmp);
 		myFlux<<"Node ID: "<<index<<" x: "<<x_tmp <<" y: "<<y_tmp<<std::endl;
 	}
 	myFlux<<std::endl<<"South real nodes: "<<std::endl;
 	for(int i=0;i<SolidIdRNodeS.size();i++)
 	{
 		index=SolidIdRNodeS[i];
 		NodeArrays.Get_coordinate(index,x_tmp, y_tmp);
 		myFlux<<"Node ID: "<<index<<" x: "<<x_tmp <<" y: "<<y_tmp<<std::endl;
 	}
 	myFlux<<std::endl<<"South Ghost nodes: "<<std::endl;
 	for(int i=0;i<SolidIdGNodeS.size();i++)
 	{
 		index=SolidIdGNodeS[i];
 		NodeArrays.Get_coordinate(index,x_tmp, y_tmp);
 		myFlux<<"Node ID: "<<index<<" x: "<<x_tmp <<" y: "<<y_tmp<<std::endl;
 	}
 	myFlux<<std::endl<<"East real nodes: "<<std::endl;
 	for(int i=0;i<SolidIdRNodeE.size();i++)
 	{
 		index=SolidIdRNodeE[i];
 		NodeArrays.Get_coordinate(index,x_tmp, y_tmp);
 		myFlux<<"Node ID: "<<index<<" x: "<<x_tmp <<" y: "<<y_tmp<<std::endl;
 	}
 	myFlux<<std::endl<<"East Ghost nodes: "<<std::endl;
 	for(int i=0;i<SolidIdGNodeE.size();i++)
 	{
 		index=SolidIdGNodeE[i];
 		NodeArrays.Get_coordinate(index,x_tmp, y_tmp);
 		myFlux<<"Node ID: "<<index<<" x: "<<x_tmp <<" y: "<<y_tmp<<std::endl;
 	}

 	myFlux<<std::endl<<std::endl;
 	myFlux<<"South West real nodes: "<<std::endl;
 	for(int i=0;i<SolidIdRNodeSW.size();i++)
 	{
 		index=SolidIdRNodeSW[i];
 		NodeArrays.Get_coordinate(index,x_tmp, y_tmp);
 		myFlux<<"Node ID: "<<index<<" x: "<<x_tmp <<" y: "<<y_tmp<<std::endl;
 	}
 	myFlux<<std::endl<<"South West Ghost nodes: "<<std::endl;
 	for(int i=0;i<SolidIdGNodeSW.size();i++)
 	{
 		index=SolidIdGNodeSW[i];
 		NodeArrays.Get_coordinate(index,x_tmp, y_tmp);
 		myFlux<<"Node ID: "<<index<<" x: "<<x_tmp <<" y: "<<y_tmp<<std::endl;
 	}
 	myFlux<<std::endl<<"North West real nodes: "<<std::endl;
 	for(int i=0;i<SolidIdRNodeNW.size();i++)
 	{
 		index=SolidIdRNodeNW[i];
 		NodeArrays.Get_coordinate(index,x_tmp, y_tmp);
 		myFlux<<"Node ID: "<<index<<" x: "<<x_tmp <<" y: "<<y_tmp<<std::endl;
 	}
 	myFlux<<std::endl<<"North West Ghost nodes: "<<std::endl;
 	for(int i=0;i<SolidIdGNodeNW.size();i++)
 	{
 		index=SolidIdGNodeNW[i];
 		NodeArrays.Get_coordinate(index,x_tmp, y_tmp);
 		myFlux<<"Node ID: "<<index<<" x: "<<x_tmp <<" y: "<<y_tmp<<std::endl;
 	}
 	myFlux<<std::endl<<"South East real nodes: "<<std::endl;
 	for(int i=0;i<SolidIdRNodeSE.size();i++)
 	{
 		index=SolidIdRNodeSE[i];
 		NodeArrays.Get_coordinate(index,x_tmp, y_tmp);
 		myFlux<<"Node ID: "<<index<<" x: "<<x_tmp <<" y: "<<y_tmp<<std::endl;
 	}
 	myFlux<<std::endl<<"South East Ghost nodes: "<<std::endl;
 	for(int i=0;i<SolidIdGNodeSE.size();i++)
 	{
 		index=SolidIdGNodeSE[i];
 		NodeArrays.Get_coordinate(index,x_tmp, y_tmp);
 		myFlux<<"Node ID: "<<index<<" x: "<<x_tmp <<" y: "<<y_tmp<<std::endl;
 	}
 	myFlux<<std::endl<<"North East real nodes: "<<std::endl;
 	for(int i=0;i<SolidIdRNodeNE.size();i++)
 	{
 		index=SolidIdRNodeNE[i];
 		NodeArrays.Get_coordinate(index,x_tmp, y_tmp);
 		myFlux<<"Node ID: "<<index<<" x: "<<x_tmp <<" y: "<<y_tmp<<std::endl;
 	}
 	myFlux<<std::endl<<"North East Ghost nodes: "<<std::endl;
 	for(int i=0;i<SolidIdGNodeNE.size();i++)
 	{
 		index=SolidIdGNodeNE[i];
 		NodeArrays.Get_coordinate(index,x_tmp, y_tmp);
 		myFlux<<"Node ID: "<<index<<" x: "<<x_tmp <<" y: "<<y_tmp<<std::endl;
 	}

}
