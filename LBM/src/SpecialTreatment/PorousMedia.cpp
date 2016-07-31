/*
 * PorousMedia.cpp
 *
 *  Created on: 17 Oct 2015
 *      Author: thomas
 */

#include "PorousMedia.h"

PorousMedia::PorousMedia() {
image=0;
image_nx=0;
image_ny=0;
image_nz=0;
xstartmedia=0;
ystartmedia=0;
zstartmedia=0;
}

PorousMedia::~PorousMedia() {
	for (int k=0;k<image_nz;k++)
	{
		for (int j=0;j<image_ny;j++)
			delete [] image[k][j];
		delete [] image[k];
	}
	delete[] image;
}

void PorousMedia::ReadBinaryImage(std::string filename,int nx, int ny, int nz)
{
	image_nx=nx;
	image_ny=ny;
	image_nz=nz;
	std::ifstream in(filename.c_str(), std::ios::in|std::ios::binary|std::ios::ate);
	std::streampos size;
	char * memblock;
	if (in.is_open())
	  {
	    size = in.tellg();
	    memblock = new char [size];
	    in.seekg (0, std::ios::beg);
	    in.read (memblock, size);
	    in.close();

	    std::cout << "the entire file content is in memory with the size of: "<<size<<std::endl;
		if(image!=0)
		{
			std::cout<<"deleting the previous image save"<<std::endl;
			delete[] image;
		}
		image =new bool**[image_nz];
		for (int k=0;k<image_nz;k++)
		{
			image[k]=new bool* [image_ny];
			for (int j=0;j<image_ny;j++)
				image[k][j]=new bool [image_nx];
		}
		int idx=0;
		for (int k=0;k<image_nz;k++)
			for (int j=0;j<image_ny;j++)
				for (int i=0;i<image_nx;i++)
				{
					if((int)memblock[idx]==1)
						image[k][j][i]=true;
					else
						image[k][j][i]=false;
					idx++;
					//if(k==0)
					//std::cout<<"x: "<<i<<" y: "<<j<<" type: "<<image[k][j][i]<<std::endl;;
				}
		delete[] memblock;
	  }
	//AdaptSolid();
}
void PorousMedia::AdaptSolid(){
	if(image_nz<=1)//2D
	{
		DomainCornerTreatment2D();
		//DomainBoundaryTreatment2D();
		DomainInteriorTreatment2D();
	}
	else
	{

	}
}
void PorousMedia::DomainCornerTreatment2D()
{
	//Corner (0,0)
	if(image[0][0][1]==true && image[0][1][0]==true)
		image[0][0][0]=true;
	if(image[0][0][1]==false && image[0][1][0]==false)
		image[0][0][0]=false;
	//Corner (nx,0)
	if(image[0][0][image_nx-2]==true && image[0][1][image_nx-1]==true)
		image[0][0][image_nx-1]=true;
	if(image[0][0][image_nx-2]==false && image[0][1][image_nx-1]==false)
		image[0][0][image_nx-1]=false;
	//Corner (0,ny)
	if(image[0][image_ny-1][1]==true && image[0][image_ny-2][0]==true)
		image[0][image_ny-1][0]=true;
	if(image[0][image_ny-1][1]==false && image[0][image_ny-2][0]==false)
		image[0][image_ny-1][0]=false;
	//Corner (nx,ny)
	if(image[0][image_ny-1][image_nx-2]==true && image[0][image_ny-2][image_nx-1]==true)
		image[0][image_ny-1][image_nx-1]=true;
	if(image[0][image_ny-1][image_nx-2]==false && image[0][image_ny-2][image_nx-1]==false)
		image[0][image_ny-1][image_nx-1]=false;
}
void PorousMedia::DomainBoundaryTreatment2D()
{
	for (int i=1;i<image_nx-2;i++)
	{
		//Bottom
		if(image[0][0][i]==true)
		{
			if((image[0][0][i+1]==true && image[0][1][i+1]==true)||(image[0][0][i-1]==true && image[0][1][i-1]==true))
			{
				image[0][1][i]=true;
			}
			if(image[0][0][i+1]==false && image[0][0][i-1]==false )
			{
				image[0][0][i]=false;
			}
		}
		//Top
		if(image[0][image_ny-1][i]==true)
		{
			if((image[0][image_ny-1][i+1]==true && image[0][image_ny-2][i+1]==true)||(image[0][image_ny-1][i-1]==true && image[0][image_ny-2][i-1]==true))
			{
				image[0][image_ny-2][i]=true;
			}
			if(image[0][image_ny-1][i+1]==false && image[0][image_ny-1][i-1]==false )
			{
				image[0][image_ny-1][i]=false;
			}
		}
	}
	for (int i=1;i<image_ny-2;i++)
	{
		//West
		if(image[0][i][0]==true)
		{
			if((image[0][i+1][0]==true && image[0][i+1][0]==true)||(image[0][i-1][0]==true && image[0][i-1][0]==true))
			{
				image[0][i][1]=true;
			}
			if(image[0][i+1][0]==false &&image[0][0][i-1]==false )
			{
				image[0][i][0]=false;
			}
		}
		//East
		if(image[0][i][image_nx-1]==true)
		{
			if((image[0][i+1][image_nx-1]==true && image[0][i+1][image_nx-2]==true)||(image[0][i-1][image_nx-1]==true && image[0][i-1][image_nx-2]==true))
			{
				image[0][i][image_nx-2]=true;
			}
			if(image[0][i+1][image_nx-1]==false && image[0][i-1][image_nx-1]==false )
			{
				image[0][i][image_nx-1]=false;
			}
		}
	}
}
void PorousMedia::DomainInteriorTreatment2D()
{
	bool allneighbours=false; //variable to remove calculation cost
	for (int j=1;j<image_ny-1;j++)
		for (int i=1;i<image_nx-1;i++)
		{
			allneighbours=image[0][j][i+1]==true && image[0][j][i-1]==true && image[0][j+1][i]==true && image[0][j-1][i]==true;
			if(image[0][j][i])// Check if it is Solid
			{
				CorrectPoreConnections2D(i, j,allneighbours);
			}
			else
			{
				//Remove lonely Fluid
				if(allneighbours)
					image[0][j][i]=true;
			}
//			if(i==57&&j==107)
//			{
//				std::cout<<" type node at x: "<<i-1<<" y: "<<j+2<< " is: "<<image[0][j+2][i-1]<<std::endl;
//				std::cout<<" type node at x: "<<i-1<<" y: "<<j+1<< " is: "<<image[0][j+1][i-1]<<std::endl;
//				std::cout<<" type node at x: "<<i-1<<" y: "<<j<< " is: "<<image[0][j][i-1]<<std::endl;
//				std::cout<<" type node at x: "<<i-1<<" y: "<<j-1<< " is: "<<image[0][j-1][i-1]<<std::endl;
//				std::cout<<" type node at x: "<<i-1<<" y: "<<j-2<< " is: "<<image[0][j-2][i-1]<<std::endl;
//				std::cout<<" type node at x: "<<i<<" y: "<<j+2<< " is: "<<image[0][j+2][i]<<std::endl;
//				std::cout<<" type node at x: "<<i<<" y: "<<j+1<< " is: "<<image[0][j+1][i]<<std::endl;
//				std::cout<<" type node at x: "<<i<<" y: "<<j<< " is: "<<image[0][j][i]<<std::endl;
//				std::cout<<" type node at x: "<<i<<" y: "<<j-1<< " is: "<<image[0][j-1][i]<<std::endl;
//				std::cout<<" type node at x: "<<i<<" y: "<<j-2<< " is: "<<image[0][j-2][i]<<std::endl;
//				std::cout<<" type node at x: "<<i+1<<" y: "<<j+2<< " is: "<<image[0][j+2][i+1]<<std::endl;
//				std::cout<<" type node at x: "<<i+1<<" y: "<<j+1<< " is: "<<image[0][j+1][i+1]<<std::endl;
//				std::cout<<" type node at x: "<<i+1<<" y: "<<j<< " is: "<<image[0][j][i+1]<<std::endl;
//				std::cout<<" type node at x: "<<i+1<<" y: "<<j-1<< " is: "<<image[0][j-1][i+1]<<std::endl;
//				std::cout<<" type node at x: "<<i+1<<" y: "<<j-2<< " is: "<<image[0][j-2][i+1]<<std::endl;
//				std::cout<<" type node at x: "<<i+2<<" y: "<<j+2<< " is: "<<image[0][j+2][i+2]<<std::endl;
//				std::cout<<" type node at x: "<<i+2<<" y: "<<j+1<< " is: "<<image[0][j+1][i+2]<<std::endl;
//				std::cout<<" type node at x: "<<i+2<<" y: "<<j<< " is: "<<image[0][j][i+2]<<std::endl;
//				std::cout<<" type node at x: "<<i+2<<" y: "<<j-1<< " is: "<<image[0][j-1][i+2]<<std::endl;
//				std::cout<<" type node at x: "<<i+2<<" y: "<<j-2<< " is: "<<image[0][j-2][i+2]<<std::endl;
//			}
		}
}
void PorousMedia::CorrectPoreConnections2D(int & i, int & j,bool & allneighbours)
{
	///Check if all neighbours are fluid
	if(image[0][j][i+1]==false  && image[0][j][i-1]==false  && image[0][j+1][i]==false  && image[0][j-1][i]==false )
		image[0][j][i]=false; //Remove lonely solid... It feels lonely, it is better to remove it

	/*///Remove two corners at the same points by adding two solids
	//Diagonal  (j+1,i-1) (j-1,i+1)
	if(allneighbours && image[0][j+1][i+1] && image[0][j-1][i-1])
		{image[0][j-1][i+1]=true; image[0][j+1][i-1]=true;}
	//Diagonal  (j-1,i-1) (j+1,i+1)
	if(allneighbours && image[0][j-1][i+1] && image[0][j+1][i-1])
		{image[0][j+1][i+1]=true; image[0][j-1][i-1]=true;}*/
	if(allneighbours)
		RemoveWrongCorners2D(i, j,allneighbours);

	///Detect 3 interiors connected and add 2 solids to remove thin wall
	//Right solid and other fluids
	if(image[0][j][i+1] && image[0][j-1][i]==false && image[0][j+1][i]==false && image[0][j][i-1]==false)
		{image[0][j+1][i]=true; image[0][j-1][i]=true;}
	//Bottom solid and other fluids
	if(image[0][j-1][i] && image[0][j][i-1]==false && image[0][j][i+1]==false && image[0][j+1][i]==false)
		{image[0][j][i-1]=true; image[0][j][i+1]=true;}
	//Left solid and other fluids
	if(image[0][j][i-1] && image[0][j-1][i]==false && image[0][j+1][i]==false && image[0][j][i+1]==false)
		{image[0][j+1][i]=true; image[0][j-1][i]=true;}
	//Top solid and other fluids
	if(image[0][j+1][i] && image[0][j-1][i]==false && image[0][j][i+1]==false && image[0][j][i-1]==false)
		{image[0][j][i-1]=true; image[0][j][i+1]=true;}
}

void PorousMedia::RemoveWrongCorners2D(int i, int j,bool allneighbours)
{
	int diagonal=0;int count=0;
	bool testWrongCorner=CheckWrongCorners2D(i, j,allneighbours, diagonal);
	if(testWrongCorner)
	{
		while(testWrongCorner && count<50)
		{
			switch(diagonal)
			{
			default :
				std::cerr<<" Error during removing wrong corners"<<std::endl;
				break;
			case 1:
				image[0][j-1][i+1]=true; image[0][j+1][i-1]=true;
				/*if(j>2)
				{
					j--;
					i++;
					allneighbours=image[0][j][i+1]==true && image[0][j][i-1]==true && image[0][j+1][i]==true && image[0][j-1][i]==true;
					if(allneighbours)
						testWrongCorner=CheckWrongCorners2D(i, j,allneighbours, diagonal);
					else
						testWrongCorner=false;
				}
				else*/
					testWrongCorner=false;
				break;
			case 2:
				image[0][j+1][i+1]=true; image[0][j-1][i-1]=true;
				/*if(j>2)
				{
					j--;
					i--;
					allneighbours=image[0][j][i+1]==true && image[0][j][i-1]==true && image[0][j+1][i]==true && image[0][j-1][i]==true;
					if(allneighbours)
					testWrongCorner=CheckWrongCorners2D(i, j,allneighbours, diagonal);
					else
						testWrongCorner=false;
				}
				else*/
					testWrongCorner=false;
				break;
			}
			count++;
		}
	}
}
bool PorousMedia::CheckWrongCorners2D(int i, int j,bool & allneighbours, int & diagonal)
{
	bool diag1S=image[0][j+1][i+1]==true && image[0][j-1][i-1]==true;
	bool diag2S=image[0][j-1][i+1]==true && image[0][j+1][i-1]==true;
	bool diag1F=image[0][j+1][i+1]==false && image[0][j-1][i-1]==false;
	bool diag2F=image[0][j-1][i+1]==false && image[0][j+1][i-1]==false;
	if(diag1S==true&&diag2S==true) //solid
	{diagonal=0; return false;}
	else
		//Diagonal  (j+1,i-1) (j-1,i+1)
		if(diag1S==true && diag2F==true)
			{diagonal=1; return true;}
		else
			//Diagonal  (j-1,i-1) (j+1,i+1)
			if(diag2S==true && diag1F==true)
				{diagonal=2; return true;}
			//No wrong Corners
			else
				{diagonal=0; return false;}
}
/*
for (int j=1;j<image_ny-1;j++)
	for (int i=1;i<image_nx-1;i++)
		DomainInteriorTreatment();
for (int j=1;j<image_ny-1;j++)
	DomainBoundaryTreatment();
	*/
