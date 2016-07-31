/*
 * UserMesh.cpp
 *
 *  Created on: 17 Oct 2015
 *      Author: thomas
 */

#include "UserMesh.h"

UserMesh::UserMesh() {

	//std::string filename("../5RockSamples/BeadPack300^3.raw");
	//std::string filename("../5RockSamples/2DPorous.raw");
	//ReadBinaryImage(filename,533, 428, 1);// if 2D put 1 for z
//	std::string filename("../5RockSamples/cheese.raw");
//	ReadBinaryImage(filename,800, 600, 1);// if 2D put 1 for z
	//	ReadBinaryImage(filename,533, 428, 1);// if 2D put 1 for z

//	std::string filename("../5RockSamples/Processed_2x_2D_Tomography.raw");
//	ReadBinaryImage(filename,1065, 855, 1);// double resolution
//	std::cout<<"**** check read file***"<<std::endl;
}


UserMesh::~UserMesh() {
	// TODO Auto-generated destructor stub
}

void UserMesh::ChangeNode(Node2D &Node, bool &solid )
{
//Image
	/*
	if(Node.get_x()>-1 && Node.get_x()<1066 && Node.get_y()>-1 && Node.get_y()<855 )
		solid=image[0][(int)Node.get_y()][(int)Node.get_x()];
	else
		solid=false;
*/

//Circle
/*	if(sqrt(pow((Node.get_x()-150),2.0)+pow((Node.get_y()-150),2.0))<=10)
		solid=true;
	else
		solid=false;
		*/


//Square
/*	if(Node.get_y()>=-1 && Node.get_y()<=250 && Node.get_x()>=250 && Node.get_x()<=750)
		{
			solid=true;
		}
		else
		{
			solid=false;
		}
*/

//Cross
/*	double xmiddle=99.5;
	double ymiddle=99.5;
	double factor=6;
	double widthsquare=(factor-1)*2+1;
	double xfirstcorner=xmiddle-0.5-factor;
	double yfirstcorner=ymiddle-0.5-factor;

	if(Node.get_x()>=xfirstcorner-widthsquare && Node.get_x()<=xfirstcorner+2*widthsquare && Node.get_y()>=yfirstcorner-widthsquare && Node.get_y()<=yfirstcorner+2*widthsquare)
		if((Node.get_x()>=xfirstcorner && Node.get_x()<=xfirstcorner+widthsquare)||( Node.get_y()>=yfirstcorner && Node.get_y()<=yfirstcorner+widthsquare))
		{
			solid=true;
		}
		else
		{
			solid=false;
		}
	else
		solid=false;*/

/*	if(Node.get_x()>=140 && Node.get_x()<=260 && Node.get_y()>=140 && Node.get_y()<=260)
		if((Node.get_y()>=170 && Node.get_y()<=230)||(Node.get_x()>=170 && Node.get_x()<=230 ))
		{
			//std::cout<<" NodeID: "<<Node.Get_index()<<" x: "<<Node.get_x()<< " y: "<<Node.get_y()<<std::endl;
			solid=true;
		}
		else
		{
			solid=false;
		}*/

/*if (Node.get_x()==70&&Node.get_y()==213)
	if(solid==true)
		std::cout<<"******Solid*****"<<std::endl;
	else
		std::cout<<"------Fluid------"<<std::endl;*/
}
void UserMesh::SetSymmetryType(SymmetryType &Type, double x, double y)
{
//	if(x==0||x>=400)
//		Type=SymmetryPressureOnNode;

}
