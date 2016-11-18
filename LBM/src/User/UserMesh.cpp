/*
 * UserMesh.cpp
 *
 *  Created on: 17 Oct 2015
 *      Author: thomas
 */

#include "UserMesh.h"

UserMesh::UserMesh() {
	epsilon=1e-10;
	pi=std::atan(1.0)*4.0 ;

	H=200;
	L=200;
	R=H/4.0;
	//Contact angle
	teta=50.0*pi/180.0;
	//tilt angle of the plate
	beta=30.0*pi/180.0;

	a=std::tan(beta);
	b=(H-L*a)/2.0;

//	std::string filename("../P0.8G1.raw");
//	ReadBinaryImage(filename,204, 204, 1);// double resolution
/*	std::string filename("../Processed_2x_2D_Tomography.raw");
	ReadBinaryImage(filename,1065, 855, 1);// double resolution
	std::cout<<"**** check read file***"<<std::endl;
*/
}


UserMesh::~UserMesh() {
	// TODO Auto-generated destructor stub
}

void UserMesh::ChangeNode(Node2D &Node, bool &solid )
{
/*	if(Node.get_x()>10 && Node.get_x()<1075 && Node.get_y()>10 && Node.get_y()<865 )
		solid=image[0][(int)Node.get_y()-10][(int)Node.get_x()-10];
	else
		solid=false;*/
	if(Node.get_x()>=20 &&Node.get_x()<=180 && Node.get_y()>=20 &&Node.get_y()<=180 )
		solid=true;
	else
		solid=false;

/*	if(Node.get_y()<=a*Node.get_x()+b+epsilon)
		solid=true;
	else
		solid=false;*/
//	if(Node.get_y()<=H/2.0+std::cos(teta)*R && Node.get_x()>=10 &&Node.get_x()<=390 && Node.get_y()>=10 )
/*	if(Node.get_y()<=a*(Node.get_x()+std::sin(beta)*R)+b+std::cos(beta)*std::cos(teta)*R && Node.get_x()>=10 &&Node.get_x()<=390 && Node.get_y()>=10 )

			solid=true;
		else
			solid=false;*/
/*	if(pow(Node.get_x()-L/2.0,2.0)+pow(Node.get_y()-H/2.0,2.0)<=R*R+1e-8)
		solid=true;
	else
		solid=false;*/
}
void UserMesh::SetSymmetryType(SymmetryType &Type, double x, double y)
{

}
