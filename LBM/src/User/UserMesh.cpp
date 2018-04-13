/*
 * UserMesh.cpp
 *
 *  Created on: 17 Oct 2015
 *      Author: thomas
 */

#include "UserMesh.h"

UserMesh::UserMesh() {
	PtrParametersUserMesh=0;
	epsilon=1e-10;
	pi=std::atan(1.0)*4.0 ;
	H=10;
	L=10;
	R=H/4.0;
	//Contact angle
	teta=pi/2;
	//tilt angle of the plate
	beta=0;
	a=std::tan(beta);
	b=(H-L*a)/2.0;

//	std::string filename("../serpentine2_58_19.raw");	
//	std::string filename("../serpentine2_119_38.raw");
	std::string filename("../serpentine2_244_74.raw");
//	std::string filename("../Serpentine_432_132.raw");

	L_Image=244;//58;
	H_Image=74;//19;
	ReadBinaryImage(filename,L_Image, H_Image, 1);
	//std::cout<<"**** check read file***"<<std::endl;
}
void UserMesh::SetUserMeshVariables(){
	H=PtrParametersUserMesh->Get_UserH();
	L=PtrParametersUserMesh->Get_UserL();
	R=H/4.0;
	//Contact angle
	teta=50.0*pi/180.0;
	//tilt angle of the plate
	beta=30.0*pi/180.0;

	a=std::tan(beta);
	b=(H-L*a)/2.0;
}

UserMesh::~UserMesh() {
	// TODO Auto-generated destructor stub
}

void UserMesh::ChangeNode(Node2D &Node, bool &solid )
{


	if(Node.get_x()>-1&& Node.get_x()<L_Image && Node.get_y()>-1 && Node.get_y()<H_Image )
		solid=image[0][(int)Node.get_y()][(int)Node.get_x()];
	else if(Node.get_y()<=0 ||Node.get_y()>=H_Image)
		solid=true;
	else
		solid=false;

}
void UserMesh::SetSymmetryType(SymmetryType &Type, double x, double y)
{

}
