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


//	std::string filename("../P0.8G1.raw");
//	ReadBinaryImage(filename,204, 204, 1);// double resolution
/*	std::string filename("../Processed_2x_2D_Tomography.raw");
	ReadBinaryImage(filename,1065, 855, 1);// double resolution
	std::cout<<"**** check read file***"<<std::endl;
*/
//	std::string filename("../Pipe_1290_455.raw");

//	std::string filename("../Perpentine_1290_371.raw");
//		ReadBinaryImage(filename,1290, 371, 1);// double resolution
//		std::cout<<"**** check read file***"<<std::endl;

//	std::string filename("../rooster_363_470.raw");
//		ReadBinaryImage(filename,363, 470, 1);// double resolution
	std::string filename("../Image_1162_filtered_filled_404_201.raw");
		ReadBinaryImage(filename,404, 201, 1);// double resolution
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
/*	if(Node.get_x()>10 && Node.get_x()<1075 && Node.get_y()>10 && Node.get_y()<865 )
		solid=image[0][(int)Node.get_y()-10][(int)Node.get_x()-10];
	else
		solid=false;*/

/*
	if(Node.get_x()>=(L/2.0-R) &&Node.get_x()<=(L/2.0+R) && Node.get_y()>=(-R) &&Node.get_y()<=(+R) )
		solid=true;
	else
		solid=false;
*/

/*	if(Node.get_y()<=a*Node.get_x()+b+epsilon)
		solid=true;
	else
		solid=false;*/
//	if(Node.get_y()<=H/2.0+std::cos(teta)*R && Node.get_x()>=10 &&Node.get_x()<=390 && Node.get_y()>=10 )
/*	if(Node.get_y()<=a*(Node.get_x()+std::sin(beta)*R)+b+std::cos(beta)*std::cos(teta)*R && Node.get_x()>=10 &&Node.get_x()<=390 && Node.get_y()>=10 )

			solid=true;
		else
			solid=false;*/
/*
	if(pow(Node.get_x()-L/2.0,2.0)+pow(Node.get_y(),2.0)<=R*R+1e-8)
		solid=true;
	else
		solid=false;
*/
	//Cross
/*	double xmiddle=(L+1)/2.0+0.5;
	double ymiddle=(H+1)/2.0+0.5;
	double factor=6;
	double widthsquare=(factor-1.0)*2+1.0;
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
	/*if(Node.get_x()==44 && Node.get_y()==44)
		std::cout<<"x: "<<Node.get_x()<<" y: "<<Node.get_y()<<" Solid type: "<<solid<<std::endl;
	if(Node.get_x()==43 && Node.get_y()==44)
		std::cout<<"x: "<<Node.get_x()<<" y: "<<Node.get_y()<<" Solid type: "<<solid<<std::endl;
	if(Node.get_x()==44 && Node.get_y()==43)
		std::cout<<"x: "<<Node.get_x()<<" y: "<<Node.get_y()<<" Solid type: "<<solid<<std::endl;
	if(Node.get_x()==43 && Node.get_y()==43)
		std::cout<<"x: "<<Node.get_x()<<" y: "<<Node.get_y()<<" Solid type: "<<solid<<std::endl;*/
/*
	if(Node.get_y()<=(H/2.0-R) || Node.get_y()>=(H/2.0+R) )
		solid=true;
	else
		solid=false;*/


/*
	if(Node.get_x()>-1&& Node.get_x()<1290 && Node.get_y()>-1 && Node.get_y()<371 )
			solid=image[0][(int)Node.get_y()][(int)Node.get_x()];
		else
			solid=false;
	*/



	//if((Node.get_x()<10 || Node.get_x()>90) &&(Node.get_y()<10 || Node.get_y()>40 ))
/*	if(Node.get_y()<20 || Node.get_y()>40)
		if((Node.get_x()>140 && Node.get_x()<160)||(Node.get_x()>190 && Node.get_x()<210)||(Node.get_x()>240 && Node.get_x()<260)||(Node.get_x()>290 && Node.get_x()<310)||(Node.get_x()>340 && Node.get_x()<360))
			solid=true;
		else
			solid=false;
	else
		solid=false;
	//force walls
	if(Node.get_x()>10 && Node.get_x()<490)
	if(Node.get_y()<=1 ||Node.get_y()>=59)
		solid=true;
//	if(Node.get_y()<2 || Node.get_y()>18)
//		solid=true;
*/
	if(Node.get_x()>100 && Node.get_x()<504)
		if(Node.get_y()<=0 || Node.get_y()>=194)
			solid=true;
		else
			solid=image[0][(int)Node.get_y()+1][(int)Node.get_x()-100];
		else
			solid=false;
}
void UserMesh::SetSymmetryType(SymmetryType &Type, double x, double y)
{

}
