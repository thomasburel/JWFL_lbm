/*
 * UserMesh.cpp
 *
 *  Created on: 17 Oct 2015
 *      Author: thomas
 */

#include "UserMesh.h"

UserMesh::UserMesh() {

}


UserMesh::~UserMesh() {
	// TODO Auto-generated destructor stub
}

void UserMesh::ChangeNode(Node2D &Node, bool &solid )
{
	if(Node.get_x()>=80 &&Node.get_x()<=90 && Node.get_y()>=15 &&Node.get_y()<=26 )
		solid=true;
	else
		solid=false;

}
void UserMesh::SetSymmetryType(SymmetryType &Type, double x, double y)
{

}
