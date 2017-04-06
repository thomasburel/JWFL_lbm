/*
 * ============================================================================
 * UserContactAngle.cpp
 *
 *  Created on: 12 Oct 2016
 *      Author: Thomas Burel
 *      Version     :
 *      Copyright   : Your copyright notice
 *      Description : LBM Code Developing by JWFL at University of Strathclyde
 * ============================================================================
 */

#include "UserContactAngle.h"

UserContactAngle::UserContactAngle() {
	// TODO Auto-generated constructor stub
	pi=atan(1.0)*4.0 ;
}

UserContactAngle::~UserContactAngle() {
	// TODO Auto-generated destructor stub
}

void UserContactAngle::Set_ContactAngle(Parameters& PtrParameters, int elem, int nodenumber, double* pos ,double& teta){
	if(pos[0]<50)
		teta=160.0*pi/180.0;
	else
		teta=30.0*pi/180.0;
}
