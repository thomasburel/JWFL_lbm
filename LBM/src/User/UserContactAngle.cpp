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
	desireCangle=90;
	beginCangle=10;
	lengthsmooth=20;
	convDegreeToRadian=pi/180.0;
}

UserContactAngle::~UserContactAngle() {
	// TODO Auto-generated destructor stub
}

void UserContactAngle::Set_ContactAngle(Parameters& PtrParameters, int elem, int nodenumber, double* pos ,double& teta){
	if(pos[0]<lengthsmooth)
		teta=(-pos[0]*(beginCangle-desireCangle)/lengthsmooth+beginCangle)*convDegreeToRadian;
	else
		teta=desireCangle*convDegreeToRadian;
}
