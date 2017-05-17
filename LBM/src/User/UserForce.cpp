/*
 * UserForce.cpp
 *
 *  Created on: 26 Apr 2016
 *      Author: thomas
 */

#include "UserForce.h"

UserForce::UserForce() {
	// TODO Auto-generated constructor stub
	UserBodyForce[0]=0.0000001;
	UserBodyForce[1]=0;
}

UserForce::~UserForce() {
	// TODO Auto-generated destructor stub
}
double UserForce::LocalForce(int const direction_i, double const Rho, double const U, double const V, double const W)
{
	return 0;
}

double UserForce::BodyForce(int const direction_xyz, double const Rho, double const U, double const  V, double const W)
{
	return UserBodyForce[direction_xyz];
}
