/*
 * UserForce.cpp
 *
 *  Created on: 26 Apr 2016
 *      Author: thomas
 */

#include "UserForce.h"

UserForce::UserForce() {
	// TODO Auto-generated constructor stub
	UserBodyForce[0]=0.0000001953125;
	UserBodyForce[1]=0;
	for (int i=0;i<9;i++)
		UserLocalForceD9[i]=0;
}

UserForce::~UserForce() {
	// TODO Auto-generated destructor stub
}
void UserForce::Set_UserForce(Parameters* PtrParameters){
	double Fx=8.0*PtrParameters->Get_ReUser()*pow(PtrParameters->Get_ViscoUser(),2.0)
				/pow(PtrParameters->Get_HchannelUser(),3.0);
	UserBodyForce[0]= Fx;
	UserBodyForce[1]=0;
	UserLocalForceD9[0]= 0.0;
	UserLocalForceD9[1]= Fx*3.0*(1.0/9.0);
	UserLocalForceD9[2]= 0.0;
	UserLocalForceD9[3]= -Fx*3.0*(1.0/9.0);
	UserLocalForceD9[4]= 0.0;
	UserLocalForceD9[5]= Fx*3.0*(1.0/36.0);
	UserLocalForceD9[6]= -Fx*3.0*(1.0/36.0);
	UserLocalForceD9[7]= -Fx*3.0*(1.0/36.0);
	UserLocalForceD9[8]= Fx*3.0*(1.0/36.0);
}
double UserForce::LocalForce(int const direction_i, double const Rho, double const U, double const V, double const W)
{
	return 0;
}

double UserForce::BodyForce(int const direction_xyz, double const Rho, double const U, double const  V, double const W)
{
	return 0;//UserBodyForce[direction_xyz];
}
