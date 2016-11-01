/*
 * UserParameters.cpp
 *
 *  Created on: 10 Jun 2015
 *      Author: thomas
 */

#include "UserParameters.h"

UserParameters::UserParameters() {

	Umax=0.01;
	H=100;
	L=100;
	Pmax=1.05;
	Pmin=1.0;
	ReUser=0;
	CaUser=0;
	DiameterUser=0;
	sigmaUser=0;
	ContactAngleUser=0;
}

UserParameters::~UserParameters() {
	// TODO Auto-generated destructor stub
}

void UserParameters::Set_UserParameters(double  Umax_,double  H_,double  L_,double  Pmax_,double  Pmin_){
	Umax=Umax_;
	H=H_;
	L=L_;
	Pmax=Pmax_;
	Pmin=Pmin_;
}
void UserParameters::Get_UserParameters(double & Umax_,double & H_,double & L_,double & Pmax_,double & Pmin_){
	Umax_=Umax;
	H_=H;
	L_=L;
	Pmax_=Pmax;
	Pmin_=Pmin;

}
void UserParameters::Set_UserDroplets(double angle, double sigma, double diameter, double Re, double Ca){
	ContactAngleUser=angle;
	sigmaUser=sigma;
	DiameterUser=diameter;
	ReUser=Re;
	CaUser=Ca;
}
void UserParameters::Get_UserDroplets(double & angle, double & sigma, double & diameter, double & Re, double & Ca){
	angle=ContactAngleUser;
	sigma=sigmaUser;
	diameter=DiameterUser;
	Re=ReUser;
	Ca=CaUser;
}
