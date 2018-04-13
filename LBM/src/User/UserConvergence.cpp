/*
 * UserConvergence.cpp
 *
 *  Created on: 8 Mar 2018
 *      Author: Thomas Burel
 *

 *
 *  Notes: element number is not used at this time and set to 0.
 */

#include "UserConvergence.h"

UserConvergence::UserConvergence() {

	Hchannel=1;
	visco=0.1;
	Re=0.001;
}

UserConvergence::~UserConvergence() {
	// TODO Auto-generated destructor stub
}
void UserConvergence::Set_UserConvergence(Parameters* PtrParameters){
	Hchannel=PtrParameters->Get_HchannelUser();
	visco=PtrParameters->Get_ViscoUser();
	Re=PtrParameters->Get_ReUser();
}
//Give the error as abs(U(lbm) - U(exact))
double UserConvergence::UserError(int nodenumber, double* pos,double& Rho, double* U, double& alpha){
	//return 0;
	return std::abs(U[0]-Re*(visco/Hchannel)*(1.0-pow(2.0*(pos[1]-Hchannel/2.0-1.0)/(Hchannel),2.0)));



}
//Give the exact result as U(exact)
double UserConvergence::UserExact(int nodenumber, double* pos,double& Rho, double* U, double& alpha){
	//return 1;
	return std::abs(Re*(visco/Hchannel)*(1.0-pow(2.0*(pos[1]-Hchannel/2.0-1.0)/(Hchannel),2.0)));

}
