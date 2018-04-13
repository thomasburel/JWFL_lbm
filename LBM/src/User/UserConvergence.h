/*
 * UserConvergence.h
 *
 *  Created on: 8 Mar 2018
 *      Author: thomas
 */

#ifndef USER_USERCONVERGENCE_H_
#define USER_USERCONVERGENCE_H_
#include <cmath>
#include <iostream>
#include "../Core/Parameters.h"
class UserConvergence{
private:

public:
	UserConvergence();
	virtual ~UserConvergence();
	void Set_UserConvergence(Parameters* PtrParameters);

protected:
	/// Call the exact solution
	double UserExact(int nodenumber, double* pos,double& Rho, double* U,double& alpha);
	double UserError(int nodenumber, double* pos,double& Rho, double* U,double& alpha);
private:
	double Re;
	double Hchannel;
	double visco;
};

#endif /* USER_USERCONVERGENCE_H_ */
