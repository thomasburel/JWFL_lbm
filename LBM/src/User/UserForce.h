/*
 * UserForce.h
 *
 *  Created on: 26 Apr 2016
 *      Author: thomas
 */

#ifndef USER_USERFORCE_H_
#define USER_USERFORCE_H_
#include "../Core/Parameters.h"
class UserForce {
public:
	UserForce();
	virtual ~UserForce();
	void Set_UserForce(Parameters* PtrParameters);
//protected:
	double LocalForce(int const direction_i, double const Rho, double const U, double const V, double const W=0);
	double BodyForce(int const direction_xyz, double const Rho, double const U, double const  V, double const W=0);

private:
	double UserBodyForce[2];
	double UserLocalForceD9[9];
};

#endif /* USER_USERFORCE_H_ */
