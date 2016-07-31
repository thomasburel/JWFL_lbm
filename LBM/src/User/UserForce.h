/*
 * UserForce.h
 *
 *  Created on: 26 Apr 2016
 *      Author: thomas
 */

#ifndef USER_USERFORCE_H_
#define USER_USERFORCE_H_

class UserForce {
public:
	UserForce();
	virtual ~UserForce();
protected:
	double LocalForce(int & direction, double & Rho, double & U, double & V, double & W);
	double BodyForce(double & Rho, double * U);
};

#endif /* USER_USERFORCE_H_ */
