/*
 * ============================================================================
 * UserContactAngle.h
 *
 *  Created on: 12 Oct 2016
 *      Author: Thomas Burel
 *      Version     :
 *      Copyright   : Your copyright notice
 *      Description : LBM Code Developing by JWFL at University of Strathclyde
 * ============================================================================
 */

#ifndef SRC_USER_USERCONTACTANGLE_H_
#define SRC_USER_USERCONTACTANGLE_H_
#include "../Core/Parameters.h"
class UserContactAngle {
public:
	UserContactAngle();
	virtual ~UserContactAngle();

protected:
	void Set_ContactAngle(Parameters& PtrParameters, int elem, int nodenumber, double* pos ,double& teta);
	double desireCangle;
private:
	double pi;
	double beginCangle;
	double lengthsmooth;
	double convDegreeToRadian;
};

#endif /* SRC_USER_USERCONTACTANGLE_H_ */
