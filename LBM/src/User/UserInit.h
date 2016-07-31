/*
 * UserInitLBM.h
 *
 *  Created on: 12 Jun 2015
 *      Author: thomas
 */

#ifndef USER_USERINIT_H_
#define USER_USERINIT_H_
#include "../Mesh/SingleBlock.h"
#include "UserParameters.h"
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

class UserInit {
private:
	friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
        ar & Block_;
    }
public:
	UserInit();
	virtual ~UserInit();
	///Call for initialise with specific to a case
	void Set_UserInit();

protected:
	/// Call for initialise boundary conditions
	void UserBc(UserParameters& PtrUserParameters,int elem, int nodenumber, double* pos,double& Rho, double* U,double& alpha);
	/// Call for initialise internal conditions
	void UserIc (UserParameters& PtrUserParameters, int elem, int nodenumber, double* pos ,double& Rho, double* U,double& alpha);



private:
	Block *Block_;
};

#endif /* USER_USERINIT_H_ */
