/*
 * ============================================================================
 * D2Q9GlobalCorner.h
 *
 *  Created on: 2 Sep 2016
 *      Author: Thomas Burel
 *      Version     :
 *      Copyright   : Your copyright notice
 *      Description : LBM Code Developing by JWFL at University of Strathclyde
 * ============================================================================
 */

#ifndef SRC_ALGORITHM_LOWORDER_BOUNDARIES_D2Q9GLOBALCORNER_H_
#define SRC_ALGORITHM_LOWORDER_BOUNDARIES_D2Q9GLOBALCORNER_H_

#include "D2Q9BcVar.h"
#include "D2Q9GenericBc.h"


class D2Q9GlobalCorner: public D2Q9BcVar {
public:
	D2Q9GlobalCorner();
	virtual ~D2Q9GlobalCorner();

	void Set_GlobalCorner(Parameters *Param,D2Q9GenericBc* D2Q9GenericBc);
	void ApplyGlobalCorner(NodeCorner2D& Node, std::map<int,NodeType> TypeOfNode_, DistriFunct* f_in);
	//Specify in the solver: set values (Rho, U) and pointers on macroscopic variables
	void ApplyGlobalCorner(NodeCorner2D& Node, double const Rho_def, double const *UDef, std::map<int,NodeType> TypeOfNode_, DistriFunct* f_in,double * & Rho, double * &U, double * &V);

private:
	void FunctionGlobalCorner(NodeCorner2D& Node, double const Rho_def, double const *UDef, std::map<int,NodeType> &TypeOfNode_, DistriFunct* &f_in,double * & Rho, double * &U, double * &V);
	D2Q9GenericBc* BcMethods;
};

#endif /* SRC_ALGORITHM_LOWORDER_BOUNDARIES_D2Q9GLOBALCORNER_H_ */
