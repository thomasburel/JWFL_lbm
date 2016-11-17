/*
 * ============================================================================
 * D2Q9Bc.h
 *
 *  Created on: 29 Aug 2016
 *      Author: Thomas Burel
 *      Version     :
 *      Copyright   : Your copyright notice
 *      Description : LBM Code Developing by JWFL at University of Strathclyde
 * ============================================================================
 */

#ifndef SRC_ALGORITHM_LOWORDER_BOUNDARIES_D2Q9BC_H_
#define SRC_ALGORITHM_LOWORDER_BOUNDARIES_D2Q9BC_H_

#include "D2Q9GenericBc.h"
#include "D2Q9GlobalCorner.h"
#include "D2Q9SpecialWall.h"

class D2Q9Bc: public D2Q9GenericBc, public D2Q9GlobalCorner, public D2Q9SpecialWall {
public:
	D2Q9Bc();
	virtual ~D2Q9Bc();
	void InitD2Q9Bc(Dictionary *PtrDic_, Parameters *Param);

};

#endif /* SRC_ALGORITHM_LOWORDER_BOUNDARIES_D2Q9BC_H_ */
