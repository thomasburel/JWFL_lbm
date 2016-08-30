/*
 * ============================================================================
 * D2Q9Boundaries.h
 *
 *  Created on: 29 Aug 2016
 *      Author: Thomas Burel
 *      Version     :
 *      Copyright   : Your copyright notice
 *      Description : LBM Code Developing by JWFL at University of Strathclyde
 * ============================================================================
 */

#ifndef SRC_ALGORITHM_LOWORDER_BOUNDARIES_D2Q9BOUNDARIES_H_
#define SRC_ALGORITHM_LOWORDER_BOUNDARIES_D2Q9BOUNDARIES_H_

#include "D2Q9Pressure.h"
#include "D2Q9Velocity.h"
#include "D2Q9Wall.h"
#include "D2Q9Corner.h"
#include "../../../Core/Dictionary.h"

class D2Q9Boundaries {
public:
	D2Q9Boundaries();
	virtual ~D2Q9Boundaries();

protected:
	void ApplyPressure(NodePressure2D& Node, double* fi);
	void ApplyVelocity(NodeVelocity2D& Node, double* fi);
	void ApplyCorner(NodeCorner2D& Node, double* fi);
	void ApplyWall(NodeWall2D& Node, double *fi);

private:
	double *Rho, *U, *V;
	Dictionary *PtrDic;
	D2Q9Pressure D2Q9PressureBc;
	D2Q9Velocity D2Q9VelocityBc;
	D2Q9Corner D2Q9CornerBc;
	D2Q9Wall D2Q9WallBc;
};

#endif /* SRC_ALGORITHM_LOWORDER_BOUNDARIES_D2Q9BOUNDARIES_H_ */
