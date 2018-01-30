/*
 * Extrapolation.h
 *
 *  Created on: 31 October 2016
 *      Author: Thomas Burel
 *
 *  Files created to manage Extrapolation types.
 *  To add a new Extrapolation type, you need to:
 *  	 Add a class in ExtrapolationDEF.h
 *  	 Defined all virtual methods in ExtrapolationDEF.cpp
 *  	 Add an enumtype
 *  	 Modified the SelectExtrapolationType function
 */

#ifndef ALGORITHM_LOWORDER_ExtrapolationS_H_
#define ALGORITHM_LOWORDER_ExtrapolationS_H_

#include "ExtrapolationDEF.h"

class Extrapolation {
public:
	Extrapolation();
	void initExtrapolation(int dimension, int nb_vel,ModelEnum::ExtrapolationType Type_);
	virtual ~Extrapolation();
	void SelectExtrapolationType(ModelEnum::ExtrapolationType Type_);

//Scalar Extrapolations
	void ExtrapolationOnWall (double * & Var, int * Connect, int & normal);
	void ExtrapolationOnCornerConcave (double * & Var, int * Connect, int & normal);
	void ExtrapolationOnCornerConvex (double * & Var, int * Connect, int & normal);

	void ExtrapolationWallToSolid (double * & Var, int * Connect, int & normal);
	void ExtrapolationCornerConcaveToSolid (double * & Var, int * Connect, int & normal);
	void ExtrapolationCornerConvexToSolid (double * & Var, int * Connect, int & normal);

private:
 // Extrapolation object to be able to select different kind of Extrapolation automatically
 ExtrapolationDEF* Extrapol;
 ModelEnum::ExtrapolationType Type;
 int dimension, nb_Vel;
};
#endif /* ALGORITHM_LOWORDER_ExtrapolationS_H_ */
