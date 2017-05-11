/*
 * Interpolation.h
 *
 *  Created on: 31 October 2016
 *      Author: Thomas Burel
 *
 *  Files created to manage Interpolation types.
 *  To add a new Interpolation type, you need to:
 *  	 Add a class in InterpolationDEF.h
 *  	 Defined all virtual methods in InterpolationDEF.cpp
 *  	 Add an enumtype
 *  	 Modified the SelectInterpolationType function
 */

#ifndef ALGORITHM_LOWORDER_INTERPOLATIONS_H_
#define ALGORITHM_LOWORDER_INTERPOLATIONS_H_

#include "InterpolationDEF.h"

class Interpolation {
public:
	Interpolation();
	void initInterpolation(int dimension, int nb_vel,ModelEnum::InterpolationType Type_,unsigned int *PtrOppositeCa_, NodeArrays2D *PtrNodes, Parameters *PtrParam);
	virtual ~Interpolation();
	void SelectInterpolationType(ModelEnum::InterpolationType Type_, NodeArrays2D *PtrNodes, Parameters *PtrParam);


//Scalar Interpolations
	void InterpolationOnWall (double * & Var, int * Connect, int & normal);
	void InterpolationOnCornerConcave (double * & Var, int * Connect, int & normal);
	void InterpolationOnCornerConvex (double * & Var, int * Connect, int & normal);

	void InterpolationOnWall (double * & Var1, double * & Var2, int * Connect, int & normal);
	void InterpolationOnCornerConcave (double * & Var1, double * & Var2, int * Connect, int & normal);
	void InterpolationOnCornerConvex (double * & Var1, double * & Var2, int * Connect, int & normal);
private:
 // Interpolation object to be able to select different kind of Interpolation automatically
 InterpolationDEF* Interpol;
 ModelEnum::InterpolationType Type;
 int dimension, nb_Vel;
 unsigned int *PtrOppositeCa; //opposite direction in the distribution function
};
#endif /* ALGORITHM_LOWORDER_INTERPOLATIONS_H_ */
