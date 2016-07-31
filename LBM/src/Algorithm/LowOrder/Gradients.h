/*
 * Gradients.h
 *
 *  Created on: 30 Jul 2016
 *      Author: Thomas Burel
 *
 *  Files created to manage gradient types.
 *  To add a new gradient type, you need to:
 *  	 Add a class in GradientsDEF.h
 *  	 Defined all virtual methods in GradientsDEF.cpp
 *  	 Add an enumtype
 *  	 Modified the SelectGradientType function
 */

#ifndef ALGORITHM_LOWORDER_GRADIENTS_H_
#define ALGORITHM_LOWORDER_GRADIENTS_H_

#include "GradientsDEF.h"

enum GradientType{FD};

class Gradients {
public:
	Gradients();
	Gradients(int dimension);
	virtual ~Gradients();
	void SelectGradientType(GradientType Type_);

//Scalar gradients
	double* Grad (double *Var, NodeInterior2D& Node);
	double* Grad (double *Var, NodeWall2D& Node);
	double* Grad (double *Var, NodeCorner2D& Node);
	double* Grad (double *Var, NodeVelocity2D& Node);
	double* Grad (double *Var, NodePressure2D& Node);
	double* Grad (double *Var, NodeSymmetry2D& Node);
//Vector gradients
	double** Grad (double *Var_x, double *Var_y, NodeInterior2D& Node);
	double** Grad (double *Var_x, double *Var_y, NodeWall2D& Node);
	double** Grad (double *Var_x, double *Var_y, NodeCorner2D& Node);
	double** Grad (double *Var_x, double *Var_y, NodeVelocity2D& Node);
	double** Grad (double *Var_x, double *Var_y, NodePressure2D& Node);
	double** Grad (double *Var_x, double *Var_y, NodeSymmetry2D& Node);

private:
 // Gradient object to be able to select different kind of gradient automatically
 GradientsDEF* grad;
 GradientType Type;
 int dimension;
};
#endif /* ALGORITHM_LOWORDER_GRADIENTS_H_ */
