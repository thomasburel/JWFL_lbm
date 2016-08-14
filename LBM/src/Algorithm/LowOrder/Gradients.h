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

class Gradients {
public:
	Gradients();
	void initGradients(int dimension, int nb_vel,GradientType Type_);
	virtual ~Gradients();
	void SelectGradientType(GradientType Type_);

//Scalar gradients
	double* Grad (double *Var, int * Connect, int & normal);
	double* GradBc (double *Var, int * Connect, int & normal);
	double* GradCorner (double *Var, int * Connect, int & normal);

//Vector gradients
	double** Grad (double *Var_x, double *Var_y, int * Connect, int & normal);
	double** GradBc (double *Var_x, double *Var_y, int * Connect, int & normal);
	double** GradCorner (double *Var_x, double *Var_y, int * Connect, int & normal);


private:
 // Gradient object to be able to select different kind of gradient automatically
 GradientsDEF* grad;
 GradientType Type;
 int dimension, nb_Vel;
};
#endif /* ALGORITHM_LOWORDER_GRADIENTS_H_ */
