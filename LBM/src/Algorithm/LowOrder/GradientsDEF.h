/*
 * GradientsDEF.h
 *
 *  Created on: 28 Jul 2016
 *      Author: Thomas Burel
 */

#ifndef ALGORITHM_LOWORDER_GRADIENTSDEF_H_
#define ALGORITHM_LOWORDER_GRADIENTSDEF_H_

#include "../../Core/Parameters.h"

///Abstract class for common function for all kind of gradients
class GradientsDEF {
public:
	GradientsDEF();
//	virtual GradientsDEF(int dimension, int nb_vel)=0;
	virtual ~GradientsDEF();
//Scalar gradient
	virtual double* Grad (double *Var, int * Connect, int & normal)=0;
	virtual double* GradBc (double *Var, int * Connect, int & normal)=0;
	virtual double* GradCorner (double *Var, int * Connect, int & normal)=0;

//Vector gradient
	virtual double** Grad (double *Var_x, double *Var_y, int * Connect, int & normal)=0;
	virtual double** GradBc (double *Var_x, double *Var_y, int * Connect, int & normal)=0;
	virtual double** GradCorner (double *Var_x, double *Var_y, int * Connect, int & normal)=0;


protected:
	double* gradient_scalar;
	double** gradient_vector;
	int dimension,nb_Vel;
};
class GradientsFD: public GradientsDEF {
public:
	GradientsFD();
	GradientsFD(int dimension, int nb_vel);
	virtual ~GradientsFD();
//Scalar gradient
	double* Grad (double *Var, int * Connect, int & normal);
	double* GradBc (double *Var, int * Connect, int & normal);
	double* GradCorner (double *Var, int * Connect, int & normal);

//Vector gradient
	double** Grad (double *Var_x, double *Var_y, int * Connect, int & normal);
	double** GradBc (double *Var_x, double *Var_y, int * Connect, int & normal);
	double** GradCorner (double *Var_x, double *Var_y, int * Connect, int & normal);
};
class GradientsLBMStencil: public GradientsDEF {
public:
	GradientsLBMStencil();
	GradientsLBMStencil(int dimension, int nb_vel);
	virtual ~GradientsLBMStencil();
//Scalar gradient
	double* Grad (double *Var, int * Connect, int & normal);
	double* GradBc (double *Var, int * Connect, int & normal);
	double* GradCorner (double *Var, int * Connect, int & normal);

//Vector gradient
	double** Grad (double *Var_x, double *Var_y, int * Connect, int & normal);
	double** GradBc (double *Var_x, double *Var_y, int * Connect, int & normal);
	double** GradCorner (double *Var_x, double *Var_y, int * Connect, int & normal);

private:
	double **Ei;//Ei[dimension][nb_velocity]
	double *omega;//omega[nb_velocity]
};
#endif /* ALGORITHM_LOWORDER_GRADIENTSDEF_H_ */
