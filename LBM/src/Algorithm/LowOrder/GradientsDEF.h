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
	GradientsDEF(int dimension, int nb_vel);
	virtual ~GradientsDEF();
//Scalar gradient
	virtual void Grad (double* grad_,double *Var, int * Connect, int & normal)=0;
	virtual void GradBc (double* grad_,double *Var, int * Connect, int & normal)=0;
	virtual void GradCorner (double* grad_,double *Var, int * Connect, int & normal)=0;

//Vector gradient
	virtual void Grad (double** grad_, double *Var_x, double *Var_y, int * Connect, int & normal)=0;
	virtual void GradBc (double** grad_, double *Var_x, double *Var_y, int * Connect, int & normal)=0;
	virtual void GradCorner (double** grad_, double *Var_x, double *Var_y, int * Connect, int & normal)=0;


protected:
//	double* gradient_scalar;
//	double** gradient_vector;
	int dimension,nb_Vel;
};
class GradientsFD: public GradientsDEF {
public:
	GradientsFD();
	GradientsFD(int dimension, int nb_vel);
	virtual ~GradientsFD();
//Scalar gradient
	void Grad (double* grad_,double *Var, int * Connect, int & normal);
	void GradBc (double* grad_,double *Var, int * Connect, int & normal);
	void GradCorner (double* grad_,double *Var, int * Connect, int & normal);

//Vector gradient
	void Grad (double** grad_,double *Var_x, double *Var_y, int * Connect, int & normal);
	void GradBc (double** grad_,double *Var_x, double *Var_y, int * Connect, int & normal);
	void GradCorner (double** grad_,double *Var_x, double *Var_y, int * Connect, int & normal);
};
class GradientsLBMStencil: public GradientsDEF {
public:
	GradientsLBMStencil();
	GradientsLBMStencil(int dimension, int nb_vel);
	virtual ~GradientsLBMStencil();
//Scalar gradient
	void Grad (double* grad_,double *Var, int * Connect, int & normal);
	void GradBc (double* grad_,double *Var, int * Connect, int & normal);
	void GradCorner (double* grad_,double *Var, int * Connect, int & normal);

//Vector gradient
	void Grad (double** grad_, double *Var_x, double *Var_y, int * Connect, int & normal);
	void GradBc (double** grad_, double *Var_x, double *Var_y, int * Connect, int & normal);
	void GradCorner (double** grad_, double *Var_x, double *Var_y, int * Connect, int & normal);

private:
	double **Ei;//Ei[dimension][nb_velocity]
	double *omega;//omega[nb_velocity]
};
#endif /* ALGORITHM_LOWORDER_GRADIENTSDEF_H_ */
