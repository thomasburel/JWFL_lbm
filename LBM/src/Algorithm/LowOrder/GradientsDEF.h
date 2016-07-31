/*
 * GradientsDEF.h
 *
 *  Created on: 28 Jul 2016
 *      Author: Thomas Burel
 */

#ifndef ALGORITHM_LOWORDER_GRADIENTSDEF_H_
#define ALGORITHM_LOWORDER_GRADIENTSDEF_H_

#include "../../Mesh/SingleBlock.h"

///Abstract class for common function for all kind of gradients
class GradientsDEF {
public:
	GradientsDEF();
	GradientsDEF(int dimension);
	virtual ~GradientsDEF();
//Scalar gradient
	virtual double* Grad (double *Var, NodeInterior2D& Node)=0;
	virtual double* Grad (double *Var, NodeWall2D& Node)=0;
	virtual double* Grad (double *Var, NodeCorner2D& Node)=0;
	virtual double* Grad (double *Var, NodeVelocity2D& Node)=0;
	virtual double* Grad (double *Var, NodePressure2D& Node)=0;
	virtual double* Grad (double *Var, NodeSymmetry2D& Node)=0;
//Vector gradient
	virtual double** Grad (double *Var_x, double *Var_y, NodeInterior2D& Node)=0;
	virtual double** Grad (double *Var_x, double *Var_y, NodeWall2D& Node)=0;
	virtual double** Grad (double *Var_x, double *Var_y, NodeCorner2D& Node)=0;
	virtual double** Grad (double *Var_x, double *Var_y, NodeVelocity2D& Node)=0;
	virtual double** Grad (double *Var_x, double *Var_y, NodePressure2D& Node)=0;
	virtual double** Grad (double *Var_x, double *Var_y, NodeSymmetry2D& Node)=0;

protected:
	double* gradient_scalar;
	double** gradient_vector;
	int dimension;
};
class GradientsFD: public GradientsDEF {
public:
	GradientsFD();
	GradientsFD(int dimension);
	virtual ~GradientsFD();
//Scalar gradient
	double* Grad (double *Var, NodeInterior2D& Node);
	double* Grad (double *Var, NodeWall2D& Node);
	double* Grad (double *Var, NodeCorner2D& Node);
	double* Grad (double *Var, NodeVelocity2D& Node);
	double* Grad (double *Var, NodePressure2D& Node);
	double* Grad (double *Var, NodeSymmetry2D& Node);
//Vector gradient
	double** Grad (double *Var_x, double *Var_y, NodeInterior2D& Node);
	double** Grad (double *Var_x, double *Var_y, NodeWall2D& Node);
	double** Grad (double *Var_x, double *Var_y, NodeCorner2D& Node);
	double** Grad (double *Var_x, double *Var_y, NodeVelocity2D& Node);
	double** Grad (double *Var_x, double *Var_y, NodePressure2D& Node);
	double** Grad (double *Var_x, double *Var_y, NodeSymmetry2D& Node);

};

#endif /* ALGORITHM_LOWORDER_GRADIENTSDEF_H_ */
