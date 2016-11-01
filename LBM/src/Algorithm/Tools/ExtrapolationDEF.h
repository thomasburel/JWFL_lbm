/*
 * ExtrapolationDEF.h
 *
 *  Created on: 31 Oct 2016
 *      Author: Thomas Burel
 */

#ifndef ALGORITHM_LOWORDER_EXTRAPOLATIONSDEF_H_
#define ALGORITHM_LOWORDER_EXTRAPOLATIONSDEF_H_

#include "../../Core/Parameters.h"

///Abstract class for common function for all kind of extrapolation
class ExtrapolationDEF {
public:
	ExtrapolationDEF();
	ExtrapolationDEF(int dimension, int nb_vel);
	virtual ~ExtrapolationDEF();
	//Extrapolation
	virtual void ExtrapolationWall (double *Var, int * Connect, int & normal)=0;
	virtual void ExtrapolationCornerConcave (double *Var, int * Connect, int & normal)=0;
	virtual void ExtrapolationCornerConvex (double *Var, int * Connect, int & normal)=0;

protected:
//	double* gradient_scalar;
//	double** gradient_vector;
	int dimension,nb_Vel;
};
class ExtrapolationTailor: public ExtrapolationDEF {
public:
	ExtrapolationTailor();
	ExtrapolationTailor(int dimension, int nb_vel);
	virtual ~ExtrapolationTailor();
//Extrapolation
	void ExtrapolationWall (double *Var, int * Connect, int & normal);
	void ExtrapolationCornerConcave (double *Var, int * Connect, int & normal);
	void ExtrapolationCornerConvex (double *Var, int * Connect, int & normal);

};
class ExtrapolationWeightDistance: public ExtrapolationDEF {
public:
	ExtrapolationWeightDistance();
	ExtrapolationWeightDistance(int dimension, int nb_vel);
	virtual ~ExtrapolationWeightDistance();
//Extrapolation
	void ExtrapolationWall (double *Var, int * Connect, int & normal);
	void ExtrapolationCornerConcave (double *Var, int * Connect, int & normal);
	void ExtrapolationCornerConvex (double *Var, int * Connect, int & normal);


private:
	double InvSqrt2;
	double InvSqrt2_5;
	double InvSumWeightWall;
	double InvSumWeightCornerConvex;
	double InvSumWeightCornerConcave;
	double InvSumWeightCornerWall;
};
#endif /* ALGORITHM_LOWORDER_EXTRAPOLATIONSDEF_H_ */
