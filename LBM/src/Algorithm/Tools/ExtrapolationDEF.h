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
	virtual void ExtrapolationOnWall (double *Var, int * Connect, int & normal)=0;
	virtual void ExtrapolationOnCornerConcave (double *Var, int * Connect, int & normal)=0;
	virtual void ExtrapolationOnCornerConvex (double *Var, int * Connect, int & normal)=0;

	virtual void ExtrapolationWallToSolid (double *Var, int * Connect, int & normal)=0;
	virtual void ExtrapolationCornerConcaveToSolid (double *Var, int * Connect, int & normal)=0;
	virtual void ExtrapolationCornerConvexToSolid (double *Var, int * Connect, int & normal)=0;

protected:
	int dimension,nb_Vel;
};
class NoExtrapolation: public ExtrapolationDEF {
public:
	NoExtrapolation();
	NoExtrapolation(int dimension, int nb_vel);
	virtual ~NoExtrapolation();
//Extrapolation
	void ExtrapolationOnWall (double *Var, int * Connect, int & normal){};
	void ExtrapolationOnCornerConcave (double *Var, int * Connect, int & normal){};
	void ExtrapolationOnCornerConvex (double *Var, int * Connect, int & normal){};

	void ExtrapolationWallToSolid (double *Var, int * Connect, int & normal){};
	void ExtrapolationCornerConcaveToSolid (double *Var, int * Connect, int & normal){};
	void ExtrapolationCornerConvexToSolid (double *Var, int * Connect, int & normal){};

};
class ExtrapolationTailor: public ExtrapolationDEF {
public:
	ExtrapolationTailor();
	ExtrapolationTailor(int dimension, int nb_vel);
	virtual ~ExtrapolationTailor();
//Extrapolation
	void ExtrapolationOnWall (double *Var, int * Connect, int & normal);
	void ExtrapolationOnCornerConcave (double *Var, int * Connect, int & normal);
	void ExtrapolationOnCornerConvex (double *Var, int * Connect, int & normal);

	void ExtrapolationWallToSolid (double *Var, int * Connect, int & normal);
	void ExtrapolationCornerConcaveToSolid (double *Var, int * Connect, int & normal);
	void ExtrapolationCornerConvexToSolid (double *Var, int * Connect, int & normal);
};
class ExtrapolationWeightDistance: public ExtrapolationDEF {
public:
	ExtrapolationWeightDistance();
	ExtrapolationWeightDistance(int dimension, int nb_vel);
	virtual ~ExtrapolationWeightDistance();
//Extrapolation
	void ExtrapolationOnWall (double *Var, int * Connect, int & normal);
	void ExtrapolationOnCornerConcave (double *Var, int * Connect, int & normal);
	void ExtrapolationOnCornerConvex (double *Var, int * Connect, int & normal);

	void ExtrapolationWallToSolid (double *Var, int * Connect, int & normal);
	void ExtrapolationCornerConcaveToSolid (double *Var, int * Connect, int & normal);
	void ExtrapolationCornerConvexToSolid (double *Var, int * Connect, int & normal);


private:
	double InvSqrt2;
	double InvSqrt2_5;
	double InvSumWeightWallToSolid;
	double InvSumWeightCornerConvexToSolid;
	double InvSumWeightCornerConcaveToSolid;
	double InvSumWeightCornerWallToSolid;
	double InvSumWeightOnWall;
	double InvSumWeightOnCornerConvex;
	double InvSumWeightOnCornerConcave;

};
#endif /* ALGORITHM_LOWORDER_EXTRAPOLATIONSDEF_H_ */
