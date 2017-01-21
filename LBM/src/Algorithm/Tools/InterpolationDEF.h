/*
 * InterpolationDEF.h
 *
 *  Created on: 31 Oct 2016
 *      Author: Thomas Burel
 */

#ifndef ALGORITHM_LOWORDER_INTERPOLATIONSDEF_H_
#define ALGORITHM_LOWORDER_INTERPOLATIONSDEF_H_

#include "../../Core/Parameters.h"
#include "../../Mesh/SingleBlock/NodeArrays.h"
///Abstract class for common function for all kind of extrapolation
class InterpolationDEF {
public:
	InterpolationDEF();
	InterpolationDEF(int dimension, int nb_vel,unsigned int *PtrOppositeInterpol_);
	virtual ~InterpolationDEF();
	virtual void InitInterpol(NodeArrays2D *PtrNodes, Parameters *PtrParam)=0;
	//Interpolation
	virtual void InterpolationOnWall (double *Var, int * Connect, int & normal)=0;
	virtual void InterpolationOnCornerConcave (double *Var, int * Connect, int & normal)=0;
	virtual void InterpolationOnCornerConvex (double *Var, int * Connect, int & normal)=0;
	virtual void InterpolationOnWall (double *Var1, double *Var2, int * Connect, int & normal)=0;
	virtual void InterpolationOnCornerConcave (double *Var1, double *Var2, int * Connect, int & normal)=0;
	virtual void InterpolationOnCornerConvex (double *Var1, double *Var2, int * Connect, int & normal)=0;
protected:
	int dimension,nb_Vel;
	unsigned int *PtrOppositeInterpol;

};
class NoInterpolation: public InterpolationDEF {
public:
	NoInterpolation();
	NoInterpolation(int dimension, int nb_vel,unsigned int *PtrOppositeInterpol_);
	virtual ~NoInterpolation();
	void InitInterpol(NodeArrays2D *PtrNodes, Parameters *PtrParam){};
//Interpolation
	void InterpolationOnWall (double *Var, int * Connect, int & normal){};
	void InterpolationOnCornerConcave (double *Var, int * Connect, int & normal){};
	void InterpolationOnCornerConvex (double *Var, int * Connect, int & normal){};
	void InterpolationOnWall (double *Var1, double *Var2, int * Connect, int & normal){};
	void InterpolationOnCornerConcave (double *Var1, double *Var2, int * Connect, int & normal){};
	void InterpolationOnCornerConvex (double *Var1, double *Var2, int * Connect, int & normal){};
};
class InterpolationLinear: public InterpolationDEF {
public:
	InterpolationLinear();
	InterpolationLinear(int dimension, int nb_vel,unsigned int *PtrOppositeInterpol_);
	virtual ~InterpolationLinear();
	void InitInterpol(NodeArrays2D *PtrNodes, Parameters *PtrParam){};
//Interpolation
	void InterpolationOnWall (double *Var, int * Connect, int & normal);
	void InterpolationOnCornerConcave (double *Var, int * Connect, int & normal);
	void InterpolationOnCornerConvex (double *Var, int * Connect, int & normal);
	void InterpolationOnWall (double *Var1, double *Var2, int * Connect, int & normal);
	void InterpolationOnCornerConcave (double *Var1, double *Var2, int * Connect, int & normal);
	void InterpolationOnCornerConvex (double *Var1, double *Var2, int * Connect, int & normal);
};
class InterpolationLinearLeastSquare: public InterpolationDEF {
public:
	InterpolationLinearLeastSquare();
	InterpolationLinearLeastSquare(int dimension, int nb_vel,unsigned int *PtrOppositeInterpol_);
	virtual ~InterpolationLinearLeastSquare();
	void InitInterpol(NodeArrays2D *PtrNodes, Parameters *PtrParam);
	void Next_WallId(NodeArrays2D *PtrNodes,int nodeId,std::vector<int>& next);
	void Mark_FluidSolidIds(NodeArrays2D *PtrNodes,int nodeId, int &countnwalls, int &countfluid, int &countsolid);
	double DistToWall(double x, double y){return std::sqrt((x-x0)*(x-x0)+(y-y0)*(y-y0));};
	void CalculLeastSquareMethod(double *Var,int wallId);
	void CalculLeastSquareMethod(double *Var1,double *Var2,int wallId);
//Interpolation
	void InterpolationOnWall (double *Var, int * Connect, int & normal);
	void InterpolationOnCornerConcave (double *Var, int * Connect, int & normal);
	void InterpolationOnCornerConvex (double *Var, int * Connect, int & normal);
	void InterpolationOnWall (double *Var1, double *Var2, int * Connect, int & normal);
	void InterpolationOnCornerConcave (double *Var1, double *Var2, int * Connect, int & normal);
	void InterpolationOnCornerConvex (double *Var1, double *Var2, int * Connect, int & normal);

private:
	double sumX,sumY,sumXY,sumX2,result1,result2;
	double x0,y0;//save position of the wall node
	short int **solidId,**fluidId;//mark nodes
	double **solidDist,**fluidDist;//keep distance
	short int *MapWallId;//Map WallId To LocalId
	short int nNodes,nNodes2, nWalls;


};
#endif /* ALGORITHM_LOWORDER_INTERPOLATIONSDEF_H_ */
