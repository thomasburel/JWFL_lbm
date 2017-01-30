/*
 * InterpolationDEF.h
 *
 *  Created on: 31 Oct 2016
 *      Author: Thomas Burel
 */

#ifndef ALGORITHM_LOWORDER_INTERPOLATIONSDEF_H_
#define ALGORITHM_LOWORDER_INTERPOLATIONSDEF_H_
#include <algorithm>
#include "../../Core/Parameters.h"
#include "../../Mesh/SingleBlock/NodeArrays.h"
//structures used to find the closest nodes
struct NextNode {
public:
	int index;
	int rankMarkNodes;
	double distance;
	~NextNode(){}
};
struct MatchNextNode{
public:
	MatchNextNode(const int& idx) : idx_(idx) {}
	 bool operator()(const NextNode& obj) const
	 {
	   return obj.index == idx_;
	 }
	 ~MatchNextNode() { }
private:
	 int idx_;
};
struct by_distance {
	bool operator()(NextNode const &a, NextNode const &b) {
		return a.distance < b.distance;
	}
};
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
class InterpolationInverseWeightDistance: public InterpolationDEF {
public:
	InterpolationInverseWeightDistance();
	InterpolationInverseWeightDistance(int dimension, int nb_vel,unsigned int *PtrOppositeInterpol_);
	virtual ~InterpolationInverseWeightDistance();
};
class InterpolationLinearLeastSquare: public InterpolationDEF {
public:
	InterpolationLinearLeastSquare();
	InterpolationLinearLeastSquare(int dimension, int nb_vel,unsigned int *PtrOppositeInterpol_);
	virtual ~InterpolationLinearLeastSquare();
	void InitInterpol(NodeArrays2D *PtrNodes, Parameters *PtrParam);
	void Next_WallId(NodeArrays2D *PtrNodes,std::vector<int> nextback,std::vector<int>& next);
	void Next_FluidId(NodeArrays2D *PtrNodes,std::vector<NextNode> nextbak,std::vector<NextNode>& next);
	void Mark_FluidSolidIds(NodeArrays2D *PtrNodes,int nodeId, int &countnwalls, int &countfluid, int &countsolid);
	void Mark_FluidIds(NodeArrays2D *PtrNodes,std::vector<NextNode>& next, int &countnwalls, int &countfluid);
	void Mark_SolidIds(NodeArrays2D *PtrNodes,std::vector<int>& next, int &countnwalls, int &countsolid);
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
	std::vector<std::vector<int> > solidId,fluidId;//mark nodes
	std::vector<int>  solidIdtmp,fluidIdtmp;//mark nodes
	std::vector<std::vector<double> > solidDist,fluidDist;//keep distance
	std::vector<double>  solidDisttmp,fluidDisttmp;//keep distance
//	double **solidDist,**fluidDist;//keep distance
	int *MapWallId;//Map WallId To LocalId
	int nNodes,nNodes2, nWalls;
	std::vector<NextNode> SolidChecked;
	std::vector<NextNode> FluidChecked;
	std::vector<int> nextwall,nextwallprevious;
	std::vector<NextNode> nextinterior,nextinteriorprevious;
	std::vector<NextNode>::iterator itnextFluid;
	std::vector<NextNode>::iterator itnextWall;
	//std::vector<int>::iterator itnextWall;

};
#endif /* ALGORITHM_LOWORDER_INTERPOLATIONSDEF_H_ */
