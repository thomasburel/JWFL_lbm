/*
 * ============================================================================
 * D2Q9Wall.h
 *
 *  Created on: 29 Aug 2016
 *      Author: Thomas Burel
 *      Version     :
 *      Copyright   : Your copyright notice
 *      Description : LBM Code Developing by JWFL at University of Strathclyde
 * ============================================================================
 */

#ifndef SRC_ALGORITHM_LOWORDER_BOUNDARIES_D2Q9WALL_H_
#define SRC_ALGORITHM_LOWORDER_BOUNDARIES_D2Q9WALL_H_

#include "D2Q9BcVar.h"

class D2Q9Wall: public D2Q9BcVar {
public:
	D2Q9Wall();
	virtual ~D2Q9Wall();

	void SetWall(Parameters *Param);
	// Apply wall treatment depending of the parameters
	void ApplyWall(int const &BcNormal,int const *Connect, DistriFunct* f_in, double const *Rho, double const *U, double const *V);
	//Special Walls are the wall on the boundary of the domain. They need two treatments: The boundary treatment and the wall treatment
	void ApplySpecialWall(NodeWall2D & Node, DistriFunct* f_in, std::map<int,NodeType> TypeOfNode_, double const *Rho, double const *U, double const *V);

private:
//Wall methods
	void ApplyDiffuseWall(int const &BcNormal,int const *Connect, DistriFunct* f_in);
	void ApplyBounceBackWall(int const &BcNormal,int const *Connect, DistriFunct* f_in);
	void ApplyHeZouWall(int const &BcNormal,int const *Connect, DistriFunct* f_in);
	void FUNC_HeZou_NoU (double & a,double & b,double & c,double & d,double & e,double & f,double & g,double & h,double & i);

//Special wall methods
	void ApplyBounceBackSymmetry(NodeWall2D & Node, DistriFunct* f_in, std::map<int,NodeType> & TypeOfNode_);
	void ApplyDiffuseWallSymmetry(NodeWall2D & Node, DistriFunct* f_in, std::map<int,NodeType> & TypeOfNode_);

// Pointers on function
///Simplify notation for pointer on a member function of D2Q9Pressure class for Pressure model used
	typedef void(D2Q9Wall::*WallMethod)(int const  &BcNormal,int const *Connect, DistriFunct* f_in);
	typedef void(D2Q9Wall::*SpecialWallMethod)(NodeWall2D & Node, DistriFunct* f_in, std::map<int,NodeType> & TypeOfNode_);
//Define name for pointers on functions
	WallMethod PtrWallMethod;
	SpecialWallMethod PtrSpecialWallMethod;

	double feq,Rho;
};

#endif /* SRC_ALGORITHM_LOWORDER_BOUNDARIES_D2Q9WALL_H_ */
