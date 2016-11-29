/*
 * ============================================================================
 * D2Q9Corner.h
 *
 *  Created on: 29 Aug 2016
 *      Author: Thomas Burel
 *      Version     :
 *      Copyright   : Your copyright notice
 *      Description : LBM Code Developing by JWFL at University of Strathclyde
 * ============================================================================
 */

#ifndef SRC_ALGORITHM_LOWORDER_BOUNDARIES_D2Q9CORNER_H_
#define SRC_ALGORITHM_LOWORDER_BOUNDARIES_D2Q9CORNER_H_

#include "D2Q9BcVar.h"


class D2Q9Corner: public D2Q9BcVar {
public:
	D2Q9Corner();
	virtual ~D2Q9Corner();

	void Set_Corner(Parameters *Param);
	void ApplyCorner(NodeCorner2D& Node, DistriFunct* f_in,double const & RhoDef,double const & UDef,double const & VDef, double *Rho, double *U, double *V);
	void ApplyCornerWall(NodeCorner2D& Node, DistriFunct* f_in, double *Rho, double *U, double *V);
	void ApplyCornerSpecialWall(NodeWall2D& Node, DistriFunct* f_in, double *Rho, double *U, double *V);
	void ApplyPreVelSpecialWall(NodeWall2D& Node, DistriFunct* f_in,double const & RhoDef,double const & UDef,double const & VDef);

private:

	void ApplyBounceBack(int const &BcNormal,int const *Connect, bool concave, DistriFunct* f_in, double Rho, double U, double V);
	void ApplyDiffuseWall(int const &BcNormal,int const *Connect, bool concave, DistriFunct* f_in, double Rho, double U, double V);

	/// Corner treat by Chih-Fung Ho, Cheng Chang, Kuen-Hau Lin and Chao-An Lin
	/// Consistent Boundary Conditions for 2D and 3D Lattice Boltzmann Simulations
	void ApplyHoChan(int const &BcNormal,int const *Connect, bool concave, DistriFunct* f_in, double Rho, double U, double V);
	void ApplyHoChanNoVel(int const &BcNormal,int const *Connect, bool concave, DistriFunct* f_in, double Rho, double U, double V);


	void FUNC_corner (double & a,double & b,double & c,double & d,double & e,double & f,double & g,double & h,double & i,double & U,double & V,double & Rho);
	void FUNC_corner_no_vel (double & a,double & b,double & c,double & d,double & e,double & f,double & g,double & h,double & i,double & Rho);

//Function for calculating the density and setting in the global variable
	//Get the density from the global variable
	double GetRho(NodeCorner2D& Node, double *Rho);
	void FixRho(NodeCorner2D& Node, double *Rho);
	//Get the density by using the two direct neighbours
	void ExtrapolationAvgRho(NodeCorner2D& Node, double *Rho);

// Pointers on function
///Simplify notation for pointer on a member function of D2Q9Pressure class for Pressure model used
	//Corner inside the domain
	typedef void(D2Q9Corner::*CornerWallMethod)(int const &BcNormal,int const *Connect, bool concave, DistriFunct* f_in, double Rho, double U, double V);
	typedef void(D2Q9Corner::*CalculRhoCornerWall)(NodeCorner2D& Node, double *Rho);
	//Corner in the corner of the domain
	typedef void(D2Q9Corner::*CornerMethod)(int const &BcNormal,int const *Connect, bool concave, DistriFunct* f_in, double Rho, double U, double V);


//Define name for pointers on functions
	CornerMethod PtrCornerMethod;
	CornerWallMethod PtrCornerWallMethod;
	CalculRhoCornerWall PtrCalculRhoCornerWall;
	Extrapolation Extrapol;
	double InvRho,InvU,InvV;
	double doubleTmpReturn;
	short int direction1,direction2;

};

#endif /* SRC_ALGORITHM_LOWORDER_BOUNDARIES_D2Q9CORNER_H_ */
