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

	void Set_Corner(Dictionary *PtrDic,NodeArrays2D* NodeArrays, Parameters *Param, double ** &Ei,unsigned int nbDistributions=1);
	void ApplyCorner(NodeCorner2D& Node, DistriFunct* f_in, unsigned int idxDistribution,double const & RhoDef,double const & UDef,double const & VDef, double *Rho, double *U, double *V);
	void ApplyCornerWall(NodeCorner2D& Node, DistriFunct* f_in, unsigned int idxDistribution, double *Rho, double *U, double *V);
	void ApplyCornerSpecialWall(NodeWall2D& Node, DistriFunct* f_in, unsigned int idxDistribution, double *Rho, double *U, double *V);
	void ApplyPreVelSpecialWall(NodeWall2D& Node, DistriFunct* f_in, unsigned int idxDistribution,double const & RhoDef,double const & UDef,double const & VDef);

	void ApplyCornerPreStream(NodeCorner2D& Node, DistriFunct* f_in, unsigned int idxDistribution,double const & RhoDef,double const & UDef,double const & VDef, double *Rho, double *U, double *V);
	void ApplyCornerWallPreStream(NodeCorner2D& Node, DistriFunct* f_in, unsigned int idxDistribution, double *Rho, double *U, double *V);
	void ApplyCornerSpecialWallPreStream(NodeWall2D& Node, DistriFunct* f_in, unsigned int idxDistribution, double *Rho, double *U, double *V);
	void ApplyPreVelSpecialWallPreStream(NodeWall2D& Node, DistriFunct* f_in, unsigned int idxDistribution,double const & RhoDef,double const & UDef,double const & VDef);
private:
// initialise specific method
	void SetHalfWayBounceBack(NodeArrays2D* NodeArrays,unsigned int nbDistributions=1);
//Wall methods
	template <class T>
	void ApplyBounceBack(T& Node,int const &BcNormal,int const *Connect, bool concave, DistriFunct* f_in, double Rho, double U, double V, unsigned int idxDistribution);
	template <class T>
	void ApplyHalfWayBounceBack(T& Node,int const &BcNormal,int const *Connect, bool concave, DistriFunct* f_in, double Rho, double U, double V, unsigned int idxDistribution);
	template <class T>
	void ApplyDiffuseWall(T& Node,int const &BcNormal,int const *Connect, bool concave, DistriFunct* f_in, double Rho, double U, double V, unsigned int idxDistribution);

	/// Corner treat by Chih-Fung Ho, Cheng Chang, Kuen-Hau Lin and Chao-An Lin
	/// Consistent Boundary Conditions for 2D and 3D Lattice Boltzmann Simulations
	template <class T>
	void ApplyHoChan(T& Node,int const &BcNormal,int const *Connect, bool concave, DistriFunct* f_in, double Rho, double U, double V, unsigned int idxDistribution);
	template <class T>
	void ApplyHoChanNoVel(T& Node,int const &BcNormal,int const *Connect, bool concave, DistriFunct* f_in, double Rho, double U, double V, unsigned int idxDistribution);

//Before streaming
	template <class T>
	void ApplyHalfWayBounceBackPreStream(T& Node,int const &BcNormal,int const *Connect, bool concave, DistriFunct* f_in, double Rho, double U, double V, unsigned int idxDistribution);

	void FUNC_corner (double & a,double & b,double & c,double & d,double & e,double & f,double & g,double & h,double & i,double & U,double & V,double & Rho);
	void FUNC_corner_no_vel (double & a,double & b,double & c,double & d,double & e,double & f,double & g,double & h,double & i,double & Rho);

//Function for calculating the density and setting in the global variable
	//Get the density from the global variable
	double GetRho(NodeCorner2D& Node, double *Rho,DistriFunct* f_in);
	void FixRho(NodeCorner2D& Node, double *Rho,DistriFunct* f_in);
	//Get the density by using the two direct neighbours
	void ExtrapolationAvgRho(NodeCorner2D& Node, double *Rho,DistriFunct* f_in);
	//Get the density by using the known distribution
	void LocalRho(NodeCorner2D& Node, double *Rho,DistriFunct* f_in);

// Pointers on function
///Simplify notation for pointer on a member function of D2Q9Pressure class for Pressure model used
	//Corner inside the domain
	typedef void(D2Q9Corner::*CornerWallMethod)(NodeCorner2D& Node,int const &BcNormal,int const *Connect, bool concave, DistriFunct* f_in, double Rho, double U, double V, unsigned int idxDistribution);
	typedef void(D2Q9Corner::*CornerWallSpecialWallMethod)(NodeWall2D& Node,int const &BcNormal,int const *Connect, bool concave, DistriFunct* f_in, double Rho, double U, double V, unsigned int idxDistribution);
	typedef void(D2Q9Corner::*CalculRhoCornerWall)(NodeCorner2D& Node, double *Rho,DistriFunct* f_in);
	//Corner in the corner of the domain
	typedef void(D2Q9Corner::*CornerMethod)(NodeCorner2D& Node,int const &BcNormal,int const *Connect, bool concave, DistriFunct* f_in, double Rho, double U, double V, unsigned int idxDistribution);
	typedef void(D2Q9Corner::*CornerSpecialWallMethod)(NodeWall2D& Node,int const &BcNormal,int const *Connect, bool concave, DistriFunct* f_in, double Rho, double U, double V, unsigned int idxDistribution);

//Define name for pointers on functions
	CornerMethod PtrCornerMethod;
	CornerWallMethod PtrCornerWallMethod;
	CornerSpecialWallMethod PtrCornerSpecialWallMethod;
	CornerMethod PtrPreStreamCornerMethod;
	CornerWallMethod PtrPreStreamCornerWallMethod;
	CornerSpecialWallMethod PtrPreStreamCornerSpecialWallMethod;
	CornerWallSpecialWallMethod PtrCornerWallSpecialWallMethod;
	CalculRhoCornerWall PtrCalculRhoCornerWall;
	Extrapolation Extrapol;
	double InvRho,InvU,InvV;
	double doubleTmpReturn;
	short int direction1,direction2;

};

#endif /* SRC_ALGORITHM_LOWORDER_BOUNDARIES_D2Q9CORNER_H_ */
