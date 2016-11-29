/*
 * D2Q9ColourFluid.h
 *
 *  Created on: 31 Jul 2016
 *      Author: Thomas Burel
 *
 *  Solver for the colour fluid model
 */

#ifndef SRC_ALGORITHM_LOWORDER_D2Q9COLOURFLUID_H_
#define SRC_ALGORITHM_LOWORDER_D2Q9COLOURFLUID_H_

#include "D2Q9_TwoPhases.h"
#include <cmath>

#define Callfunction(pointer)  (this->*(pointer))

class D2Q9ColourFluid: public D2Q9TwoPhases {
public:
	D2Q9ColourFluid();
	virtual ~D2Q9ColourFluid();
	D2Q9ColourFluid(MultiBlock* MultiBlock__,ParallelManager* parallel__,WriterManager* Writer__, Parameters* Parameters__,InitLBM& ini);
	virtual void run();
	virtual void run(Parameters* UpdatedParam);
	virtual void UpdateAllDomain(Parameters* UpdatedParam,InitLBM& ini);
	virtual void UpdateDomainBc(Parameters* UpdatedParam,InitLBM& ini);
	virtual void UpdateWall(Parameters* UpdatedParam,InitLBM& ini);
	virtual void UpdateInterior(Parameters* UpdatedParam,InitLBM& ini);
private:
	//! Initialise the colour fluid model.
    /*!
      \param ini : initialisation class (generic initialised methods).
    */
	void InitColourFluid(InitLBM& ini);
	void InitColourFluidAllDomain(InitLBM& ini);
	void InitColourFluidDomainBc(InitLBM& ini);
	void InitColourFluidWall(InitLBM& ini);
	void InitColourFluidInterior(InitLBM& ini);
    //! Set Pointers On Functions for selecting the right model dynamically.
    /*!
     *
      \sa Set_Collide(), Set_Colour_gradient(), Set_Recolouring(), Set_Macro().
    */
	void Set_PointersOnFunctions();

	double Calcul_Error();
	double Error_RhoN(int &idx);

	void UpdateMacroVariables();

// Macroscopic calculation
	void Set_Macro();
	void MacroVariables(int& idx);
	void MacroVariablesWithForce(int& idx);
	void MacroVariablesWithNormalDensity(int& idx);
	void MacroVariablesWithNormalDensityAndForce(int& idx);
// Calculate Force for macro

// Multiphase member functions
	void ColourFluid_Collision();

// Colour gradient definition
	void Set_Colour_gradient();
	void Colour_gradient();


	void Colour_gradient_Gunstensen(int & nodenumber, int* connect,int & normal);
	void Colour_gradient_DensityGrad(int & nodenumber, int* connect,int & normal);
	void Colour_gradient_DensityGradBc(int & nodenumber, int* connect,int & normal);
	void Colour_gradient_DensityGradCorner(int & nodenumber, int* connect,int & normal);
	void Colour_gradient_DensityNormalGrad(int & nodenumber, int* connect,int & normal);
	void Colour_gradient_DensityNormalGradBc(int & nodenumber, int* connect,int & normal);
	void Colour_gradient_DensityNormalGradCorner(int & nodenumber, int* connect,int & normal);

	void Set_NormalWall();
	void CalculNormal(int & nodenumber, int* connect,int & normal); //Normalise the colour gradient
	void Not_CalculNormal(int & nodenumber, int* connect,int & normal){}; //Do nothing
	void CalculNormal_NoTeta(int & nodenumber, int* connect,int & normal);
	void CalculNormal_FixTeta(int & nodenumber, int* connect,int & normal);
	void CalculNormal_NoCstTeta(int & nodenumber, int* connect,int & normal);

	void Extrapolate_NormalInSolid();
	void Impose_ContactAngleInSolid();
	void Select_ContactAngle(int & normal,double & Nx, double & Ny);
	void ContactAngleConcaveCornerInSolid(int* connect, int & normal);
	void ContactAngleConvexCornerInSolid(int* connect, int & normal);
	void ContactAngleWallInSolid(int* connect, int & normal);
	void ContactAngleConcaveCornerOnBc(int* connect, int & normal);
	void ContactAngleConvexCornerOnBc(int* connect, int & normal);
	void ContactAngleWallOnBc(int* connect, int & normal);
	void Synchronise_Colour_gradient();

// Recolouring definition
	void Set_Recolouring();
	double CosPhi(int nodenumber, int & direction,double & F_Norm);
	void Recolouring_Latva(int & nodenumber, double * fi_tmp);
	void Recolouring_Wall(int & nodenumber, double * fi_tmp);
//Collide definition
//	void CollideD2Q9ColourFluid(int & direction, double & fi,double &rho,double*  F,double & F_Norm, double & InvTau_, double &u, double &v);
	void Set_Collide();
	void Select_Colour_Operator(ColourOperatorType OperatorType_);
	double TwoPhase_Collision_operator(int & i, double* F);
	void SurfaceForce(int & nodenumber, int* connect,int & normal,double & Fx,double & Fy);
	void SurfaceForceWall(int & nodenumber, int* connect,int & normal,double & Fx,double & Fy);
	double& Collision_operator_Grunau(int & i, int & nodenumber, double Ak);
	double& Collision_operator_Reis(int & i, int & nodenumber, double Ak);
	void Collision_Grunau(int & nodenumber, int* connect,int & normal,double* fi);
	void Collision_Reis(int & nodenumber, int* connect,int & normal,double* fi);
	void Collision_SurfaceForce(int & nodenumber, int* connect,int & normal,double* fi);
	void Collision_SurfaceForceWall(int & nodenumber, int* connect,int & normal,double* fi);

	double Curvature(int & nodenumber, int* connect,int & normal);
	double CurvatureWall(int & nodenumber, int* connect,int & normal);

	double TwoPhase_Collision_operator(int & nodenumber, int & direction, double & Ak, double* F, double & F_Norm);
	double Convert_Alpha_To_Rho(double alpha);
	double Convert_Rho_To_Alpha(double Rho);
	double Cal_RhoR_Corner(NodeCorner2D& Node);
	double Cal_RhoB_Corner(NodeCorner2D& Node);
	void Extrapol_Density_Corner();

	double ExtrapolationSpacial2ndOrder(double* Var,int index1,int index2){return 2.0*Var[index1]-Var[index2];};
	void NormalDensityExtrapolationSpacial2ndOrder(int const & idxNodeArray, int & nodenumber, int* connect,int & normal);
	void NormalDensityExtrapolationWeight(int const & idxNodeArray, int & nodenumber, int* connect,int & normal);
	void NormalDensityNoExtrapolation(int const & idxNodeArray, int & nodenumber, int* connect,int & normal){};
	void NormalDensityNoTeta(int const & idxNodeArray, int & nodenumber, int* connect,int & normal);
// Boundary conditions depend of the model. Some functions has to be rewritten.
	void ApplyBc();
/*	void ApplyGlobalCorner(NodeCorner2D& NodeIn);
	void ApplyBounceBack(NodeWall2D& Node);
	void ApplyCorner(NodeCorner2D& Node);
	void ApplyBounceBack(NodeCorner2D& Node);
	void ApplySymmetryPressureOnNode(NodeSymmetry2D& NodeIn);*/


private:
//Multiphase variables
	double *testVar;
	double *RhoN;// Normal density
	double *Rhor;// Density red fluid
	double *Rhob;// DensityBlue fluid
	double **F, **G;///< Surface Force and Colour gradient/density gradient
	double **Normal;///< normal of the interface
	double *G_Norm;///< Norm of the colour gradient
	double *Curv;//Store curvature
	double tension;///< Surface tension
	double *teta;///< Contact angle
	double beta;///< Separation coefficient for the recolouring method
	double A1;///< Guntensen parameter for the red fluid
	double A2;///< Guntensen parameter for the blue fluid
	double Bi[9];///< Reis correction
	double Rho_limiter;///< Approximation to null density
	int I_tmp;// Temporary integer
	int& IntRef(int I_input){I_tmp=I_input;return I_tmp;};
	double D_tmp;// Temporary double
	double& doubleRef(int D_input){D_tmp=D_input;return D_tmp;};
	double* PtrD_tmp;// Temporary pointer for a double
	double DVec_2D_tmp[2];// Temporary vector 2D for a double
	double DArray_2D_tmp[2][2];// Temporary vector 2D for a double

	double n1[9][2],n2[9][2];
	double D1,D2,r,rMinus1;
	double costeta,sinteta;

// Pointers on function

///Simplify notation for pointer on a member function of D2Q9ColourFluid class for Colour Gradient methods
	typedef void(D2Q9ColourFluid::*ColourGrad)(int & nodenumber, int* connect,int & normal);
///Simplify notation for pointer on a member function of D2Q9ColourFluid class for recolouring methods
	typedef void(D2Q9ColourFluid::*Recolour)(int & nodenumber, double * fi_tmp);
///Simplify notation for pointer on a member function of D2Q9ColourFluid class for macroscopic variables calculation methods
	typedef void(D2Q9ColourFluid::*Macro)(int & nodenumber);
///Simplify notation for pointer on a member function of D2Q9ColourFluid class for extrapolation density in solid by considering a mixture fluid
	typedef void(D2Q9ColourFluid::*ExtrapolDensity)(int const & idxNodeArray, int & nodenumber, int* connect,int & normal);
///Simplify notation for pointer on a member function of D2Q9ColourFluid class for collision models
	typedef void(D2Q9ColourFluid::*Collision)(int & nodenumber, int* connect,int & normal,double* fi);
///Simplify notation for pointer on a member function of D2Q9ColourFluid class for the normalisation of the Colour Gradient
	typedef void(D2Q9ColourFluid::*CalNormal)(int & nodenumber, int* connect,int & normal);
//Define name for pointers on functions
	ColourGrad PtrColourGrad;///< Colour gradient pointer for interior nodes
	ColourGrad PtrColourGradBc;///< Colour gradient pointer for boundary condition nodes
	ColourGrad PtrColourGradWall;///< Colour gradient pointer for wall nodes (force the contact angle)
	ColourGrad PtrColourGradCorner;///< Colour gradient pointer for corner nodes
	Recolour PtrRecolour;///< Recolouring pointer
	Macro PtrMacro;///< Macroscopic pointer
	ExtrapolDensity PtrExtrapolDensity;///< Pointer on Density extrapolation in Solid
	Collision PtrCollision;///< Collision pointer
	Collision PtrCollisionWall;///< Collision pointer
	CalNormal PtrCalNormal;///< Calcul Normal
	CalNormal PtrCalNormalWall;///< Calcul or fix normal at the wall
};

#endif /* SRC_ALGORITHM_LOWORDER_D2Q9COLOURFLUID_H_ */
