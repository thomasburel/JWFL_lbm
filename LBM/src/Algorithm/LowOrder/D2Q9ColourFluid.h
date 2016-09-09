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


#define Callfunction(pointer)  (this->*(pointer))

class D2Q9ColourFluid: public D2Q9TwoPhases {
public:
	D2Q9ColourFluid();
	virtual ~D2Q9ColourFluid();
	D2Q9ColourFluid(MultiBlock* MultiBlock__,ParallelManager* parallel__,WriterManager* Writer__, Parameters* Parameters__,InitLBM& ini);
	virtual void run();

private:
    //! Initialise the colour fluid model.
    /*!
      \param ini : initialisation class (generic initialised methods).
    */
	void InitColourFluid(InitLBM& ini);
    //! Set Pointers On Functions for selecting the right model dynamically.
    /*!
     *
      \sa Set_Collide(), Set_Colour_gradient(), Set_Recolouring(), Set_Macro().
    */
	void Set_PointersOnFunctions();

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

// Recolouring definition
	void Set_Recolouring();
	double CosPhi(int nodenumber, int & direction,double & F_Norm);
	void Recolouring_Latva(int & nodenumber, double * fi_tmp);

//Collide definition
//	void CollideD2Q9ColourFluid(int & direction, double & fi,double &rho,double*  F,double & F_Norm, double & InvTau_, double &u, double &v);
	void Set_Collide();
	void Select_Colour_Operator(ColourOperatorType OperatorType_);
	double TwoPhase_Collision_operator(int & i, double* F);
	void SurfaceForce(int & nodenumber, int* connect,int & normal,double & Fx,double & Fy);
	double Curvature(int & nodenumber, int* connect,int & normal);
	double& Collision_operator_Grunau(int & i, int & nodenumber, double Ak);
	double& Collision_operator_Reis(int & i, int & nodenumber, double Ak);
	void Collision_Grunau(int & nodenumber, int* connect,int & normal,double* fi);
	void Collision_Reis(int & nodenumber, int* connect,int & normal,double* fi);
	void Collision_SurfaceForce(int & nodenumber, int* connect,int & normal,double* fi);

	double TwoPhase_Collision_operator(int & nodenumber, int & direction, double & Ak, double* F, double & F_Norm);
	double Convert_Alpha_To_Rho(double alpha);
	double Convert_Rho_To_Alpha(double Rho);
	double Cal_RhoR_Corner(NodeCorner2D& Node);
	double Cal_RhoB_Corner(NodeCorner2D& Node);

// Boundary conditions depend of the model. Some functions has to be rewritten.
	void ApplyBc();
/*	void ApplyGlobalCorner(NodeCorner2D& NodeIn);
	void ApplyBounceBack(NodeWall2D& Node);
	void ApplyCorner(NodeCorner2D& Node);
	void ApplyBounceBack(NodeCorner2D& Node);
	void ApplySymmetryPressureOnNode(NodeSymmetry2D& NodeIn);*/


private:
//Multiphase variables
	double *RhoN;// Normal density
	double *Rhor;// Density red fluid
	double *Rhob;// DensityBlue fluid
	double **F, **G;///< Surface Force and Colour gradient/density gradient
	double *G_Norm;///< Norm of the colour gradient
	double tension;///< Surface tension
	double beta;///< Separation coefficient for the recolouring method
	double A1;///< Guntensen parameter for the red fluid
	double A2;///< Guntensen parameter for the blue fluid
	double Bi[9];///< Reis correction
	double Rho_limiter;///< Approximation to null density
	double D_tmp;// Temporary double
	double* PtrD_tmp;// Temporary pointer for a double
	double DVec_2D_tmp[2];// Temporary vector 2D for a double
	double DArray_2D_tmp[2][2];// Temporary vector 2D for a double
// Pointers on function

///Simplify notation for pointer on a member function of D2Q9ColourFluid class for Colour Gradient methods
	typedef void(D2Q9ColourFluid::*ColourGrad)(int & nodenumber, int* connect,int & normal);
///Simplify notation for pointer on a member function of D2Q9ColourFluid class for recolouring methods
	typedef void(D2Q9ColourFluid::*Recolour)(int & nodenumber, double * fi_tmp);
///Simplify notation for pointer on a member function of D2Q9ColourFluid class for macroscopic variables calculation methods
	typedef void(D2Q9ColourFluid::*Macro)(int & nodenumber);
///Simplify notation for pointer on a member function of D2Q9ColourFluid class for collision models
	typedef void(D2Q9ColourFluid::*Collision)(int & nodenumber, int* connect,int & normal,double* fi);
//Define name for pointers on functions
	ColourGrad PtrColourGrad;///< Colour gradient pointer for interior nodes
	ColourGrad PtrColourGradBc;///< Colour gradient pointer for boundary condition nodes
	ColourGrad PtrColourGradCorner;///< Colour gradient pointer for corner nodes
	Recolour PtrRecolour;///< Recolouring pointer
	Macro PtrMacro;///< Macroscopic pointer
	Collision PtrCollision;///< Collision pointer

};

#endif /* SRC_ALGORITHM_LOWORDER_D2Q9COLOURFLUID_H_ */
