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
	void InitMultiphase(InitLBM& ini);
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
	double CosPhi(int & direction, double* F,double & F_Norm);
// Colour gradient definition
	void Set_Colour_gradient();
	void Colour_gradient();
	void Colour_gradient_Gunstensen(int & nodenumber, int* connect,int & normal, double* F);
	void Colour_gradient_DensityGrad(int & nodenumber, int* connect,int & normal, double* F);
	void Colour_gradient_DensityGradBc(int & nodenumber, int* connect,int & normal, double* F);
	void Colour_gradient_DensityGradCorner(int & nodenumber, int* connect,int & normal, double* F);
	void Colour_gradient_DensityNormalGrad(int & nodenumber, int* connect,int & normal, double* F);
	void Colour_gradient_DensityNormalGradBc(int & nodenumber, int* connect,int & normal, double* F);
	void Colour_gradient_DensityNormalGradCorner(int & nodenumber, int* connect,int & normal, double* F);

// Recolouring definition
	void Set_Recolouring();
	void Recolouring_Latva(int & nodenumber,int & i, double & ftmp, double* F,double & F_Norm);

//Collide definition
	void CollideD2Q9ColourFluid(int & direction, double & fi,double &rho,double*  F,double & F_Norm, double & InvTau_, double &u, double &v);
	void Set_Collide();
	double TwoPhase_Collision_operator(int & i, double* F, double & F_Norm);
	double* SurfaceForce(int & nodenumber, int* connect,int & normal, double & F_Norm);
	double Curvature(int & nodenumber, int* connect,int & normal);
	double* TwoPhase_Collision_BodyForce(int & nodenumber, int & i, double & InvTau_tmp, double* F_tmp, double & F_Norm);

	double TwoPhase_Collision_operator(int & nodenumber, int & direction, double & Ak, double* F, double & F_Norm);
	double Convert_Alpha_To_Rho(double alpha);
	double Convert_Rho_To_Alpha(double Rho);
	double Cal_RhoR_Corner(NodeCorner2D& Node);
	double Cal_RhoB_Corner(NodeCorner2D& Node);

// Boundary conditions depend of the model. Some functions has to be rewritten.
	void ApplyBc();
	void ApplyGlobalCorner(NodeCorner2D& NodeIn);
	void ApplyBounceBack(NodeWall2D& Node);
	void ApplyCorner(NodeCorner2D& Node);
	void ApplyBounceBack(NodeCorner2D& Node);
	void ApplySymmetryPressureOnNode(NodeSymmetry2D& NodeIn);


private:
//Multiphase variables
	double **F, **G;//Surface Force and Colour gradient/density gradient
	double tension;
	double beta,A1,A2;
	double Rho_limiter;
	double D_tmp;// Temporary double
	double* PtrD_tmp;// Temporary pointer for a double
	double DVec_2D_tmp[2];// Temporary vector 2D for a double
// Pointers on function
	//Simplify notation for pointer on member functions
	typedef void(D2Q9ColourFluid::*ColourGrad)(int & nodenumber, int* connect,int & normal, double* F);
	typedef void(D2Q9ColourFluid::*Recolour)(int & nodenumber,int & i, double & ftmp, double* F,double & F_Norm);
	typedef void(D2Q9ColourFluid::*Macro)(int & nodenumber);
	typedef void(D2Q9ColourFluid::*TwoPhaseOperator)(int & nodenumber, int & direction, double & Ak, double* F, double & F_Norm);
	//Define name for pointers on functions
	ColourGrad PtrColourGrad;
	ColourGrad PtrColourGradBc;
	ColourGrad PtrColourGradCorner;
	Recolour PtrRecolour;
	Macro PtrMacro;
	TwoPhaseOperator PtrTwoPhaseOperator;

};

#endif /* SRC_ALGORITHM_LOWORDER_D2Q9COLOURFLUID_H_ */
