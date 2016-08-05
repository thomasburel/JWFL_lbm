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


//#define CallRecoloring(objet,pointer)  ((objet).*(pointer))

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

	void MacroVariables(int& idx);
	void MacroVariablesWithForce(int& idx);
	void MacroVariablesWithNormalDensity(int& idx);
	void MacroVariablesWithNormalDensityAndForce(int& idx);

// Multiphase member functions
	void ColourFluid_Collision();
// Colour gradient definition
	void Set_Colour_gradient();
	void Colour_gradient(int & nodenumber, double* F);
	void Colour_gradient_Gunstensen(int & nodenumber, double* F);
	void Colour_gradient_DensityGrad(int & nodenumber, double* F);
	void Colour_gradient_DensityNormalGrad(int & nodenumber, double* F);

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
	double *Rhor, *Rhob;

// Pointers on function
	//Simplify notation for pointer on member functions
	typedef void(D2Q9ColourFluid::*ColourGrad)(int & nodenumber, double* F);
	typedef void(D2Q9ColourFluid::*Recolour)(double & f, double & fr, double & fb, double & Rho, double & Rho_r, double & Rho_b);

	ColourGrad PtrColourGrad;
	Recolour PtrRecolour;

};

#endif /* SRC_ALGORITHM_LOWORDER_D2Q9COLOURFLUID_H_ */
