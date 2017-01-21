/*
 * ============================================================================
 * ContactAngle.h
 *
 *  Created on: 10 Jan 2017
 *      Author: Thomas Burel
 *      Version     :
 *      Copyright   : Your copyright notice
 *      Description : LBM Code Developing by JWFL at University of Strathclyde
 * ============================================================================
 */

#ifndef SRC_ALGORITHM_TOOLS_CONTACTANGLE_H_
#define SRC_ALGORITHM_TOOLS_CONTACTANGLE_H_

#include "../../Core/Parameters.h"
#include "../../Mesh/SingleBlock/NodeArrays.h"
#include "../../Core/GlobalDef.h"

#include "Extrapolation.h"
#include "Interpolation.h"

class ContactAngle {
public:
	ContactAngle();
	virtual ~ContactAngle();
	void InitContactAngle(NodeArrays2D *PtrNode,Parameters *PtrParam,unsigned int *PtrOpposite, double * &tetaIn);
	void AllocateTeta(NodeArrays2D *PtrNode,Parameters *PtrParam,double * &tetaIn);
	void ApplyContactAngle2D(double **&Normal);
	double*& Get_PtrTeta(){return teta;};

private:
	void Select_ContactAngle2D(int &nodeWallIdx, int & normal,double & Nx, double & Ny);
	void ContactAngleonWalls2D(int &nodeWallIdx, double **&Normal,int* Connect, int & normal);
	void ContactAngleConcaveCornerInSolid2D(int &nodeWallIdx, double **&Normal,int* connect, int & normal);
	void ContactAngleConvexCornerInSolid2D(int &nodeWallIdx, double **&Normal,int* connect, int & normal);
	void ContactAngleWallInSolid2D(int &nodeWallIdx, double **&Normal,int* Connect, int & normal);

	void Set_TwoChoiceOfContactAngle2D(int &nodeWallIdx);
	void LinearSwitchContactAngle2D(int &nodeWallIdx, double & r,int & normal,double & Nx, double & Ny);
	void BinarySwitchContactAngle2D(int &nodeWallIdx, double & r,int & normal,double & Nx, double & Ny);
	typedef void(ContactAngle::*SwitchContactAngle)(int &nodeWallIdx, double & r,int & normal,double & Nx, double & Ny);
	SwitchContactAngle Switch2D;

	void ApplyNoTetaMethod(double **&Normal);
	void ApplyStandardMethod(double **&Normal);
	void ApplyInterpolMethod(double **&Normal);
	typedef void(ContactAngle::*SwitchModel)(double **&Normal);
	SwitchModel Model;
	void NoExtrapolate_NormalInSolid(double **&Normal){};
	void Extrapolate_NormalInSolid2D(double **&Normal);
	typedef void(ContactAngle::*ExtrapolNormal)(double **&Normal);
	ExtrapolNormal ExtrapolNormalInSolid;
	void Impose_ContactAngleOnWall2DFixTeta(double **&Normal);
	void Impose_ContactAngleOnWall2DNonCstTeta(double **&Normal);
	void Impose_ContactAngleInSolidAndInterpol2DFixTeta(double **&Normal);
	void Impose_ContactAngleInSolidAndInterpol2DNonCstTeta(double **&Normal);
	typedef void(ContactAngle::*ImposeNormal)(double **&Normal);
	ImposeNormal ImposeNormalOnWall;
private:
	NodeArrays2D *PtrNodeCa;
	Parameters *PtrParamCa;
	Extrapolation ExtraPolNormalCa;
	Interpolation InterPolNormalCa;
	int I_tmp;
	double ***n1,***n2;
	double D1,D2,r,rMinus1;
	double costeta,sinteta;
	double *teta;///< Contact angle
	double D_tmp,epsilon;
	double InvSqrt2;
	unsigned int *PtrOppositeCa;
	short int *MapWallId;//Map WallId To LocalId
	short int *MapCornerConcaveId;//Map CornerCancaveId To LocalId
	short int *MapCornerConvexId;//Map CornerConvexId To LocalId

//Tools
inline	int& IntRef(int I_input){I_tmp=I_input;return I_tmp;};
inline	void Normalise(double* &Var_x,double* &Var_y, int nodenumber){
	// Normalise
	D_tmp=sqrt(Var_x[nodenumber]*Var_x[nodenumber]+Var_y[nodenumber]*Var_y[nodenumber]);
	if(D_tmp>0)
		{Var_x[nodenumber]/=D_tmp;Var_y[nodenumber]/=D_tmp;}
	else
		{Var_x[nodenumber]=0.0; Var_y[nodenumber]=0.0;}
}
inline	void Normalise(double* &Var_x,double* &Var_y,double* &Var_z, int nodenumber){
	// Normalise
	D_tmp=sqrt(Var_x[nodenumber]*Var_x[nodenumber]+Var_y[nodenumber]*Var_y[nodenumber]+Var_z[nodenumber]*Var_z[nodenumber]);
	if(D_tmp>0)
		{Var_x[nodenumber]/=D_tmp;Var_y[nodenumber]/=D_tmp;Var_z[nodenumber]/=D_tmp;}
	else
		{Var_x[nodenumber]=0.0; Var_y[nodenumber]=0.0;Var_z[nodenumber]=0.0;}
}
};

#endif /* SRC_ALGORITHM_TOOLS_CONTACTANGLE_H_ */
