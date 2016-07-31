/*
 * GradientsFD.cpp
 *
 *  Created on: 28 Jul 2016
 *      Author: thomas
 */

#include "GradientsDEF.h"

GradientsDEF::GradientsDEF() {
	gradient_scalar=0;
	gradient_vector=0;
	dimension=0;
}

GradientsDEF::~GradientsDEF() {
	if(dimension!=0)
	{
		delete[]gradient_scalar;
		for(int i=0;i<dimension;i++)
			delete[] gradient_vector[i];
		delete[] gradient_vector;
	}
}

GradientsDEF::GradientsDEF(int dimension_){
	dimension=dimension_;
	gradient_scalar=new double[dimension];
	gradient_vector=new double*[dimension];
	for(int i=0;i<2;i++)
		gradient_vector[i]=new double[dimension];
}
GradientsFD::GradientsFD() {
	gradient_scalar=0;
	gradient_vector=0;
	dimension=0;
}

GradientsFD::~GradientsFD() {
}

GradientsFD::GradientsFD(int dimension_){
	dimension=dimension_;
	gradient_scalar=new double[dimension];
	gradient_vector=new double*[dimension];
	for(int i=0;i<2;i++)
		gradient_vector[i]=new double[dimension];
}
///Finite difference scheme: centre and second order
double* GradientsFD::Grad (double *Var, NodeInterior2D& Node){
	gradient_scalar[0]=Var[Node.Get_connect()[1]]-Var[Node.Get_connect()[3]];
	gradient_scalar[1]=Var[Node.Get_connect()[2]]-Var[Node.Get_connect()[4]];
	return gradient_scalar;
}
///Finite difference scheme: First order in the direction of the normal and second order (centre) in the tangent direction of the boundary
double* GradientsFD::Grad (double *Var, NodeWall2D& Node){
	switch(Node.Get_BcNormal())
	{
	case 1:
		gradient_scalar[0]=Var[Node.Get_connect()[1]]-Var[Node.Get_index()];
		gradient_scalar[1]=Var[Node.Get_connect()[2]]-Var[Node.Get_connect()[4]];
		break;
	case 2:
		gradient_scalar[0]=Var[Node.Get_connect()[1]]-Var[Node.Get_connect()[3]];
		gradient_scalar[1]=Var[Node.Get_connect()[2]]-Var[Node.Get_index()];
		break;
	case 3:
		gradient_scalar[0]=Var[Node.Get_index()]-Var[Node.Get_connect()[3]];
		gradient_scalar[1]=Var[Node.Get_connect()[2]]-Var[Node.Get_connect()[4]];
		break;
	case 4:
		gradient_scalar[0]=Var[Node.Get_connect()[1]]-Var[Node.Get_connect()[3]];
		gradient_scalar[1]=Var[Node.Get_index()]-Var[Node.Get_connect()[4]];
		 break;
	}
	return gradient_scalar;
}
///Finite difference scheme: centre and second order
double* GradientsFD::Grad (double *Var, NodeCorner2D& Node){
		gradient_scalar[0]=Var[Node.Get_connect()[1]]-Var[Node.Get_connect()[3]];
		gradient_scalar[1]=Var[Node.Get_connect()[2]]-Var[Node.Get_connect()[4]];
	return gradient_scalar;
}
///Finite difference scheme: First order in the direction of the normal and second order (centre) in the tangent direction of the boundary
double* GradientsFD::Grad (double *Var, NodeVelocity2D& Node){
	switch(Node.Get_BcNormal())
	{
	case 1:
		gradient_scalar[0]=Var[Node.Get_connect()[1]]-Var[Node.Get_index()];
		gradient_scalar[1]=Var[Node.Get_connect()[2]]-Var[Node.Get_connect()[4]];
		break;
	case 2:
		gradient_scalar[0]=Var[Node.Get_connect()[1]]-Var[Node.Get_connect()[3]];
		gradient_scalar[1]=Var[Node.Get_connect()[2]]-Var[Node.Get_index()];
		break;
	case 3:
		gradient_scalar[0]=Var[Node.Get_index()]-Var[Node.Get_connect()[3]];
		gradient_scalar[1]=Var[Node.Get_connect()[2]]-Var[Node.Get_connect()[4]];
		break;
	case 4:
		gradient_scalar[0]=Var[Node.Get_connect()[1]]-Var[Node.Get_connect()[3]];
		gradient_scalar[1]=Var[Node.Get_index()]-Var[Node.Get_connect()[4]];
		 break;
	}
	return gradient_scalar;
}
///Finite difference scheme: First order in the direction of the normal and second order (centre) in the tangent direction of the boundary
double* GradientsFD::Grad (double *Var, NodePressure2D& Node){
	switch(Node.Get_BcNormal())
	{
	case 1:
		gradient_scalar[0]=Var[Node.Get_connect()[1]]-Var[Node.Get_index()];
		gradient_scalar[1]=Var[Node.Get_connect()[2]]-Var[Node.Get_connect()[4]];
		break;
	case 2:
		gradient_scalar[0]=Var[Node.Get_connect()[1]]-Var[Node.Get_connect()[3]];
		gradient_scalar[1]=Var[Node.Get_connect()[2]]-Var[Node.Get_index()];
		break;
	case 3:
		gradient_scalar[0]=Var[Node.Get_index()]-Var[Node.Get_connect()[3]];
		gradient_scalar[1]=Var[Node.Get_connect()[2]]-Var[Node.Get_connect()[4]];
		break;
	case 4:
		gradient_scalar[0]=Var[Node.Get_connect()[1]]-Var[Node.Get_connect()[3]];
		gradient_scalar[1]=Var[Node.Get_index()]-Var[Node.Get_connect()[4]];
		 break;
	}
	return gradient_scalar;
}
///Finite difference scheme: First order in the direction of the normal and second order (centre) in the tangent direction of the boundary
double* GradientsFD::Grad (double *Var, NodeSymmetry2D& Node){
	switch(Node.Get_BcNormal())
	{
	case 1:
		gradient_scalar[0]=Var[Node.Get_connect()[1]]-Var[Node.Get_index()];
		gradient_scalar[1]=Var[Node.Get_connect()[2]]-Var[Node.Get_connect()[4]];
		break;
	case 2:
		gradient_scalar[0]=Var[Node.Get_connect()[1]]-Var[Node.Get_connect()[3]];
		gradient_scalar[1]=Var[Node.Get_connect()[2]]-Var[Node.Get_index()];
		break;
	case 3:
		gradient_scalar[0]=Var[Node.Get_index()]-Var[Node.Get_connect()[3]];
		gradient_scalar[1]=Var[Node.Get_connect()[2]]-Var[Node.Get_connect()[4]];
		break;
	case 4:
		gradient_scalar[0]=Var[Node.Get_connect()[1]]-Var[Node.Get_connect()[3]];
		gradient_scalar[1]=Var[Node.Get_index()]-Var[Node.Get_connect()[4]];
		 break;
	}
	return gradient_scalar;
}
///Finite difference scheme: centre and second order
double** GradientsFD::Grad (double *Var_x, double *Var_y, NodeInterior2D& Node){
	// du/dx
	gradient_vector[0][0]=Var_x[Node.Get_connect()[1]]-Var_x[Node.Get_connect()[3]];
	// du/dy
	gradient_vector[0][1]=Var_x[Node.Get_connect()[2]]-Var_x[Node.Get_connect()[4]];
	// dv/dx
	gradient_vector[1][0]=Var_y[Node.Get_connect()[1]]-Var_y[Node.Get_connect()[3]];
	// dv/dy
	gradient_vector[1][1]=Var_y[Node.Get_connect()[2]]-Var_y[Node.Get_connect()[4]];
	return gradient_vector;
}
double** GradientsFD::Grad (double *Var_x, double *Var_y, NodeWall2D& Node){
	return gradient_vector;
}
double** GradientsFD::Grad (double *Var_x, double *Var_y, NodeCorner2D& Node){
	return gradient_vector;
}
double** GradientsFD::Grad (double *Var_x, double *Var_y, NodeVelocity2D& Node){
	return gradient_vector;
}
double** GradientsFD::Grad (double *Var_x, double *Var_y, NodePressure2D& Node){
	return gradient_vector;
}
double** GradientsFD::Grad (double *Var_x, double *Var_y, NodeSymmetry2D& Node){
	return gradient_vector;
}
