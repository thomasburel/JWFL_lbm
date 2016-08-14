/*
 * GradientsFD.cpp
 *
 *  Created on: 28 Jul 2016
 *      Author: thomas
 */

#include "GradientsDEF.h"
#include <iostream>
GradientsDEF::GradientsDEF() {
	gradient_scalar=0;
	gradient_vector=0;
	dimension=0;
	nb_Vel=0;
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

/*GradientsDEF::GradientsDEF(int dimension_, int nb_vel){
	dimension=dimension_;
	nb_Vel=nb_vel;
	gradient_scalar=new double[dimension];
	gradient_vector=new double*[dimension];
	for(int i=0;i<2;i++)
		gradient_vector[i]=new double[dimension];
}*/
GradientsFD::GradientsFD() {
	gradient_scalar=0;
	gradient_vector=0;
	dimension=0;
}

GradientsFD::~GradientsFD() {
}

GradientsFD::GradientsFD(int dimension_, int nb_vel){
	dimension=dimension_;
	nb_Vel=nb_vel;
	gradient_scalar=new double[dimension];
	gradient_vector=new double*[dimension];
	for(int i=0;i<2;i++)
		gradient_vector[i]=new double[dimension];
}
///Finite difference scheme: centre and second order
double* GradientsFD::Grad (double *Var, int * Connect, int & normal){
	gradient_scalar[0]=0.5*(Var[Connect[1]]-Var[Connect[3]]);//second order center
	gradient_scalar[1]=0.5*(Var[Connect[2]]-Var[Connect[4]]);//second order center
	return gradient_scalar;
}
///Finite difference scheme: First order in the direction of the normal and second order (centre) in the tangent direction of the boundary
double* GradientsFD::GradBc (double *Var, int * Connect, int & normal){
	switch(normal)
	{
	case 1:
		gradient_scalar[0]=Var[Connect[1]]-Var[Connect[0]];//first order decenter
		gradient_scalar[1]=0.5*(Var[Connect[2]]-Var[Connect[4]]);//second order center
		break;
	case 2:
		gradient_scalar[0]=0.5*(Var[Connect[1]]-Var[Connect[3]]);//second order center
		gradient_scalar[1]=Var[Connect[2]]-Var[Connect[0]];//first order decenter
		break;
	case 3:
		gradient_scalar[0]=Var[Connect[0]]-Var[Connect[3]];//first order decenter
		gradient_scalar[1]=0.5*(Var[Connect[2]]-Var[Connect[4]]);//second order center
		break;
	case 4:
		gradient_scalar[0]=0.5*(Var[Connect[1]]-Var[Connect[3]]);//second order center
		gradient_scalar[1]=Var[Connect[0]]-Var[Connect[4]];//first order decenter
		 break;
//Global corner gradients
	case 5:
		gradient_scalar[0]=Var[Connect[1]]-Var[Connect[0]];//first order decenter
		gradient_scalar[1]=Var[Connect[2]]-Var[Connect[0]];//first order decenter
		break;
	case 6:
		gradient_scalar[0]=Var[Connect[0]]-Var[Connect[3]];//first order decenter
		gradient_scalar[1]=Var[Connect[2]]-Var[Connect[0]];//first order decenter
		break;
	case 7:
		gradient_scalar[0]=Var[Connect[0]]-Var[Connect[3]];//first order decenter
		gradient_scalar[1]=Var[Connect[0]]-Var[Connect[4]];//first order decenter
		break;
	case 8:
		gradient_scalar[0]=Var[Connect[1]]-Var[Connect[0]];//first order decenter
		gradient_scalar[1]=Var[Connect[0]]-Var[Connect[4]];//first order decenter
		break;
	}
	return gradient_scalar;
}
///Finite difference scheme: centre and second order
double* GradientsFD::GradCorner (double *Var, int * Connect, int & normal){
		gradient_scalar[0]=0.5*(Var[Connect[1]]-Var[Connect[3]]);//second order center
		gradient_scalar[1]=0.5*(Var[Connect[2]]-Var[Connect[4]]);//second order center
	return gradient_scalar;
}

///Finite difference scheme: centre and second order
double** GradientsFD::Grad (double *Var_x, double *Var_y, int * Connect, int & normal){
	// du/dx
	gradient_vector[0][0]=0.5*(Var_x[Connect[1]]-Var_x[Connect[3]]);//second order center
	// du/dy
	gradient_vector[0][1]=0.5*(Var_x[Connect[2]]-Var_x[Connect[4]]);//second order center
	// dv/dx
	gradient_vector[1][0]=0.5*(Var_y[Connect[1]]-Var_y[Connect[3]]);//second order center
	// dv/dy
	gradient_vector[1][1]=0.5*(Var_y[Connect[2]]-Var_y[Connect[4]]);//second order center
	return gradient_vector;
}
double** GradientsFD::GradBc (double *Var_x, double *Var_y, int * Connect, int & normal){
	switch(normal)
	{
	case 1:
		// du/dx
		gradient_vector[0][0]=Var_x[Connect[1]]-Var_x[Connect[0]];//first order decenter
		// du/dy
		gradient_vector[0][1]=0.5*(Var_x[Connect[2]]-Var_x[Connect[4]]);//second order center
		// dv/dx
		gradient_vector[1][0]=Var_y[Connect[1]]-Var_y[Connect[0]];//first order decenter
		// dv/dy
		gradient_vector[1][1]=0.5*(Var_y[Connect[2]]-Var_y[Connect[4]]);//second order center
		break;
	case 2:
		// du/dx
		gradient_vector[0][0]=0.5*(Var_x[Connect[1]]-Var_x[Connect[3]]);//second order center
		// du/dy
		gradient_vector[0][1]=Var_x[Connect[2]]-Var_x[Connect[0]];//first order decenter
		// dv/dx
		gradient_vector[1][0]=0.5*(Var_y[Connect[1]]-Var_y[Connect[3]]);//second order center
		// dv/dy
		gradient_vector[1][1]=Var_y[Connect[2]]-Var_y[Connect[0]];//first order decenter
		break;
	case 3:
		// du/dx
		gradient_vector[0][0]=Var_x[Connect[0]]-Var_x[Connect[3]];//first order decenter
		// du/dy
		gradient_vector[0][1]=0.5*(Var_x[Connect[2]]-Var_x[Connect[4]]);//second order center
		// dv/dx
		gradient_vector[1][0]=Var_y[Connect[0]]-Var_y[Connect[3]];//first order decenter
		// dv/dy
		gradient_vector[1][1]=0.5*(Var_y[Connect[2]]-Var_y[Connect[4]]);//second order center
		break;
	case 4:
		// du/dx
		gradient_vector[0][0]=0.5*(Var_x[Connect[1]]-Var_x[Connect[3]]);//second order center
		// du/dy
		gradient_vector[0][1]=Var_x[Connect[0]]-Var_x[Connect[4]];//first order decenter
		// dv/dx
		gradient_vector[1][0]=0.5*(Var_y[Connect[1]]-Var_y[Connect[3]]);//second order center
		// dv/dy
		gradient_vector[1][1]=Var_y[Connect[0]]-Var_y[Connect[4]];//first order decenter
		 break;
	}
	return gradient_vector;
}
double** GradientsFD::GradCorner (double *Var_x, double *Var_y, int * Connect, int & normal){
	// du/dx
	gradient_vector[0][0]=0.5*(Var_x[Connect[1]]-Var_x[Connect[3]]);//second order center
	// du/dy
	gradient_vector[0][1]=0.5*(Var_x[Connect[2]]-Var_x[Connect[4]]);//second order center
	// dv/dx
	gradient_vector[1][0]=0.5*(Var_y[Connect[1]]-Var_y[Connect[3]]);//second order center
	// dv/dy
	gradient_vector[1][1]=0.5*(Var_y[Connect[2]]-Var_y[Connect[4]]);//second order center
	return gradient_vector;
}
GradientsLBMStencil::GradientsLBMStencil(){
	Ei=0;
	omega=0;
}
GradientsLBMStencil::GradientsLBMStencil(int dimension_,int nb_vel){
	dimension=dimension_;
	nb_Vel=nb_vel;
	gradient_scalar=new double[dimension];
	gradient_vector=new double*[dimension];
	for(int i=0;i<2;i++)
		gradient_vector[i]=new double[dimension];
	Ei=new double*[dimension];
	for(int i=0;i<dimension;i++)
		Ei[i]=new double[nb_vel];
	omega=new double[nb_vel];

	//Distribution velocity
	Ei[0][0]= 0.0;
	Ei[1][0]= 0.0;
	Ei[0][1]= 1.0;
	Ei[1][1]= 0.0;
	Ei[0][2]= 0.0;
	Ei[1][2]= 1.0;
	Ei[0][3]= -1.0;
	Ei[1][3]= 0.0;
	Ei[0][4]= 0.0;
	Ei[1][4]= -1.0;
	Ei[0][5]= 1.0;
	Ei[1][5]= 1.0;
	Ei[0][6]= -1.0;
	Ei[1][6]= 1.0;
	Ei[0][7]= -1.0;
	Ei[1][7]= -1.0;
	Ei[0][8]= 1.0;
	Ei[1][8]= -1.0;
	//Weigh depending of the direction
	omega[0]=4.0/9.0;
	omega[1]=1.0/9.0;
	omega[2]=1.0/9.0;
	omega[3]=1.0/9.0;
	omega[4]=1.0/9.0;
	omega[5]=1.0/36.0;
	omega[6]=1.0/36.0;
	omega[7]=1.0/36.0;
	omega[8]=1.0/36.0;
}
GradientsLBMStencil::~GradientsLBMStencil(){
	for(int i=0;i<dimension;i++)
		delete [] Ei[i];
	delete [] omega;
	delete[] Ei;
}
//Scalar gradient
double* GradientsLBMStencil::Grad (double *Var, int * Connect, int & normal){
	gradient_scalar[0]=0;
	gradient_scalar[1]=0;
	for(int i=0;i<nb_Vel;i++)
	{
		gradient_scalar[0]+=omega[i]*Var[Connect[i]]*Ei[0][i];//fourth order compact scheme
		gradient_scalar[1]+=omega[i]*Var[Connect[i]]*Ei[1][i];//fourth order compact scheme
	}
	gradient_scalar[0]*=3;
	gradient_scalar[1]*=3;
return gradient_scalar;
}
double* GradientsLBMStencil::GradBc (double *Var, int * Connect, int & normal){
/*	gradient_scalar[0]=0;
	gradient_scalar[1]=0;
	for(int i=0;i<nb_Vel;i++)
	{
		gradient_scalar[0]+=omega[i]*Var[Connect[i]]*Ei[0][i];//fourth order compact scheme
		gradient_scalar[1]+=omega[i]*Var[Connect[i]]*Ei[1][i];//fourth order compact scheme
	}
	gradient_scalar[0]*=3;
	gradient_scalar[1]*=3;*/
	switch(normal)
	{
	case 1:
		gradient_scalar[0]=Var[Connect[1]]-Var[Connect[0]];//first order decenter
		gradient_scalar[1]=0.5*(Var[Connect[2]]-Var[Connect[4]]);//second order center
		break;
	case 2:
		gradient_scalar[0]=0.5*(Var[Connect[1]]-Var[Connect[3]]);//second order center
		gradient_scalar[1]=Var[Connect[2]]-Var[Connect[0]];//first order decenter
		break;
	case 3:
		gradient_scalar[0]=Var[Connect[0]]-Var[Connect[3]];//first order decenter
		gradient_scalar[1]=0.5*(Var[Connect[2]]-Var[Connect[4]]);//second order center
		break;
	case 4:
		gradient_scalar[0]=0.5*(Var[Connect[1]]-Var[Connect[3]]);//second order center
		gradient_scalar[1]=Var[Connect[0]]-Var[Connect[4]];//first order decenter
		 break;
//Global corner gradients
	case 5:
		gradient_scalar[0]=Var[Connect[1]]-Var[Connect[0]];//first order decenter
		gradient_scalar[1]=Var[Connect[2]]-Var[Connect[0]];//first order decenter
		break;
	case 6:
		gradient_scalar[0]=Var[Connect[0]]-Var[Connect[3]];//first order decenter
		gradient_scalar[1]=Var[Connect[2]]-Var[Connect[0]];//first order decenter
		break;
	case 7:
		gradient_scalar[0]=Var[Connect[0]]-Var[Connect[3]];//first order decenter
		gradient_scalar[1]=Var[Connect[0]]-Var[Connect[4]];//first order decenter
		break;
	case 8:
		gradient_scalar[0]=Var[Connect[1]]-Var[Connect[0]];//first order decenter
		gradient_scalar[1]=Var[Connect[0]]-Var[Connect[4]];//first order decenter
		break;
	}
	return gradient_scalar;
return gradient_scalar;
}
double* GradientsLBMStencil::GradCorner (double *Var, int * Connect, int & normal){
/*	gradient_scalar[0]=0;
	gradient_scalar[1]=0;
	for(int i=0;i<nb_Vel;i++)
	{
		gradient_scalar[0]+=omega[i]*Var[Connect[i]]*Ei[0][i];//fourth order compact scheme
		gradient_scalar[1]+=omega[i]*Var[Connect[i]]*Ei[1][i];//fourth order compact scheme
	}
	gradient_scalar[0]*=3;
	gradient_scalar[1]*=3;*/
	gradient_scalar[0]=0.5*(Var[Connect[1]]-Var[Connect[3]]);//second order center
	gradient_scalar[1]=0.5*(Var[Connect[2]]-Var[Connect[4]]);//second order center
return gradient_scalar;
}

//Vector gradient
double** GradientsLBMStencil::Grad (double *Var_x, double *Var_y, int * Connect, int & normal){
	// du/dx
	gradient_vector[0][0]=0.5*(Var_x[Connect[1]]-Var_x[Connect[3]]);//second order center
	// du/dy
	gradient_vector[0][1]=0.5*(Var_x[Connect[2]]-Var_x[Connect[4]]);//second order center
	// dv/dx
	gradient_vector[1][0]=0.5*(Var_y[Connect[1]]-Var_y[Connect[3]]);//second order center
	// dv/dy
	gradient_vector[1][1]=0.5*(Var_y[Connect[2]]-Var_y[Connect[4]]);//second order center
	return gradient_vector;
}
double** GradientsLBMStencil::GradBc (double *Var_x, double *Var_y, int * Connect, int & normal){
	// du/dx
	gradient_vector[0][0]=0.5*(Var_x[Connect[1]]-Var_x[Connect[3]]);//second order center
	// du/dy
	gradient_vector[0][1]=0.5*(Var_x[Connect[2]]-Var_x[Connect[4]]);//second order center
	// dv/dx
	gradient_vector[1][0]=0.5*(Var_y[Connect[1]]-Var_y[Connect[3]]);//second order center
	// dv/dy
	gradient_vector[1][1]=0.5*(Var_y[Connect[2]]-Var_y[Connect[4]]);//second order center
	return gradient_vector;
}
double** GradientsLBMStencil::GradCorner (double *Var_x, double *Var_y, int * Connect, int & normal){
	// du/dx
	gradient_vector[0][0]=0.5*(Var_x[Connect[1]]-Var_x[Connect[3]]);//second order center
	// du/dy
	gradient_vector[0][1]=0.5*(Var_x[Connect[2]]-Var_x[Connect[4]]);//second order center
	// dv/dx
	gradient_vector[1][0]=0.5*(Var_y[Connect[1]]-Var_y[Connect[3]]);//second order center
	// dv/dy
	gradient_vector[1][1]=0.5*(Var_y[Connect[2]]-Var_y[Connect[4]]);//second order center
	return gradient_vector;
}
