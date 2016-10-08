/*
 * GradientsFD.cpp
 *
 *  Created on: 28 Jul 2016
 *      Author: thomas
 */

#include "GradientsDEF.h"
#include <iostream>
GradientsDEF::GradientsDEF() {
/*	gradient_scalar=0;
	gradient_vector=0;*/
	dimension=0;
	nb_Vel=0;
}
GradientsDEF::GradientsDEF(int dimension_, int nb_vel) {
/*	gradient_scalar=0;
	gradient_vector=0;*/
	dimension=dimension_;
	nb_Vel=nb_vel;
}
GradientsDEF::~GradientsDEF() {
/*	if(dimension!=0)
	{
		delete[]gradient_scalar;
		for(int i=0;i<dimension;i++)
			delete[] gradient_vector[i];
		delete[] gradient_vector;
	}*/
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
/*	gradient_scalar=0;
	gradient_vector=0;*/
	dimension=0;
	nb_Vel=0;
}

GradientsFD::~GradientsFD() {
}

GradientsFD::GradientsFD(int dimension_, int nb_vel){
	dimension=dimension_;
	nb_Vel=nb_vel;
/*	gradient_scalar=new double[dimension];
	gradient_vector=new double*[dimension];
	for(int i=0;i<2;i++)
		gradient_vector[i]=new double[dimension];*/
}
///Finite difference scheme: centre and second order
void GradientsFD::Grad (double* grad_, double *Var, int * Connect, int & normal){
	grad_[0]=0.5*(Var[Connect[1]]-Var[Connect[3]]);//second order center
	grad_[1]=0.5*(Var[Connect[2]]-Var[Connect[4]]);//second order center
}
///Finite difference scheme: First order in the direction of the normal and second order (centre) in the tangent direction of the boundary
void GradientsFD::GradBc (double* grad_, double *Var, int * Connect, int & normal){
	switch(normal)
	{
	case 1:
		grad_[0]=Var[Connect[1]]-Var[Connect[0]];//first order decenter
		grad_[1]=0.5*(Var[Connect[2]]-Var[Connect[4]]);//second order center
		break;
	case 2:
		grad_[0]=0.5*(Var[Connect[1]]-Var[Connect[3]]);//second order center
		grad_[1]=Var[Connect[2]]-Var[Connect[0]];//first order decenter
		break;
	case 3:
		grad_[0]=Var[Connect[0]]-Var[Connect[3]];//first order decenter
		grad_[1]=0.5*(Var[Connect[2]]-Var[Connect[4]]);//second order center
		break;
	case 4:
		grad_[0]=0.5*(Var[Connect[1]]-Var[Connect[3]]);//second order center
		grad_[1]=Var[Connect[0]]-Var[Connect[4]];//first order decenter
		 break;
//Global corner gradients
	case 5:
		grad_[0]=Var[Connect[1]]-Var[Connect[0]];//first order decenter
		grad_[1]=Var[Connect[2]]-Var[Connect[0]];//first order decenter
		break;
	case 6:
		grad_[0]=Var[Connect[0]]-Var[Connect[3]];//first order decenter
		grad_[1]=Var[Connect[2]]-Var[Connect[0]];//first order decenter
		break;
	case 7:
		grad_[0]=Var[Connect[0]]-Var[Connect[3]];//first order decenter
		grad_[1]=Var[Connect[0]]-Var[Connect[4]];//first order decenter
		break;
	case 8:
		grad_[0]=Var[Connect[1]]-Var[Connect[0]];//first order decenter
		grad_[1]=Var[Connect[0]]-Var[Connect[4]];//first order decenter
		break;
	}
}
///Finite difference scheme: centre and second order
void GradientsFD::GradCorner (double* grad_, double *Var, int * Connect, int & normal){
		grad_[0]=0.5*(Var[Connect[1]]-Var[Connect[3]]);//second order center
		grad_[1]=0.5*(Var[Connect[2]]-Var[Connect[4]]);//second order center
}

///Finite difference scheme: centre and second order
void GradientsFD::Grad (double** grad_, double *Var_x, double *Var_y, int * Connect, int & normal){
	// du/dx
	grad_[0][0]=0.5*(Var_x[Connect[1]]-Var_x[Connect[3]]);//second order center
	// du/dy
	grad_[0][1]=0.5*(Var_x[Connect[2]]-Var_x[Connect[4]]);//second order center
	// dv/dx
	grad_[1][0]=0.5*(Var_y[Connect[1]]-Var_y[Connect[3]]);//second order center
	// dv/dy
	grad_[1][1]=0.5*(Var_y[Connect[2]]-Var_y[Connect[4]]);//second order center

}
void GradientsFD::GradBc (double** grad_, double *Var_x, double *Var_y, int * Connect, int & normal){
	switch(normal)
	{
	case 1:
		// du/dx
		grad_[0][0]=Var_x[Connect[1]]-Var_x[Connect[0]];//first order decenter
		// du/dy
		grad_[0][1]=0.5*(Var_x[Connect[2]]-Var_x[Connect[4]]);//second order center
		// dv/dx
		grad_[1][0]=Var_y[Connect[1]]-Var_y[Connect[0]];//first order decenter
		// dv/dy
		grad_[1][1]=0.5*(Var_y[Connect[2]]-Var_y[Connect[4]]);//second order center
		break;
	case 2:
		// du/dx
		grad_[0][0]=0.5*(Var_x[Connect[1]]-Var_x[Connect[3]]);//second order center
		// du/dy
		grad_[0][1]=Var_x[Connect[2]]-Var_x[Connect[0]];//first order decenter
		// dv/dx
		grad_[1][0]=0.5*(Var_y[Connect[1]]-Var_y[Connect[3]]);//second order center
		// dv/dy
		grad_[1][1]=Var_y[Connect[2]]-Var_y[Connect[0]];//first order decenter
		break;
	case 3:
		// du/dx
		grad_[0][0]=Var_x[Connect[0]]-Var_x[Connect[3]];//first order decenter
		// du/dy
		grad_[0][1]=0.5*(Var_x[Connect[2]]-Var_x[Connect[4]]);//second order center
		// dv/dx
		grad_[1][0]=Var_y[Connect[0]]-Var_y[Connect[3]];//first order decenter
		// dv/dy
		grad_[1][1]=0.5*(Var_y[Connect[2]]-Var_y[Connect[4]]);//second order center
		break;
	case 4:
		// du/dx
		grad_[0][0]=0.5*(Var_x[Connect[1]]-Var_x[Connect[3]]);//second order center
		// du/dy
		grad_[0][1]=Var_x[Connect[0]]-Var_x[Connect[4]];//first order decenter
		// dv/dx
		grad_[1][0]=0.5*(Var_y[Connect[1]]-Var_y[Connect[3]]);//second order center
		// dv/dy
		grad_[1][1]=Var_y[Connect[0]]-Var_y[Connect[4]];//first order decenter
		 break;
	}

}
void GradientsFD::GradCorner (double** grad_, double *Var_x, double *Var_y, int * Connect, int & normal){
	// du/dx
	grad_[0][0]=0.5*(Var_x[Connect[1]]-Var_x[Connect[3]]);//second order center
	// du/dy
	grad_[0][1]=0.5*(Var_x[Connect[2]]-Var_x[Connect[4]]);//second order center
	// dv/dx
	grad_[1][0]=0.5*(Var_y[Connect[1]]-Var_y[Connect[3]]);//second order center
	// dv/dy
	grad_[1][1]=0.5*(Var_y[Connect[2]]-Var_y[Connect[4]]);//second order center

}
GradientsLBMStencil::GradientsLBMStencil(){
	Ei=0;
	omega=0;
}
GradientsLBMStencil::GradientsLBMStencil(int dimension_,int nb_vel){
	dimension=dimension_;
	nb_Vel=nb_vel;
/*	gradient_scalar=new double[dimension];
	gradient_vector=new double*[dimension];
	for(int i=0;i<2;i++)
		gradient_vector[i]=new double[dimension];*/
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
void GradientsLBMStencil::Grad (double* grad_, double *Var, int * Connect, int & normal){
	grad_[0]=0;
	grad_[1]=0;
	for(int i=0;i<nb_Vel;i++)
	{
		grad_[0]+=omega[i]*Var[Connect[i]]*Ei[0][i];//fourth order compact scheme
		grad_[1]+=omega[i]*Var[Connect[i]]*Ei[1][i];//fourth order compact scheme
	}
	grad_[0]*=3.0;
	grad_[1]*=3.0;
}
void GradientsLBMStencil::GradBc (double* grad_, double *Var, int * Connect, int & normal){
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
		grad_[0]=Var[Connect[1]]-Var[Connect[0]];//first order decenter
		grad_[1]=0.5*(Var[Connect[2]]-Var[Connect[4]]);//second order center
		break;
	case 2:
		grad_[0]=0.5*(Var[Connect[1]]-Var[Connect[3]]);//second order center
		grad_[1]=Var[Connect[2]]-Var[Connect[0]];//first order decenter
		break;
	case 3:
		grad_[0]=Var[Connect[0]]-Var[Connect[3]];//first order decenter
		grad_[1]=0.5*(Var[Connect[2]]-Var[Connect[4]]);//second order center
		break;
	case 4:
		grad_[0]=0.5*(Var[Connect[1]]-Var[Connect[3]]);//second order center
		grad_[1]=Var[Connect[0]]-Var[Connect[4]];//first order decenter
		 break;
//Global corner gradients
	case 5:
		grad_[0]=Var[Connect[1]]-Var[Connect[0]];//first order decenter
		grad_[1]=Var[Connect[2]]-Var[Connect[0]];//first order decenter
		break;
	case 6:
		grad_[0]=Var[Connect[0]]-Var[Connect[3]];//first order decenter
		grad_[1]=Var[Connect[2]]-Var[Connect[0]];//first order decenter
		break;
	case 7:
		grad_[0]=Var[Connect[0]]-Var[Connect[3]];//first order decenter
		grad_[1]=Var[Connect[0]]-Var[Connect[4]];//first order decenter
		break;
	case 8:
		grad_[0]=Var[Connect[1]]-Var[Connect[0]];//first order decenter
		grad_[1]=Var[Connect[0]]-Var[Connect[4]];//first order decenter
		break;
	}
}
void GradientsLBMStencil::GradCorner (double* grad_,  double *Var, int * Connect, int & normal){
/*	gradient_scalar[0]=0;
	gradient_scalar[1]=0;
	for(int i=0;i<nb_Vel;i++)
	{
		gradient_scalar[0]+=omega[i]*Var[Connect[i]]*Ei[0][i];//fourth order compact scheme
		gradient_scalar[1]+=omega[i]*Var[Connect[i]]*Ei[1][i];//fourth order compact scheme
	}
	gradient_scalar[0]*=3;
	gradient_scalar[1]*=3;*/
	grad_[0]=0.5*(Var[Connect[1]]-Var[Connect[3]]);//second order center
	grad_[1]=0.5*(Var[Connect[2]]-Var[Connect[4]]);//second order center

}

//Vector gradient
void GradientsLBMStencil::Grad (double** grad_, double *Var_x, double *Var_y, int * Connect, int & normal){
	// du/dx
	grad_[0][0]=0.5*(Var_x[Connect[1]]-Var_x[Connect[3]]);//second order center
	// du/dy
	grad_[0][1]=0.5*(Var_x[Connect[2]]-Var_x[Connect[4]]);//second order center
	// dv/dx
	grad_[1][0]=0.5*(Var_y[Connect[1]]-Var_y[Connect[3]]);//second order center
	// dv/dy
	grad_[1][1]=0.5*(Var_y[Connect[2]]-Var_y[Connect[4]]);//second order center
}
void GradientsLBMStencil::GradBc (double** grad_, double *Var_x, double *Var_y, int * Connect, int & normal){
	// du/dx
	grad_[0][0]=0.5*(Var_x[Connect[1]]-Var_x[Connect[3]]);//second order center
	// du/dy
	grad_[0][1]=0.5*(Var_x[Connect[2]]-Var_x[Connect[4]]);//second order center
	// dv/dx
	grad_[1][0]=0.5*(Var_y[Connect[1]]-Var_y[Connect[3]]);//second order center
	// dv/dy
	grad_[1][1]=0.5*(Var_y[Connect[2]]-Var_y[Connect[4]]);//second order center

}
void GradientsLBMStencil::GradCorner (double** grad_, double *Var_x, double *Var_y, int * Connect, int & normal){
	// du/dx
	grad_[0][0]=0.5*(Var_x[Connect[1]]-Var_x[Connect[3]]);//second order center
	// du/dy
	grad_[0][1]=0.5*(Var_x[Connect[2]]-Var_x[Connect[4]]);//second order center
	// dv/dx
	grad_[1][0]=0.5*(Var_y[Connect[1]]-Var_y[Connect[3]]);//second order center
	// dv/dy
	grad_[1][1]=0.5*(Var_y[Connect[2]]-Var_y[Connect[4]]);//second order center

}
