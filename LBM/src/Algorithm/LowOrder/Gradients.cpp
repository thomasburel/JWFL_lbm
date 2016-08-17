/*
 * Gradients.cpp
 *
 *  Created on: 30 Jul 2016
 *      Author: Thomas Burel
 */
#include "Gradients.h"

Gradients::Gradients(){
	grad=0;
	Type=FD;
	dimension=0;
	nb_Vel=0;
}
void Gradients::initGradients(int dimension_, int nb_vel,GradientType Type_){
	Type=FD;
	dimension=dimension_;
	nb_Vel=nb_vel;
	SelectGradientType(Type_);
}
Gradients::~Gradients(){
	delete grad;
}
void Gradients::SelectGradientType(GradientType Type_){
//	if(Type_!=Type)
//	{
		delete grad;
		Type=Type_;
	// Add new gradient type here
		switch(Type)
		{
		case FD:
			grad=new GradientsFD(dimension, nb_Vel);
			break;
		case LBMStencil:
			grad=new GradientsLBMStencil(dimension, nb_Vel);
			break;
		default:
			std::cerr<<" Gradient type not found"<<std::endl;
		}
//	}
}

//Scalar gradients
void Gradients::Grad (double* grad_, double *Var, int * Connect, int & normal){
	grad->Grad(grad_,Var,Connect,normal);
}
void Gradients::GradBc (double* grad_, double *Var, int * Connect, int & normal){
	grad->GradBc(grad_,Var,Connect,normal);
}
void Gradients::GradCorner (double* grad_, double *Var, int * Connect, int & normal){
	grad->GradCorner(grad_,Var,Connect,normal);
}

//Vector gradients
void Gradients::Grad (double** grad_, double *Var_x, double *Var_y, int * Connect, int & normal){
	grad->Grad(grad_,Var_x,Var_y,Connect,normal);
}
void Gradients::GradBc (double** grad_, double *Var_x, double *Var_y, int * Connect, int & normal){
	grad->GradBc(grad_,Var_x,Var_y,Connect,normal);
}
void Gradients::GradCorner (double** grad_, double *Var_x, double *Var_y, int * Connect, int & normal){
	grad->GradCorner(grad_,Var_x,Var_y,Connect,normal);
}



