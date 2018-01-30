/*
 * Gradients.cpp
 *
 *  Created on: 30 Jul 2016
 *      Author: Thomas Burel
 */
#include "../Tools/Gradients.h"

Gradients::Gradients(){
	grad=0;
	Type=ModelEnum::FD;
	dimension=0;
	nb_Vel=0;
}
void Gradients::initGradients(int dimension_, int nb_vel,ModelEnum::GradientType Type_){
	Type=ModelEnum::FD;
	dimension=dimension_;
	nb_Vel=nb_vel;
	SelectGradientType(Type_);
}
Gradients::~Gradients(){
	delete grad;
}
void Gradients::SelectGradientType(ModelEnum::GradientType Type_){
//	if(Type_!=Type)
//	{
		delete grad;
		Type=Type_;
	// Add new gradient type here
		switch(Type)
		{
		case ModelEnum::FD:
			grad=new GradientsFD(dimension, nb_Vel);
			break;
		case ModelEnum::LBMStencil:
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



