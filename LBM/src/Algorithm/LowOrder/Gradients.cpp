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
}
Gradients::Gradients(int dimension_){
	Type=FD;
	dimension=dimension_;
	grad=new GradientsFD(dimension);
}
Gradients::~Gradients(){
	delete grad;
}
void Gradients::SelectGradientType(GradientType Type_){
	if(Type_!=Type)
	{
		delete grad;
		Type=Type_;
	// Add new gradient type here
		switch(Type)
		{
		case FD:
			grad=new GradientsFD(dimension);
			break;
		}
	}
}

//Scalar gradients
double* Gradients::Grad (double *Var, NodeInterior2D& Node){
	return grad->Grad(Var,Node);
}
double* Gradients::Grad (double *Var, NodeWall2D& Node){
	return grad->Grad(Var,Node);
}
double* Gradients::Grad (double *Var, NodeCorner2D& Node){
	return grad->Grad(Var,Node);
}
double* Gradients::Grad (double *Var, NodeVelocity2D& Node){
	return grad->Grad(Var,Node);
}
double* Gradients::Grad (double *Var, NodePressure2D& Node){
	return grad->Grad(Var,Node);
}
double* Gradients::Grad (double *Var, NodeSymmetry2D& Node){
	return grad->Grad(Var,Node);
}
//Vector gradients
double** Gradients::Grad (double *Var_x, double *Var_y, NodeInterior2D& Node){
	return grad->Grad(Var_x,Var_y,Node);
}
double** Gradients::Grad (double *Var_x, double *Var_y, NodeWall2D& Node){
	return grad->Grad(Var_x,Var_y,Node);
}
double** Gradients::Grad (double *Var_x, double *Var_y, NodeCorner2D& Node){
	return grad->Grad(Var_x,Var_y,Node);
}
double** Gradients::Grad (double *Var_x, double *Var_y, NodeVelocity2D& Node){
	return grad->Grad(Var_x,Var_y,Node);
}
double** Gradients::Grad (double *Var_x, double *Var_y, NodePressure2D& Node){
	return grad->Grad(Var_x,Var_y,Node);
}
double** Gradients::Grad (double *Var_x, double *Var_y, NodeSymmetry2D& Node){
	return grad->Grad(Var_x,Var_y,Node);
}


