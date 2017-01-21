/*
 * Interpolation.cpp
 *
 *  Created on: 31 Oct 2016
 *      Author: Thomas Burel
 */
#include "Interpolation.h"

Interpolation::Interpolation(){
	Interpol=0;
	Type=LinearInterpol;
	dimension=0;
	nb_Vel=0;
	PtrOppositeCa=0;
}
void Interpolation::initInterpolation(int dimension_, int nb_vel,InterpolationType Type_,unsigned int *PtrOppositeCa_, NodeArrays2D *PtrNodes, Parameters *PtrParam){
	Type=LinearInterpol;
	dimension=dimension_;
	nb_Vel=nb_vel;
	PtrOppositeCa=PtrOppositeCa_;
	SelectInterpolationType(Type_,PtrNodes,PtrParam);

}
Interpolation::~Interpolation(){
	delete Interpol;
}
void Interpolation::SelectInterpolationType(InterpolationType Type_, NodeArrays2D *PtrNodes, Parameters *PtrParam){
//	if(Type_!=Type)
//	{
		delete Interpol;
		Type=Type_;
	// Add new Interpolation type here
		switch(Type)
		{
		case NoInterpol:
			Interpol=new NoInterpolation();
			break;
		case LinearInterpol:
			Interpol=new InterpolationLinear(dimension, nb_Vel,PtrOppositeCa);
			break;
		case LinearLeastSquareInterpol:
			Interpol=new InterpolationLinearLeastSquare(dimension, nb_Vel,PtrOppositeCa);
			Interpol->InitInterpol(PtrNodes,PtrParam);
			break;
		default:
			Interpol=new InterpolationLinear(dimension, nb_Vel,PtrOppositeCa);
			std::cerr<<" Interpolation type not found"<<std::endl;
		}
//	}
}

//Scalar Interpolation
void Interpolation::InterpolationOnWall (double * & Var, int * Connect, int & normal){
	Interpol->InterpolationOnWall(Var,Connect,normal);
}
void Interpolation::InterpolationOnCornerConcave (double * & Var, int * Connect, int & normal){
	Interpol->InterpolationOnCornerConcave(Var,Connect,normal);
}
void Interpolation::InterpolationOnCornerConvex (double * & Var, int * Connect, int & normal){
	Interpol->InterpolationOnCornerConvex(Var,Connect,normal);
}
void Interpolation::InterpolationOnWall (double * & Var1, double * & Var2, int * Connect, int & normal){
	Interpol->InterpolationOnWall(Var1,Var2,Connect,normal);
}
void Interpolation::InterpolationOnCornerConcave (double * & Var1, double * & Var2, int * Connect, int & normal){
	Interpol->InterpolationOnCornerConcave(Var1,Var2,Connect,normal);
}
void Interpolation::InterpolationOnCornerConvex (double * & Var1, double * & Var2, int * Connect, int & normal){
	Interpol->InterpolationOnCornerConvex(Var1,Var2,Connect,normal);
}



