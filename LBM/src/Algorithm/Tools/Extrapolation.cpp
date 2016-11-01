/*
 * Extrapolation.cpp
 *
 *  Created on: 31 Oct 2016
 *      Author: Thomas Burel
 */
#include "../Tools/Extrapolation.h"

Extrapolation::Extrapolation(){
	Extrapol=0;
	Type=TailorExtrapol;
	dimension=0;
	nb_Vel=0;
}
void Extrapolation::initExtrapolation(int dimension_, int nb_vel,ExtrapolationType Type_){
	Type=TailorExtrapol;
	dimension=dimension_;
	nb_Vel=nb_vel;
	SelectExtrapolationType(Type_);
}
Extrapolation::~Extrapolation(){
	delete Extrapol;
}
void Extrapolation::SelectExtrapolationType(ExtrapolationType Type_){
//	if(Type_!=Type)
//	{
		delete Extrapol;
		Type=Type_;
	// Add new Extrapolation type here
		switch(Type)
		{
		case TailorExtrapol:
			Extrapol=new ExtrapolationTailor(dimension, nb_Vel);
			break;
		case WeightDistanceExtrapol:
			Extrapol=new ExtrapolationWeightDistance(dimension, nb_Vel);
			break;
		default:
			std::cerr<<" Extrapolation type not found"<<std::endl;
		}
//	}
}

//Scalar Extrapolation
void Extrapolation::ExtrapolationWall (double *Var, int * Connect, int & normal){
	Extrapol->ExtrapolationWall(Var,Connect,normal);
}
void Extrapolation::ExtrapolationCornerConcave (double *Var, int * Connect, int & normal){
	Extrapol->ExtrapolationCornerConcave(Var,Connect,normal);
}
void Extrapolation::ExtrapolationCornerConvex (double *Var, int * Connect, int & normal){
	Extrapol->ExtrapolationCornerConvex(Var,Connect,normal);
}




