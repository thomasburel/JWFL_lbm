/*
 * Extrapolation.cpp
 *
 *  Created on: 31 Oct 2016
 *      Author: Thomas Burel
 */
#include "../Tools/Extrapolation.h"

Extrapolation::Extrapolation(){
	Extrapol=0;
	Type=ModelEnum::TailorExtrapol;
	dimension=0;
	nb_Vel=0;
}
void Extrapolation::initExtrapolation(int dimension_, int nb_vel,ModelEnum::ExtrapolationType Type_){
	Type=ModelEnum::TailorExtrapol;
	dimension=dimension_;
	nb_Vel=nb_vel;
	SelectExtrapolationType(Type_);
}
Extrapolation::~Extrapolation(){
	delete Extrapol;
}
void Extrapolation::SelectExtrapolationType(ModelEnum::ExtrapolationType Type_){
//	if(Type_!=Type)
//	{
		delete Extrapol;
		Type=Type_;
	// Add new Extrapolation type here
		switch(Type)
		{
		case ModelEnum::NoExtrapol:
			Extrapol=new NoExtrapolation();
			break;
		case ModelEnum::TailorExtrapol:
			Extrapol=new ExtrapolationTailor(dimension, nb_Vel);
			break;
		case ModelEnum::WeightDistanceExtrapol:
			Extrapol=new ExtrapolationWeightDistance(dimension, nb_Vel);
			break;
		default:
			Extrapol=new ExtrapolationWeightDistance(dimension, nb_Vel);
			std::cerr<<" Extrapolation type not found"<<std::endl;
		}
//	}
}

//Scalar Extrapolation
void Extrapolation::ExtrapolationOnWall (double * & Var, int * Connect, int & normal){
	Extrapol->ExtrapolationOnWall(Var,Connect,normal);
}
void Extrapolation::ExtrapolationOnCornerConcave (double * & Var, int * Connect, int & normal){
	Extrapol->ExtrapolationOnCornerConcave(Var,Connect,normal);
}
void Extrapolation::ExtrapolationOnCornerConvex (double * & Var, int * Connect, int & normal){
	Extrapol->ExtrapolationOnCornerConvex(Var,Connect,normal);
}
void Extrapolation::ExtrapolationWallToSolid (double * & Var, int * Connect, int & normal){
	Extrapol->ExtrapolationWallToSolid(Var,Connect,normal);
}
void Extrapolation::ExtrapolationCornerConcaveToSolid (double * & Var, int * Connect, int & normal){
	Extrapol->ExtrapolationCornerConcaveToSolid(Var,Connect,normal);
}
void Extrapolation::ExtrapolationCornerConvexToSolid (double * & Var, int * Connect, int & normal){
	Extrapol->ExtrapolationCornerConvexToSolid(Var,Connect,normal);
}




