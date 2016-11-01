/*
 * ExtrapolationTailor.cpp
 *
 *  Created on: 28 Jul 2016
 *      Author: thomas
 */

#include "../Tools/ExtrapolationDEF.h"

#include <iostream>
ExtrapolationDEF::ExtrapolationDEF() {
	dimension=0;
	nb_Vel=0;
}
ExtrapolationDEF::ExtrapolationDEF(int dimension_, int nb_vel) {
	dimension=dimension_;
	nb_Vel=nb_vel;
}
ExtrapolationDEF::~ExtrapolationDEF() {

}


ExtrapolationTailor::ExtrapolationTailor() {
	dimension=0;
	nb_Vel=0;
}

ExtrapolationTailor::~ExtrapolationTailor() {
}

ExtrapolationTailor::ExtrapolationTailor(int dimension_, int nb_vel){
	dimension=dimension_;
	nb_Vel=nb_vel;
}
///Second order Tailor extrapolation
void ExtrapolationTailor::ExtrapolationWall (double *Var, int * Connect, int & normal){

	switch(normal)
	{
	case 1:
		Var[Connect[3]]=2*Var[Connect[0]]-Var[Connect[1]];
		break;
	case 2:
		Var[Connect[4]]=2*Var[Connect[0]]-Var[Connect[2]];
		break;
	case 3:
		Var[Connect[1]]=2*Var[Connect[0]]-Var[Connect[3]];
		break;
	case 4:
		Var[Connect[2]]=2*Var[Connect[0]]-Var[Connect[4]];
		break;
	}
}
///Second order Tailor extrapolation
void ExtrapolationTailor::ExtrapolationCornerConcave (double *Var, int * Connect, int & normal){
	switch(normal)
	{
	case 5:
		Var[Connect[3]]=2*Var[Connect[0]]-Var[Connect[1]];
		Var[Connect[4]]=2*Var[Connect[0]]-Var[Connect[2]];
		Var[Connect[7]]=2*Var[Connect[0]]-Var[Connect[5]];
		break;
	case 6:
		Var[Connect[1]]=2*Var[Connect[0]]-Var[Connect[3]];
		Var[Connect[4]]=2*Var[Connect[0]]-Var[Connect[2]];
		Var[Connect[8]]=2*Var[Connect[0]]-Var[Connect[6]];
		break;
	case 7:
		Var[Connect[1]]=2*Var[Connect[0]]-Var[Connect[3]];
		Var[Connect[2]]=2*Var[Connect[0]]-Var[Connect[4]];
		Var[Connect[5]]=2*Var[Connect[0]]-Var[Connect[7]];
		break;
	case 8:
		Var[Connect[2]]=2*Var[Connect[0]]-Var[Connect[4]];
		Var[Connect[3]]=2*Var[Connect[0]]-Var[Connect[1]];
		Var[Connect[6]]=2*Var[Connect[0]]-Var[Connect[8]];
		break;
	}
}
///Second order Tailor extrapolation
void ExtrapolationTailor::ExtrapolationCornerConvex (double *Var, int * Connect, int & normal){
	switch(normal)
	{
	case 5:
		Var[Connect[7]]=2*Var[Connect[0]]-Var[Connect[5]];
		break;
	case 6:
		Var[Connect[8]]=2*Var[Connect[0]]-Var[Connect[6]];
		break;
	case 7:
		Var[Connect[5]]=2*Var[Connect[0]]-Var[Connect[7]];
		break;
	case 8:
		Var[Connect[6]]=2*Var[Connect[0]]-Var[Connect[8]];
		break;
	}
}


ExtrapolationWeightDistance::ExtrapolationWeightDistance(){
	InvSqrt2=std::sqrt(0.5);
	InvSqrt2_5=std::sqrt(0.4);//0.4=2/5
	InvSumWeightWall=1.0/(1.0+2.0*InvSqrt2);
	InvSumWeightCornerConvex=1.0/(2.0+1.0*InvSqrt2);
	InvSumWeightCornerConcave=1.0/(1.0+2.0*InvSqrt2_5);//0.4=2/5
	InvSumWeightCornerWall=1.0/(1.0+1.0*InvSqrt2);
}
ExtrapolationWeightDistance::ExtrapolationWeightDistance(int dimension_,int nb_vel){
	dimension=dimension_;
	nb_Vel=nb_vel;
	InvSqrt2=std::sqrt(0.5);
	InvSqrt2_5=std::sqrt(0.4);//0.4=2/5
	InvSumWeightWall=1.0/(1.0+2.0*InvSqrt2);
	InvSumWeightCornerConvex=1.0/(2.0+1.0*InvSqrt2);
	InvSumWeightCornerConcave=1.0/(1.0+2.0*InvSqrt2_5);//0.4=2/5
	InvSumWeightCornerWall=1.0/(1.0+1.0*InvSqrt2);
}
ExtrapolationWeightDistance::~ExtrapolationWeightDistance(){

}

void ExtrapolationWeightDistance::ExtrapolationWall (double *Var, int * Connect, int & normal){
	switch(normal)
	{
	case 1:
		Var[Connect[3]]=InvSumWeightWall*(InvSqrt2*(Var[Connect[2]]+Var[Connect[4]])+Var[Connect[0]]);
		break;
	case 2:
		Var[Connect[4]]=InvSumWeightWall*(InvSqrt2*(Var[Connect[1]]+Var[Connect[3]])+Var[Connect[0]]);
		break;
	case 3:
		Var[Connect[1]]=InvSumWeightWall*(InvSqrt2*(Var[Connect[2]]+Var[Connect[4]])+Var[Connect[0]]);
		break;
	case 4:
		Var[Connect[2]]=InvSumWeightWall*(InvSqrt2*(Var[Connect[1]]+Var[Connect[3]])+Var[Connect[0]]);
		break;
	}
}
void ExtrapolationWeightDistance::ExtrapolationCornerConcave (double *Var, int * Connect, int & normal){
	switch(normal)
	{
	case 5:
		Var[Connect[3]]=InvSumWeightCornerWall*(InvSqrt2*Var[Connect[2]]+Var[Connect[0]]);
		Var[Connect[4]]=InvSumWeightCornerWall*(InvSqrt2*Var[Connect[1]]+Var[Connect[0]]);
		Var[Connect[7]]=InvSumWeightCornerConcave*(InvSqrt2_5*(Var[Connect[3]]+Var[Connect[4]])+Var[Connect[0]]);
		break;
	case 6:
		Var[Connect[1]]=InvSumWeightCornerWall*(InvSqrt2*Var[Connect[2]]+Var[Connect[0]]);
		Var[Connect[4]]=InvSumWeightCornerWall*(InvSqrt2*Var[Connect[3]]+Var[Connect[0]]);
		Var[Connect[8]]=InvSumWeightCornerConcave*(InvSqrt2_5*(Var[Connect[3]]+Var[Connect[2]])+Var[Connect[0]]);
		break;
	case 7:
		Var[Connect[1]]=InvSumWeightCornerWall*(InvSqrt2*Var[Connect[4]]+Var[Connect[0]]);
		Var[Connect[2]]=InvSumWeightCornerWall*(InvSqrt2*Var[Connect[3]]+Var[Connect[0]]);
		Var[Connect[5]]=InvSumWeightCornerConcave*(InvSqrt2_5*(Var[Connect[3]]+Var[Connect[4]])+Var[Connect[0]]);
		break;
	case 8:
		Var[Connect[3]]=InvSumWeightCornerWall*(InvSqrt2*Var[Connect[4]]+Var[Connect[0]]);
		Var[Connect[2]]=InvSumWeightCornerWall*(InvSqrt2*Var[Connect[1]]+Var[Connect[0]]);
		Var[Connect[6]]=InvSumWeightCornerConcave*(InvSqrt2_5*(Var[Connect[1]]+Var[Connect[4]])+Var[Connect[0]]);

		break;
	}
}
void ExtrapolationWeightDistance::ExtrapolationCornerConvex (double *Var, int * Connect, int & normal){
	switch(normal)
	{
	case 5:
		Var[Connect[7]]=InvSumWeightCornerConvex*(Var[Connect[3]]+Var[Connect[4]]+InvSqrt2*Var[Connect[0]]);
		break;
	case 6:
		Var[Connect[8]]=InvSumWeightCornerConvex*(Var[Connect[1]]+Var[Connect[4]]+InvSqrt2*Var[Connect[0]]);
		break;
	case 7:
		Var[Connect[5]]=InvSumWeightCornerConvex*(Var[Connect[1]]+Var[Connect[2]]+InvSqrt2*Var[Connect[0]]);

		break;
	case 8:
		Var[Connect[6]]=InvSumWeightCornerConvex*(Var[Connect[3]]+Var[Connect[2]]+InvSqrt2*Var[Connect[0]]);
		break;
	}

}
