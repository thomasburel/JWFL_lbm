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
NoExtrapolation::NoExtrapolation(){}
NoExtrapolation::NoExtrapolation(int dimension, int nb_vel){}
NoExtrapolation::~NoExtrapolation(){}

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
///First order Tailor extrapolation
void ExtrapolationTailor::ExtrapolationOnWall (double *Var, int * Connect, int & normal){

	switch(normal)
	{
	case 1:
		Var[Connect[0]]=Var[Connect[1]];
		break;
	case 2:
		Var[Connect[0]]=Var[Connect[2]];
		break;
	case 3:
		Var[Connect[0]]=Var[Connect[3]];
		break;
	case 4:
		Var[Connect[0]]=Var[Connect[4]];
		break;
	}
}
///First order Tailor extrapolation
void ExtrapolationTailor::ExtrapolationOnCornerConcave (double *Var, int * Connect, int & normal){
	switch(normal)
	{
	case 5:
		Var[Connect[0]]=Var[Connect[5]];
		break;
	case 6:
		Var[Connect[0]]=Var[Connect[6]];
		break;
	case 7:
		Var[Connect[0]]=Var[Connect[7]];
		break;
	case 8:
		Var[Connect[0]]=Var[Connect[8]];
		break;
	}
}
///First order Tailor extrapolation
void ExtrapolationTailor::ExtrapolationOnCornerConvex(double *Var, int * Connect, int & normal){
	switch(normal)
	{
	case 5:
		Var[Connect[0]]=Var[Connect[5]];
		break;
	case 6:
		Var[Connect[0]]=Var[Connect[6]];
		break;
	case 7:
		Var[Connect[0]]=Var[Connect[7]];
		break;
	case 8:
		Var[Connect[0]]=Var[Connect[8]];
		break;
	}
}
///Second order Tailor extrapolation
void ExtrapolationTailor::ExtrapolationWallToSolid (double *Var, int * Connect, int & normal){

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
void ExtrapolationTailor::ExtrapolationCornerConcaveToSolid (double *Var, int * Connect, int & normal){
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
void ExtrapolationTailor::ExtrapolationCornerConvexToSolid (double *Var, int * Connect, int & normal){
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
	InvSumWeightWallToSolid=1.0/(1.0+2.0*InvSqrt2);
	InvSumWeightCornerConvexToSolid=1.0/(2.0+1.0*InvSqrt2);
	InvSumWeightCornerConcaveToSolid=1.0/(1.0+2.0*InvSqrt2_5);//0.4=2/5
	InvSumWeightCornerWallToSolid=1.0/(1.0+1.0*InvSqrt2);
	InvSumWeightOnCornerConvex=1.0/(4.0+3.0*InvSqrt2);
	InvSumWeightOnCornerConcave=1.0/(2.0+1.0*InvSqrt2);
	InvSumWeightOnWall=1.0/(1.0+2.0*InvSqrt2);
}
ExtrapolationWeightDistance::ExtrapolationWeightDistance(int dimension_,int nb_vel){
	dimension=dimension_;
	nb_Vel=nb_vel;
	InvSqrt2=std::sqrt(0.5);
	InvSqrt2_5=std::sqrt(0.4);//0.4=2/5
	InvSumWeightWallToSolid=1.0/(1.0+2.0*InvSqrt2);
	InvSumWeightCornerConvexToSolid=1.0/(2.0+1.0*InvSqrt2);
	InvSumWeightCornerConcaveToSolid=1.0/(1.0+2.0*InvSqrt2_5);
	InvSumWeightCornerWallToSolid=1.0/(1.0+1.0*InvSqrt2);
	InvSumWeightOnCornerConvex=1.0/(4.0+3.0*InvSqrt2);
	InvSumWeightOnCornerConcave=1.0/(2.0+1.0*InvSqrt2);
	InvSumWeightOnWall=1.0/(1.0+2.0*InvSqrt2);
}
ExtrapolationWeightDistance::~ExtrapolationWeightDistance(){

}
void ExtrapolationWeightDistance::ExtrapolationOnWall(double *Var, int * Connect, int & normal){
	switch(normal)
	{
	case 1:
		Var[Connect[0]]=InvSumWeightOnWall*(InvSqrt2*(Var[Connect[5]]+Var[Connect[8]])+Var[Connect[1]]);
		break;
	case 2:
		Var[Connect[0]]=InvSumWeightOnWall*(InvSqrt2*(Var[Connect[5]]+Var[Connect[6]])+Var[Connect[2]]);
		break;
	case 3:
		Var[Connect[0]]=InvSumWeightOnWall*(InvSqrt2*(Var[Connect[6]]+Var[Connect[7]])+Var[Connect[3]]);
		break;
	case 4:
		Var[Connect[0]]=InvSumWeightOnWall*(InvSqrt2*(Var[Connect[7]]+Var[Connect[8]])+Var[Connect[4]]);
		break;
	}
}
void ExtrapolationWeightDistance::ExtrapolationOnCornerConcave (double *Var, int * Connect, int & normal){
	switch(normal)
	{
	case 5:
		Var[Connect[0]]=InvSumWeightOnCornerConcave*(Var[Connect[1]]+Var[Connect[2]]+InvSqrt2*Var[Connect[5]]);
		break;
	case 6:
		Var[Connect[0]]=InvSumWeightOnCornerConcave*(Var[Connect[3]]+Var[Connect[2]]+InvSqrt2*Var[Connect[6]]);
		break;
	case 7:
		Var[Connect[0]]=InvSumWeightOnCornerConcave*(Var[Connect[3]]+Var[Connect[4]]+InvSqrt2*Var[Connect[7]]);
		break;
	case 8:
		Var[Connect[0]]=InvSumWeightOnCornerConcave*(Var[Connect[1]]+Var[Connect[4]]+InvSqrt2*Var[Connect[8]]);
		break;
	}
}
void ExtrapolationWeightDistance::ExtrapolationOnCornerConvex (double *Var, int * Connect, int & normal){
	switch(normal)
	{
	case 5:
		Var[Connect[0]]=InvSumWeightOnCornerConvex*(Var[Connect[1]]+Var[Connect[2]]+Var[Connect[3]]+Var[Connect[4]]+InvSqrt2*(Var[Connect[5]]+Var[Connect[6]]+Var[Connect[8]]));
		break;
	case 6:
		Var[Connect[0]]=InvSumWeightOnCornerConvex*(Var[Connect[1]]+Var[Connect[2]]+Var[Connect[3]]+Var[Connect[4]]+InvSqrt2*(Var[Connect[5]]+Var[Connect[6]]+Var[Connect[7]]));
		break;
	case 7:
		Var[Connect[0]]=InvSumWeightOnCornerConvex*(Var[Connect[1]]+Var[Connect[2]]+Var[Connect[3]]+Var[Connect[4]]+InvSqrt2*(Var[Connect[6]]+Var[Connect[7]]+Var[Connect[8]]));
		break;
	case 8:
		Var[Connect[0]]=InvSumWeightOnCornerConvex*(Var[Connect[1]]+Var[Connect[2]]+Var[Connect[3]]+Var[Connect[4]]+InvSqrt2*(Var[Connect[5]]+Var[Connect[7]]+Var[Connect[8]]));
		break;
	}

}
void ExtrapolationWeightDistance::ExtrapolationWallToSolid (double *Var, int * Connect, int & normal){
	switch(normal)
	{
	case 1:
		Var[Connect[3]]=InvSumWeightWallToSolid*(InvSqrt2*(Var[Connect[2]]+Var[Connect[4]])+Var[Connect[0]]);
		break;
	case 2:
		Var[Connect[4]]=InvSumWeightWallToSolid*(InvSqrt2*(Var[Connect[1]]+Var[Connect[3]])+Var[Connect[0]]);
		break;
	case 3:
		Var[Connect[1]]=InvSumWeightWallToSolid*(InvSqrt2*(Var[Connect[2]]+Var[Connect[4]])+Var[Connect[0]]);
		break;
	case 4:
		Var[Connect[2]]=InvSumWeightWallToSolid*(InvSqrt2*(Var[Connect[1]]+Var[Connect[3]])+Var[Connect[0]]);
		break;
	}
}
void ExtrapolationWeightDistance::ExtrapolationCornerConcaveToSolid (double *Var, int * Connect, int & normal){
	switch(normal)
	{
	case 5:
		Var[Connect[3]]=InvSumWeightCornerWallToSolid*(InvSqrt2*Var[Connect[2]]+Var[Connect[0]]);
		Var[Connect[4]]=InvSumWeightCornerWallToSolid*(InvSqrt2*Var[Connect[1]]+Var[Connect[0]]);
		Var[Connect[7]]=InvSumWeightCornerConcaveToSolid*(InvSqrt2_5*(Var[Connect[1]]+Var[Connect[2]])+Var[Connect[0]]);
		break;
	case 6:
		Var[Connect[1]]=InvSumWeightCornerWallToSolid*(InvSqrt2*Var[Connect[2]]+Var[Connect[0]]);
		Var[Connect[4]]=InvSumWeightCornerWallToSolid*(InvSqrt2*Var[Connect[3]]+Var[Connect[0]]);
		Var[Connect[8]]=InvSumWeightCornerConcaveToSolid*(InvSqrt2_5*(Var[Connect[3]]+Var[Connect[2]])+Var[Connect[0]]);
		break;
	case 7:
		Var[Connect[1]]=InvSumWeightCornerWallToSolid*(InvSqrt2*Var[Connect[4]]+Var[Connect[0]]);
		Var[Connect[2]]=InvSumWeightCornerWallToSolid*(InvSqrt2*Var[Connect[3]]+Var[Connect[0]]);
		Var[Connect[5]]=InvSumWeightCornerConcaveToSolid*(InvSqrt2_5*(Var[Connect[3]]+Var[Connect[4]])+Var[Connect[0]]);
		break;
	case 8:
		Var[Connect[3]]=InvSumWeightCornerWallToSolid*(InvSqrt2*Var[Connect[4]]+Var[Connect[0]]);
		Var[Connect[2]]=InvSumWeightCornerWallToSolid*(InvSqrt2*Var[Connect[1]]+Var[Connect[0]]);
		Var[Connect[6]]=InvSumWeightCornerConcaveToSolid*(InvSqrt2_5*(Var[Connect[1]]+Var[Connect[4]])+Var[Connect[0]]);

		break;
	}
}
void ExtrapolationWeightDistance::ExtrapolationCornerConvexToSolid (double *Var, int * Connect, int & normal){
	switch(normal)
	{
	case 5:
		Var[Connect[7]]=InvSumWeightCornerConvexToSolid*(Var[Connect[3]]+Var[Connect[4]]+InvSqrt2*Var[Connect[0]]);
		break;
	case 6:
		Var[Connect[8]]=InvSumWeightCornerConvexToSolid*(Var[Connect[1]]+Var[Connect[4]]+InvSqrt2*Var[Connect[0]]);
		break;
	case 7:
		Var[Connect[5]]=InvSumWeightCornerConvexToSolid*(Var[Connect[1]]+Var[Connect[2]]+InvSqrt2*Var[Connect[0]]);
		break;
	case 8:
		Var[Connect[6]]=InvSumWeightCornerConvexToSolid*(Var[Connect[3]]+Var[Connect[2]]+InvSqrt2*Var[Connect[0]]);
		break;
	}

}
