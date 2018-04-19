/*
 * ============================================================================
 * D2Q9BcVar.cpp
 *
 *  Created on: 29 Aug 2016
 *      Author: Thomas Burel
 *      Version     :
 *      Copyright   : Your copyright notice
 *      Description : LBM Code Developing by JWFL at University of Strathclyde
 * ============================================================================
 */

#include "D2Q9BcVar.h"

D2Q9BcVar::D2Q9BcVar() {
	ru1=0;
	ru2=0;
	rhodiff=0;
// He Zhou parameters
	q23=2.0/3.0;
	q16=1.0/6.0;
	q13=1.0/3.0;

	//Distribution velocity
	EiBc=0;
	invCs2=3.0;
/*	EiBc[0][0]= 0.0;
	EiBc[0][1]= 0.0;
	EiBc[1][0]= 1.0;
	EiBc[1][1]= 0.0;
	EiBc[2][0]= 0.0;
	EiBc[2][1]= 1.0;
	EiBc[3][0]= -1.0;
	EiBc[3][1]= 0.0;
	EiBc[4][0]= 0.0;
	EiBc[4][1]= -1.0;
	EiBc[5][0]= 1.0;
	EiBc[5][1]= 1.0;
	EiBc[6][0]= -1.0;
	EiBc[6][1]= 1.0;
	EiBc[7][0]= -1.0;
	EiBc[7][1]= -1.0;
	EiBc[8][0]= 1.0;*/
// OppositeBc direction
	OppositeBc[0]=0;
	OppositeBc[1]=3;
	OppositeBc[2]=4;
	OppositeBc[3]=1;
	OppositeBc[4]=2;
	OppositeBc[5]=7;
	OppositeBc[6]=8;
	OppositeBc[7]=5;
	OppositeBc[8]=6;

	//Weigh depending of the direction
	omegaBc[0]=4.0/9.0;
	omegaBc[1]=1.0/9.0;
	omegaBc[2]=1.0/9.0;
	omegaBc[3]=1.0/9.0;
	omegaBc[4]=1.0/9.0;
	omegaBc[5]=1.0/36.0;
	omegaBc[6]=1.0/36.0;
	omegaBc[7]=1.0/36.0;
	omegaBc[8]=1.0/36.0;

	//Diffuse variables
	SumWeightS=omegaBc[4]+omegaBc[7]+omegaBc[8];
	SumWeightE=omegaBc[1]+omegaBc[5]+omegaBc[8];
	SumWeightN=omegaBc[2]+omegaBc[5]+omegaBc[6];
	SumWeightW=omegaBc[3]+omegaBc[6]+omegaBc[7];
	SumWeightConcaveSE=omegaBc[1]+omegaBc[4]+omegaBc[8];
	SumWeightConcaveNE=omegaBc[1]+omegaBc[2]+omegaBc[5];
	SumWeightConcaveNW=omegaBc[2]+omegaBc[3]+omegaBc[6];
	SumWeightConcaveSW=omegaBc[3]+omegaBc[4]+omegaBc[7];
	SumWeightConvexSE=omegaBc[2]+omegaBc[3]+omegaBc[6];
	SumWeightConvexNE=omegaBc[3]+omegaBc[4]+omegaBc[7];
	SumWeightConvexNW=omegaBc[1]+omegaBc[4]+omegaBc[8];
	SumWeightConvexSW=omegaBc[1]+omegaBc[2]+omegaBc[5];

	for(int i=0;i<2;i++)
		DiffusekWallConnect[0][i]=0;

	for(int i=1;i<5;i++)
		BounceBackWallConnect[i]=new int [2];
	for(int i=5;i<9;i++)
		BounceBackWallConnect[i]=new int [4];

	BounceBackWallConnect[0]=0;
	BounceBackWallConnect[1][0]=5;
	BounceBackWallConnect[1][1]=8;
	BounceBackWallConnect[2][0]=5;
	BounceBackWallConnect[2][1]=6;
	BounceBackWallConnect[3][0]=6;
	BounceBackWallConnect[3][1]=7;
	BounceBackWallConnect[4][0]=7;
	BounceBackWallConnect[4][1]=8;

	BounceBackWallConnect[5][0]=1;
	BounceBackWallConnect[5][1]=2;
	BounceBackWallConnect[5][2]=6;
	BounceBackWallConnect[5][3]=8;

	BounceBackWallConnect[6][0]=2;
	BounceBackWallConnect[6][1]=3;
	BounceBackWallConnect[6][2]=5;
	BounceBackWallConnect[6][3]=7;

	BounceBackWallConnect[7][0]=3;
	BounceBackWallConnect[7][1]=4;
	BounceBackWallConnect[7][2]=6;
	BounceBackWallConnect[7][3]=8;

	BounceBackWallConnect[8][0]=1;
	BounceBackWallConnect[8][1]=4;
	BounceBackWallConnect[8][2]=5;
	BounceBackWallConnect[8][3]=7;
}

D2Q9BcVar::~D2Q9BcVar() {
	// TODO Auto-generated destructor stub
	for(int i=1;i<5;i++)
		delete [] BounceBackWallConnect[i];
	for(int i=5;i<9;i++)
		delete [] BounceBackWallConnect[i];
}

