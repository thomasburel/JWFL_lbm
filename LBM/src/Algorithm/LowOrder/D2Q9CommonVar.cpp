/*
 * ============================================================================
 * D2Q9CommonVar.cpp
 *
 *  Created on: 6 Sep 2016
 *      Author: Thomas Burel
 *      Version     :
 *      Copyright   : Your copyright notice
 *      Description : LBM Code Developing by JWFL at University of Strathclyde
 * ============================================================================
 */

#include "D2Q9CommonVar.h"

D2Q9CommonVar::D2Q9CommonVar() {
	Ei=0;
	omega=0;
	Opposite=0;
	Ei=new double* [9]; for(int i=0;i<9;i++) Ei[i]=new double [2];
	omega=new double [9];
	Opposite=new short int [9];
	//Distribution velocity
	Ei[0][0]= 0.0;
	Ei[0][1]= 0.0;
	Ei[1][0]= 1.0;
	Ei[1][1]= 0.0;
	Ei[2][0]= 0.0;
	Ei[2][1]= 1.0;
	Ei[3][0]= -1.0;
	Ei[3][1]= 0.0;
	Ei[4][0]= 0.0;
	Ei[4][1]= -1.0;
	Ei[5][0]= 1.0;
	Ei[5][1]= 1.0;
	Ei[6][0]= -1.0;
	Ei[6][1]= 1.0;
	Ei[7][0]= -1.0;
	Ei[7][1]= -1.0;
	Ei[8][0]= 1.0;
	Ei[8][1]= -1.0;
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
	// Opposite direction
	Opposite[0]=0;
	Opposite[1]=3;
	Opposite[2]=4;
	Opposite[3]=1;
	Opposite[4]=2;
	Opposite[5]=7;
	Opposite[6]=8;
	Opposite[7]=5;
	Opposite[8]=6;

}

D2Q9CommonVar::~D2Q9CommonVar() {
	// TODO Auto-generated destructor stub
}

