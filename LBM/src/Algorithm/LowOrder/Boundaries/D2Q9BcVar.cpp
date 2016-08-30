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
// He Zhou parameters
	q23=2.0/3.0;
	q16=1.0/6.0;
	q13=1.0/3.0;

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

D2Q9BcVar::~D2Q9BcVar() {
	// TODO Auto-generated destructor stub
}

