/*
 * ============================================================================
 * D2Q9Bc.cpp
 *
 *  Created on: 29 Aug 2016
 *      Author: Thomas Burel
 *      Version     :
 *      Copyright   : Your copyright notice
 *      Description : LBM Code Developing by JWFL at University of Strathclyde
 * ============================================================================
 */


#include "D2Q9Bc.h"

D2Q9Bc::D2Q9Bc() {

}
void D2Q9Bc::InitD2Q9Bc(Dictionary *PtrDic_, Parameters *Param, double ** &Ei){
//Initialise Variables shared between boundary conditions
	PtrD2Q9BcDic=	PtrDic_;
	//pointers for generic boundary conditions
	bool Var_found;
	PtrD2Q9BcDic->Get_PtrVar("Density",PtrD2Q9BcRho,Var_found);
	PtrD2Q9BcDic->Get_PtrVar("VelocityX",PtrD2Q9BcU,Var_found);
	PtrD2Q9BcDic->Get_PtrVar("VelocityY",PtrD2Q9BcV,Var_found);

//Initialise the Boundaries conditions classes
	SetBcObjects(Param,Ei);
	Set_GlobalCorner(Param,this);
	Set_SpecialWall(Param,this);
}
D2Q9Bc::~D2Q9Bc() {
	// TODO Auto-generated destructor stub
}

