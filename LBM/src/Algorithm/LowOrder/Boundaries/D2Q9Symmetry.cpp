/*
 * ============================================================================
 * D2Q9Symmetry.cpp
 *
 *  Created on: 2 Sep 2016
 *      Author: Thomas Burel
 *      Version     :
 *      Copyright   : Your copyright notice
 *      Description : LBM Code Developing by JWFL at University of Strathclyde
 * ============================================================================
 */

#include "D2Q9Symmetry.h"

D2Q9Symmetry::D2Q9Symmetry() {
	PtrSymmetryMethod=0;
	LocalForce=0;
}

D2Q9Symmetry::~D2Q9Symmetry() {
	// TODO Auto-generated destructor stub
}
void D2Q9Symmetry::Set_Symmetry(Dictionary *PtrDic,NodeArrays2D* NodeArrays, Parameters *Param, double ** &Ei){
	SetSymmetry(Param->Get_SymmetryType());
	EiBc=Ei;
}
void D2Q9Symmetry::SetSymmetry(SymmetryType SymmetryType_){
	SetSymmetryType(SymmetryType_);
}
void D2Q9Symmetry::SetSymmetryType(SymmetryType SymmetryType_){
	PtrSymmetryMethod=&D2Q9Symmetry::ApplySymmetryOnNode;
}
void D2Q9Symmetry::ApplySymmetry(int const &BcNormal,int const *Connect, double const &Rho_def, double const *UDef, DistriFunct* f_in, double *Rho, double *U, double *V){
	(this->*PtrSymmetryMethod)(BcNormal,Connect,Rho_def,UDef,LocalForce,f_in);
}
//Symmetry treatment on node
void D2Q9Symmetry::ApplySymmetryOnNode(int const &BcNormal,int const *Connect, double const &Rho_def, double const *UDef, double *LocalForce, DistriFunct* f_in){
	switch(BcNormal)
			{
			case 2:
				f_in->f[2][Connect[0]]=f_in->f[4][Connect[0]];
				f_in->f[5][Connect[0]]=f_in->f[8][Connect[0]];
				f_in->f[6][Connect[0]]=f_in->f[7][Connect[0]];
				break;
			case 4:
				f_in->f[4][Connect[0]]=f_in->f[2][Connect[0]];
				f_in->f[7][Connect[0]]=f_in->f[6][Connect[0]];
				f_in->f[8][Connect[0]]=f_in->f[5][Connect[0]];
				break;
			case 1:
				f_in->f[1][Connect[0]]=f_in->f[3][Connect[0]];
				f_in->f[5][Connect[0]]=f_in->f[6][Connect[0]];
				f_in->f[8][Connect[0]]=f_in->f[7][Connect[0]];
				break;
			case 3:
				f_in->f[3][Connect[0]]=f_in->f[1][Connect[0]];
				f_in->f[6][Connect[0]]=f_in->f[5][Connect[0]];
				f_in->f[7][Connect[0]]=f_in->f[8][Connect[0]];
				break;
			default:
				std::cerr<<"Direction symmetry not found"<<std::endl;
				break;
			}
}
