/*
 * PatchBcTailor.cpp
 *
 *  Created on: 28 Jul 2016
 *      Author: thomas
 */

#include "PatchBcDEF.h"


PatchBcDEF::PatchBcDEF(){
	Type=SolverEnum::Symmetry;
	extrapolationAlpha=false;
}

PatchBcDEF::PatchBcDEF(std::string PatchName_){
	PatchName=PatchName_;
	Type=SolverEnum::Symmetry;
	extrapolationAlpha=false;
}
PatchBcDEF::~PatchBcDEF(){

}





