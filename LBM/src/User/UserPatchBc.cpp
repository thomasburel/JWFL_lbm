/*
 * ============================================================================
 * UserPatchBc.cpp
 *
 *  Created on: 7 Feb 2017
 *      Author: Thomas Burel
 *      Version     :
 *      Copyright   : Your copyright notice
 *      Description : LBM Code Developing by JWFL at University of Strathclyde
 * ============================================================================
 */

#include "UserPatchBc.h"

UserPatchBc::UserPatchBc() {
	// TODO Auto-generated constructor stub

}

UserPatchBc::~UserPatchBc() {
	// TODO Auto-generated destructor stub
}

void UserPatchBc::IntitialiseUserPatchBc(Parameters& PtrParameters,int &NumberOfPatchBc){
	//Set Number of Patch needed (User)
	NumberOfPatchBc=3;
	//Define each patch
		//Patch 0
		//AddPressurePatch("Inlet",HeZouP, FixP);
		AddVelocityPatch("Inlet",HeZouV, FixV);
		//Patch 1
		AddPressurePatch("Outlet",HeZouP, FixP);
		//Patch 2
		AddSymmetryPatch("Symmetry",OnNode);
		//Patch 3
	//	AddWallPatch("Wall",BounceBack);

}
void UserPatchBc::SetUserPatchBc(Parameters& PtrParameters,int elem, int nodenumber, double* pos, int &IdPatchBc){
	if(pos[0]==0)
		IdPatchBc=0;
	if(pos[0]==PtrParameters.Get_Nx())
		IdPatchBc=1;
	if(pos[1]==PtrParameters.Get_Ny())
		IdPatchBc=2;
	if(pos[1]==0)
		IdPatchBc=2;
}


void UserPatchBc::AddPressurePatch(string PatchBcNames,PressureModel PModel_, PressureType PType_){
	PressurePatch.push_back(PressurePatchBc(PatchBcNames,PModel_,PType_));
	PatchIdInType.push_back(PressurePatch.size()-1);
	PatchTypeInType.push_back(SolverEnum::Pressure);
}
void UserPatchBc::AddVelocityPatch(std::string PatchName_,VelocityModel VModel_, VelocityType VType_){
	VelocityPatch.push_back(VelocityPatchBc(PatchName_,VModel_,VType_));
	PatchIdInType.push_back(VelocityPatch.size()-1);
	PatchTypeInType.push_back(SolverEnum::Velocity);
}
void UserPatchBc::AddSymmetryPatch(std::string PatchName_,SymmetryType SType_){
	SymmetryPatch.push_back(SymmetryPatchBc(PatchName_,SType_));
	PatchIdInType.push_back(SymmetryPatch.size()-1);
	PatchTypeInType.push_back(SolverEnum::Symmetry);
}
void UserPatchBc::AddPeriodicPatch(std::string PatchName_,PeriodicType PType_){
	PeriodicPatch.push_back(PeriodicPatchBc(PatchName_,PType_));
	PatchIdInType.push_back(PeriodicPatch.size()-1);
	PatchTypeInType.push_back(SolverEnum::Periodic);
}
void UserPatchBc::AddWallPatch(std::string PatchName_,WallType WType_){
	WallPatch.push_back(WallPatchBc(PatchName_,WType_));
	PatchIdInType.push_back(WallPatch.size()-1);
	PatchTypeInType.push_back(SolverEnum::Wall);
}
