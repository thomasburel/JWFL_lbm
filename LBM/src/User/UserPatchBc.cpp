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
	NumberOfPatchBc=4;
	//Define each patch
		//Patch 0
	//	AddWallPatch("Wall",BounceBack);

		AddPressurePatch("Inlet",HeZouP, FixP);
		PressurePatch.back().Set_Inlet(true);
		//AddVelocityPatch("Inlet",HeZouV, FixV);
		//Patch 1
		AddPressurePatch("Outlet",HeZouP, FixP);
		//AddVelocityPatch("Outlet",HeZouV, zeroVGrad1st); //FixV,zeroVGrad1st
		PressurePatch.back().Set_extrapolationAlpha(true);
		PressurePatch.back().Set_Outlet(true);
		//PressurePatch.back().Set_extrapolationNormal(true);
		//AddPressurePatch("Outlet",HeZouP, FixP);
		//Patch 2
		//AddWallPatch("Wall",BounceBack);
		AddSymmetryPatch("Symmetry",OnNode);
		//Patch 3
		//AddWallPatch("Wall",BounceBack);
		AddSymmetryPatch("Symmetry",OnNode);

}
void UserPatchBc::SetUserPatchBc(Parameters& PtrParameters,int elem, int nodenumber, double* pos, int &IdPatchBc1, int &IdPatchBc2){
//left side
	if(pos[0]==0)
		// left bottom corner
		if (pos[1]==0)
		{
			IdPatchBc1=0;
			IdPatchBc2=2;
		}
	// Left top corner
		else if(pos[1]==PtrParameters.Get_Ny())
		{
			IdPatchBc1=0;
			IdPatchBc2=3;
		}
	//left side
		else
		{
			IdPatchBc1=0;
		}
//right side
	else if(pos[0]==PtrParameters.Get_Nx())
		// right bottom corner
		if (pos[1]==0)
		{
			IdPatchBc1=1;
			IdPatchBc2=2;
		}
	// right top corner
		else if(pos[1]==PtrParameters.Get_Ny())
		{
			IdPatchBc1=1;
			IdPatchBc2=3;
		}
	//right side
		else
		{
			IdPatchBc1=1;
		}
	else
// Bottom
		if (pos[1]==0)
		{
			IdPatchBc1=2;
		}
// Top
		else if(pos[1]==PtrParameters.Get_Ny())
		{
			IdPatchBc1=3;
		}



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
