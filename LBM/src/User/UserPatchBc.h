/*
 * ============================================================================
 * UserPatchBc.h
 *
 *  Created on: 7 Feb 2017
 *      Author: Thomas Burel
 *      Version     :
 *      Copyright   : Your copyright notice
 *      Description : LBM Code Developing by JWFL at University of Strathclyde
 * ============================================================================
 */

#ifndef SRC_USER_USERPATCHBC_H_
#define SRC_USER_USERPATCHBC_H_
#include "../Mesh/SingleBlock/Patch/PatchBcDEF.h"
#include "../Core/Parameters.h"

class UserPatchBc {


// Do not modify under this line
public:
	UserPatchBc();
	virtual ~UserPatchBc();
	void IntitialiseUserPatchBc(Parameters& PtrParameters,int &NumberOfPatchBc);
	void SetUserPatchBc(Parameters& PtrParameters,int elem, int nodenumber, double* pos, int &IdPatchBc1, int &IdPatchBc2);

protected:
	void AddPressurePatch(string PatchBcNames,PressureModel PModel_, PressureType PType_);
	void AddVelocityPatch(std::string PatchName_,VelocityModel VModel_, VelocityType VType_);
	void AddSymmetryPatch(std::string PatchName_,SymmetryType SType_);
	void AddPeriodicPatch(std::string PatchName_,PeriodicType PType_);
	void AddWallPatch(std::string PatchName_,WallType WType_);

protected:
	std::vector<PressurePatchBc> PressurePatch;
	std::vector<VelocityPatchBc> VelocityPatch;
	std::vector<SymmetryPatchBc> SymmetryPatch;
	std::vector<PeriodicPatchBc> PeriodicPatch;
	std::vector<WallPatchBc> WallPatch;
	std::vector<int> PatchIdInType;
	std::vector<SolverEnum::PatchType> PatchTypeInType;
};

#endif /* SRC_USER_USERPATCHBC_H_ */
