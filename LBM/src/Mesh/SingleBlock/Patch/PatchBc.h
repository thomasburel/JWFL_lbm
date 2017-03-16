/*
 * Interpolation.h
 *
 *  Created on: 31 October 2016
 *      Author: Thomas Burel
 *
 *  Files created to manage Interpolation types.
 *  To add a new Interpolation type, you need to:
 *  	 Add a class in InterpolationDEF.h
 *  	 Defined all virtual methods in InterpolationDEF.cpp
 *  	 Add an enumtype
 *  	 Modified the SelectInterpolationType function
 */

#ifndef SRC_MESH_SINGLEBLOCK_PATCH_PATCHBC_H_
#define SRC_MESH_SINGLEBLOCK_PATCH_PATCHBC_H_

#include "../../../User/UserPatchBc.h"
#include "PatchBcDEF.h"
#include "mpi.h"
class PatchBc : public UserPatchBc {
public:
	PatchBc();
	virtual ~PatchBc();
	void initPatchBc(std::vector<Node2D*> Node, Parameters *PtrParam);
	int Get_NumberOfPatchBc(){return NumberOfPatchBc;};
	std::vector<int> Get_PatchIdInType(){return PatchIdInType;};
	std::vector<SolverEnum::PatchType> Get_PatchTypeInType(){return PatchTypeInType;};
	std::vector<PressurePatchBc>& Get_PressurePatch(){return PressurePatch;};
	std::vector<VelocityPatchBc>& Get_VelocityPatch(){return VelocityPatch;};
	std::vector<SymmetryPatchBc>& Get_SymmetryPatch(){return SymmetryPatch;};
	std::vector<PeriodicPatchBc>& Get_PeriodicPatch(){return  PeriodicPatch;};
	std::vector<WallPatchBc>& Get_WallPatch(){return WallPatch;};
	void Set_NodeIndexByType(SolverEnum::PatchType Type_, int IdIntype, std::vector<int>& IdxNode_);
	void Set_NodeIndexSpecialWalls(SolverEnum::PatchType Type_, int IdIntype, std::vector<int>& IdxNode_);
	void Set_NodeIndexGlobalCorner(SolverEnum::PatchType Type_, int IdIntype, std::vector<int>& IdxNode_);
private:
	void SetNodeIdxForPatchBc(SolverEnum::PatchType Type_, int IdIntype,int PatchId,std::vector<int> &nodeIdx);
	void SelectPatchBc(SolverEnum::PatchType Type_, Parameters *PtrParam, string PatchBcNames, int PatchId);
	void SelectPatchBc(NodeType Type_, Parameters *PtrParam, string PatchBcNames, int PatchId);
	void RemovePatchBc(int PatchId);
	void ReorderingPatches();
private:
	int NumberOfPatchBc;

};
#endif /* SRC_MESH_SINGLEBLOCK_PATCH_PATCHBC_H_ */
