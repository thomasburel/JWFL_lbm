/*
 * PatchBcDEF.h
 *
 *  Created on: 31 Oct 2016
 *      Author: Thomas Burel
 */

#ifndef SRC_MESH_SINGLEBLOCK_PATCH_PATCHBCDEF_H_
#define SRC_MESH_SINGLEBLOCK_PATCH_PATCHBCDEF_H_

#include "../../../Core/Parameters.h"
#include "../NodeArrays.h"

///Abstract class for common function for all kind of extrapolation
class PatchBcDEF {
public:
	PatchBcDEF();
	PatchBcDEF(std::string PatchName_);
	virtual ~PatchBcDEF();
	void Set_PatchName(std::string PatchName_){PatchName=PatchName_;};
	std::string Get_PatchName(){return PatchName;};
	void Set_PatchType(SolverEnum::PatchType PatchType_){Type=PatchType_;};
	SolverEnum::PatchType Get_PatchType(){return Type;};
	void addNodeIndex(int NodeIndex){IdxNode.push_back(NodeIndex);};
	void Set_NodeIndex(std::vector<int>& IdxNode_){IdxNode=IdxNode_;};
	std::vector<int>& Get_NodeIndex(){return IdxNode;};
	void Set_NodeIndexSpecialWalls(std::vector<int>& IdxNode_){NodeIndexSpecialWalls=IdxNode_;};
	std::vector<int>& Get_NodeIndexSpecialWalls(){return NodeIndexSpecialWalls;};
	void Set_NodeIndexGlobalCorner(std::vector<int>& IdxNode_){NodeIndexGlobalCorner=IdxNode_;};
	std::vector<int>& Get_NodeIndexGlobalCorner(){return NodeIndexGlobalCorner;};
	void Set_NodeIndexByType(std::vector<int>& IdxNode_){NodeIndexByType=IdxNode_;};
	std::vector<int>& Get_NodeIndexByType(){return NodeIndexByType;};
	void Set_NodeIndexByTypeSpecialWalls(std::vector<int>& IdxNode_){NodeIndexByTypeSpecialWalls=IdxNode_;};
	std::vector<int>& Get_NodeIndexByTypeSpecialWalls(){return NodeIndexByTypeSpecialWalls;};
	void Set_NodeIndexByTypeGlobalCorner(std::vector<int>& IdxNode_){NodeIndexByTypeGlobalCorner=IdxNode_;};
	std::vector<int>& Get_NodeIndexByTypeGlobalCorner(){return NodeIndexByTypeGlobalCorner;};
	void Set_extrapolationAlpha(bool extrapol){extrapolationAlpha=extrapol;};
	bool Get_extrapolationAlpha(){return extrapolationAlpha;};
	void Set_extrapolationNormal(bool extrapol){extrapolationNormal=extrapol;};
	bool Get_extrapolationNormal(){return extrapolationNormal;};
	void Set_Inlet(bool inlet){Inlet=inlet;if(Inlet) Outlet=false;};
	bool Get_Inlet(){return Inlet;};
	void Set_Outlet(bool outlet){Outlet=outlet;if(Outlet) Inlet=false;};
	bool Get_Outlet(){return Outlet;};


protected:
	std::vector<int> IdxNode,NodeIndexSpecialWalls,NodeIndexGlobalCorner, NodeIndexByType,NodeIndexByTypeSpecialWalls,NodeIndexByTypeGlobalCorner;
	std::string PatchName;
	SolverEnum::PatchType Type;
	bool extrapolationAlpha,extrapolationNormal,Inlet,Outlet;
	int orientation;

};
class PressurePatchBc: public PatchBcDEF {
public:
	PressurePatchBc(std::string PatchName_,PressureModel PModel_, PressureType PType_)
		{PatchName=PatchName_;Type=SolverEnum::Pressure;PModel=PModel_;PType=PType_;};
	virtual ~PressurePatchBc(){};
	PressureModel& Get_PressureModel(){return PModel;};
	PressureType& Get_PressureType(){return PType;};
private:
	PressureModel PModel;
	PressureType PType;
};
class VelocityPatchBc: public PatchBcDEF {
public:
	VelocityPatchBc();
	VelocityPatchBc(std::string PatchName_,VelocityModel VModel_, VelocityType VType_)
		{PatchName=PatchName_;Type=SolverEnum::Velocity;VModel=VModel_;VType=VType_;};
	virtual ~VelocityPatchBc(){};
	VelocityModel& Get_VelocityModel(){return VModel;};
	VelocityType& Get_VelocityType(){return VType;};
//PatchBc
private:
	VelocityModel VModel;
	VelocityType VType;
};
class SymmetryPatchBc: public PatchBcDEF {
public:
	SymmetryPatchBc(std::string PatchName_,SymmetryType SType_)
		{PatchName=PatchName_;Type=SolverEnum::Symmetry;SType=SType_;};
	virtual ~SymmetryPatchBc(){};
	SymmetryType& Get_SymmetryType(){return SType;};
private:
	SymmetryType SType;
};
class PeriodicPatchBc: public PatchBcDEF {
public:
	PeriodicPatchBc(std::string PatchName_,PeriodicType PType_)
		{PatchName=PatchName_;Type=SolverEnum::Periodic;PType=PType_;};
	virtual ~PeriodicPatchBc(){};
	PeriodicType& Get_PeriodicType(){return PType;};

private:
	PeriodicType PType;
};
class WallPatchBc: public PatchBcDEF {
public:
	WallPatchBc(std::string PatchName_,WallType WType_)
		{PatchName=PatchName_;Type=SolverEnum::Wall;WType=WType_;};
	virtual ~WallPatchBc(){};
	WallType& Get_WallType(){return WType;};
private:
	WallType WType;
};

#endif /* SRC_MESH_SINGLEBLOCK_PATCH_PATCHBCDEF_H_ */
