/*
 * PatchBc.cpp
 *
 *  Created on: 31 Oct 2016
 *      Author: Thomas Burel
 */
#include "PatchBc.h"

PatchBc::PatchBc(){

}
void PatchBc::initPatchBc(std::vector<Node2D*> Node, Parameters *PtrParam){
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	//Initialise user patches
	NumberOfPatchBc=0;
	IntitialiseUserPatchBc(*PtrParam,NumberOfPatchBc);
	//Add Default Patches for the boundary of the Cartesian domain
	SelectPatchBc(PtrParam->Get_GlobalBcType(3), PtrParam, string("Left"), NumberOfPatchBc+0);
	SelectPatchBc(PtrParam->Get_GlobalBcType(1), PtrParam, string("Right"), NumberOfPatchBc+1);
	SelectPatchBc(PtrParam->Get_GlobalBcType(0), PtrParam, string("Bottom"), NumberOfPatchBc+2);
	SelectPatchBc(PtrParam->Get_GlobalBcType(2), PtrParam, string("Top"), NumberOfPatchBc+3);
	//Define temporary variables
	double *pos;int IdPatchBc=-1;
	//std::vector<int>  *IdxNode;
	//IdxNode=new std::vector<int>[NumberOfPatchBc+4];
	std::vector<std::vector<int> >  IdxNode;
	for (int i=0;i<NumberOfPatchBc+4;i++)
		IdxNode.push_back(std::vector<int> (0));
	int nx=PtrParam->Get_Nx();
	int ny=PtrParam->Get_Ny();

	pos=new double[2];
	int Bcbottom=0;int BcRight=0;int BcTop=0;int BcLeft=0;

	//Get index for patches
	for(int i=0;i<Node.size();i++)
	{
		pos[0]=Node[i]->get_x();pos[1]=Node[i]->get_y();
		if(Node[i]->get_NodeType()!=Ghost && Node[i]->get_NodeType()!=SolidGhost &&
				(pos[0]==0 || pos[0]==nx || pos[1]==0 || pos[1]==ny))
		{
			IdPatchBc=-1;
			SetUserPatchBc(*PtrParam,0, i, pos, IdPatchBc);
			if(IdPatchBc>=0 && IdPatchBc<NumberOfPatchBc)
			{
				if(IdPatchBc<NumberOfPatchBc)
					IdxNode[IdPatchBc].push_back(i);
			}
			else
			{
				if(pos[0]==0)
				{
					if(pos[1]==0)
					{
						IdxNode[NumberOfPatchBc+0].push_back(i);
						IdxNode[NumberOfPatchBc+2].push_back(i);
					}
					else if(pos[1]==ny)
					{
						IdxNode[NumberOfPatchBc+0].push_back(i);
						IdxNode[NumberOfPatchBc+3].push_back(i);
					}
					else
					{
						IdxNode[NumberOfPatchBc+0].push_back(i);
					}
				}
				else if(pos[0]==nx)
				{
					if(pos[1]==0)
					{
						IdxNode[NumberOfPatchBc+1].push_back(i);
						IdxNode[NumberOfPatchBc+2].push_back(i);
					}
					else if(pos[1]==ny)
					{
						IdxNode[NumberOfPatchBc+1].push_back(i);
						IdxNode[NumberOfPatchBc+3].push_back(i);
					}
					else
					{
						IdxNode[NumberOfPatchBc+1].push_back(i);
					}
				}
				else
				{
					if(pos[1]==0)
					{
						IdxNode[NumberOfPatchBc+2].push_back(i);
					}
					else if(pos[1]==ny)
					{
						IdxNode[NumberOfPatchBc+3].push_back(i);
					}
				}
			}
		}
	}
	//Remove empty patches
	std::vector<int> IdToRemove;
	//mark patches need to be removed
	for(int i=NumberOfPatchBc;i<NumberOfPatchBc+4;i++)
	{
		if(IdxNode[i].empty())
			IdToRemove.push_back(i);
	}
	//remove patches
	for(int i=0;i<IdToRemove.size();i++)
		RemovePatchBc(IdToRemove[i]);;
	//Adjust arrays
	for(int i=IdToRemove.size()-1;i>=0;i--)
	{
		PatchIdInType.erase(PatchIdInType.begin()+IdToRemove[i]);
		PatchTypeInType.erase(PatchTypeInType.begin()+IdToRemove[i]);
		IdxNode.erase(IdxNode.begin()+IdToRemove[i]);
	}
	NumberOfPatchBc=NumberOfPatchBc+4-IdToRemove.size();
	// add index in the patches
	for (int i=0;i<NumberOfPatchBc;i++)
	{
		SetNodeIdxForPatchBc(PatchTypeInType[i], PatchIdInType[i],i,IdxNode[i]);
	}
	//Reordering the variables PatchTypeInType and PatchIdInType to have first periodic then symmetry then others for treatment of Special walls and Global corners
	ReorderingPatches();
}
PatchBc::~PatchBc(){

}
void PatchBc::ReorderingPatches(){
	std::vector<SolverEnum::PatchType> PatchTypeInTypetmp;
	std::vector<int> PatchIdInTypetmp;

	for (int i=0;i<NumberOfPatchBc;i++)
	{
		if(PatchTypeInType[i]==SolverEnum::Periodic)
		{
			PatchTypeInTypetmp.push_back(PatchTypeInType[i]);
			PatchIdInTypetmp.push_back(PatchIdInType[i]);
		}
	}
	for (int i=0;i<NumberOfPatchBc;i++)
	{
		if(PatchTypeInType[i]==SolverEnum::Symmetry)
		{
			PatchTypeInTypetmp.push_back(PatchTypeInType[i]);
			PatchIdInTypetmp.push_back(PatchIdInType[i]);
		}
	}
	for (int i=0;i<NumberOfPatchBc;i++)
	{
		if(PatchTypeInType[i]==SolverEnum::Velocity)
		{
			PatchTypeInTypetmp.push_back(PatchTypeInType[i]);
			PatchIdInTypetmp.push_back(PatchIdInType[i]);
		}
	}
	for (int i=0;i<NumberOfPatchBc;i++)
	{
		if(PatchTypeInType[i]==SolverEnum::Pressure)
		{
			PatchTypeInTypetmp.push_back(PatchTypeInType[i]);
			PatchIdInTypetmp.push_back(PatchIdInType[i]);
		}
	}
	for (int i=0;i<NumberOfPatchBc;i++)
	{
		if(PatchTypeInType[i]==SolverEnum::Wall)
		{
			PatchTypeInTypetmp.push_back(PatchTypeInType[i]);
			PatchIdInTypetmp.push_back(PatchIdInType[i]);
		}
	}
	PatchTypeInType=PatchTypeInTypetmp;
	PatchIdInType=PatchIdInTypetmp;
}
void PatchBc::SetNodeIdxForPatchBc(SolverEnum::PatchType Type_, int IdIntype,int PatchId,std::vector<int> & nodeIdx){
	// Add new PatchBc type here
	switch(Type_)
	{
	case SolverEnum::Pressure:
		PressurePatch[IdIntype].Set_NodeIndex(nodeIdx);
		break;
	case SolverEnum::Velocity:
		VelocityPatch[IdIntype].Set_NodeIndex(nodeIdx);
		break;
	case SolverEnum::Symmetry:
		SymmetryPatch[IdIntype].Set_NodeIndex(nodeIdx);
		break;
	case SolverEnum::Periodic:
		PeriodicPatch[IdIntype].Set_NodeIndex(nodeIdx);
		break;
/*	case SolverEnum::Wall:
		WallPatch[IdIntype].Set_NodeIndex(nodeIdx);
		break;*/
	}
}
void PatchBc::Set_NodeIndexByType(SolverEnum::PatchType Type_, int IdIntype,std::vector<int>& nodeIdx){
	// Add new PatchBc type here
	switch(Type_)
	{
	case SolverEnum::Pressure:
		PressurePatch[IdIntype].Set_NodeIndexByType(nodeIdx);
		break;
	case SolverEnum::Velocity:
		VelocityPatch[IdIntype].Set_NodeIndexByType(nodeIdx);
		break;
	case SolverEnum::Symmetry:
		SymmetryPatch[IdIntype].Set_NodeIndexByType(nodeIdx);
		break;
	case SolverEnum::Periodic:
		PeriodicPatch[IdIntype].Set_NodeIndexByType(nodeIdx);
		break;
//	case SolverEnum::Wall:
//		WallPatch[IdIntype].Set_NodeIndexByType(nodeIdx);
//		break;
	}
}

void PatchBc::Set_NodeIndexSpecialWalls(SolverEnum::PatchType Type_, int IdIntype,std::vector<int>& nodeIdx){
	// Add new PatchBc type here
	switch(Type_)
	{
	case SolverEnum::Pressure:
		PressurePatch[IdIntype].Set_NodeIndexSpecialWalls(nodeIdx);
		break;
	case SolverEnum::Velocity:
		VelocityPatch[IdIntype].Set_NodeIndexSpecialWalls(nodeIdx);
		break;
	case SolverEnum::Symmetry:
		SymmetryPatch[IdIntype].Set_NodeIndexSpecialWalls(nodeIdx);
		break;
	case SolverEnum::Periodic:
		PeriodicPatch[IdIntype].Set_NodeIndexSpecialWalls(nodeIdx);
		break;
//	case SolverEnum::Wall:
//		WallPatch[IdIntype].Set_NodeIndexByType(nodeIdx);
//		break;
	}
}

void PatchBc::Set_NodeIndexGlobalCorner(SolverEnum::PatchType Type_, int IdIntype,std::vector<int>& nodeIdx){
	// Add new PatchBc type here
	switch(Type_)
	{
	case SolverEnum::Pressure:
		PressurePatch[IdIntype].Set_NodeIndexGlobalCorner(nodeIdx);
		break;
	case SolverEnum::Velocity:
		VelocityPatch[IdIntype].Set_NodeIndexGlobalCorner(nodeIdx);
		break;
	case SolverEnum::Symmetry:
		SymmetryPatch[IdIntype].Set_NodeIndexGlobalCorner(nodeIdx);
		break;
	case SolverEnum::Periodic:
		PeriodicPatch[IdIntype].Set_NodeIndexGlobalCorner(nodeIdx);
		break;
//	case SolverEnum::Wall:
//		WallPatch[IdIntype].Set_NodeIndexByType(nodeIdx);
//		break;
	}
}
void PatchBc::SelectPatchBc(SolverEnum::PatchType Type_, Parameters *PtrParam, string PatchBcNames, int PatchId){
	// Add new PatchBc type here
	switch(Type_)
	{
	case SolverEnum::Pressure:
		AddPressurePatch(PatchBcNames, PtrParam->Get_PressureModel(), PtrParam->Get_PressureType());

		break;
	case SolverEnum::Velocity:
		AddVelocityPatch(PatchBcNames, PtrParam->Get_VelocityModel(), PtrParam->Get_VelocityType());
		break;
	case SolverEnum::Symmetry:
		AddSymmetryPatch(PatchBcNames, PtrParam->Get_SymmetryType());
		break;
	case SolverEnum::Periodic:
		AddPeriodicPatch(PatchBcNames, PtrParam->Get_PeriodicType());
		break;
//	case SolverEnum::Wall:
//		AddWallPatch(PatchBcNames, PtrParam->Get_WallType());
//		break;
	}
}
void PatchBc::SelectPatchBc(NodeType Type_, Parameters *PtrParam, string PatchBcNames, int PatchId){
	// Add new PatchBc type here
	switch(Type_)
	{
	case Pressure:
		AddPressurePatch(PatchBcNames, PtrParam->Get_PressureModel(), PtrParam->Get_PressureType());
		break;
	case Velocity:
		AddVelocityPatch(PatchBcNames, PtrParam->Get_VelocityModel(), PtrParam->Get_VelocityType());
		break;
	case Symmetry:
		AddSymmetryPatch(PatchBcNames, PtrParam->Get_SymmetryType());
		break;
	case Periodic:
		AddPeriodicPatch(PatchBcNames, PtrParam->Get_PeriodicType());
		break;
//	case Wall:
//		AddWallPatch(PatchBcNames, PtrParam->Get_WallType());
//		break;
	}
}
void PatchBc::RemovePatchBc(int PatchId){
	switch(PatchTypeInType[PatchId])
		{
		case SolverEnum::Pressure:
			if (PatchIdInType[PatchId]<=PressurePatch.size()-1)
				for (int i=PatchId;i<NumberOfPatchBc+4;i++)
					if (PatchTypeInType[i]==Pressure)
						PatchIdInType[i]=PatchIdInType[i]-1;
			if(!PressurePatch.empty())
				PressurePatch.erase(PressurePatch.begin()+PatchIdInType[PatchId]);
			break;
		case SolverEnum::Velocity:
			if (PatchIdInType[PatchId]<=VelocityPatch.size()-1)
				for (int i=PatchId;i<NumberOfPatchBc+4;i++)
					if (PatchTypeInType[i]==Velocity)
						PatchIdInType[i]=PatchIdInType[i]-1;
			if(!VelocityPatch.empty())
				VelocityPatch.erase(VelocityPatch.begin()+PatchIdInType[PatchId]);
			break;
		case SolverEnum::Symmetry:
			if (PatchIdInType[PatchId]<=SymmetryPatch.size()-1)
				for (int i=PatchId;i<NumberOfPatchBc+4;i++)
					if (PatchTypeInType[i]==Symmetry)
						PatchIdInType[i]=PatchIdInType[i]-1;
			if(!SymmetryPatch.empty())
				SymmetryPatch.erase(SymmetryPatch.begin()+PatchIdInType[PatchId]);
			break;
		case SolverEnum::Periodic:
			if (PatchIdInType[PatchId]<=PeriodicPatch.size()-1)
				for (int i=PatchId;i<NumberOfPatchBc+4;i++)
					if (PatchTypeInType[i]==Periodic)
						PatchIdInType[i]=PatchIdInType[i]-1;
			if(!PeriodicPatch.empty())
				PeriodicPatch.erase(PeriodicPatch.begin()+PatchIdInType[PatchId]);
			break;
/*		case SolverEnum::Wall:
			if (PatchIdInType[PatchId]<=WallPatch.size()-1)
				for (int i=PatchId;i<NumberOfPatchBc+4;i++)
					if (PatchTypeInType[i]==Wall)
						PatchIdInType[i]=PatchIdInType[i]-1;
			if(!WallPatch.empty())
				WallPatch.erase(WallPatch.begin()+PatchIdInType[PatchId]);
			break;
			*/
		}
}
