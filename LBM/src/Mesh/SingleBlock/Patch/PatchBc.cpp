/*
 * PatchBc.cpp
 *
 *  Created on: 31 Oct 2016
 *      Author: Thomas Burel
 */
#include "PatchBc.h"

PatchBc::PatchBc(){
	NumberOfPatchBc=0;
}
void PatchBc::initPatchBc(std::vector<Node2D*> Node, Parameters *PtrParam){
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	//Initialise user patches
	NumberOfPatchBc=0;
	IntitialiseUserPatchBc(*PtrParam,NumberOfPatchBc);
	//check if the number of patch is less or equal to the number of patch and less than 5
	if(NumberOfPatchBc>NumberOfPatches()||NumberOfPatchBc>4)
	{
		if(rank==0)
			std::cout<<"Number of Patches set: "<<NumberOfPatchBc<<std::endl;
		if(NumberOfPatchBc>NumberOfPatches())
			NumberOfPatchBc=NumberOfPatches();
		else
			NumberOfPatchBc=4;
		if(rank==0)
			std::cout<<"Number of Patches corrected: "<<NumberOfPatchBc<<std::endl;
	}
	//Add Default Patches for the boundary of the Cartesian domain
	SelectPatchBc(PtrParam->Get_GlobalBcType(3), PtrParam, string("Left"), NumberOfPatchBc+0);
	SelectPatchBc(PtrParam->Get_GlobalBcType(1), PtrParam, string("Right"), NumberOfPatchBc+1);
	SelectPatchBc(PtrParam->Get_GlobalBcType(0), PtrParam, string("Bottom"), NumberOfPatchBc+2);
	SelectPatchBc(PtrParam->Get_GlobalBcType(2), PtrParam, string("Top"), NumberOfPatchBc+3);
	//Define temporary variables
	double *pos;int IdPatchBc1=-1;int IdPatchBc2=-1;
//	std::vector<int>  *IdxNode;
//	IdxNode=new std::vector<int>[NumberOfPatchBc+4];

	std::vector<std::vector<int> >  IdxNode;
	for (int i=0;i<NumberOfPatchBc+4;i++)
		IdxNode.push_back(std::vector<int> (0));
	//Initialise vector to empty vectors
	for (int i=0;i<NumberOfPatchBc+4;i++)
		IdxNode[i].clear();

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
			IdPatchBc1=-1;IdPatchBc2=-1;
			SetUserPatchBc(*PtrParam,0, i, pos, IdPatchBc1, IdPatchBc2);
			if((IdPatchBc1>=0 && IdPatchBc1<NumberOfPatchBc)||(IdPatchBc2>=0 && IdPatchBc2<NumberOfPatchBc && IdPatchBc1!=IdPatchBc2))
			{
				if(IdPatchBc1>=0 && IdPatchBc1<NumberOfPatchBc)
				{
					if(IdPatchBc1<NumberOfPatchBc)
						IdxNode[IdPatchBc1].push_back(i);
				}
				if(IdPatchBc2>=0 && IdPatchBc2<NumberOfPatchBc && IdPatchBc1!=IdPatchBc2)
				{
					if(IdPatchBc2<NumberOfPatchBc)
						IdxNode[IdPatchBc2].push_back(i);
				}
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
	//Set Inlet or Outlet
	Set_Inlet_outletPatch();
	//Set orientation
	Set_orientation(Node);
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
	case SolverEnum::Wall:
		WallPatch[IdIntype].Set_NodeIndex(nodeIdx);
		break;
	}
}

void PatchBc::Set_NodeIndex(SolverEnum::PatchType Type_, int IdIntype,std::vector<int>& nodeIdx){
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
	case SolverEnum::Wall:
		WallPatch[IdIntype].Set_NodeIndex(nodeIdx);
		break;
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
	case SolverEnum::Wall:
		//WallPatch[IdIntype].Set_NodeIndexByType(nodeIdx);
		break;
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
	case SolverEnum::Wall:
		//WallPatch[IdIntype].Set_NodeIndexByType(nodeIdx);
		break;
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
	case SolverEnum::Wall:
		WallPatch[IdIntype].Set_NodeIndexByType(nodeIdx);
		break;
	}
}

void PatchBc::Set_NodeIndexByTypeSpecialWalls(SolverEnum::PatchType Type_, int IdIntype,std::vector<int>& nodeIdx){
	// Add new PatchBc type here
	switch(Type_)
	{
	case SolverEnum::Pressure:
		PressurePatch[IdIntype].Set_NodeIndexByTypeSpecialWalls(nodeIdx);
		break;
	case SolverEnum::Velocity:
		VelocityPatch[IdIntype].Set_NodeIndexByTypeSpecialWalls(nodeIdx);
		break;
	case SolverEnum::Symmetry:
		SymmetryPatch[IdIntype].Set_NodeIndexByTypeSpecialWalls(nodeIdx);
		break;
	case SolverEnum::Periodic:
		PeriodicPatch[IdIntype].Set_NodeIndexByTypeSpecialWalls(nodeIdx);
		break;
	case SolverEnum::Wall:
		//WallPatch[IdIntype].Set_NodeIndexByType(nodeIdx);
		break;
	}
}

void PatchBc::Set_NodeIndexByTypeGlobalCorner(SolverEnum::PatchType Type_, int IdIntype,std::vector<int>& nodeIdx){
	// Add new PatchBc type here
	switch(Type_)
	{
	case SolverEnum::Pressure:
		PressurePatch[IdIntype].Set_NodeIndexByTypeGlobalCorner(nodeIdx);
		break;
	case SolverEnum::Velocity:
		VelocityPatch[IdIntype].Set_NodeIndexByTypeGlobalCorner(nodeIdx);
		break;
	case SolverEnum::Symmetry:
		SymmetryPatch[IdIntype].Set_NodeIndexByTypeGlobalCorner(nodeIdx);
		break;
	case SolverEnum::Periodic:
		PeriodicPatch[IdIntype].Set_NodeIndexByTypeGlobalCorner(nodeIdx);
		break;
	case SolverEnum::Wall:
		//WallPatch[IdIntype].Set_NodeIndexByType(nodeIdx);
		break;
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
	case SolverEnum::Wall:
		AddWallPatch(PatchBcNames, PtrParam->Get_WallType());
		break;
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
	case Wall:
		AddWallPatch(PatchBcNames, PtrParam->Get_WallType());
		break;
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
		case SolverEnum::Wall:
			if (PatchIdInType[PatchId]<=WallPatch.size()-1)
				for (int i=PatchId;i<NumberOfPatchBc+4;i++)
					if (PatchTypeInType[i]==Wall)
						PatchIdInType[i]=PatchIdInType[i]-1;
			if(!WallPatch.empty())
				WallPatch.erase(WallPatch.begin()+PatchIdInType[PatchId]);
			break;

		}
}
bool PatchBc::Get_NodeIndex(int idPatch,std::vector<int>* &NodeIndex,std::vector<int>* &NodeIndexSpecialWalls,std::vector<int>* &NodeIndexGlobalCorner){
	bool NonEmptyPatch;
	switch(PatchTypeInType[idPatch])
			{
			case SolverEnum::Periodic:
				if(PeriodicPatch[PatchIdInType[idPatch]].Get_NodeIndex().empty()&&
						PeriodicPatch[PatchIdInType[idPatch]].Get_NodeIndexSpecialWalls().empty()&&
						PeriodicPatch[PatchIdInType[idPatch]].Get_NodeIndexGlobalCorner().empty())
				{
					NodeIndex=0;NodeIndexSpecialWalls=0;NodeIndexGlobalCorner=0;
					NonEmptyPatch= false;
				}
				else
				{
					NodeIndex=&PeriodicPatch[PatchIdInType[idPatch]].Get_NodeIndex();
					NodeIndexSpecialWalls=&PeriodicPatch[PatchIdInType[idPatch]].Get_NodeIndexSpecialWalls();
					NodeIndexGlobalCorner=&PeriodicPatch[PatchIdInType[idPatch]].Get_NodeIndexGlobalCorner();
					NonEmptyPatch= true;
				}
				break;
			case SolverEnum::Symmetry:
				if(SymmetryPatch[PatchIdInType[idPatch]].Get_NodeIndex().empty()&&
						SymmetryPatch[PatchIdInType[idPatch]].Get_NodeIndexSpecialWalls().empty()&&
						SymmetryPatch[PatchIdInType[idPatch]].Get_NodeIndexGlobalCorner().empty())
				{
					NodeIndex=0;NodeIndexSpecialWalls=0;NodeIndexGlobalCorner=0;
					NonEmptyPatch= false;
				}
				else
				{
					NodeIndex=&SymmetryPatch[PatchIdInType[idPatch]].Get_NodeIndex();
					NodeIndexSpecialWalls=&SymmetryPatch[PatchIdInType[idPatch]].Get_NodeIndexSpecialWalls();
					NodeIndexGlobalCorner=&SymmetryPatch[PatchIdInType[idPatch]].Get_NodeIndexGlobalCorner();
					NonEmptyPatch= true;
				}
				break;
			case SolverEnum::Pressure:
				if(PressurePatch[PatchIdInType[idPatch]].Get_NodeIndex().empty()&&
						PressurePatch[PatchIdInType[idPatch]].Get_NodeIndexSpecialWalls().empty()&&
						PressurePatch[PatchIdInType[idPatch]].Get_NodeIndexGlobalCorner().empty())
				{
					NodeIndex=0;NodeIndexSpecialWalls=0;NodeIndexGlobalCorner=0;
					NonEmptyPatch= false;
				}
				else
				{
					NodeIndex=&PressurePatch[PatchIdInType[idPatch]].Get_NodeIndex();
					NodeIndexSpecialWalls=&PressurePatch[PatchIdInType[idPatch]].Get_NodeIndexSpecialWalls();
					NodeIndexGlobalCorner=&PressurePatch[PatchIdInType[idPatch]].Get_NodeIndexGlobalCorner();
					NonEmptyPatch= true;
				}
				break;
			case SolverEnum::Velocity:
				if(VelocityPatch[PatchIdInType[idPatch]].Get_NodeIndex().empty()&&
						VelocityPatch[PatchIdInType[idPatch]].Get_NodeIndexSpecialWalls().empty()&&
						VelocityPatch[PatchIdInType[idPatch]].Get_NodeIndexGlobalCorner().empty())
				{
					NodeIndex=0;NodeIndexSpecialWalls=0;NodeIndexGlobalCorner=0;
					NonEmptyPatch= false;
				}
				else
				{
					NodeIndex=&VelocityPatch[PatchIdInType[idPatch]].Get_NodeIndex();
					NodeIndexSpecialWalls=&VelocityPatch[PatchIdInType[idPatch]].Get_NodeIndexSpecialWalls();
					NodeIndexGlobalCorner=&VelocityPatch[PatchIdInType[idPatch]].Get_NodeIndexGlobalCorner();
					NonEmptyPatch= true;
				}
				break;
			case SolverEnum::Wall:
				NodeIndex=&WallPatch[PatchIdInType[idPatch]].Get_NodeIndex();
				NodeIndexSpecialWalls=&WallPatch[PatchIdInType[idPatch]].Get_NodeIndexSpecialWalls();
				NodeIndexGlobalCorner=&WallPatch[PatchIdInType[idPatch]].Get_NodeIndexGlobalCorner();
				NonEmptyPatch= false;
				break;
			}
	return NonEmptyPatch;
}
void PatchBc::Set_Inlet_outletPatch(){

	for (int i=0;i<NumberOfPatchBc;i++)
	{
		switch(PatchTypeInType[i])
				{
				case SolverEnum::Periodic:
						inlet.push_back(PeriodicPatch[PatchIdInType[i]].Get_Inlet());
						outlet.push_back(PeriodicPatch[PatchIdInType[i]].Get_Outlet());
					break;
				case SolverEnum::Symmetry:
					inlet.push_back(SymmetryPatch[PatchIdInType[i]].Get_Inlet());
					outlet.push_back(SymmetryPatch[PatchIdInType[i]].Get_Outlet());
					break;
				case SolverEnum::Pressure:
					inlet.push_back(PressurePatch[PatchIdInType[i]].Get_Inlet());
					outlet.push_back(PressurePatch[PatchIdInType[i]].Get_Outlet());
					break;
				case SolverEnum::Velocity:

						inlet.push_back(VelocityPatch[PatchIdInType[i]].Get_Inlet());
						outlet.push_back(VelocityPatch[PatchIdInType[i]].Get_Outlet());
					break;
				case SolverEnum::Wall:
					inlet.push_back(false);
					outlet.push_back(false);
					break;
				}
	}

}
void PatchBc::Set_orientation(std::vector<Node2D*> Node){
	for (int i=0;i<NumberOfPatchBc;i++)
	{
		switch(PatchTypeInType[i])
				{
				case SolverEnum::Periodic:
					if(PeriodicPatch[PatchIdInType[i]].Get_NodeIndex().size()>1)
					{
						if(Node[PeriodicPatch[PatchIdInType[i]].Get_NodeIndex()[0]]->get_x()==Node[PeriodicPatch[PatchIdInType[i]].Get_NodeIndex()[1]]->get_x())
							Orientation.push_back(0);
						else
							Orientation.push_back(1);
					}
					else
					{
						//std::cout<<"Orientation of the Patch cannot be determined."<<std::endl;
						Orientation.push_back(-1);//Set to error value
					}
					break;
				case SolverEnum::Symmetry:
					if(SymmetryPatch[PatchIdInType[i]].Get_NodeIndex().size()>1)
					{
						if(Node[SymmetryPatch[PatchIdInType[i]].Get_NodeIndex()[0]]->get_x()==Node[SymmetryPatch[PatchIdInType[i]].Get_NodeIndex()[1]]->get_x())
							Orientation.push_back(0);
						else
							Orientation.push_back(1);
					}
					else
					{
						//std::cout<<"Orientation of the Patch cannot be determined."<<std::endl;
						Orientation.push_back(-1);//Set to error value
					}
					break;
				case SolverEnum::Pressure:
					if(PressurePatch[PatchIdInType[i]].Get_NodeIndex().size()>1)
					{
						if(Node[PressurePatch[PatchIdInType[i]].Get_NodeIndex()[0]]->get_x()==Node[PressurePatch[PatchIdInType[i]].Get_NodeIndex()[1]]->get_x())
							Orientation.push_back(0);
						else
							Orientation.push_back(1);
					}
					else
					{
						//std::cout<<"Orientation of the Patch cannot be determined."<<std::endl;
						Orientation.push_back(-1);//Set to error value
					}
					break;
				case SolverEnum::Velocity:
					if(VelocityPatch[PatchIdInType[i]].Get_NodeIndex().size()>1)
					{
						if(Node[VelocityPatch[PatchIdInType[i]].Get_NodeIndex()[0]]->get_x()==Node[VelocityPatch[PatchIdInType[i]].Get_NodeIndex()[1]]->get_x())
							Orientation.push_back(0);
						else
							Orientation.push_back(1);
					}
					else
					{
						//std::cout<<"Orientation of the Patch cannot be determined."<<std::endl;
						Orientation.push_back(-1);//Set to error value
					}
					break;
				case SolverEnum::Wall:
					if(WallPatch[PatchIdInType[i]].Get_NodeIndex().size()>1)
					{
						if(Node[WallPatch[PatchIdInType[i]].Get_NodeIndex()[0]]->get_x()==Node[WallPatch[PatchIdInType[i]].Get_NodeIndex()[1]]->get_x())
							Orientation.push_back(0);
						else
							Orientation.push_back(1);
					}
					else
					{
						//std::cout<<"Orientation of the Patch cannot be determined."<<std::endl;
						Orientation.push_back(-1);//Set to error value
					}
					break;
				}
	}
}
