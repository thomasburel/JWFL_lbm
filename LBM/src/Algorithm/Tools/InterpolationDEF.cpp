/*
 * InterpolationTailor.cpp
 *
 *  Created on: 28 Jul 2016
 *      Author: thomas
 */

#include "../Tools/InterpolationDEF.h"

#include <iostream>
InterpolationDEF::InterpolationDEF() {
	dimension=0;
	nb_Vel=0;
	PtrOppositeInterpol=0;
}
InterpolationDEF::InterpolationDEF(int dimension_, int nb_vel,unsigned int *PtrOppositeInterpol_) {
	dimension=dimension_;
	nb_Vel=nb_vel;
	PtrOppositeInterpol=PtrOppositeInterpol_;
}
InterpolationDEF::~InterpolationDEF() {

}
NoInterpolation::NoInterpolation(){}
NoInterpolation::NoInterpolation(int dimension, int nb_vel,unsigned int *PtrOppositeInterpol_){}
NoInterpolation::~NoInterpolation(){}

InterpolationLinear::InterpolationLinear() {
	dimension=0;
	nb_Vel=0;
}

InterpolationLinear::~InterpolationLinear() {
}

InterpolationLinear::InterpolationLinear(int dimension_, int nb_vel,unsigned int *PtrOppositeInterpol_){
	dimension=dimension_;
	nb_Vel=nb_vel;
	PtrOppositeInterpol=PtrOppositeInterpol_;
}
///First order Tailor extrapolation
void InterpolationLinear::InterpolationOnWall (double *Var, int * Connect, int & normal){
	Var[Connect[0]]=0.5*(Var[Connect[PtrOppositeInterpol[normal]]]+Var[Connect[normal]]);
}
void InterpolationLinear::InterpolationOnWall (double *Var1, double *Var2, int * Connect, int & normal){
	Var1[Connect[0]]=0.5*(Var1[Connect[PtrOppositeInterpol[normal]]]+Var1[Connect[normal]]);
	Var2[Connect[0]]=0.5*(Var2[Connect[PtrOppositeInterpol[normal]]]+Var2[Connect[normal]]);
}
///First order Tailor extrapolation
void InterpolationLinear::InterpolationOnCornerConcave (double *Var, int * Connect, int & normal){
	Var[Connect[0]]=0.5*(Var[Connect[PtrOppositeInterpol[normal]]]+Var[Connect[normal]]);
}
void InterpolationLinear::InterpolationOnCornerConcave (double *Var1, double *Var2, int * Connect, int & normal){
	Var1[Connect[0]]=0.5*(Var1[Connect[PtrOppositeInterpol[normal]]]+Var1[Connect[normal]]);
	Var2[Connect[0]]=0.5*(Var2[Connect[PtrOppositeInterpol[normal]]]+Var2[Connect[normal]]);}
///First order Tailor extrapolation
void InterpolationLinear::InterpolationOnCornerConvex(double *Var, int * Connect, int & normal){
	Var[Connect[0]]=0.5*(Var[Connect[PtrOppositeInterpol[normal]]]+Var[Connect[normal]]);
}
void InterpolationLinear::InterpolationOnCornerConvex(double *Var1, double *Var2, int * Connect, int & normal){
	Var1[Connect[0]]=0.5*(Var1[Connect[PtrOppositeInterpol[normal]]]+Var1[Connect[normal]]);
	Var2[Connect[0]]=0.5*(Var2[Connect[PtrOppositeInterpol[normal]]]+Var2[Connect[normal]]);
}
InterpolationInverseWeightDistance::InterpolationInverseWeightDistance(){

}
InterpolationInverseWeightDistance::InterpolationInverseWeightDistance(int dimension, int nb_vel,unsigned int *PtrOppositeInterpol_){

}
InterpolationInverseWeightDistance::~InterpolationInverseWeightDistance(){

}
InterpolationLinearLeastSquare::InterpolationLinearLeastSquare(){
	sumX=0;	sumY=0;sumXY=0;sumX2=0;
	MapWallId=0;
	nWalls=0;
	x0=0;y0=0;
	result1=0;result2=0;
}
InterpolationLinearLeastSquare::InterpolationLinearLeastSquare(int dimension_,int nb_vel,unsigned int *PtrOppositeInterpol_){
	dimension=dimension_;
	nb_Vel=nb_vel;
	PtrOppositeInterpol=PtrOppositeInterpol_;
	sumX=0;	sumY=0;sumXY=0;sumX2=0;
	MapWallId=0;
	nWalls=0;
	x0=0;y0=0;
	result1=0;result2=0;
}
InterpolationLinearLeastSquare::~InterpolationLinearLeastSquare(){

	delete []MapWallId;
}
void InterpolationLinearLeastSquare::InitInterpol(NodeArrays2D *PtrNodes, Parameters *PtrParam){
	MapWallId=new int [PtrNodes->TypeOfNode.size()];
	nWalls=PtrNodes->NodeWall.size()+PtrNodes->CornerConcave.size()+PtrNodes->CornerConvex.size();

	nSolidNodes=PtrParam->Get_NumberOfInterpolNodeInSolid();
	nFluidNodes=PtrParam->Get_NumberOfInterpolNodeInFluid();

	int countfluid=0;int countsolid=0; int maxwhile=2*max(nSolidNodes,nFluidNodes); int countwhile=0; int countnWalls=0;
	NextNode Nodetmp;

	if(PtrNodes->CornerConcave.size()>0)
	for (int i=0;i<PtrNodes->CornerConcave.size();i++)
	{
		//Reset variables
		SolidChecked.clear();FluidChecked.clear();
		nextwall.clear(); nextinterior.clear();
		countfluid=0;countsolid=0;countwhile=0;
		//Save the mapping
		MapWallId[PtrNodes->NodeCorner[PtrNodes->CornerConcave[i]].Get_index()]=countnWalls;
		//Keep the origin node position to calculate the distance
		x0=PtrNodes->NodeCorner[PtrNodes->CornerConcave[i]].get_x();y0=PtrNodes->NodeCorner[PtrNodes->CornerConcave[i]].get_y();
		//Set the first wall
		nextwall.push_back(PtrNodes->NodeCorner[PtrNodes->CornerConcave[i]].Get_index());
		//set the first interior
		Nodetmp.index=(PtrNodes->NodeCorner[PtrNodes->CornerConcave[i]].Get_connect()[PtrNodes->NodeCorner[PtrNodes->CornerConcave[i]].Get_BcNormal()]);
		Nodetmp.distance=DistToWall(PtrNodes->NodeInterior[PtrNodes->NodeIndexByType[Nodetmp.index]].get_x(),PtrNodes->NodeInterior[PtrNodes->NodeIndexByType[Nodetmp.index]].get_y());
		Nodetmp.rankMarkNodes=0;
		nextinterior.push_back(Nodetmp);

		while((countfluid<nFluidNodes || countsolid<nSolidNodes) && countwhile < maxwhile)
		{
			if(nextwall.size()>0 && countsolid<nSolidNodes)
			{
				Mark_SolidIds(PtrNodes,nextwall,countnWalls,countsolid);
				nextwallprevious=nextwall;
				nextwall.clear();
				Next_WallId(PtrNodes,nextwallprevious,nextwall);
			}
			if(nextinterior.size()>0 && countfluid<nFluidNodes)
			{
				Mark_FluidIds(PtrNodes,nextinterior,countnWalls,countfluid);
				nextinteriorprevious=nextinterior;
				nextinterior.clear();
				Next_FluidId(PtrNodes,nextinteriorprevious,nextinterior);
			}
			countwhile++;
		}
		solidId.push_back(solidIdtmp);
		solidIdtmp.clear();
		fluidId.push_back(fluidIdtmp);
		fluidIdtmp.clear();
		solidDist.push_back(solidDisttmp);
		solidDisttmp.clear();
		fluidDist.push_back(fluidDisttmp);
		fluidDisttmp.clear();
		countnWalls++;
	}
		if(PtrNodes->NodeWall.size()>0)
	for (int i=0;i<PtrNodes->NodeWall.size();i++)
	{
		//Reset variables
		SolidChecked.clear();FluidChecked.clear();
		nextwall.clear(); nextinterior.clear();
		countfluid=0;countsolid=0;countwhile=0;
		//Save the mapping
		MapWallId[PtrNodes->NodeWall[i].Get_index()]=countnWalls;
		//Keep the origin node position to calculate the distance
		x0=PtrNodes->NodeWall[i].get_x();y0=PtrNodes->NodeWall[i].get_y();
		//Set the first wall
		nextwall.push_back(PtrNodes->NodeWall[i].Get_index());
		//set the first interior
		Nodetmp.index=(PtrNodes->NodeWall[i].Get_connect()[PtrNodes->NodeWall[i].Get_BcNormal()]);
		Nodetmp.distance=DistToWall(PtrNodes->NodeInterior[PtrNodes->NodeIndexByType[Nodetmp.index]].get_x(),PtrNodes->NodeInterior[PtrNodes->NodeIndexByType[Nodetmp.index]].get_y());
		Nodetmp.rankMarkNodes=0;
		nextinterior.push_back(Nodetmp);

		while((countfluid<nFluidNodes || countsolid<nSolidNodes) && countwhile < maxwhile)
		{
			if(nextwall.size()>0 && countsolid<nSolidNodes)
			{
				Mark_SolidIds(PtrNodes,nextwall,countnWalls,countsolid);
				nextwallprevious=nextwall;
				nextwall.clear();
				Next_WallId(PtrNodes,nextwallprevious,nextwall);
			}
			if(nextinterior.size()>0 && countfluid<nFluidNodes)
			{
				Mark_FluidIds(PtrNodes,nextinterior,countnWalls,countfluid);
				nextinteriorprevious=nextinterior;
				nextinterior.clear();
				Next_FluidId(PtrNodes,nextinteriorprevious,nextinterior);
			}
			countwhile++;
		}
		solidId.push_back(solidIdtmp);
		solidIdtmp.clear();
		fluidId.push_back(fluidIdtmp);
		fluidIdtmp.clear();
		solidDist.push_back(solidDisttmp);
		solidDisttmp.clear();
		fluidDist.push_back(fluidDisttmp);
		fluidDisttmp.clear();
		countnWalls++;
	}
	if(PtrNodes->CornerConvex.size()>0)
	for (int i=0;i<PtrNodes->CornerConvex.size();i++)
	{
		//Reset variables
		SolidChecked.clear();FluidChecked.clear();
		nextwall.clear(); nextinterior.clear();
		countfluid=0;countsolid=0;countwhile=0;
		//Save the mapping
		MapWallId[PtrNodes->NodeCorner[PtrNodes->CornerConvex[i]].Get_index()]=countnWalls;
		//Keep the origin node position to calculate the distance
		x0=PtrNodes->NodeCorner[PtrNodes->CornerConvex[i]].get_x();y0=PtrNodes->NodeCorner[PtrNodes->CornerConvex[i]].get_y();
		//Set the first wall
		nextwall.push_back(PtrNodes->NodeCorner[PtrNodes->CornerConvex[i]].Get_index());
		//set the first interior
		Nodetmp.index=(PtrNodes->NodeCorner[PtrNodes->CornerConvex[i]].Get_connect()[PtrNodes->NodeCorner[PtrNodes->CornerConvex[i]].Get_BcNormal()]);
		Nodetmp.distance=DistToWall(PtrNodes->NodeInterior[PtrNodes->NodeIndexByType[Nodetmp.index]].get_x(),PtrNodes->NodeInterior[PtrNodes->NodeIndexByType[Nodetmp.index]].get_y());
		Nodetmp.rankMarkNodes=0;
		nextinterior.push_back(Nodetmp);
		while((countfluid<nFluidNodes || countsolid<nSolidNodes) && countwhile < maxwhile)
		{
			if(nextwall.size()>0 && countsolid<nSolidNodes)
			{
				Mark_SolidIds(PtrNodes,nextwall,countnWalls,countsolid);
				nextwallprevious=nextwall;
				nextwall.clear();
				Next_WallId(PtrNodes,nextwallprevious,nextwall);
			}
			if(nextinterior.size()>0 && countfluid<nFluidNodes)
			{
				Mark_FluidIds(PtrNodes,nextinterior,countnWalls,countfluid);
				nextinteriorprevious=nextinterior;
				nextinterior.clear();
				Next_FluidId(PtrNodes,nextinteriorprevious,nextinterior);
			}
			countwhile++;
		}
		solidId.push_back(solidIdtmp);
		solidIdtmp.clear();
		fluidId.push_back(fluidIdtmp);
		fluidIdtmp.clear();
		solidDist.push_back(solidDisttmp);
		solidDisttmp.clear();
		fluidDist.push_back(fluidDisttmp);
		fluidDisttmp.clear();
		countnWalls++;
	}

}
void InterpolationLinearLeastSquare::Mark_FluidIds(NodeArrays2D *PtrNodes,std::vector<NextNode>& next, int &countnwalls, int &countfluid){
	int nodeIdlocal=0;
	for(int i=0;i<next.size();i++)
	{
		if(countfluid<nFluidNodes)
		{
			countfluid++;
			nodeIdlocal=PtrNodes->NodeIndexByType[next[i].index];
			//fluidId[countnwalls][countfluid]=nodeIdlocal;
			fluidIdtmp.push_back(next[i].index);
			//fluidDist[countnwalls][countfluid]=next[i].distance;
			fluidDisttmp.push_back(next[i].distance);
			next[i].rankMarkNodes=countfluid;
			FluidChecked.push_back(next[i]);

		}
	}
}
void InterpolationLinearLeastSquare::Mark_SolidIds(NodeArrays2D *PtrNodes,std::vector<int>& next, int &countnwalls, int &countsolid){
	int nodeIdlocal=0;
	int nodeIdlocaltmp=0;
	NextNode NextNodetmp;
	for(int i=0;i<next.size();i++)
	{
		nodeIdlocal=PtrNodes->NodeIndexByType[next[i]];
		switch (PtrNodes->TypeOfNode[next[i]])
		{
		case Wall:
			if(countsolid<nSolidNodes)
			{
				countsolid++;
				solidIdtmp.push_back(PtrNodes->NodeWall[nodeIdlocal].Get_connect()[PtrOppositeInterpol[PtrNodes->NodeWall[nodeIdlocal].Get_BcNormal()]]);
				NextNodetmp.index=solidIdtmp.back();
				//solidId[countnwalls][countsolid]=(PtrNodes->NodeWall[nodeIdlocal].Get_connect()[PtrOppositeInterpol[PtrNodes->NodeWall[nodeIdlocal].Get_BcNormal()]]);
				//NextNodetmp.index=solidId[countnwalls][countsolid];
				nodeIdlocaltmp=PtrNodes->NodeIndexByType[solidIdtmp.back()];
				//nodeIdlocaltmp=PtrNodes->NodeIndexByType[solidId[countnwalls][countsolid]];
	//			solidDist[countnwalls][countsolid]=-DistToWall(PtrNodes->NodeSolid[nodeIdlocal].get_x(),PtrNodes->NodeSolid[nodeIdlocal].get_y());
				solidDisttmp.push_back(-DistToWall(PtrNodes->NodeSolid[nodeIdlocaltmp].get_x(),PtrNodes->NodeSolid[nodeIdlocaltmp].get_y()));
			//	NextNodetmp.distance=-solidDist[countnwalls][countsolid];
				NextNodetmp.distance=-solidDisttmp.back();
				NextNodetmp.rankMarkNodes=countsolid;
				SolidChecked.push_back(NextNodetmp);

			}
			break;
		case Corner:
			if(countsolid<nSolidNodes)
			{
				countsolid++;
				solidIdtmp.push_back(PtrNodes->NodeCorner[nodeIdlocal].Get_connect()[PtrOppositeInterpol[PtrNodes->NodeCorner[nodeIdlocal].Get_BcNormal()]]);
				NextNodetmp.index=solidIdtmp.back();
				//solidId[countnwalls][countsolid]=(PtrNodes->NodeCorner[nodeIdlocal].Get_connect()[PtrOppositeInterpol[PtrNodes->NodeCorner[nodeIdlocal].Get_BcNormal()]]);
				//NextNodetmp.index=solidId[countnwalls][countsolid];
				nodeIdlocaltmp=PtrNodes->NodeIndexByType[solidIdtmp.back()];
				//nodeIdlocaltmp=PtrNodes->NodeIndexByType[solidId[countnwalls][countsolid]];
				//solidDist[countnwalls][countsolid]=-DistToWall(PtrNodes->NodeSolid[nodeIdlocal].get_x(),PtrNodes->NodeSolid[nodeIdlocal].get_y());
				solidDisttmp.push_back(-DistToWall(PtrNodes->NodeSolid[nodeIdlocaltmp].get_x(),PtrNodes->NodeSolid[nodeIdlocaltmp].get_y()));
				//NextNodetmp.distance=-solidDist[countnwalls][countsolid];
				NextNodetmp.distance=-solidDisttmp.back();
				NextNodetmp.rankMarkNodes=countsolid;
				SolidChecked.push_back(NextNodetmp);
			}
			break;
		case ConcaveCorner:
			if(countsolid<nSolidNodes)
			{
				countsolid++;
				solidIdtmp.push_back(PtrNodes->NodeCorner[nodeIdlocal].Get_connect()[PtrOppositeInterpol[PtrNodes->NodeCorner[nodeIdlocal].Get_BcNormal()]]);
				NextNodetmp.index=solidIdtmp.back();
	//			solidId[countnwalls][countsolid]=(PtrNodes->NodeCorner[nodeIdlocal].Get_connect()[PtrOppositeInterpol[PtrNodes->NodeCorner[nodeIdlocal].Get_BcNormal()]]);
	//			NextNodetmp.index=solidId[countnwalls][countsolid];
				nodeIdlocaltmp=PtrNodes->NodeIndexByType[solidIdtmp.back()];
				//nodeIdlocaltmp=PtrNodes->NodeIndexByType[solidId[countnwalls][countsolid]];
				//solidDist[countnwalls][countsolid]=-DistToWall(PtrNodes->NodeSolid[nodeIdlocal].get_x(),PtrNodes->NodeSolid[nodeIdlocal].get_y());
				solidDisttmp.push_back(-DistToWall(PtrNodes->NodeSolid[nodeIdlocaltmp].get_x(),PtrNodes->NodeSolid[nodeIdlocaltmp].get_y()));
				//NextNodetmp.distance=-solidDist[countnwalls][countsolid];
				NextNodetmp.distance=-solidDisttmp.back();
				NextNodetmp.rankMarkNodes=countsolid;
				SolidChecked.push_back(NextNodetmp);
			}
			break;
		case ConvexCorner:
			if(countsolid<nSolidNodes)
			{
				countsolid++;
				solidIdtmp.push_back(PtrNodes->NodeCorner[nodeIdlocal].Get_connect()[PtrOppositeInterpol[PtrNodes->NodeCorner[nodeIdlocal].Get_BcNormal()]]);
				NextNodetmp.index=solidIdtmp.back();
				//solidId[countnwalls][countsolid]=(PtrNodes->NodeCorner[nodeIdlocal].Get_connect()[PtrOppositeInterpol[PtrNodes->NodeCorner[nodeIdlocal].Get_BcNormal()]]);
				//NextNodetmp.index=solidId[countnwalls][countsolid];
				nodeIdlocaltmp=PtrNodes->NodeIndexByType[solidIdtmp.back()];
				//nodeIdlocaltmp=PtrNodes->NodeIndexByType[solidId[countnwalls][countsolid]];
				//solidDist[countnwalls][countsolid]=-DistToWall(PtrNodes->NodeSolid[nodeIdlocal].get_x(),PtrNodes->NodeSolid[nodeIdlocal].get_y());
				solidDisttmp.push_back(-DistToWall(PtrNodes->NodeSolid[nodeIdlocaltmp].get_x(),PtrNodes->NodeSolid[nodeIdlocaltmp].get_y()));
				//NextNodetmp.distance=-solidDist[countnwalls][countsolid];
				NextNodetmp.distance=-solidDisttmp.back();
				NextNodetmp.rankMarkNodes=countsolid;
				SolidChecked.push_back(NextNodetmp);
			}
			break;
		case SolidGhost:

			break;
		case GlobalCorner:

			break;
		case SpecialWall:

		/*	if(countsolid<nNodes)
			{
				countsolid++;
				solidIdtmp.push_back(PtrNodes->NodeSpecialWall[nodeIdlocal].Get_connect()[PtrOppositeInterpol[PtrNodes->NodeSpecialWall[nodeIdlocal].Get_BcNormal()]]);
				NextNodetmp.index=solidIdtmp.back();
//				solidId[countnwalls][countsolid]=(PtrNodes->NodeSpecialWall[nodeIdlocal].Get_connect()[PtrOppositeInterpol[PtrNodes->NodeSpecialWall[nodeIdlocal].Get_BcNormal()]]);
//				NextNodetmp.index=solidId[countnwalls][countsolid];
				nodeIdlocaltmp=PtrNodes->NodeIndexByType[solidIdtmp.back()];
				//nodeIdlocaltmp=PtrNodes->NodeIndexByType[solidId[countnwalls][countsolid]];
				solidDist[countnwalls][countsolid]=-DistToWall(PtrNodes->NodeSolid[nodeIdlocal].get_x(),PtrNodes->NodeSolid[nodeIdlocal].get_y());
				NextNodetmp.distance=-solidDist[countnwalls][countsolid];
				NextNodetmp.rankMarkNodes=countsolid;
				SolidChecked.push_back(NextNodetmp);
			}*/
			break;
		default:
			std::cerr<<"Node type not found for interpolation. Type of Node :"<<PtrNodes->TypeOfNode[next[i]]<<std::endl;
			break;
		}

	}
}

void InterpolationLinearLeastSquare::Next_WallId(NodeArrays2D *PtrNodes,std::vector<int> nextback,std::vector<int>& next){
	int nodeIdlocal=0;
	next.clear();
	std::vector<NextNode> VectNextNodetmp;
	NextNode NextNodetmp;
	for(int i=0;i<nextback.size();i++)
	{
		nodeIdlocal=PtrNodes->NodeIndexByType[nextback[i]];

		switch (PtrNodes->TypeOfNode[nextback[i]])
		{
		case Wall:
			switch(PtrNodes->NodeWall[nodeIdlocal].Get_BcNormal())
			{
			case 1:
				NextNodetmp.index=PtrNodes->NodeWall[nodeIdlocal].Get_connect()[2];
				NextNodetmp.distance=DistToWall(PtrNodes->NodeWall[PtrNodes->NodeIndexByType[NextNodetmp.index]].get_x(),PtrNodes->NodeWall[PtrNodes->NodeIndexByType[NextNodetmp.index]].get_y());
				VectNextNodetmp.push_back(NextNodetmp);
				NextNodetmp.index=PtrNodes->NodeWall[nodeIdlocal].Get_connect()[4];
				NextNodetmp.distance=DistToWall(PtrNodes->NodeWall[PtrNodes->NodeIndexByType[NextNodetmp.index]].get_x(),PtrNodes->NodeWall[PtrNodes->NodeIndexByType[NextNodetmp.index]].get_y());
				VectNextNodetmp.push_back(NextNodetmp);
				break;
			case 2:
				NextNodetmp.index=PtrNodes->NodeWall[nodeIdlocal].Get_connect()[1];
				NextNodetmp.distance=DistToWall(PtrNodes->NodeWall[PtrNodes->NodeIndexByType[NextNodetmp.index]].get_x(),PtrNodes->NodeWall[PtrNodes->NodeIndexByType[NextNodetmp.index]].get_y());
				VectNextNodetmp.push_back(NextNodetmp);
				NextNodetmp.index=PtrNodes->NodeWall[nodeIdlocal].Get_connect()[3];
				NextNodetmp.distance=DistToWall(PtrNodes->NodeWall[PtrNodes->NodeIndexByType[NextNodetmp.index]].get_x(),PtrNodes->NodeWall[PtrNodes->NodeIndexByType[NextNodetmp.index]].get_y());
				VectNextNodetmp.push_back(NextNodetmp);
				break;
			case 3:
				NextNodetmp.index=PtrNodes->NodeWall[nodeIdlocal].Get_connect()[2];
				NextNodetmp.distance=DistToWall(PtrNodes->NodeWall[PtrNodes->NodeIndexByType[NextNodetmp.index]].get_x(),PtrNodes->NodeWall[PtrNodes->NodeIndexByType[NextNodetmp.index]].get_y());
				VectNextNodetmp.push_back(NextNodetmp);
				NextNodetmp.index=PtrNodes->NodeWall[nodeIdlocal].Get_connect()[4];
				NextNodetmp.distance=DistToWall(PtrNodes->NodeWall[PtrNodes->NodeIndexByType[NextNodetmp.index]].get_x(),PtrNodes->NodeWall[PtrNodes->NodeIndexByType[NextNodetmp.index]].get_y());
				VectNextNodetmp.push_back(NextNodetmp);
				break;
			case 4:
				NextNodetmp.index=PtrNodes->NodeWall[nodeIdlocal].Get_connect()[1];
				NextNodetmp.distance=DistToWall(PtrNodes->NodeWall[PtrNodes->NodeIndexByType[NextNodetmp.index]].get_x(),PtrNodes->NodeWall[PtrNodes->NodeIndexByType[NextNodetmp.index]].get_y());
				VectNextNodetmp.push_back(NextNodetmp);
				NextNodetmp.index=PtrNodes->NodeWall[nodeIdlocal].Get_connect()[3];
				NextNodetmp.distance=DistToWall(PtrNodes->NodeWall[PtrNodes->NodeIndexByType[NextNodetmp.index]].get_x(),PtrNodes->NodeWall[PtrNodes->NodeIndexByType[NextNodetmp.index]].get_y());
				VectNextNodetmp.push_back(NextNodetmp);
				break;
			}
			break;
			case SpecialWall:
			/*	switch(PtrNodes->NodeSpecialWall[nodeIdlocal].Get_BcNormal())
				{
				case 1:
					next.push_back(PtrNodes->NodeSpecialWall[nodeIdlocal].Get_connect()[2]);
					next.push_back(PtrNodes->NodeSpecialWall[nodeIdlocal].Get_connect()[4]);
					break;
				case 2:
					next.push_back(PtrNodes->NodeSpecialWall[nodeIdlocal].Get_connect()[1]);
					next.push_back(PtrNodes->NodeSpecialWall[nodeIdlocal].Get_connect()[3]);
					break;
				case 3:
					next.push_back(PtrNodes->NodeSpecialWall[nodeIdlocal].Get_connect()[2]);
					next.push_back(PtrNodes->NodeSpecialWall[nodeIdlocal].Get_connect()[4]);
					break;
				case 4:
					next.push_back(PtrNodes->NodeSpecialWall[nodeIdlocal].Get_connect()[1]);
					next.push_back(PtrNodes->NodeSpecialWall[nodeIdlocal].Get_connect()[3]);
					break;
				}*/
				break;
		case Corner:
			switch(PtrNodes->NodeCorner[nodeIdlocal].Get_BcNormal())
			{
			case 5:
				if(PtrNodes->NodeCorner[nodeIdlocal].Get_CornerType()==Concave)
				{
					NextNodetmp.index=PtrNodes->NodeCorner[nodeIdlocal].Get_connect()[1];
					NextNodetmp.distance=DistToWall(PtrNodes->NodeCorner[PtrNodes->NodeIndexByType[NextNodetmp.index]].get_x(),PtrNodes->NodeCorner[PtrNodes->NodeIndexByType[NextNodetmp.index]].get_y());
					VectNextNodetmp.push_back(NextNodetmp);
					NextNodetmp.index=PtrNodes->NodeCorner[nodeIdlocal].Get_connect()[2];
					NextNodetmp.distance=DistToWall(PtrNodes->NodeCorner[PtrNodes->NodeIndexByType[NextNodetmp.index]].get_x(),PtrNodes->NodeCorner[PtrNodes->NodeIndexByType[NextNodetmp.index]].get_y());
					VectNextNodetmp.push_back(NextNodetmp);
				}
				else
				{
					NextNodetmp.index=PtrNodes->NodeCorner[nodeIdlocal].Get_connect()[3];
					NextNodetmp.distance=DistToWall(PtrNodes->NodeCorner[PtrNodes->NodeIndexByType[NextNodetmp.index]].get_x(),PtrNodes->NodeCorner[PtrNodes->NodeIndexByType[NextNodetmp.index]].get_y());
					VectNextNodetmp.push_back(NextNodetmp);
					NextNodetmp.index=PtrNodes->NodeCorner[nodeIdlocal].Get_connect()[4];
					NextNodetmp.distance=DistToWall(PtrNodes->NodeCorner[PtrNodes->NodeIndexByType[NextNodetmp.index]].get_x(),PtrNodes->NodeCorner[PtrNodes->NodeIndexByType[NextNodetmp.index]].get_y());
					VectNextNodetmp.push_back(NextNodetmp);
				}
				break;
			case 6:
				if(PtrNodes->NodeCorner[nodeIdlocal].Get_CornerType()==Concave)
				{
					NextNodetmp.index=PtrNodes->NodeCorner[nodeIdlocal].Get_connect()[2];
					NextNodetmp.distance=DistToWall(PtrNodes->NodeCorner[PtrNodes->NodeIndexByType[NextNodetmp.index]].get_x(),PtrNodes->NodeCorner[PtrNodes->NodeIndexByType[NextNodetmp.index]].get_y());
					VectNextNodetmp.push_back(NextNodetmp);
					NextNodetmp.index=PtrNodes->NodeCorner[nodeIdlocal].Get_connect()[3];
					NextNodetmp.distance=DistToWall(PtrNodes->NodeCorner[PtrNodes->NodeIndexByType[NextNodetmp.index]].get_x(),PtrNodes->NodeCorner[PtrNodes->NodeIndexByType[NextNodetmp.index]].get_y());
					VectNextNodetmp.push_back(NextNodetmp);
				}
				else
				{
					NextNodetmp.index=PtrNodes->NodeCorner[nodeIdlocal].Get_connect()[1];
					NextNodetmp.distance=DistToWall(PtrNodes->NodeCorner[PtrNodes->NodeIndexByType[NextNodetmp.index]].get_x(),PtrNodes->NodeCorner[PtrNodes->NodeIndexByType[NextNodetmp.index]].get_y());
					VectNextNodetmp.push_back(NextNodetmp);
					NextNodetmp.index=PtrNodes->NodeCorner[nodeIdlocal].Get_connect()[4];
					NextNodetmp.distance=DistToWall(PtrNodes->NodeCorner[PtrNodes->NodeIndexByType[NextNodetmp.index]].get_x(),PtrNodes->NodeCorner[PtrNodes->NodeIndexByType[NextNodetmp.index]].get_y());
					VectNextNodetmp.push_back(NextNodetmp);
				}
				break;
			case 7:
				if(PtrNodes->NodeCorner[nodeIdlocal].Get_CornerType()==Concave)
				{
					NextNodetmp.index=PtrNodes->NodeCorner[nodeIdlocal].Get_connect()[3];
					NextNodetmp.distance=DistToWall(PtrNodes->NodeCorner[PtrNodes->NodeIndexByType[NextNodetmp.index]].get_x(),PtrNodes->NodeCorner[PtrNodes->NodeIndexByType[NextNodetmp.index]].get_y());
					VectNextNodetmp.push_back(NextNodetmp);
					NextNodetmp.index=PtrNodes->NodeCorner[nodeIdlocal].Get_connect()[4];
					NextNodetmp.distance=DistToWall(PtrNodes->NodeCorner[PtrNodes->NodeIndexByType[NextNodetmp.index]].get_x(),PtrNodes->NodeCorner[PtrNodes->NodeIndexByType[NextNodetmp.index]].get_y());
					VectNextNodetmp.push_back(NextNodetmp);
				}
				else
				{
					NextNodetmp.index=PtrNodes->NodeCorner[nodeIdlocal].Get_connect()[1];
					NextNodetmp.distance=DistToWall(PtrNodes->NodeCorner[PtrNodes->NodeIndexByType[NextNodetmp.index]].get_x(),PtrNodes->NodeCorner[PtrNodes->NodeIndexByType[NextNodetmp.index]].get_y());
					VectNextNodetmp.push_back(NextNodetmp);
					NextNodetmp.index=PtrNodes->NodeCorner[nodeIdlocal].Get_connect()[2];
					NextNodetmp.distance=DistToWall(PtrNodes->NodeCorner[PtrNodes->NodeIndexByType[NextNodetmp.index]].get_x(),PtrNodes->NodeCorner[PtrNodes->NodeIndexByType[NextNodetmp.index]].get_y());
					VectNextNodetmp.push_back(NextNodetmp);
				}
				break;
			case 8:
				if(PtrNodes->NodeCorner[nodeIdlocal].Get_CornerType()==Concave)
				{
					NextNodetmp.index=PtrNodes->NodeCorner[nodeIdlocal].Get_connect()[1];
					NextNodetmp.distance=DistToWall(PtrNodes->NodeCorner[PtrNodes->NodeIndexByType[NextNodetmp.index]].get_x(),PtrNodes->NodeCorner[PtrNodes->NodeIndexByType[NextNodetmp.index]].get_y());
					VectNextNodetmp.push_back(NextNodetmp);
					NextNodetmp.index=PtrNodes->NodeCorner[nodeIdlocal].Get_connect()[4];
					NextNodetmp.distance=DistToWall(PtrNodes->NodeCorner[PtrNodes->NodeIndexByType[NextNodetmp.index]].get_x(),PtrNodes->NodeCorner[PtrNodes->NodeIndexByType[NextNodetmp.index]].get_y());
					VectNextNodetmp.push_back(NextNodetmp);
				}
				else
				{
					NextNodetmp.index=PtrNodes->NodeCorner[nodeIdlocal].Get_connect()[2];
					NextNodetmp.distance=DistToWall(PtrNodes->NodeCorner[PtrNodes->NodeIndexByType[NextNodetmp.index]].get_x(),PtrNodes->NodeCorner[PtrNodes->NodeIndexByType[NextNodetmp.index]].get_y());
					VectNextNodetmp.push_back(NextNodetmp);
					NextNodetmp.index=PtrNodes->NodeCorner[nodeIdlocal].Get_connect()[3];
					NextNodetmp.distance=DistToWall(PtrNodes->NodeCorner[PtrNodes->NodeIndexByType[NextNodetmp.index]].get_x(),PtrNodes->NodeCorner[PtrNodes->NodeIndexByType[NextNodetmp.index]].get_y());
					VectNextNodetmp.push_back(NextNodetmp);
				}
				break;
			}
			break;
			case ConcaveCorner:
				switch(PtrNodes->NodeCorner[nodeIdlocal].Get_BcNormal())
				{
				case 5:

					NextNodetmp.index=PtrNodes->NodeCorner[nodeIdlocal].Get_connect()[1];
					NextNodetmp.distance=DistToWall(PtrNodes->NodeCorner[PtrNodes->NodeIndexByType[NextNodetmp.index]].get_x(),PtrNodes->NodeCorner[PtrNodes->NodeIndexByType[NextNodetmp.index]].get_y());
					VectNextNodetmp.push_back(NextNodetmp);
					NextNodetmp.index=PtrNodes->NodeCorner[nodeIdlocal].Get_connect()[2];
					NextNodetmp.distance=DistToWall(PtrNodes->NodeCorner[PtrNodes->NodeIndexByType[NextNodetmp.index]].get_x(),PtrNodes->NodeCorner[PtrNodes->NodeIndexByType[NextNodetmp.index]].get_y());
					VectNextNodetmp.push_back(NextNodetmp);
					break;
				case 6:
					NextNodetmp.index=PtrNodes->NodeCorner[nodeIdlocal].Get_connect()[2];
					NextNodetmp.distance=DistToWall(PtrNodes->NodeCorner[PtrNodes->NodeIndexByType[NextNodetmp.index]].get_x(),PtrNodes->NodeCorner[PtrNodes->NodeIndexByType[NextNodetmp.index]].get_y());
					VectNextNodetmp.push_back(NextNodetmp);
					NextNodetmp.index=PtrNodes->NodeCorner[nodeIdlocal].Get_connect()[3];
					NextNodetmp.distance=DistToWall(PtrNodes->NodeCorner[PtrNodes->NodeIndexByType[NextNodetmp.index]].get_x(),PtrNodes->NodeCorner[PtrNodes->NodeIndexByType[NextNodetmp.index]].get_y());
					VectNextNodetmp.push_back(NextNodetmp);
					break;
				case 7:
					NextNodetmp.index=PtrNodes->NodeCorner[nodeIdlocal].Get_connect()[3];
					NextNodetmp.distance=DistToWall(PtrNodes->NodeCorner[PtrNodes->NodeIndexByType[NextNodetmp.index]].get_x(),PtrNodes->NodeCorner[PtrNodes->NodeIndexByType[NextNodetmp.index]].get_y());
					VectNextNodetmp.push_back(NextNodetmp);
					NextNodetmp.index=PtrNodes->NodeCorner[nodeIdlocal].Get_connect()[4];
					NextNodetmp.distance=DistToWall(PtrNodes->NodeCorner[PtrNodes->NodeIndexByType[NextNodetmp.index]].get_x(),PtrNodes->NodeCorner[PtrNodes->NodeIndexByType[NextNodetmp.index]].get_y());
					VectNextNodetmp.push_back(NextNodetmp);
					break;
				case 8:
					NextNodetmp.index=PtrNodes->NodeCorner[nodeIdlocal].Get_connect()[1];
					NextNodetmp.distance=DistToWall(PtrNodes->NodeCorner[PtrNodes->NodeIndexByType[NextNodetmp.index]].get_x(),PtrNodes->NodeCorner[PtrNodes->NodeIndexByType[NextNodetmp.index]].get_y());
					VectNextNodetmp.push_back(NextNodetmp);
					NextNodetmp.index=PtrNodes->NodeCorner[nodeIdlocal].Get_connect()[4];
					NextNodetmp.distance=DistToWall(PtrNodes->NodeCorner[PtrNodes->NodeIndexByType[NextNodetmp.index]].get_x(),PtrNodes->NodeCorner[PtrNodes->NodeIndexByType[NextNodetmp.index]].get_y());
					VectNextNodetmp.push_back(NextNodetmp);
					break;
				}
				break;
			case ConvexCorner:
				switch(PtrNodes->NodeCorner[nodeIdlocal].Get_BcNormal())
				{
				case 5:
					NextNodetmp.index=PtrNodes->NodeCorner[nodeIdlocal].Get_connect()[3];
					NextNodetmp.distance=DistToWall(PtrNodes->NodeCorner[PtrNodes->NodeIndexByType[NextNodetmp.index]].get_x(),PtrNodes->NodeCorner[PtrNodes->NodeIndexByType[NextNodetmp.index]].get_y());
					VectNextNodetmp.push_back(NextNodetmp);
					NextNodetmp.index=PtrNodes->NodeCorner[nodeIdlocal].Get_connect()[4];
					NextNodetmp.distance=DistToWall(PtrNodes->NodeCorner[PtrNodes->NodeIndexByType[NextNodetmp.index]].get_x(),PtrNodes->NodeCorner[PtrNodes->NodeIndexByType[NextNodetmp.index]].get_y());
					VectNextNodetmp.push_back(NextNodetmp);
					break;
				case 6:
					NextNodetmp.index=PtrNodes->NodeCorner[nodeIdlocal].Get_connect()[1];
					NextNodetmp.distance=DistToWall(PtrNodes->NodeCorner[PtrNodes->NodeIndexByType[NextNodetmp.index]].get_x(),PtrNodes->NodeCorner[PtrNodes->NodeIndexByType[NextNodetmp.index]].get_y());
					VectNextNodetmp.push_back(NextNodetmp);
					NextNodetmp.index=PtrNodes->NodeCorner[nodeIdlocal].Get_connect()[4];
					NextNodetmp.distance=DistToWall(PtrNodes->NodeCorner[PtrNodes->NodeIndexByType[NextNodetmp.index]].get_x(),PtrNodes->NodeCorner[PtrNodes->NodeIndexByType[NextNodetmp.index]].get_y());
					VectNextNodetmp.push_back(NextNodetmp);
					break;
				case 7:
					NextNodetmp.index=PtrNodes->NodeCorner[nodeIdlocal].Get_connect()[1];
					NextNodetmp.distance=DistToWall(PtrNodes->NodeCorner[PtrNodes->NodeIndexByType[NextNodetmp.index]].get_x(),PtrNodes->NodeCorner[PtrNodes->NodeIndexByType[NextNodetmp.index]].get_y());
					VectNextNodetmp.push_back(NextNodetmp);
					NextNodetmp.index=PtrNodes->NodeCorner[nodeIdlocal].Get_connect()[2];
					NextNodetmp.distance=DistToWall(PtrNodes->NodeCorner[PtrNodes->NodeIndexByType[NextNodetmp.index]].get_x(),PtrNodes->NodeCorner[PtrNodes->NodeIndexByType[NextNodetmp.index]].get_y());
					VectNextNodetmp.push_back(NextNodetmp);
					break;
				case 8:
					NextNodetmp.index=PtrNodes->NodeCorner[nodeIdlocal].Get_connect()[2];
					NextNodetmp.distance=DistToWall(PtrNodes->NodeCorner[PtrNodes->NodeIndexByType[NextNodetmp.index]].get_x(),PtrNodes->NodeCorner[PtrNodes->NodeIndexByType[NextNodetmp.index]].get_y());
					VectNextNodetmp.push_back(NextNodetmp);
					NextNodetmp.index=PtrNodes->NodeCorner[nodeIdlocal].Get_connect()[3];
					NextNodetmp.distance=DistToWall(PtrNodes->NodeCorner[PtrNodes->NodeIndexByType[NextNodetmp.index]].get_x(),PtrNodes->NodeCorner[PtrNodes->NodeIndexByType[NextNodetmp.index]].get_y());
					VectNextNodetmp.push_back(NextNodetmp);
					break;
				}
				break;
			case SolidGhost:
				break;
			case GlobalCorner:
				break;
		default:
			std::cerr<<"Node type not found for interpolation for next wall in interpolation. Node type is: "<<PtrNodes->TypeOfNode[nextback[i]]<<std::endl;
			break;
		}
	}
/*	if(next.size()>0)
	{
		for(int i=0;i<next.size();i++)
			if(PtrNodes->TypeOfNode[next[i]]!= Wall||PtrNodes->TypeOfNode[next[i]]!= Corner||PtrNodes->TypeOfNode[next[i]]!= ConcaveCorner||PtrNodes->TypeOfNode[next[i]]!= ConvexCorner)
				next.erase(next.begin()+i);
		if(SolidChecked.size()>0)
			for(int i=0;i<SolidChecked.size();i++)
			{
				itnextWall = find(next.begin(),next.end(), SolidChecked[i].index);

				if ( itnextWall != next.end() )
				  next.erase(itnextWall);
			}

	}*/
	if(VectNextNodetmp.size()>0)
	{
		for(int i=0;i<VectNextNodetmp.size();i++)
			if(PtrNodes->TypeOfNode[VectNextNodetmp[i].index]!= Wall||PtrNodes->TypeOfNode[VectNextNodetmp[i].index]!= Corner||PtrNodes->TypeOfNode[VectNextNodetmp[i].index]!= ConcaveCorner||PtrNodes->TypeOfNode[VectNextNodetmp[i].index]!= ConvexCorner)
				VectNextNodetmp.erase(VectNextNodetmp.begin()+i);
		if(SolidChecked.size()>0)
			for(int i=0;i<SolidChecked.size();i++)
			{
				itnextWall = find_if(VectNextNodetmp.begin(),VectNextNodetmp.end(), MatchNextNode(SolidChecked[i].index));

				if ( itnextWall != VectNextNodetmp.end() )
					VectNextNodetmp.erase(itnextWall);
					//VectNextNodetmp.erase(VectNextNodetmp.begin()+i);
			}

		std::sort(VectNextNodetmp.begin(), VectNextNodetmp.end(), by_distance());
		for(int i=0;i<VectNextNodetmp.size();i++)
			next.push_back(VectNextNodetmp[i].index);
	}
}
void InterpolationLinearLeastSquare::Next_FluidId(NodeArrays2D *PtrNodes,std::vector<NextNode> nextbak,std::vector<NextNode>& next){
	next.clear();int nodeIdlocal=0;
	NextNode NextNodetmp;
	for(int i=0;i<nextbak.size();i++)
	{
		nodeIdlocal=PtrNodes->NodeIndexByType[nextbak[i].index];
		switch (PtrNodes->TypeOfNode[nextbak[i].index])
		{
		case Interior:
			for(int i=1;i<9;i++)
				if(PtrNodes->TypeOfNode[PtrNodes->NodeInterior[nodeIdlocal].Get_connect()[i]]==Interior)
				{
					NextNodetmp.index=PtrNodes->NodeInterior[nodeIdlocal].Get_connect()[i];
					NextNodetmp.distance=DistToWall(PtrNodes->NodeInterior[PtrNodes->NodeIndexByType[NextNodetmp.index]].get_x(),PtrNodes->NodeInterior[PtrNodes->NodeIndexByType[NextNodetmp.index]].get_y());
					next.push_back(NextNodetmp);
				}

			break;
		}
	}
	//Remove node already checked. To speed up, it is split in two cases.
	if(next.size()>0)
	{
		if(FluidChecked.size()<next.size())
			for(int i=0;i<FluidChecked.size();i++)
			{
				itnextFluid = find_if(next.begin(),next.end(), MatchNextNode(FluidChecked[i].index));

				if ( itnextFluid != next.end() )
				  next.erase(itnextFluid);
			}
		else
			for(int i=0;i<next.size();i++)
			{
				itnextFluid = find_if(FluidChecked.begin(),FluidChecked.end(), MatchNextNode(next[i].index));

				if ( itnextFluid != FluidChecked.end() )
					next.erase(next.begin()+i);
			}

		std::sort(next.begin(), next.end(), by_distance());
	}

}
void InterpolationLinearLeastSquare::InterpolationOnWall(double *Var, int * Connect, int & normal){
	CalculLeastSquareMethod(Var,MapWallId[Connect[0]]);
	Var[Connect[0]]=result1;
}
void InterpolationLinearLeastSquare::InterpolationOnCornerConcave (double *Var, int * Connect, int & normal){
	CalculLeastSquareMethod(Var,MapWallId[Connect[0]]);
	Var[Connect[0]]=result1;
}
void InterpolationLinearLeastSquare::InterpolationOnCornerConvex (double *Var, int * Connect, int & normal){
	CalculLeastSquareMethod(Var,MapWallId[Connect[0]]);
	Var[Connect[0]]=result1;
}
void InterpolationLinearLeastSquare::InterpolationOnWall(double *Var1, double *Var2, int * Connect, int & normal){
	CalculLeastSquareMethod(Var1,Var2,MapWallId[Connect[0]]);
	Var1[Connect[0]]=result1;
	Var2[Connect[0]]=result2;
}
void InterpolationLinearLeastSquare::InterpolationOnCornerConcave (double *Var1, double *Var2, int * Connect, int & normal){
	CalculLeastSquareMethod(Var1,Var2,MapWallId[Connect[0]]);
	Var1[Connect[0]]=result1;
	Var2[Connect[0]]=result2;
}
void InterpolationLinearLeastSquare::InterpolationOnCornerConvex (double *Var1, double *Var2, int * Connect, int & normal){
	CalculLeastSquareMethod(Var1,Var2,MapWallId[Connect[0]]);
	Var1[Connect[0]]=result1;
	Var2[Connect[0]]=result2;
}
void InterpolationLinearLeastSquare::CalculLeastSquareMethod(double *Var,int wallId){
	sumX=0;	sumY=0;sumXY=0;sumX2=0;
	for(int i=0;i< fluidId[wallId].size();i++)
	{
		sumX+=fluidDist[wallId][i];
		sumX2+=fluidDist[wallId][i]*fluidDist[wallId][i];
		sumY+=Var[fluidId[wallId][i]];
		sumXY+=fluidDist[wallId][i]*Var[fluidId[wallId][i]];
	}
	for(int i=0;i< solidId[wallId].size();i++)
	{
		sumX+=solidDist[wallId][i];
		sumX2+=solidDist[wallId][i]*solidDist[wallId][i];
		sumY+=Var[solidId[wallId][i]];
		sumXY+=solidDist[wallId][i]*Var[solidId[wallId][i]];
	}
	result1=(sumY*sumX2-sumX*sumXY)/((fluidId[wallId].size()+solidId[wallId].size())*sumX2-sumX*sumX);

}
void InterpolationLinearLeastSquare::CalculLeastSquareMethod(double *Var1,double *Var2,int wallId){
	sumX=0;	sumY=0;sumXY=0;sumX2=0;
	for(int i=0;i< fluidId[wallId].size();i++)
	{
		sumX+=fluidDist[wallId][i];
		sumX2+=fluidDist[wallId][i]*fluidDist[wallId][i];
		sumY+=Var1[fluidId[wallId][i]];
		sumXY+=fluidDist[wallId][i]*Var1[fluidId[wallId][i]];
	}
	for(int i=0;i< solidId[wallId].size();i++)
	{
		sumX+=solidDist[wallId][i];
		sumX2+=solidDist[wallId][i]*solidDist[wallId][i];
		sumY+=Var1[solidId[wallId][i]];
		sumXY+=solidDist[wallId][i]*Var1[solidId[wallId][i]];
	}
	result1=(sumY*sumX2-sumX*sumXY)/((fluidId[wallId].size()+solidId[wallId].size())*sumX2-sumX*sumX);
	sumY=0;sumXY=0;

	for(int i=0;i< fluidId[wallId].size();i++)
	{
		sumY+=Var2[fluidId[wallId][i]];
		sumXY+=fluidDist[wallId][i]*Var2[fluidId[wallId][i]];
	}
	for(int i=0;i< solidId[wallId].size();i++)
	{
		sumY+=Var2[solidId[wallId][i]];
		sumXY+=solidDist[wallId][i]*Var2[solidId[wallId][i]];
	}
	result2=(sumY*sumX2-sumX*sumXY)/((fluidId[wallId].size()+solidId[wallId].size())*sumX2-sumX*sumX);
}
