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
	Var2[Connect[0]]=0.5*(Var2[Connect[PtrOppositeInterpol[normal]]]+Var2[Connect[normal]]);}
InterpolationLinearLeastSquare::InterpolationLinearLeastSquare(){
	sumX=0;	sumY=0;sumXY=0;sumX2=0;
	solidId=0;fluidId=0;MapWallId=0;
	nNodes=0;nNodes2=0;nWalls=0;
	solidDist=0;fluidDist=0;
	x0=0;y0=0;
	result1=0;result2=0;
}
InterpolationLinearLeastSquare::InterpolationLinearLeastSquare(int dimension_,int nb_vel,unsigned int *PtrOppositeInterpol_){
	dimension=dimension_;
	nb_Vel=nb_vel;
	PtrOppositeInterpol=PtrOppositeInterpol_;
	sumX=0;	sumY=0;sumXY=0;sumX2=0;
	solidId=0;fluidId=0;MapWallId=0;
	nNodes=0;nNodes2=0;	nWalls=0;
	solidDist=0;fluidDist=0;
	x0=0;y0=0;
	result1=0;result2=0;
}
InterpolationLinearLeastSquare::~InterpolationLinearLeastSquare(){
	for (int i=0;i<nWalls;i++)
	{
		delete []solidId[i];delete []fluidId[i];
		delete []solidDist[i];delete []fluidDist[i];
	}
	delete []solidId;delete []fluidId;delete []solidDist;delete []fluidDist;delete []MapWallId;
}
void InterpolationLinearLeastSquare::InitInterpol(NodeArrays2D *PtrNodes, Parameters *PtrParam){
	MapWallId=new short int [PtrNodes->TypeOfNode.size()];
	nWalls=PtrNodes->NodeWall.size()+PtrNodes->CornerConcave.size()+PtrNodes->CornerConvex.size();
	nNodes=4;
	nNodes2=2*nNodes;
	std::vector<int> next;
	fluidId=new short int* [nWalls];solidId=new short int* [nWalls];
	fluidDist=new double* [nWalls];solidDist=new double* [nWalls];
	for (int i=0;i<nWalls;i++)
	{
		fluidId[i]=new short int[nNodes];solidId[i]=new short int[nNodes];
		fluidDist[i]=new double[nNodes];solidDist[i]=new double[nNodes];
	}
	int countfluid=0;int countsolid=0;
	int count=0;
	for (int i=0;i<PtrNodes->NodeWall.size();i++)
	{
		countfluid=0;countsolid=0;
		next.push_back(PtrNodes->NodeWall[i].Get_index());
		MapWallId[next[0]]=i;
		x0=PtrNodes->NodeWall[i].get_x();y0=PtrNodes->NodeWall[i].get_y();
		while(countfluid<nNodes || countsolid<nNodes)
		{
			Mark_FluidSolidIds(PtrNodes,next[count],i,countfluid,countsolid);
			count++;
			if(count>=next.size())
			{
				Next_WallId(PtrNodes,PtrNodes->NodeWall[i].Get_index(),next);
				count=0;
			}
		}
	}
}
void InterpolationLinearLeastSquare::Mark_FluidSolidIds(NodeArrays2D *PtrNodes,int nodeId, int &countnwalls, int &countfluid, int &countsolid){

		int nodeIdlocal=PtrNodes->NodeIndexByType[nodeId];
		int nodeIdlocaltmp=0;
		switch (PtrNodes->TypeOfNode[nodeId])
		{
		case Wall:
			if(countfluid<nNodes)
			{
				fluidId[countnwalls][countfluid]=(PtrNodes->NodeWall[nodeIdlocal].Get_connect()[PtrNodes->NodeWall[nodeIdlocal].Get_BcNormal()]);
				nodeIdlocaltmp=PtrNodes->NodeIndexByType[fluidId[countnwalls][countfluid]];
				fluidDist[countnwalls][countfluid]=DistToWall(PtrNodes->NodeInterior[nodeIdlocaltmp].get_x(),PtrNodes->NodeInterior[nodeIdlocaltmp].get_y());
			}
			if(countsolid<nNodes)
			{
				solidId[countnwalls][countsolid]=(PtrNodes->NodeWall[nodeIdlocal].Get_connect()[PtrOppositeInterpol[PtrNodes->NodeWall[nodeIdlocal].Get_BcNormal()]]);
				nodeIdlocaltmp=PtrNodes->NodeIndexByType[solidId[countnwalls][countsolid]];
				solidDist[countnwalls][countsolid]=-DistToWall(PtrNodes->NodeSolid[nodeIdlocal].get_x(),PtrNodes->NodeSolid[nodeIdlocal].get_y());
			}
			break;
		case Corner:
			if(countfluid<nNodes)
			{
				fluidId[countnwalls][countfluid]=(PtrNodes->NodeCorner[nodeIdlocal].Get_connect()[PtrNodes->NodeCorner[nodeIdlocal].Get_BcNormal()]);
				nodeIdlocaltmp=PtrNodes->NodeIndexByType[fluidId[countnwalls][countfluid]];
				fluidDist[countnwalls][countfluid]=DistToWall(PtrNodes->NodeInterior[nodeIdlocaltmp].get_x(),PtrNodes->NodeInterior[nodeIdlocaltmp].get_y());
			}
			if(countsolid<nNodes)
			{
				solidId[countnwalls][countsolid]=(PtrNodes->NodeCorner[nodeIdlocal].Get_connect()[PtrOppositeInterpol[PtrNodes->NodeCorner[nodeIdlocal].Get_BcNormal()]]);
				nodeIdlocaltmp=PtrNodes->NodeIndexByType[solidId[countnwalls][countsolid]];
				solidDist[countnwalls][countsolid]=-DistToWall(PtrNodes->NodeSolid[nodeIdlocal].get_x(),PtrNodes->NodeSolid[nodeIdlocal].get_y());
			}
			break;
		case SolidGhost:
			break;
		case GlobalCorner:
			break;
		case SpecialWall:
			if(countfluid<nNodes)
			{
				fluidId[countnwalls][countfluid]=(PtrNodes->NodeSpecialWall[nodeIdlocal].Get_connect()[PtrNodes->NodeSpecialWall[nodeIdlocal].Get_BcNormal()]);
				nodeIdlocaltmp=PtrNodes->NodeIndexByType[fluidId[countnwalls][countfluid]];
				fluidDist[countnwalls][countfluid]=DistToWall(PtrNodes->NodeInterior[nodeIdlocaltmp].get_x(),PtrNodes->NodeInterior[nodeIdlocaltmp].get_y());
			}
			if(countsolid<nNodes)
			{
				solidId[countnwalls][countsolid]=(PtrNodes->NodeSpecialWall[nodeIdlocal].Get_connect()[PtrOppositeInterpol[PtrNodes->NodeSpecialWall[nodeIdlocal].Get_BcNormal()]]);
				nodeIdlocaltmp=PtrNodes->NodeIndexByType[solidId[countnwalls][countsolid]];
				solidDist[countnwalls][countsolid]=-DistToWall(PtrNodes->NodeSolid[nodeIdlocal].get_x(),PtrNodes->NodeSolid[nodeIdlocal].get_y());
			}
			break;
		default:
			std::cerr<<"Node type not found for interpolation. Type of Node :"<<PtrNodes->TypeOfNode[nodeId]<<std::endl;
			break;
		}
		if(countfluid<nNodes)
			countfluid++;
		if(countsolid<nNodes)
			countsolid++;

}
void InterpolationLinearLeastSquare::Next_WallId(NodeArrays2D *PtrNodes,int nodeId,std::vector<int>& next){

	int nodeIdlocal=PtrNodes->NodeIndexByType[nodeId];
	next.clear();
	switch (PtrNodes->TypeOfNode[nodeId])
	{
	case Wall:
		switch(PtrNodes->NodeWall[nodeIdlocal].Get_BcNormal())
		{
		case 1:
			next.push_back(PtrNodes->NodeWall[nodeIdlocal].Get_connect()[2]);
			next.push_back(PtrNodes->NodeWall[nodeIdlocal].Get_connect()[4]);
			break;
		case 2:
			next.push_back(PtrNodes->NodeWall[nodeIdlocal].Get_connect()[1]);
			next.push_back(PtrNodes->NodeWall[nodeIdlocal].Get_connect()[3]);
			break;
		case 3:
			next.push_back(PtrNodes->NodeWall[nodeIdlocal].Get_connect()[2]);
			next.push_back(PtrNodes->NodeWall[nodeIdlocal].Get_connect()[4]);
			break;
		case 4:
			next.push_back(PtrNodes->NodeWall[nodeIdlocal].Get_connect()[1]);
			next.push_back(PtrNodes->NodeWall[nodeIdlocal].Get_connect()[3]);
			break;
		}
		break;
	case Corner:
		switch(PtrNodes->NodeCorner[nodeIdlocal].Get_BcNormal())
		{
		case 5:
			if(PtrNodes->NodeCorner[nodeIdlocal].Get_CornerType()==Concave)
			{
				next.push_back(PtrNodes->NodeCorner[nodeIdlocal].Get_connect()[1]);
				next.push_back(PtrNodes->NodeCorner[nodeIdlocal].Get_connect()[2]);
			}
			else
			{
				next.push_back(PtrNodes->NodeCorner[nodeIdlocal].Get_connect()[3]);
				next.push_back(PtrNodes->NodeCorner[nodeIdlocal].Get_connect()[4]);
			}
			break;
		case 6:
			if(PtrNodes->NodeCorner[nodeIdlocal].Get_CornerType()==Concave)
			{
				next.push_back(PtrNodes->NodeCorner[nodeIdlocal].Get_connect()[2]);
				next.push_back(PtrNodes->NodeCorner[nodeIdlocal].Get_connect()[3]);
			}
			else
			{
				next.push_back(PtrNodes->NodeCorner[nodeIdlocal].Get_connect()[1]);
				next.push_back(PtrNodes->NodeCorner[nodeIdlocal].Get_connect()[4]);
			}
			break;
		case 7:
			if(PtrNodes->NodeCorner[nodeIdlocal].Get_CornerType()==Concave)
			{
				next.push_back(PtrNodes->NodeCorner[nodeIdlocal].Get_connect()[3]);
				next.push_back(PtrNodes->NodeCorner[nodeIdlocal].Get_connect()[4]);
			}
			else
			{
				next.push_back(PtrNodes->NodeCorner[nodeIdlocal].Get_connect()[1]);
				next.push_back(PtrNodes->NodeCorner[nodeIdlocal].Get_connect()[2]);
			}
			break;
		case 8:
			if(PtrNodes->NodeCorner[nodeIdlocal].Get_CornerType()==Concave)
			{
				next.push_back(PtrNodes->NodeCorner[nodeIdlocal].Get_connect()[1]);
				next.push_back(PtrNodes->NodeCorner[nodeIdlocal].Get_connect()[4]);
			}
			else
			{
				next.push_back(PtrNodes->NodeCorner[nodeIdlocal].Get_connect()[2]);
				next.push_back(PtrNodes->NodeCorner[nodeIdlocal].Get_connect()[3]);
			}
			break;
		}
		break;
	default:
		std::cerr<<"Node type not found for interpolation"<<std::endl;
		break;
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
	for(int i=0;i<nNodes;i++)
	{
		sumX+=fluidDist[wallId][i]+solidDist[wallId][i];
		sumX2+=fluidDist[wallId][i]*fluidDist[wallId][i]+solidDist[wallId][i]*solidDist[wallId][i];
		sumY+=Var[fluidId[wallId][i]]+Var[solidId[wallId][i]];
		sumXY+=fluidDist[wallId][i]*Var[fluidId[wallId][i]]+solidDist[wallId][i]*Var[solidId[wallId][i]];
	}
	result1=(sumY*sumX2-sumX*sumXY)/(nNodes2*sumX2-sumX*sumX);
}
void InterpolationLinearLeastSquare::CalculLeastSquareMethod(double *Var1,double *Var2,int wallId){
	sumX=0;	sumY=0;sumXY=0;sumX2=0;
	for(int i=0;i<nNodes;i++)
	{
		sumX+=fluidDist[wallId][i]+solidDist[wallId][i];
		sumX2+=fluidDist[wallId][i]*fluidDist[wallId][i]+solidDist[wallId][i]*solidDist[wallId][i];
		sumY+=Var1[fluidId[wallId][i]]+Var1[solidId[wallId][i]];
		sumXY+=fluidDist[wallId][i]*Var1[fluidId[wallId][i]]+solidDist[wallId][i]*Var1[solidId[wallId][i]];
	}
	result1=(sumY*sumX2-sumX*sumXY)/(nNodes2*sumX2-sumX*sumX);
	for(int i=0;i<nNodes;i++)
	{
		sumY+=Var2[fluidId[wallId][i]]+Var2[solidId[wallId][i]];
		sumXY+=fluidDist[wallId][i]*Var2[fluidId[wallId][i]]+solidDist[wallId][i]*Var2[solidId[wallId][i]];
	}
	result2=(sumY*sumX2-sumX*sumXY)/(nNodes2*sumX2-sumX*sumX);
}
