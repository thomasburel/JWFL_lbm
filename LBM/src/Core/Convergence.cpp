/*
 * ============================================================================
 * Convergence.cpp
 *
 *  Created on: 30 Sep 2016
 *      Author: Thomas Burel
 *      Version     :
 *      Copyright   : Your copyright notice
 *      Description : LBM Code Developing by JWFL at University of Strathclyde
 * ============================================================================
 */

#include "Convergence.h"

Convergence::Convergence() {
	PtrMultiBlockConv=0;
	PtrDicConv=0;
	PtrParmConv=0;
	PtrNodeArraysConv=0;
	PtrViscosityConv=0;
	PtrPatchBcConv=0;
	Error=0;
	Error_tmp=0;
	Error_sum=0;
	Error_avg=0;
	Sum_Current=1;
	avg=0;porosity=0;

	Scalar_last=0;
	Vector_last=0;
	NbNodes=0;
	Scalar_CurrentTime=0;
	Vector_CurrentTime=0;


	RhoNProductionRate=0;
	UPerm=0;RhoPerm=0;RhoNPerm=0;
	NodeId=0;
	NodeIdSpeWall=0;
	NodeIdGloCorner=0;

	LengthMedium=100;
	SectionMedium=40;

}

Convergence::~Convergence() {
	// TODO Auto-generated destructor stub
}
void Convergence::Set_Convergence(){
	//Temporary
	switch(PtrParmConv->Get_ErrorVariable())
	{
	case SolverEnum::Density:
		PtrDicConv->Get_PtrVar("Density",Scalar_CurrentTime);
		break;
	case SolverEnum::RhoN:
		PtrDicConv->Get_PtrVar("RhoN",Scalar_CurrentTime);
		break;
	case SolverEnum::VelocityX:
		PtrDicConv->Get_PtrVar("VelocityX",Scalar_CurrentTime);
		break;
	case SolverEnum::VelocityY:
		PtrDicConv->Get_PtrVar("VelocityY",Scalar_CurrentTime);
		break;
	default:
		if(PtrMultiBlockConv->IsMainProcessor())
			std::cerr<<"Error variable type not found. Density will be used."<<std::endl;
		PtrDicConv->Get_PtrVar("Density",Scalar_CurrentTime);
	}
	NbNodes=PtrDicConv->Get_NbNodes();
	PtrDicConv->AddVar(Scalar,"ScalarError_last",false,true,false,Scalar_last);
	for(int i=0;i<NbNodes;i++)
	{
		Scalar_last[i]=1.0;
	}


//Porous media cases
	if(PtrParmConv->IsPorousMediaCase())
	{
		if(PtrParmConv->IsCalculatePorosity())
		{
			porosity=Porosity();
				if(PtrMultiBlockConv->IsMainProcessor())
					std::cout<<"Porosity of the domain is: "<<porosity<<std::endl;
		}
		if(PtrParmConv->IsCalculateProductionRate())
		{
			if(OutletPatchId.empty())
			{
				for(int i=0;i<PtrPatchBcConv->Get_NumberOfPatchBc();i++)
					if(PtrPatchBcConv->IsOutlet(i))
						OutletPatchId.push_back(i);
			}
			if(OutletPatchId.size()>0)
			{
				PtrDicConv->Get_PtrVar("RhoN",RhoNProductionRate);
				if(PtrMultiBlockConv->IsMainProcessor())
				{
					ofstream ProductionRatefile;
					ProductionRatefile.open("ProductionRate.txt",ios::out | ios::trunc);
					ProductionRatefile.close();
				}
			}
			else
			{
				if(PtrMultiBlockConv->IsMainProcessor())
					std::cerr<<"No outlet Patch is set. No calculation of Production Rate will be done."<<std::endl;
				PtrParmConv->CalculateProductionRate(false);
			}
		}
		if(PtrParmConv->IsCalculatePermeability())
		{
			if(OutletPatchId.empty())
			{
				for(int i=0;i<PtrPatchBcConv->Get_NumberOfPatchBc();i++)
					if(PtrPatchBcConv->IsOutlet(i))
						OutletPatchId.push_back(i);
			}
			if(InletPatchId.empty())
			{
				for(int i=0;i<PtrPatchBcConv->Get_NumberOfPatchBc();i++)
					if(PtrPatchBcConv->IsInlet(i))
						InletPatchId.push_back(i);
			}
			if(OutletPatchId.size()>0 && InletPatchId.size()>0)
			{
				PtrDicConv->Get_PtrVar("RhoN",RhoNPerm);
				PtrDicConv->Get_PtrVar("Density",RhoPerm);
				UPerm=new double*[2];
				PtrDicConv->Get_PtrVar("VelocityX",UPerm[0]);
				PtrDicConv->Get_PtrVar("VelocityY",UPerm[1]);

				if(PtrMultiBlockConv->IsMainProcessor())
				{
					ofstream Permeabilityfile;
					Permeabilityfile.open("Permeability.txt",ios::out | ios::trunc);
					Permeabilityfile.close();
				}
			}
			else
			{
				if(PtrMultiBlockConv->IsMainProcessor())
					std::cerr<<"No inlet and/or outlet Patch is/are set. No calculation of Permeability will be done."<<std::endl;
				PtrParmConv->CalculatePermeability(false);
			}
		}
		if(!(PtrParmConv->IsCalculatePermeability() || PtrParmConv->IsCalculateProductionRate()))
		{
			std::cout<<"Porous media case detected but no Production rate or Permeability requested."<<std::endl;
			PtrParmConv->PorousMediaCase(false);
		}
	}
}
void Convergence::Calcul_Error(int &Time){

	Calcul_Error_ScalarField();
	if(PtrParmConv->IsPorousMediaCase())
		Calcul_PorousMediaConvergence(Time);
}
void Convergence::Calcul_Error_ScalarField(){
	Calcul_Error_ScalarFieldInOneProc();
	Error_tmp=PtrMultiBlockConv->SumAllProcessors(&Error_tmp);
	Sum_Current=PtrMultiBlockConv->SumAllProcessors(&Sum_Current);
	Error=Error_tmp/(Sum_Current+1e-15);
}
double Convergence::Calcul_Error_ScalarFieldInOneProc(){
	Sum_Current=0;Error_tmp=0;
	for(int i=0;i<NbNodes;i++)
	{
		Error_tmp+=std::abs(Scalar_CurrentTime[i]-Scalar_last[i]);
		Sum_Current+=std::abs(Scalar_CurrentTime[i]);
		Scalar_last[i]=Scalar_CurrentTime[i];//Save the new value
	}
//	Error_tmp=Error_tmp/(Sum_Current+1e-15);//std::abs(Scalar_CurrentTime[idx]-Scalar_last[idx])/(std::abs(Scalar_CurrentTime[idx])+1e-15);

	return Error_tmp;
}
void Convergence::Calcul_PorousMediaConvergence(int &Time){
	if(PtrParmConv->IsCalculateProductionRate())
		Calcul_ProductionRate(Time);
	if(PtrParmConv->IsCalculatePermeability())
		Calcul_DarcyPermeability_SinglePhase(Time);
}
///Calcul Production Rate and write it in a file
double Convergence::Calcul_ProductionRate(int &Time){

	double ProductionRate=0;
	double sum=0,sumtmp=0,nbsumnodes=0;
	for (int i=0;i<OutletPatchId.size();i++)
	{
		if(PtrPatchBcConv->Get_NodeIndex(OutletPatchId[i],NodeId,NodeIdSpeWall,NodeIdGloCorner))
		{
			Sum_ScalarNodeArray(*NodeId,*NodeIdSpeWall,*NodeIdGloCorner,RhoNProductionRate,sumtmp);
			sum+=sumtmp;
			nbsumnodes+=NodeId->size()+0.5*(NodeIdSpeWall->size()+NodeIdGloCorner->size());
		}
	}
	sum=PtrMultiBlockConv->SumAllProcessors(&sum);
	nbsumnodes=PtrMultiBlockConv->SumAllProcessors(&nbsumnodes);
	ProductionRate=(sum/nbsumnodes+1.0)*0.5;
	if(PtrMultiBlockConv->IsMainProcessor())
	{
		ofstream ProductionRatefile;
		ProductionRatefile.open("ProductionRate.txt",ios::out | ios::app);
		ProductionRatefile<<Time<<'\t'<<ProductionRate<<std::endl;
		ProductionRatefile.close();
		std::cout<<"Production Rate is: "<<ProductionRate<<std::endl;
	}
	return ProductionRate;
}
double Convergence::Calcul_DarcyPermeability_SinglePhase(int &Time){
	double Permeability=0;
	double sumtmp=0,nbsumnodesIn,nbsumnodesOut=0;
	double sumUIn=0,sumUOut=0,sumRhoIn=0,sumRhoOut=0,sumRhoNIn=0,sumRhoNOut=0;
	for (int i=0;i<InletPatchId.size();i++)
	{
		if(PtrPatchBcConv->Get_NodeIndex(InletPatchId[i],NodeId,NodeIdSpeWall,NodeIdGloCorner))
		{
			Sum_ScalarNodeArray(*NodeId,*NodeIdSpeWall,*NodeIdGloCorner,RhoPerm,sumtmp);
			sumRhoIn=sumtmp;
			Sum_ScalarNodeArray(*NodeId,*NodeIdSpeWall,*NodeIdGloCorner,UPerm[PtrPatchBcConv->Get_Orientation(InletPatchId[i])],sumtmp);
			sumUIn+=sumtmp;
			Sum_ScalarNodeArray(*NodeId,*NodeIdSpeWall,*NodeIdGloCorner,RhoNPerm,sumtmp);
			sumRhoNIn+=sumtmp;
			nbsumnodesIn+=NodeId->size()+0.5*(NodeIdSpeWall->size()+NodeIdGloCorner->size());
		}
	}
	sumRhoIn=PtrMultiBlockConv->SumAllProcessors(&sumRhoIn);
	sumRhoNIn=PtrMultiBlockConv->SumAllProcessors(&sumRhoNIn);
	sumUIn=PtrMultiBlockConv->SumAllProcessors(&sumUIn);
	nbsumnodesIn=PtrMultiBlockConv->SumAllProcessors(&nbsumnodesIn);
	for (int i=0;i<OutletPatchId.size();i++)
	{
		if(PtrPatchBcConv->Get_NodeIndex(OutletPatchId[i],NodeId,NodeIdSpeWall,NodeIdGloCorner))
		{
			Sum_ScalarNodeArray(*NodeId,*NodeIdSpeWall,*NodeIdGloCorner,RhoPerm,sumtmp);
			sumRhoOut+=sumtmp;
			Sum_ScalarNodeArray(*NodeId,*NodeIdSpeWall,*NodeIdGloCorner,UPerm[PtrPatchBcConv->Get_Orientation(OutletPatchId[i])],sumtmp);
			sumUOut+=sumtmp;
			Sum_ScalarNodeArray(*NodeId,*NodeIdSpeWall,*NodeIdGloCorner,RhoNPerm,sumtmp);
			sumRhoNOut+=sumtmp;
			nbsumnodesOut+=NodeId->size()+0.5*(NodeIdSpeWall->size()+NodeIdGloCorner->size());
		}
	}
	sumRhoOut=PtrMultiBlockConv->SumAllProcessors(&sumRhoOut);
	sumRhoNOut=PtrMultiBlockConv->SumAllProcessors(&sumRhoNOut);
	sumUOut=PtrMultiBlockConv->SumAllProcessors(&sumUOut);
	nbsumnodesOut=PtrMultiBlockConv->SumAllProcessors(&nbsumnodesOut);

	double avgRhoIn=sumRhoIn/nbsumnodesIn;
	double avgRhoOut=sumRhoOut/nbsumnodesOut;
	double avgU=(sumUIn/nbsumnodesIn+sumUOut/nbsumnodesOut)*0.5;
	double avgRhoN=(sumRhoNIn/nbsumnodesIn+sumRhoNOut/nbsumnodesOut)*0.5;
	double avgRho=(avgRhoIn+avgRhoOut)*0.5;
	double deltaP=avgRhoOut-avgRhoIn;
	double avgViscoNu=PtrViscosityConv->Get_Nu(avgRho,avgRhoN);
	if(std::abs(deltaP)>0)
		Permeability=LengthMedium*avgViscoNu*avgRho*avgU/(SectionMedium*std::abs(deltaP));
	else
		Permeability=1e10;

	if(PtrMultiBlockConv->IsMainProcessor())
	{
		ofstream Permeabilityfile;
		Permeabilityfile.open("Permeability.txt",ios::out | ios::app);
		Permeabilityfile<<Time<<'\t'<<Permeability<<'\t'<<avgU<<'\t'<<deltaP<<'\t'<<avgViscoNu<<'\t'<<avgRho<<'\t'<<avgRhoIn<<'\t'<<avgRhoOut<<std::endl;
		Permeabilityfile.close();
		std::cout<<"Permeability is: "<<Permeability<<
				" Average Velocity: "<<avgU<<" Delta P: "<<std::abs(deltaP)
		<<" Viscosity: "<<avgViscoNu<<" Rho: "<<avgRho
		<<" Rho In: "<<avgRhoIn<<" Rho Out: "<<avgRhoOut<<std::endl;
	}
	return Permeability;

}
void Convergence::Avg_ScalarNodeArray(std::vector<int>& NodeArray,std::vector<int>& NodeArraySpecialWall,std::vector<int>& NodeArrayGlobalCorner,double *&Var1,double &sum){
	Sum_ScalarNodeArray(NodeArray,NodeArraySpecialWall,NodeArrayGlobalCorner,Var1,sum);
	sum/=(NodeArray.size()+NodeArraySpecialWall.size()+NodeArrayGlobalCorner.size());
}
void Convergence::Avg_VectorNodeArray(std::vector<int>& NodeArray,std::vector<int>& NodeArraySpecialWall,std::vector<int>& NodeArrayGlobalCorner,double *&Var1,double *&Var2,double &sum1,double &sum2){
	Sum_VectorNodeArray(NodeArray,NodeArraySpecialWall,NodeArrayGlobalCorner,Var1,Var2,sum1,sum2);
	int nbNodes=NodeArray.size()+NodeArraySpecialWall.size()+NodeArrayGlobalCorner.size();
	sum1/=nbNodes;sum2/=nbNodes;
}
void Convergence::Sum_ScalarNodeArray(std::vector<int>& NodeArray,std::vector<int>& NodeArraySpecialWall,std::vector<int>& NodeArrayGlobalCorner,double *&Var1,double &sum){
	sum=0;
	double sum_tmp=0;
	for(int i=0;i<NodeArray.size();i++)
		sum+=Var1[NodeArray[i]];
	for(int i=0;i<NodeArraySpecialWall.size();i++)
		sum_tmp+=Var1[NodeArraySpecialWall[i]];
	for(int i=0;i<NodeArrayGlobalCorner.size();i++)
		sum_tmp+=Var1[NodeArrayGlobalCorner[i]];
	sum_tmp/=2.0;
	sum+=sum_tmp;
}
void Convergence::Sum_VectorNodeArray(std::vector<int>& NodeArray,std::vector<int>& NodeArraySpecialWall,std::vector<int>& NodeArrayGlobalCorner,double *&Var1,double *&Var2,double &sum1,double &sum2){
	sum1=0;sum2=0;
	for(int i=0;i<NodeArray.size();i++)
		{sum1+=Var1[NodeArray[i]];sum2+=Var2[NodeArray[i]];}
	double sum_tmp1=0;double sum_tmp2=0;
	for(int i=0;i<NodeArraySpecialWall.size();i++)
		{sum_tmp1+=Var1[NodeArraySpecialWall[i]];sum_tmp2+=Var2[NodeArraySpecialWall[i]];}
	for(int i=0;i<NodeArrayGlobalCorner.size();i++)
		{sum_tmp1+=Var1[NodeArrayGlobalCorner[i]];sum_tmp2+=Var2[NodeArrayGlobalCorner[i]];}
	sum_tmp1/=2.0;sum_tmp2/=2.0;
	sum1+=sum_tmp1;sum2+=sum_tmp2;
}
double Convergence::Porosity(){
	double fluidVolumesum=0;double solidVolumesum=0;
	SumFluidVolume(fluidVolumesum);
	SumSolidVolume(solidVolumesum);
	fluidVolumesum=PtrMultiBlockConv->SumAllProcessors(&fluidVolumesum);
	solidVolumesum=PtrMultiBlockConv->SumAllProcessors(&solidVolumesum);
	return fluidVolumesum/(fluidVolumesum+solidVolumesum);
}
void Convergence::SumFluidVolume(double & sumfluidvolume){
	sumfluidvolume=0;
	sumfluidvolume+=PtrNodeArraysConv->Get_SizeNodeIdInterior();
	sumfluidvolume+=0.5*(PtrNodeArraysConv->Get_SizeNodeIdPressure()+PtrNodeArraysConv->Get_SizeNodeIdVelocity()+PtrNodeArraysConv->Get_SizeNodeIdPeriodic()+PtrNodeArraysConv->Get_SizeNodeIdSymmetry());
	sumfluidvolume+=0.5*PtrNodeArraysConv->Get_SizeNodeIdWall();
	sumfluidvolume+=0.25*(PtrNodeArraysConv->Get_SizeNodeIdGlobalCorner()+PtrNodeArraysConv->Get_SizeNodeIdCornerConcave()+PtrNodeArraysConv->Get_SizeNodeIdSpecialWall());
	sumfluidvolume+=0.75*PtrNodeArraysConv->Get_SizeNodeIdCornerConvex();
}
void Convergence::SumSolidVolume(double & sumsolidvolume){
	sumsolidvolume=0;

	for (int i=0;i<PtrNodeArraysConv->Get_SizeNodeIdSolid();i++)
		if(IsBoundary(PtrNodeArraysConv->Get_NodeIdSolid(i)))
			sumsolidvolume+=0.5;
		else
			sumsolidvolume++;
	for (int i=0;i<PtrNodeArraysConv->Get_SizeNodeIdCornerConvex();i++)
		if(!IsBoundary(PtrNodeArraysConv->Get_NodeIdCornerConvex(i)))
			sumsolidvolume+=0.25;
	for (int i=0;i<PtrNodeArraysConv->Get_SizeNodeIdWall();i++)
		if(!IsBoundary(PtrNodeArraysConv->Get_NodeIdWall(i)))
			sumsolidvolume+=0.5;
	for (int i=0;i<PtrNodeArraysConv->Get_SizeNodeIdSpecialWall();i++)
		if(!IsGlobalCornerASpecialWall(PtrNodeArraysConv->Get_NodeIdSpecialWall(i)))
			sumsolidvolume+=0.25;
	for (int i=0;i<PtrNodeArraysConv->Get_SizeNodeIdCornerConcave();i++)
		if(IsBoundary(PtrNodeArraysConv->Get_NodeIdCornerConcave(i)))
		{
			if(IsGlobalCornerASpecialWall(!PtrNodeArraysConv->Get_NodeIdCornerConcave(i)))
				sumsolidvolume+=0.25;
		}
		else
			sumsolidvolume+=0.75;
}
bool Convergence::IsBoundary(int idx){
	double x,y;
	PtrNodeArraysConv->Get_coordinate(idx,x,y);
	if(x==0 || x==PtrParmConv->Get_Nx()||y==0 || y==PtrParmConv->Get_Ny())
		return true;
	else
		return false;
}

bool Convergence::IsGlobalCornerASpecialWall(int idx){
	double x,y;
	PtrNodeArraysConv->Get_coordinate(idx,x,y);
	if((x==0 && y==0)|| (x==0 && y==PtrParmConv->Get_Ny()) ||(x==PtrParmConv->Get_Nx() && y==0) || (x==PtrParmConv->Get_Nx() && y==PtrParmConv->Get_Ny()))
		return true;
	else
		return false;
}
