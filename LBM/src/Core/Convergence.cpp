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
	PtrPermeabilityDarcy=0;
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
		if(PtrParmConv->Get_Model()==SolverEnum::SinglePhase)
		{
			PtrDicConv->Get_PtrVar("Density",Scalar_CurrentTime);
			if(PtrMultiBlockConv->IsMainProcessor())
				std::cout<<"Single Phase case: Normal density is not used. Density will be used for the convergence."<<std::endl;
		}
		else
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
		LengthMedium=PtrParmConv->Get_LenghtMedia();
		SectionMedium=PtrParmConv->Get_SectionMedia();
		PtrParmConv->Get_PositionMedia(minXMedia,maxXMedia,minYMedia,maxYMedia,minZMedia,maxZMedia);
		if(minXMedia<0) minXMedia=0;if(maxXMedia>PtrParmConv->Get_Nx()) maxXMedia=PtrParmConv->Get_Nx();
		if(minXMedia>maxXMedia)minXMedia=maxXMedia;
		if(minYMedia<0) minYMedia=0;if(maxYMedia>PtrParmConv->Get_Ny()) maxYMedia=PtrParmConv->Get_Ny();
		if(minYMedia>maxYMedia)minYMedia=maxYMedia;
		if(minZMedia<0) minZMedia=0;if(maxZMedia>PtrParmConv->Get_Nz()) maxZMedia=PtrParmConv->Get_Nz();
		if(minZMedia>maxZMedia)minZMedia=maxZMedia;
		if(PtrParmConv->IsCalculatePorosity())
		{
			porosity=Porosity();
		}
		if(PtrParmConv->IsCalculateProductionRate())
		{
			if(PtrParmConv->Get_Model()==SolverEnum::SinglePhase)
			{
				if(PtrMultiBlockConv->IsMainProcessor())
					std::cout<<"Single phase Case: no calculation of Production Rate will be done."<<std::endl;
				PtrParmConv->CalculateProductionRate(false);
			}
			else
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
						ProductionRatefile<<"Time,Production Rate"<<std::endl;
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
				PtrDicConv->Get_PtrVar("Density",RhoPerm);
				UPerm=new double*[2];
				PtrDicConv->Get_PtrVar("VelocityX",UPerm[0]);
				PtrDicConv->Get_PtrVar("VelocityY",UPerm[1]);

				if(!(PtrParmConv->Get_Model()==SolverEnum::SinglePhase))
				{
					PtrDicConv->Get_PtrVar("RhoN",RhoNPerm);
					PtrPermeabilityDarcy=&Convergence::Calcul_DarcyPermeability_TwoPhases;
				}
				else
					PtrPermeabilityDarcy=&Convergence::Calcul_DarcyPermeability_SinglePhase;

				if(PtrMultiBlockConv->IsMainProcessor())
				{
					ofstream Permeabilityfile;
					Permeabilityfile.open("Permeability.txt",ios::out | ios::trunc);
					if(PtrParmConv->Get_Model()==SolverEnum::SinglePhase)
						Permeabilityfile<<"Time,Permeability,Average Velocity,Average Pore Velocity,Average Viscosity,Delta P,Average Pressure Inlet,Average Pressure Outlet"<<std::endl;
					else
						Permeabilityfile<<"Time,Permeability,Permeability Phase 1,Permeability Phase 2,Average Velocity,Average Velocity Phase 1,Average Velocity Phase 2,Average Pore Velocity,Average Pore Velocity Phase 1,Average Pore Velocity Phase 2,Average Viscosity,Average Viscosity Phase 1,Average Viscosity Phase 2,Delta P,Average Pressure Inlet,Average Pressure Outlet,Alpha In,Alpha Out"<<std::endl;
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
		(this->*PtrPermeabilityDarcy)(Time);
		//Calcul_DarcyPermeability_TwoPhases(Time);
}
///Calcul Production Rate and write it in a file
double Convergence::Calcul_ProductionRate(int &Time){

	double ProductionRate=0;
	double sum=0,sumtmp=0,nbsumnodes=0;
	for (int i=0;i<OutletPatchId.size();i++)
	{
		if(PtrPatchBcConv->Get_NodeIndex(OutletPatchId[i],NodeId,NodeIdSpeWall,NodeIdGloCorner))
		{
			Sum_WeightNodeArray(*NodeId,*NodeIdSpeWall,*NodeIdGloCorner,RhoNProductionRate,true,sumtmp);
			//Sum_ScalarNodeArray(*NodeId,*NodeIdSpeWall,*NodeIdGloCorner,RhoNProductionRate,sumtmp);
			sum+=sumtmp;
			nbsumnodes+=NodeId->size()+0.5*(NodeIdSpeWall->size()+NodeIdGloCorner->size());
		}
	}
	sum=PtrMultiBlockConv->SumAllProcessors(&sum);
	nbsumnodes=PtrMultiBlockConv->SumAllProcessors(&nbsumnodes);
	ProductionRate=sum/nbsumnodes;//(sum/nbsumnodes+1.0)*0.5;
	if(PtrMultiBlockConv->IsMainProcessor())
	{
		ofstream ProductionRatefile;
		ProductionRatefile.open("ProductionRate.txt",ios::out | ios::app);
		ProductionRatefile<<Time<<','<<ProductionRate<<std::endl;
		ProductionRatefile.close();
		std::cout<<"Production Rate is: "<<ProductionRate<<std::endl;
	}
	return ProductionRate;
}

double Convergence::Calcul_DarcyPermeability_SinglePhase(int &Time){

	double sumtmp=0,nbsumnodesIn=0,nbsumnodesOut=0;
	double sumRhoIn=0,sumRhoOut=0,sumUIn=0,sumUOut=0;
	double sumMuIn=0,sumMuOut=0;
	for (int i=0;i<InletPatchId.size();i++)
	{
		if(PtrPatchBcConv->Get_NodeIndex(InletPatchId[i],NodeId,NodeIdSpeWall,NodeIdGloCorner))
		{
			//sum Rho1 and Rho2 with weight (alpha*Rhor or (1-alpha)*Rhob)
			Sum_ScalarNodeArray(*NodeId,*NodeIdSpeWall,*NodeIdGloCorner,RhoPerm,sumtmp);
			sumRhoIn=sumtmp;
			//sum U1 and U2 with weight (alpha*U1 or (1-alpha)*U2)
			Sum_ScalarNodeArray(*NodeId,*NodeIdSpeWall,*NodeIdGloCorner,UPerm[PtrPatchBcConv->Get_Orientation(InletPatchId[i])],sumtmp);
			sumUIn+=sumtmp;
			//Sum Mu1 and Mu2
			Sum_ViscosityNodeArray(*NodeId,*NodeIdSpeWall,*NodeIdGloCorner,RhoPerm,sumtmp);
			sumMuIn+=sumtmp;
			//Sum number of nodes in the patch
			nbsumnodesIn+=NodeId->size()+0.5*(NodeIdSpeWall->size()+NodeIdGloCorner->size());
		}
	}
	sumRhoIn=PtrMultiBlockConv->SumAllProcessors(&sumRhoIn);
	sumUIn=PtrMultiBlockConv->SumAllProcessors(&sumUIn);;
	sumMuIn=PtrMultiBlockConv->SumAllProcessors(&sumMuIn);
	nbsumnodesIn=PtrMultiBlockConv->SumAllProcessors(&nbsumnodesIn);
	for (int i=0;i<OutletPatchId.size();i++)
	{
		if(PtrPatchBcConv->Get_NodeIndex(OutletPatchId[i],NodeId,NodeIdSpeWall,NodeIdGloCorner))
		{
			//sum Rho1 and Rho2 with weight (alpha*Rhor or (1-alpha)*Rhob)
			Sum_ScalarNodeArray(*NodeId,*NodeIdSpeWall,*NodeIdGloCorner,RhoPerm,sumtmp);
			sumRhoOut=sumtmp;
			//sum U1 and U2 with weight (alpha*U1 or (1-alpha)*U2)
			Sum_ScalarNodeArray(*NodeId,*NodeIdSpeWall,*NodeIdGloCorner,UPerm[PtrPatchBcConv->Get_Orientation(OutletPatchId[i])],sumtmp);
			sumUOut+=sumtmp;
			//Sum Mu1 and Mu2
			Sum_ViscosityNodeArray(*NodeId,*NodeIdSpeWall,*NodeIdGloCorner,RhoPerm,sumtmp);
			sumMuOut+=sumtmp;
			//Sum number of nodes in the patch
			nbsumnodesOut+=NodeId->size()+0.5*(NodeIdSpeWall->size()+NodeIdGloCorner->size());
		}
	}
	sumRhoOut=PtrMultiBlockConv->SumAllProcessors(&sumRhoOut);
	sumUOut=PtrMultiBlockConv->SumAllProcessors(&sumUOut);
	sumMuOut=PtrMultiBlockConv->SumAllProcessors(&sumMuOut);
	nbsumnodesOut=PtrMultiBlockConv->SumAllProcessors(&nbsumnodesOut);

	double avgRhoIn=sumRhoIn/nbsumnodesIn;
	double avgRhoOut=sumRhoOut/nbsumnodesOut;
	double avgU=(sumUIn+sumUOut)/(2.0*SectionMedium);
	double avgUPore=(sumUIn+sumUOut)/(nbsumnodesIn+nbsumnodesOut);
	double avgMu=(sumMuIn+sumMuOut)/(nbsumnodesIn+nbsumnodesOut);
	double deltaP=RhoToP(avgRhoIn-avgRhoOut);

	double Permeability=0;//As single Phase

	if(std::abs(deltaP)>0)
	{
		Permeability=LengthMedium*avgMu*avgU/deltaP;
	}
	else
	{
		Permeability=1;
	}


	if(PtrMultiBlockConv->IsMainProcessor())
	{
		ofstream Permeabilityfile;
		Permeabilityfile.open("Permeability.txt",ios::out | ios::app);
		Permeabilityfile<<Time<<","<<Permeability<<","<<avgU<<","<<avgUPore<<","<<avgMu<<","<<deltaP<<","<<avgRhoIn<<","<<avgRhoOut<<std::endl;
		Permeabilityfile.close();
		std::cout<<"Permeability: "<<Permeability<<
				" Average Velocity: "<<avgU<<" Average Pore Velocity: "<<avgUPore<<" Delta P: "<<std::abs(deltaP)
		<<" Viscosity: "<<avgMu
		<<" Rho In: "<<avgRhoIn<<" Rho Out: "<<avgRhoOut<<std::endl;
	}
	return Permeability;

}
double Convergence::Calcul_DarcyPermeability_TwoPhases(int &Time){

	double sumtmp=0,nbsumnodesIn=0,nbsumnodesOut=0;
	//double sumUIn=0,sumUOut=0,sumRhoIn=0,sumRhoOut=0,sumRhoNIn=0,sumRhoNOut=0;
	double sumRho1In=0,sumRho2In=0,sumRho1Out=0,sumRho2Out=0,sumU1In=0,sumU1Out=0,sumU2In=0,sumU2Out=0;
	double sumMu1In=0,sumMu2In=0,sumMu1Out=0,sumMu2Out=0,sumAlpha1In=0,sumAlpha2In=0,sumAlpha1Out=0,sumAlpha2Out=0;
	for (int i=0;i<InletPatchId.size();i++)
	{
		if(PtrPatchBcConv->Get_NodeIndex(InletPatchId[i],NodeId,NodeIdSpeWall,NodeIdGloCorner))
		{
			//sum Rho1 and Rho2 with weight (alpha*Rhor or (1-alpha)*Rhob)
			Sum_ScalarNodeArray(*NodeId,*NodeIdSpeWall,*NodeIdGloCorner,RhoPerm,RhoNPerm,true,sumtmp);
			sumRho1In=sumtmp;
			Sum_ScalarNodeArray(*NodeId,*NodeIdSpeWall,*NodeIdGloCorner,RhoPerm,RhoNPerm,false,sumtmp);
			sumRho2In=sumtmp;
			//sum U1 and U2 with weight (alpha*U1 or (1-alpha)*U2)
			int orienta=PtrPatchBcConv->Get_Orientation(InletPatchId[i]);
			Sum_ScalarNodeArray(*NodeId,*NodeIdSpeWall,*NodeIdGloCorner,UPerm[PtrPatchBcConv->Get_Orientation(InletPatchId[i])],RhoNPerm,true,sumtmp);
			sumU1In+=sumtmp;
			Sum_ScalarNodeArray(*NodeId,*NodeIdSpeWall,*NodeIdGloCorner,UPerm[PtrPatchBcConv->Get_Orientation(InletPatchId[i])],RhoNPerm,false,sumtmp);
			sumU2In+=sumtmp;
			//Sum alpha and (1-alpha)
			Sum_WeightNodeArray(*NodeId,*NodeIdSpeWall,*NodeIdGloCorner,RhoNPerm,true,sumtmp);
			sumAlpha1In+=sumtmp;
			Sum_WeightNodeArray(*NodeId,*NodeIdSpeWall,*NodeIdGloCorner,RhoNPerm,false,sumtmp);
			sumAlpha2In+=sumtmp;
			//Sum Mu1 and Mu2
			Sum_ViscosityNodeArray(*NodeId,*NodeIdSpeWall,*NodeIdGloCorner,RhoPerm,RhoNPerm,true,sumtmp);
			sumMu1In+=sumtmp;
			Sum_ViscosityNodeArray(*NodeId,*NodeIdSpeWall,*NodeIdGloCorner,RhoPerm,RhoNPerm,false,sumtmp);
			sumMu2In+=sumtmp;
			//Sum number of nodes in the patch
			nbsumnodesIn+=NodeId->size()+0.5*(NodeIdSpeWall->size()+NodeIdGloCorner->size());
		}
	}
	sumRho1In=PtrMultiBlockConv->SumAllProcessors(&sumRho1In);
	sumRho2In=PtrMultiBlockConv->SumAllProcessors(&sumRho2In);
	sumU1In=PtrMultiBlockConv->SumAllProcessors(&sumU1In);
	sumU2In=PtrMultiBlockConv->SumAllProcessors(&sumU2In);
	sumAlpha1In=PtrMultiBlockConv->SumAllProcessors(&sumAlpha1In);
	sumAlpha2In=PtrMultiBlockConv->SumAllProcessors(&sumAlpha2In);
	sumMu1In=PtrMultiBlockConv->SumAllProcessors(&sumMu1In);
	sumMu2In=PtrMultiBlockConv->SumAllProcessors(&sumMu2In);
	nbsumnodesIn=PtrMultiBlockConv->SumAllProcessors(&nbsumnodesIn);
	for (int i=0;i<OutletPatchId.size();i++)
	{
		if(PtrPatchBcConv->Get_NodeIndex(OutletPatchId[i],NodeId,NodeIdSpeWall,NodeIdGloCorner))
		{
			//sum Rho1 and Rho2 with weight (alpha*Rhor or (1-alpha)*Rhob)
			Sum_ScalarNodeArray(*NodeId,*NodeIdSpeWall,*NodeIdGloCorner,RhoPerm,RhoNPerm,true,sumtmp);
			sumRho1Out=sumtmp;
			Sum_ScalarNodeArray(*NodeId,*NodeIdSpeWall,*NodeIdGloCorner,RhoPerm,RhoNPerm,false,sumtmp);
			sumRho2Out=sumtmp;
			//sum U1 and U2 with weight (alpha*U1 or (1-alpha)*U2)
			Sum_ScalarNodeArray(*NodeId,*NodeIdSpeWall,*NodeIdGloCorner,UPerm[PtrPatchBcConv->Get_Orientation(OutletPatchId[i])],RhoNPerm,true,sumtmp);
			sumU1Out+=sumtmp;
			Sum_ScalarNodeArray(*NodeId,*NodeIdSpeWall,*NodeIdGloCorner,UPerm[PtrPatchBcConv->Get_Orientation(OutletPatchId[i])],RhoNPerm,false,sumtmp);
			sumU2Out+=sumtmp;
			//Sum alpha and (1-alpha)
			Sum_WeightNodeArray(*NodeId,*NodeIdSpeWall,*NodeIdGloCorner,RhoNPerm,true,sumtmp);
			sumAlpha1Out+=sumtmp;
			Sum_WeightNodeArray(*NodeId,*NodeIdSpeWall,*NodeIdGloCorner,RhoNPerm,false,sumtmp);
			sumAlpha2Out+=sumtmp;
			//Sum Mu1 and Mu2
			Sum_ViscosityNodeArray(*NodeId,*NodeIdSpeWall,*NodeIdGloCorner,RhoPerm,RhoNPerm,true,sumtmp);
			sumMu1Out+=sumtmp;
			Sum_ViscosityNodeArray(*NodeId,*NodeIdSpeWall,*NodeIdGloCorner,RhoPerm,RhoNPerm,false,sumtmp);
			sumMu2Out+=sumtmp;
			//Sum number of nodes in the patch
			nbsumnodesOut+=NodeId->size()+0.5*(NodeIdSpeWall->size()+NodeIdGloCorner->size());
		}
	}
	sumRho1Out=PtrMultiBlockConv->SumAllProcessors(&sumRho1Out);
	sumRho2Out=PtrMultiBlockConv->SumAllProcessors(&sumRho2Out);
	sumU1Out=PtrMultiBlockConv->SumAllProcessors(&sumU1Out);
	sumU2Out=PtrMultiBlockConv->SumAllProcessors(&sumU2Out);
	sumAlpha1Out=PtrMultiBlockConv->SumAllProcessors(&sumAlpha1Out);
	sumAlpha2Out=PtrMultiBlockConv->SumAllProcessors(&sumAlpha2Out);
	sumMu1Out=PtrMultiBlockConv->SumAllProcessors(&sumMu1Out);
	sumMu2Out=PtrMultiBlockConv->SumAllProcessors(&sumMu2Out);
	nbsumnodesOut=PtrMultiBlockConv->SumAllProcessors(&nbsumnodesOut);

	double avgRho1In=sumRho1In/sumAlpha1In;
	double avgRho2In=sumRho2In/sumAlpha2In;
	double avgRho1Out=sumRho1Out/sumAlpha1Out;
	double avgRho2Out=sumRho2Out/sumAlpha2Out;
	double avgRhoIn=(sumRho1In+sumRho2In)/(sumAlpha1In+sumAlpha2In);
	double avgRhoOut=(sumRho1Out+sumRho2Out)/(sumAlpha1Out+sumAlpha2Out);

	double avgAlphaIn=sumAlpha1In/nbsumnodesIn;
	double avgAlphaOut=sumAlpha1Out/nbsumnodesOut;

	double scale=(nbsumnodesIn+nbsumnodesOut)/(2.0*SectionMedium);
	double avgU1Pore=(sumU1In+sumU1Out)/(sumAlpha1In+sumAlpha1Out);
	double avgU2Pore=(sumU2In+sumU2Out)/(sumAlpha2In+sumAlpha2Out);
	double avgUPore=(sumU1In+sumU1Out+sumU2In+sumU2Out)/(sumAlpha1In+sumAlpha1Out+sumAlpha2In+sumAlpha2Out);

	double avgU1=scale*avgU1Pore;
	double avgU2=scale*avgU2Pore;
	double avgU=scale*avgUPore;

	double avgMu1=(sumMu1In+sumMu1Out)/(sumAlpha1In+sumAlpha1Out);
	double avgMu2=(sumMu2In+sumMu2Out)/(sumAlpha2In+sumAlpha2Out);
	double avgMu=(sumMu1In+sumMu1Out+sumMu2In+sumMu2Out)/(sumAlpha1In+sumAlpha1Out+sumAlpha2In+sumAlpha2Out);


	double deltaP1=RhoToP(avgRho1In-avgRho1Out);
	double deltaP2=RhoToP(avgRho2In-avgRho2Out);
	double deltaP=RhoToP(avgRhoIn-avgRhoOut);

	double Permeability=0;//As single Phase
	double Permeability1=0,Permeability2=0;//for each Phase
	if(std::abs(deltaP)>0)
	{
		Permeability=LengthMedium*avgMu*avgU/deltaP;
		Permeability1=LengthMedium*avgMu1*avgU1/deltaP;
		Permeability2=LengthMedium*avgMu2*avgU2/deltaP;
	}
	else
	{
		Permeability=1;
		Permeability1=1;
		Permeability2=1;
	}

//	double RelativePermeability1=Permeability1/Permeability;
//	double RelativePermeability2=Permeability2/Permeability;

	if(PtrMultiBlockConv->IsMainProcessor())
	{
		ofstream Permeabilityfile;
		Permeabilityfile.open("Permeability.txt",ios::out | ios::app);
		Permeabilityfile<<Time<<","<<Permeability<<","<<Permeability1<<","<<Permeability2<<","<<avgU<<","<<avgU1<<","<<avgU2<<","<<avgUPore<<","<<avgU1Pore<<","<<avgU2Pore<<","<<avgMu<<","<<avgMu1<<","<<avgMu2<<","<<deltaP<<","<<avgRhoIn<<","<<avgRhoOut<<","<<avgAlphaIn<<","<<avgAlphaOut<<std::endl;
		Permeabilityfile.close();
		std::cout<<"Permeability: "<<Permeability<<" Permeability 1: "<<Permeability1<<" Permeability 2: "<<Permeability2<<
				" Alpha In: "<<avgAlphaIn<<" Alpha Out: "<<avgAlphaOut<<
				" Average Velocity: "<<avgU<<" Average Pore Velocity: "<<avgUPore<<" Delta P: "<<std::abs(deltaP)
		<<" Viscosity: "<<avgMu
		<<" Rho In: "<<avgRhoIn<<" Rho Out: "<<avgRhoOut<<std::endl;
	}
	return Permeability;

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
void Convergence::Sum_ScalarNodeArray(std::vector<int>& NodeArray,std::vector<int>& NodeArraySpecialWall,std::vector<int>& NodeArrayGlobalCorner,double *&Var1,double *&weight,bool phase1,double &sum){
	sum=0;
	double sum_tmp=0;
	for(int i=0;i<NodeArray.size();i++)
		sum+=Get_Weigth(weight[NodeArray[i]],phase1)*Var1[NodeArray[i]];
	for(int i=0;i<NodeArraySpecialWall.size();i++)
		sum_tmp+=Get_Weigth(weight[NodeArraySpecialWall[i]],phase1)*Var1[NodeArraySpecialWall[i]];
	for(int i=0;i<NodeArrayGlobalCorner.size();i++)
		sum_tmp+=Get_Weigth(weight[NodeArrayGlobalCorner[i]],phase1)*Var1[NodeArrayGlobalCorner[i]];
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
void Convergence::Sum_VectorNodeArray(std::vector<int>& NodeArray,std::vector<int>& NodeArraySpecialWall,std::vector<int>& NodeArrayGlobalCorner,double *&Var1,double *&Var2,double *&weight,bool phase1,double &sum1,double &sum2){
	sum1=0;sum2=0;
	for(int i=0;i<NodeArray.size();i++)
		{sum1+=Get_Weigth(weight[NodeArray[i]],phase1)*Var1[NodeArray[i]];sum2+=Get_Weigth(weight[NodeArray[i]],phase1)*Var2[NodeArray[i]];}
	double sum_tmp1=0;double sum_tmp2=0;
	for(int i=0;i<NodeArraySpecialWall.size();i++)
		{sum_tmp1+=Get_Weigth(weight[NodeArraySpecialWall[i]],phase1)*Var1[NodeArraySpecialWall[i]];sum_tmp2+=Get_Weigth(weight[NodeArraySpecialWall[i]],phase1)*Var2[NodeArraySpecialWall[i]];}
	for(int i=0;i<NodeArrayGlobalCorner.size();i++)
		{sum_tmp1+=Get_Weigth(weight[NodeArrayGlobalCorner[i]],phase1)*Var1[NodeArrayGlobalCorner[i]];sum_tmp2+=Get_Weigth(weight[NodeArrayGlobalCorner[i]],phase1)*Var2[NodeArrayGlobalCorner[i]];}
	sum_tmp1/=2.0;sum_tmp2/=2.0;
	sum1+=sum_tmp1;sum2+=sum_tmp2;
}
void Convergence::Sum_WeightNodeArray(std::vector<int>& NodeArray,std::vector<int>& NodeArraySpecialWall,std::vector<int>& NodeArrayGlobalCorner,double *&weight,bool phase1,double &sum){
	sum=0;
	double sum_tmp=0;
	for(int i=0;i<NodeArray.size();i++)
		sum+=Get_Weigth(weight[NodeArray[i]],phase1);
	for(int i=0;i<NodeArraySpecialWall.size();i++)
		sum_tmp+=Get_Weigth(weight[NodeArraySpecialWall[i]],phase1);
	for(int i=0;i<NodeArrayGlobalCorner.size();i++)
		sum_tmp+=Get_Weigth(weight[NodeArrayGlobalCorner[i]],phase1);
	sum_tmp/=2.0;
	sum+=sum_tmp;
}
///Sum the dynamic viscosity (&mu)
void Convergence::Sum_ViscosityNodeArray(std::vector<int>& NodeArray,std::vector<int>& NodeArraySpecialWall,std::vector<int>& NodeArrayGlobalCorner,double *&Rho,double *&weight,bool phase1,double &sum){
	sum=0;
	double sum_tmp=0;
	for(int i=0;i<NodeArray.size();i++)
		sum+=Get_Weigth(weight[NodeArray[i]],phase1)*PtrViscosityConv->Get_Mu(Rho[NodeArray[i]],weight[NodeArray[i]]);
	for(int i=0;i<NodeArraySpecialWall.size();i++)
		sum_tmp+=Get_Weigth(weight[NodeArraySpecialWall[i]],phase1)*PtrViscosityConv->Get_Mu(Rho[NodeArraySpecialWall[i]],weight[NodeArraySpecialWall[i]]);
	for(int i=0;i<NodeArrayGlobalCorner.size();i++)
		sum_tmp+=Get_Weigth(weight[NodeArrayGlobalCorner[i]],phase1)*PtrViscosityConv->Get_Mu(Rho[NodeArrayGlobalCorner[i]],weight[NodeArrayGlobalCorner[i]]);
	sum_tmp/=2.0;
	sum+=sum_tmp;
}
void Convergence::Sum_ViscosityNodeArray(std::vector<int>& NodeArray,std::vector<int>& NodeArraySpecialWall,std::vector<int>& NodeArrayGlobalCorner,double *&Rho,double &sum){
	sum=0;
	double sum_tmp=0;
	for(int i=0;i<NodeArray.size();i++)
		sum+=PtrViscosityConv->Get_Mu();
	for(int i=0;i<NodeArraySpecialWall.size();i++)
		sum_tmp+=PtrViscosityConv->Get_Mu();
	for(int i=0;i<NodeArrayGlobalCorner.size();i++)
		sum_tmp+=PtrViscosityConv->Get_Mu();
	sum_tmp/=2.0;
	sum+=sum_tmp;
}
double Convergence::Porosity(){
	double fluidVolumesum=0;double solidVolumesum=0;
	SumFluidVolume(fluidVolumesum);
	SumSolidVolume(solidVolumesum);
	fluidVolumesum=PtrMultiBlockConv->SumAllProcessors(&fluidVolumesum);
	solidVolumesum=PtrMultiBlockConv->SumAllProcessors(&solidVolumesum);
	if(PtrMultiBlockConv->IsMainProcessor())
		std::cout<<"Porosity of the domain is: "<<fluidVolumesum/(fluidVolumesum+solidVolumesum)
		<<" Fluid volumes: "<<fluidVolumesum<<" Solid Volumes: "<<solidVolumesum<<std::endl;
	return fluidVolumesum/(fluidVolumesum+solidVolumesum);
}
void Convergence::SumFluidVolume(double & sumfluidvolume){
	sumfluidvolume=0;
	//sumfluidvolume+=PtrNodeArraysConv->Get_SizeNodeIdInterior();
	for (int i=0;i<PtrNodeArraysConv->Get_SizeNodeIdInterior();i++)
		if(IsGlobalCornerASpecialWall(PtrNodeArraysConv->Get_NodeIdInterior(i)))
					sumfluidvolume+=0.25;
		else if(IsInsideDomain(PtrNodeArraysConv->Get_NodeIdInterior(i)))
			sumfluidvolume+=1;
		else if(IsInDomain(PtrNodeArraysConv->Get_NodeIdInterior(i)))
			sumfluidvolume+=0.5;


	//sumfluidvolume+=0.5*(PtrNodeArraysConv->Get_SizeNodeIdPressure()+PtrNodeArraysConv->Get_SizeNodeIdVelocity()+PtrNodeArraysConv->Get_SizeNodeIdPeriodic()+PtrNodeArraysConv->Get_SizeNodeIdSymmetry());
	for (int i=0;i<PtrNodeArraysConv->Get_SizeNodeIdPressure();i++)
			if(IsGlobalCornerASpecialWall(PtrNodeArraysConv->Get_NodeIdPressure(i)))
				sumfluidvolume+=0.25;
			else if(IsInDomain(PtrNodeArraysConv->Get_NodeIdPressure(i)))
				sumfluidvolume+=0.5;
	for (int i=0;i<PtrNodeArraysConv->Get_SizeNodeIdVelocity();i++)
			if(IsGlobalCornerASpecialWall(PtrNodeArraysConv->Get_NodeIdVelocity(i)))
				sumfluidvolume+=0.25;
			else if(IsInDomain(PtrNodeArraysConv->Get_NodeIdVelocity(i)))
				sumfluidvolume+=0.5;
	for (int i=0;i<PtrNodeArraysConv->Get_SizeNodeIdPeriodic();i++)
			if(IsGlobalCornerASpecialWall(PtrNodeArraysConv->Get_NodeIdPeriodic(i)))
				sumfluidvolume+=0.25;
			else if(IsInDomain(PtrNodeArraysConv->Get_NodeIdPeriodic(i)))
				sumfluidvolume+=0.5;
	for (int i=0;i<PtrNodeArraysConv->Get_SizeNodeIdSymmetry();i++)
			if(IsGlobalCornerASpecialWall(PtrNodeArraysConv->Get_NodeIdSymmetry(i)))
				sumfluidvolume+=0.25;
			else if(IsInDomain(PtrNodeArraysConv->Get_NodeIdSymmetry(i)))
				sumfluidvolume+=0.5;

// treatment of walls
	for (int i=0;i<PtrNodeArraysConv->Get_SizeNodeIdWall();i++)
		if(IsInDomain(PtrNodeArraysConv->Get_NodeIdWall(i)))
			if(IsInsideDomain(PtrNodeArraysConv->Get_NodeIdWall(i)))
				sumfluidvolume+=0.5;
	//treatment border of the domain
			else if(IsGlobalCornerASpecialWall(PtrNodeArraysConv->Get_NodeIdWall(i)))
			{
				//need the normal toward the domain
				if(IsNormalLimitDomain(PtrNodeArraysConv->Get_NodeIdWall(i)))
					sumfluidvolume+=0.25;
			}
			else
			{

				if(IsWrongSideDomain(PtrNodeArraysConv->Get_NodeIdWall(i)))
				{
					//treat as a "special wall" of the porous media domain
					if(IsNormalLimitDomain(PtrNodeArraysConv->Get_NodeIdWall(i)))
						sumfluidvolume+=0.25;
				}
				else
					//wall on the boundary of the porous media domain but toward the domain
					sumfluidvolume+=0.5;
			}

// treatment global corner of the computational domain
	for (int i=0;i<PtrNodeArraysConv->Get_SizeNodeIdGlobalCorner();i++)
			if(IsInDomain(PtrNodeArraysConv->Get_NodeIdGlobalCorner(i)))
				sumfluidvolume+=0.25;

//treatment concave corner
	for (int i=0;i<PtrNodeArraysConv->Get_SizeNodeIdCornerConcave();i++)
		//Concave corner in the domain
		if(IsInDomain(PtrNodeArraysConv->Get_NodeIdCornerConcave(i)))
			if(IsInsideDomain(PtrNodeArraysConv->Get_NodeIdCornerConcave(i)))
				sumfluidvolume+=0.25;
			//treatment border of the domain
			//check the orientation toward or not the domain
			else if(IsWrongSideDomain(PtrNodeArraysConv->Get_NodeIdCornerConcave(i)))
			{
				//keep corner in the domain
				if(IsNormalLimitDomain(PtrNodeArraysConv->Get_NodeIdCornerConcave(i)))
					sumfluidvolume+=0.25;
			}
			else
				sumfluidvolume+=0.25;


//treatment of the special wall of the computational domain
	for (int i=0;i<PtrNodeArraysConv->Get_SizeNodeIdSpecialWall();i++)
		if(IsGlobalCornerASpecialWall(PtrNodeArraysConv->Get_NodeIdSpecialWall(i)))
		{
			//keep only the special node oriented toward the domain at the border of the porous media domain
			if(IsNormalLimitDomain(PtrNodeArraysConv->Get_NodeIdSpecialWall(i)))
				sumfluidvolume+=0.25;
		}
		else if(IsInDomain(PtrNodeArraysConv->Get_NodeIdSpecialWall(i)))
				sumfluidvolume+=0.25;


//	sumfluidvolume+=0.75*PtrNodeArraysConv->Get_SizeNodeIdCornerConvex();
	for (int i=0;i<PtrNodeArraysConv->Get_SizeNodeIdCornerConvex();i++)
			if(IsInDomain(PtrNodeArraysConv->Get_NodeIdCornerConvex(i)))
				if(IsInsideDomain(PtrNodeArraysConv->Get_NodeIdCornerConvex(i)))
					sumfluidvolume+=0.75;
	//exclude convex corner in the corner of the porous media and with normal toward the outside diagonal
				else if(IsGlobalCornerASpecialWall(PtrNodeArraysConv->Get_NodeIdCornerConvex(i)))
						{if(!IsConvexCornerNormalOutside(PtrNodeArraysConv->Get_NodeIdCornerConvex(i)))
							sumfluidvolume+=0.25;}
				else if(IsWrongSideDomain(PtrNodeArraysConv->Get_NodeIdCornerConvex(i)))
					sumfluidvolume+=0.25;
				else
					sumfluidvolume+=0.5;
}
void Convergence::SumSolidVolume(double & sumsolidvolume){
	sumsolidvolume=0;
//Solid
	for (int i=0;i<PtrNodeArraysConv->Get_SizeNodeIdSolid();i++)
		if(IsInDomain(PtrNodeArraysConv->Get_NodeIdSolid(i)))
			if(IsGlobalCornerASpecialWall(PtrNodeArraysConv->Get_NodeIdSolid(i)))
				sumsolidvolume+=0.25;
			else if(IsBoundary(PtrNodeArraysConv->Get_NodeIdSolid(i)))
				sumsolidvolume+=0.5;
			else
				sumsolidvolume++;
//Convex corner
	for (int i=0;i<PtrNodeArraysConv->Get_SizeNodeIdCornerConvex();i++)
		if(IsInDomain(PtrNodeArraysConv->Get_NodeIdCornerConvex(i)))
			if(IsInsideDomain(PtrNodeArraysConv->Get_NodeIdCornerConvex(i)))
				sumsolidvolume+=0.25;
//keep convex corner in the corner of the porous media and with normal toward the outside diagonal
			else if(IsGlobalCornerASpecialWall(PtrNodeArraysConv->Get_NodeIdCornerConvex(i)))
					{if(IsConvexCornerNormalOutside(PtrNodeArraysConv->Get_NodeIdCornerConvex(i)))
						sumsolidvolume+=0.25;}
			else if(IsWrongSideDomain(PtrNodeArraysConv->Get_NodeIdCornerConvex(i)))
				sumsolidvolume+=0.25;

// wall
	for (int i=0;i<PtrNodeArraysConv->Get_SizeNodeIdWall();i++)
		if(IsInDomain(PtrNodeArraysConv->Get_NodeIdWall(i)))
			if(IsInsideDomain(PtrNodeArraysConv->Get_NodeIdWall(i)))
				sumsolidvolume+=0.5;
	//treatment border of the domain
			else if(IsGlobalCornerASpecialWall(PtrNodeArraysConv->Get_NodeIdWall(i)))
			{
				//need the normal toward the domain
				if(!IsNormalLimitDomain(PtrNodeArraysConv->Get_NodeIdWall(i)))
					sumsolidvolume+=0.25;
			}
			else
			{
				if(IsWrongSideDomain(PtrNodeArraysConv->Get_NodeIdWall(i)))
					//wall on the boundary of the porous media domain but toward the domain
					sumsolidvolume+=0.5;
				else if(IsNormalLimitDomain(PtrNodeArraysConv->Get_NodeIdWall(i)))
					//treat as a "special wall" of the porous media domain
					sumsolidvolume+=0.25;
			}

//Special wall
	for (int i=0;i<PtrNodeArraysConv->Get_SizeNodeIdSpecialWall();i++)
		if(IsInsideDomain(PtrNodeArraysConv->Get_NodeIdSpecialWall(i)))
			if(!IsGlobalCornerASpecialWall(PtrNodeArraysConv->Get_NodeIdSpecialWall(i)))
				sumsolidvolume+=0.25;

//Corner concave
	for (int i=0;i<PtrNodeArraysConv->Get_SizeNodeIdCornerConcave();i++)

	//Concave corner in the domain
	if(IsInDomain(PtrNodeArraysConv->Get_NodeIdCornerConcave(i)))
		if(IsInsideDomain(PtrNodeArraysConv->Get_NodeIdCornerConcave(i)))
			sumsolidvolume+=0.75;
		//treatment border of the domain
		//check the orientation toward or not the domain
		else if(IsWrongSideDomain(PtrNodeArraysConv->Get_NodeIdCornerConcave(i)))
		{
			if(IsGlobalCornerASpecialWall(PtrNodeArraysConv->Get_NodeIdCornerConcave(i)))
				sumsolidvolume+=0.25;
			//keep corner in the domain
			else if(IsNormalLimitDomain(PtrNodeArraysConv->Get_NodeIdCornerConcave(i)))
				sumsolidvolume+=0.25;
			else
				sumsolidvolume+=0.5;
		}

}
bool Convergence::IsBoundary(int idx){
	double x,y;
	PtrNodeArraysConv->Get_coordinate(idx,x,y);
	if(x==minXMedia || x==maxXMedia||y==minYMedia || y==maxYMedia)
		return true;
	else
		return false;
}
bool Convergence::IsInDomain(int idx){
	double x,y;
	PtrNodeArraysConv->Get_coordinate(idx,x,y);
	if(x>=minXMedia && x<=maxXMedia && y>=minYMedia && y<=maxYMedia)
		return true;
	else
		return false;
}
bool Convergence::IsInsideDomain(int idx){
	double x,y;
	PtrNodeArraysConv->Get_coordinate(idx,x,y);
	if(x>minXMedia && x<maxXMedia && y>minYMedia && y<maxYMedia)
		return true;
	else
		return false;
}
bool Convergence::IsGlobalCornerASpecialWall(int idx){
	double x,y;
	PtrNodeArraysConv->Get_coordinate(idx,x,y);
	if((x==minXMedia && y==minYMedia)|| (x==minXMedia && y==maxYMedia) ||(x==maxXMedia && y==minYMedia) || (x==maxXMedia && y==maxYMedia))
		return true;
	else
		return false;
}
//Check if the normal of the wall is in the domain
bool Convergence::IsWrongSideDomain(int idx){
	double x,y;
	PtrNodeArraysConv->Get_CoordinateNextNodeAtNormal(idx,x,y);
	if(x>minXMedia && x<maxXMedia && y>minYMedia && y<maxYMedia)
		return false;
	else
		return true;
}
//check if the wall is in the limit of the porous media box similar of a special wall for the computational domain
bool Convergence::IsNormalLimitDomain(int idx){
	double x,y;
	PtrNodeArraysConv->Get_CoordinateNextNodeAtNormal(idx,x,y);
	if((x==minXMedia || x==maxXMedia || y==minYMedia || y==maxYMedia) &&
			(x>=minXMedia && x<=maxXMedia && y>=minYMedia && y<=maxYMedia))
		return true;
	else
		return false;
}
bool Convergence::IsConvexCornerNormalOutside(int idx){
	double x,y;
	PtrNodeArraysConv->Get_CoordinateNextNodeAtNormal(idx,x,y);
	if((x<minXMedia && y<minYMedia)|| (x<minXMedia && y>maxYMedia) ||(x>maxXMedia && y<minYMedia) || (x>maxXMedia && y>maxYMedia))
		return true;
	else
		return false;
}
