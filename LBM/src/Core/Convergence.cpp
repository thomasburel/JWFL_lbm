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
	bool Var_found;
	switch(PtrParmConv->Get_ErrorVariable())
	{
	case SolverEnum::Density:
		PtrDicConv->Get_PtrVar("Density",Scalar_CurrentTime,Var_found);
		break;
	case SolverEnum::RhoN:
		if(PtrParmConv->Get_Model()==SolverEnum::SinglePhase)
		{
			PtrDicConv->Get_PtrVar("Density",Scalar_CurrentTime,Var_found);
			if(PtrMultiBlockConv->IsMainProcessor())
				std::cout<<"Single Phase case: Normal density is not used. Density will be used for the convergence."<<std::endl;
		}
		else
		PtrDicConv->Get_PtrVar("RhoN",Scalar_CurrentTime,Var_found);
		break;
	case SolverEnum::VelocityX:
		PtrDicConv->Get_PtrVar("VelocityX",Scalar_CurrentTime,Var_found);
		break;
	case SolverEnum::VelocityY:
		PtrDicConv->Get_PtrVar("VelocityY",Scalar_CurrentTime,Var_found);
		break;
	default:
		if(PtrMultiBlockConv->IsMainProcessor())
			std::cerr<<"Error variable type not found. Density will be used."<<std::endl;
		PtrDicConv->Get_PtrVar("Density",Scalar_CurrentTime,Var_found);
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
		if(PtrParmConv->IsCalculatePorosity() || PtrParmConv->IsCalculatePermeability())
		{
			Set_MarkFluidNodes();
			if(PtrMultiBlockConv->IsMainProcessor())
				std::cout<<"Number of fluid nodes:  "<<MarkFluidNode_V1.size()+0.75*MarkFluidNode_V075.size()+0.5*MarkFluidNode_V05.size()+0.25*MarkFluidNode_V025.size()<<std::endl;;
			porosity=Porosity();
			LuToPhy2=PtrParmConv->Get_deltax()*PtrParmConv->Get_deltax();
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
					PtrDicConv->Get_PtrVar("RhoN",RhoNProductionRate,Var_found);
					if(PtrMultiBlockConv->IsMainProcessor()&&Var_found)
					{
						ofstream ProductionRatefile;
						ProductionRatefile.open("ProductionRate.txt",ios::out | ios::trunc);
						ProductionRatefile<<"Time,Production Rate"<<std::endl;
						ProductionRatefile.close();
					}
					if(!Var_found)
					{
						if(PtrMultiBlockConv->IsMainProcessor())
							std::cerr<<"Normal density is not found. No calculation of Production Rate will be done."<<std::endl;
						PtrParmConv->CalculateProductionRate(false);
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
// Average of permeability in the domain
			DensityGradient.initGradients(2, 9,ModelEnum::FD);
			PtrDicConv->Get_PtrVar("Pressure",PressureConv,Var_found);
			DeltaP=new double*[2];
			PtrDicConv->AddVar(Vector,"Delata P",false, false,false,DeltaP[0],DeltaP[1]);
			for(int i=0;i<PtrNodeArraysConv->TypeOfNode.size();i++)
			{
				DeltaP[0][i]=0;DeltaP[1][i]=0;
			}
			//PtrDicConv->AddVar(Scalar,"Norm Delta P",false, false,false,NormDeltaP);
			if(!(PtrParmConv->Get_Model()==SolverEnum::SinglePhase))
			{
				PtrDicConv->Get_PtrVar("RhoN",RhoNPerm,Var_found);
				if(!Var_found)
				{
					if(PtrMultiBlockConv->IsMainProcessor())
						std::cerr<<"Normal density not found. No calculation of Permeability will be done."<<std::endl;
					PtrParmConv->CalculatePermeability(false);
				}
				PtrPermeability=&Convergence::Calcul_Permeability_TwoPhases;
			}
			else
				PtrPermeability=&Convergence::Calcul_Permeability_SinglePhase;

			if(PtrMultiBlockConv->IsMainProcessor())
			{
				ofstream Permeabilityfile;
				Permeabilityfile.open("Permeability.txt",ios::out | ios::trunc);
				if(PtrParmConv->Get_Model()==SolverEnum::SinglePhase)
					Permeabilityfile<<"Time,Average Permeability,Global Permeability,Pore-scale Permeability,Pore-scale Global Permeability,Average Velocity,Average Viscosity,Average Delta P"<<std::endl;
				else
					Permeabilityfile<<"Time,Average Permeability,Average Permeability Phase 1,Average Permeability Phase 2,Global Permeability,Global Permeability Phase 1,Global Permeability Phase 2"
							",Pore-scale Permeability,Pore-scale Permeability Phase 1,Pore-scale Permeability Phase 2,Pore-scale Global Permeability,Pore-scale Global Permeability Phase 1,Pore-scale Global Permeability Phase 2"
							",Average Velocity,Average Velocity Phase 1,Average Velocity Phase 2,Average Pore Velocity,Average Pore Velocity Phase 1,Average Pore Velocity Phase 2,Average Viscosity,Average Viscosity Phase 1,Average Viscosity Phase 2,Average Delta P,Average Delta P Phase 1,Average Delta P Phase 2,Alpha,1-Alpha"<<std::endl;
				Permeabilityfile.close();
			}
		}

		if(PtrParmConv->IsCalculateDarcyPermeability())
		{
// Calculate the Darcy Permeability (from inlet and outlet)
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
				PtrDicConv->Get_PtrVar("Density",RhoPerm,Var_found);
				UPerm=new double*[2];
				PtrDicConv->Get_PtrVar("VelocityX",UPerm[0],Var_found);
				PtrDicConv->Get_PtrVar("VelocityY",UPerm[1],Var_found);

				if(!(PtrParmConv->Get_Model()==SolverEnum::SinglePhase))
				{
					PtrDicConv->Get_PtrVar("RhoN",RhoNPerm,Var_found);
					if(!Var_found)
					{
						if(PtrMultiBlockConv->IsMainProcessor())
							std::cerr<<"Normal density not found. No calculation of Permeability will be done."<<std::endl;
						PtrParmConv->CalculatePermeability(false);
					}
					PtrPermeabilityDarcy=&Convergence::Calcul_DarcyPermeability_TwoPhases;
				}
				else
					PtrPermeabilityDarcy=&Convergence::Calcul_DarcyPermeability_SinglePhase;

				if(PtrMultiBlockConv->IsMainProcessor())
				{
					ofstream Permeabilityfile;
					Permeabilityfile.open("Permeability_Darcy.txt",ios::out | ios::trunc);
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
					std::cerr<<"No inlet and/or outlet Patch is/are set. No calculation of Darcy permeability will be done."<<std::endl;
				PtrParmConv->CalculatePermeability(false);
			}
		}
		if(!(PtrParmConv->IsCalculatePermeability() || PtrParmConv->IsCalculateDarcyPermeability() || PtrParmConv->IsCalculateProductionRate()))
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
	if(PtrParmConv->IsCalculateDarcyPermeability())
		(this->*PtrPermeabilityDarcy)(Time);
	if(PtrParmConv->IsCalculatePermeability())
		(this->*PtrPermeability)(Time);
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
		Permeabilityfile.open("Permeability_Darcy.txt",ios::out | ios::app);
		Permeabilityfile<<Time<<","<<Permeability<<","<<avgU<<","<<avgUPore<<","<<avgMu<<","<<deltaP<<","<<avgRhoIn<<","<<avgRhoOut<<std::endl;
		Permeabilityfile.close();
		std::cout<<"Permeability: "<<Permeability<<
				" Average Velocity: "<<avgU<<" Average Pore Velocity: "<<avgUPore<<" Delta P: "<<std::abs(deltaP)
		<<" Viscosity: "<<avgMu
		<<" Rho In: "<<avgRhoIn<<" Rho Out: "<<avgRhoOut<<std::endl;
	}
	return Permeability;

}
double Convergence::Calcul_Permeability_SinglePhase(int &Time){
	double Umag=0;double DeltaPmag=0;double visco;
	double sum_Umag=0;double sum_DeltaPmag=0;double sum_visco=0;
	double sum_Umag_tmp=0;double sum_DeltaPmag_tmp=0;double sum_visco_tmp=0;
	double sum_permeability=0;double sum_tmp=0;
	//Update the locate DeltaP
	Calcul_localDeltaP();


	for (int i=0;i<MarkFluidNode_V1.size();i++){
		DeltaPmag=std::sqrt(DeltaP[0][MarkFluidNode_V1[i]]*DeltaP[0][MarkFluidNode_V1[i]]+DeltaP[1][MarkFluidNode_V1[i]]*DeltaP[1][MarkFluidNode_V1[i]]);
		sum_DeltaPmag+=DeltaPmag;
		Umag=std::sqrt(UPerm[0][MarkFluidNode_V1[i]]*UPerm[0][MarkFluidNode_V1[i]]+UPerm[1][MarkFluidNode_V1[i]]*UPerm[1][MarkFluidNode_V1[i]]);
		sum_Umag+=Umag;
		visco=PtrViscosityConv->Get_Mu();
		sum_visco+=visco;
		if(DeltaPmag>0)
			sum_permeability+=Permeability(DeltaPmag,Umag,visco);
	}
	sum_tmp=0;sum_Umag_tmp=0;sum_DeltaPmag_tmp=0;sum_visco_tmp=0;

	for (int i=0;i<MarkFluidNode_V075.size();i++){
		DeltaPmag=std::sqrt(DeltaP[0][MarkFluidNode_V075[i]]*DeltaP[0][MarkFluidNode_V075[i]]+DeltaP[1][MarkFluidNode_V075[i]]*DeltaP[1][MarkFluidNode_V075[i]]);
		sum_DeltaPmag_tmp+=DeltaPmag;
		Umag=std::sqrt(UPerm[0][MarkFluidNode_V075[i]]*UPerm[0][MarkFluidNode_V075[i]]+UPerm[1][MarkFluidNode_V075[i]]*UPerm[1][MarkFluidNode_V075[i]]);
		sum_Umag_tmp+=Umag;
		visco=PtrViscosityConv->Get_Mu();
		sum_visco_tmp+=visco;
		if(DeltaPmag>0)
			sum_tmp+=Permeability(DeltaPmag,Umag,visco);
	}
	sum_permeability+=0.75*sum_tmp;sum_Umag+=0.75*sum_Umag_tmp;sum_visco+=0.75*sum_visco_tmp;sum_DeltaPmag+=0.75*sum_DeltaPmag_tmp;
	sum_tmp=0;sum_Umag_tmp=0;sum_DeltaPmag_tmp=0;sum_visco_tmp=0;

	for (int i=0;i<MarkFluidNode_V05.size();i++){
		DeltaPmag=std::sqrt(DeltaP[0][MarkFluidNode_V05[i]]*DeltaP[0][MarkFluidNode_V05[i]]+DeltaP[1][MarkFluidNode_V05[i]]*DeltaP[1][MarkFluidNode_V05[i]]);
		sum_DeltaPmag_tmp+=DeltaPmag;
		Umag=std::sqrt(UPerm[0][MarkFluidNode_V05[i]]*UPerm[0][MarkFluidNode_V05[i]]+UPerm[1][MarkFluidNode_V05[i]]*UPerm[1][MarkFluidNode_V05[i]]);
		sum_Umag_tmp+=Umag;
		visco=PtrViscosityConv->Get_Mu();
		sum_visco_tmp+=visco;
		if(DeltaPmag>0)
			sum_tmp+=Permeability(DeltaPmag,Umag,visco);
	}
	sum_permeability+=0.5*sum_tmp;sum_Umag+=0.5*sum_Umag_tmp;sum_visco+=0.5*sum_visco_tmp;sum_DeltaPmag+=0.5*sum_DeltaPmag_tmp;
	sum_tmp=0;sum_Umag_tmp=0;sum_DeltaPmag_tmp=0;sum_visco_tmp=0;

	for (int i=0;i<MarkFluidNode_V025.size();i++){
		DeltaPmag=std::sqrt(DeltaP[0][MarkFluidNode_V025[i]]*DeltaP[0][MarkFluidNode_V025[i]]+DeltaP[1][MarkFluidNode_V025[i]]*DeltaP[1][MarkFluidNode_V025[i]]);
		sum_DeltaPmag_tmp+=DeltaPmag;
		Umag=std::sqrt(UPerm[0][MarkFluidNode_V025[i]]*UPerm[0][MarkFluidNode_V025[i]]+UPerm[1][MarkFluidNode_V025[i]]*UPerm[1][MarkFluidNode_V025[i]]);
		sum_Umag_tmp+=Umag;
		visco=PtrViscosityConv->Get_Mu();
		sum_visco_tmp+=visco;
		if(DeltaPmag>0)
			sum_tmp+=Permeability(DeltaPmag,Umag,visco);
	}
	sum_permeability+=0.25*sum_tmp;sum_Umag+=0.25*sum_Umag_tmp;sum_visco+=0.25*sum_visco_tmp;sum_DeltaPmag+=0.25*sum_DeltaPmag_tmp;
	sum_tmp=0;sum_Umag_tmp=0;sum_DeltaPmag_tmp=0;sum_visco_tmp=0;

	sum_permeability=PtrMultiBlockConv->SumAllProcessors(&sum_permeability);
	sum_Umag=PtrMultiBlockConv->SumAllProcessors(&sum_Umag);
	sum_visco=PtrMultiBlockConv->SumAllProcessors(&sum_visco);
	sum_DeltaPmag=PtrMultiBlockConv->SumAllProcessors(&sum_DeltaPmag);

	double Avg_Umag=sum_Umag/fluidVolumesum;
	double Avg_visco=sum_visco/fluidVolumesum;
	double Avg_DeltaPmag=sum_DeltaPmag/fluidVolumesum;

	double PorePermeability=sum_permeability/fluidVolumesum;
	double PoreGlobalPermeability=Avg_Umag*Avg_visco/Avg_DeltaPmag;
	double Permeability=LuToPhy2*porosity*PorePermeability;
	double GlobalPermeability=LuToPhy2*porosity*PoreGlobalPermeability;

	if(PtrMultiBlockConv->IsMainProcessor())
	{
		ofstream Permeabilityfile;
		Permeabilityfile.open("Permeability.txt",ios::out | ios::app);
		Permeabilityfile<<Time<<","<<Permeability<<","<<GlobalPermeability<<","<<PorePermeability<<","<<PoreGlobalPermeability<<","<<Avg_Umag<<","<<Avg_visco<<","<<Avg_DeltaPmag<<std::endl;
		Permeabilityfile.close();
		std::cout<<"Average Permeability [m2]: "<<Permeability<<" Global Permeability [m2]: "<<GlobalPermeability<<" Pore-scale Average Permeability: "<<PorePermeability<<" Pore-scale Global Permeability: "<<PoreGlobalPermeability<<" Average Velocity: "<<Avg_Umag<<" Average DeltaP: "<<Avg_DeltaPmag<<" Average viscosity: "<<Avg_visco<<std::endl;
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
	double sumAlpha1=sumAlpha1In+sumAlpha1Out;
	double sumAlpha2=sumAlpha2In+sumAlpha2Out;
	double avgU1Pore=(sumU1In+sumU1Out)/(sumAlpha1);
	double avgU2Pore=(sumU2In+sumU2Out)/(sumAlpha2);
	double avgUPore=(sumU1In+sumU1Out+sumU2In+sumU2Out)/(sumAlpha1+sumAlpha2);

	double avgMu1=(sumMu1In+sumMu1Out)/(sumAlpha1);
	double avgMu2=(sumMu2In+sumMu2Out)/(sumAlpha2);
	double avgMu=(sumMu1In+sumMu1Out+sumMu2In+sumMu2Out)/(sumAlpha1+sumAlpha2);

	if(sumAlpha1In==0)
		avgRho1In=0;
	if(sumAlpha1Out==0)
		avgRho1Out=0;
	if(sumAlpha2In==0)
		avgRho2In=0;
	if(sumAlpha2Out==0)
		avgRho2Out=0;
	if(sumAlpha1==0)
	{
		avgU1Pore=0;
		avgMu1=0;
	}
	if(sumAlpha2==0)
	{
		avgU2Pore=0;
		avgMu2=0;
	}

	double avgU1=scale*avgU1Pore;
	double avgU2=scale*avgU2Pore;
	double avgU=scale*avgUPore;

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
		Permeability=0;
		Permeability1=0;
		Permeability2=0;
	}

//	double RelativePermeability1=Permeability1/Permeability;
//	double RelativePermeability2=Permeability2/Permeability;

	if(PtrMultiBlockConv->IsMainProcessor())
	{
		ofstream Permeabilityfile;
		Permeabilityfile.open("Permeability_Darcy.txt",ios::out | ios::app);
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
double Convergence::Calcul_Permeability_TwoPhases(int &Time){
	double Umag=0;double DeltaPmag=0;double visco;
	//double sum_Umag=0;double sum_DeltaPmag=0;double sum_visco=0;
//	double sum_Umag_tmp=0;double sum_DeltaPmag_tmp=0;double sum_visco_tmp=0;
	double sum_permeability=0,sum_permeability1=0,sum_permeability2=0;
	double sum_permeability_tmp=0,sum_permeability1_tmp=0,sum_permeability2_tmp=0;
	double alpha=0.0;double alphaM1=0.0;

	double sum_DeltaPmag=0,sum_DeltaP1mag=0,sum_DeltaP2mag=0;
	double sum_DeltaPmag_tmp=0,sum_DeltaP1mag_tmp=0,sum_DeltaP2mag_tmp=0;
	double sum_Umag=0,sum_U1mag=0,sum_U2mag=0;
	double sum_Umag_tmp=0,sum_U1mag_tmp=0,sum_U2mag_tmp=0;
	double sum_Mu=0,sum_Mu1=0,sum_Mu2=0;
	double sum_Mu_tmp=0,sum_Mu1_tmp=0,sum_Mu2_tmp=0;
	double sum_Alpha1=0,sum_Alpha2=0;
	double sum_Alpha1_tmp=0,sum_Alpha2_tmp=0;

	//Update the locate DeltaP
	Calcul_localDeltaP();


	for (int i=0;i<MarkFluidNode_V1.size();i++){
		alpha=Convert_RhoNToAlpha(RhoNPerm[MarkFluidNode_V1[i]],true);
		alphaM1=1.0-alpha;
		sum_Alpha1_tmp+=alpha;
		sum_Alpha2_tmp+=alphaM1;
		DeltaPmag=std::sqrt(DeltaP[0][MarkFluidNode_V1[i]]*DeltaP[0][MarkFluidNode_V1[i]]+DeltaP[1][MarkFluidNode_V1[i]]*DeltaP[1][MarkFluidNode_V1[i]]);
		sum_DeltaPmag_tmp+=DeltaPmag;
		sum_DeltaP1mag_tmp+=alpha*DeltaPmag;
		sum_DeltaP2mag_tmp+=alphaM1*DeltaPmag;
		Umag=std::sqrt(UPerm[0][MarkFluidNode_V1[i]]*UPerm[0][MarkFluidNode_V1[i]]+UPerm[1][MarkFluidNode_V1[i]]*UPerm[1][MarkFluidNode_V1[i]]);
		sum_Umag_tmp+=Umag;
		sum_U1mag_tmp+=alpha*Umag;
		sum_U2mag_tmp+=alphaM1*Umag;
		visco=PtrViscosityConv->Get_Mu(1,RhoNPerm[MarkFluidNode_V1[i]]);
		sum_Mu_tmp+=visco;
		sum_Mu1_tmp+=alpha*visco;
		sum_Mu2_tmp+=alphaM1*visco;
		if(DeltaPmag>0)
		{
			sum_permeability_tmp+=Permeability(DeltaPmag,Umag,visco);
			sum_permeability1_tmp+=Permeability(DeltaPmag,alpha*Umag,alpha*visco);
			sum_permeability2_tmp+=Permeability(DeltaPmag,alphaM1*Umag,alphaM1*visco);
		}
	}
	sum_permeability+=sum_permeability_tmp;sum_permeability1+=sum_permeability1_tmp;sum_permeability2+=sum_permeability2_tmp;
	sum_Umag+=sum_Umag_tmp;sum_U1mag+=sum_U1mag_tmp;sum_U2mag+=sum_U2mag_tmp;
	sum_Mu+=sum_Mu_tmp;sum_Mu1+=sum_Mu1_tmp;sum_Mu2+=sum_Mu2_tmp;
	sum_DeltaPmag+=sum_DeltaPmag_tmp;sum_DeltaP1mag+=sum_DeltaP1mag_tmp;sum_DeltaP2mag+=sum_DeltaP2mag_tmp;
	sum_Alpha1+=sum_Alpha1_tmp;sum_Alpha2+=sum_Alpha2_tmp;

	sum_permeability_tmp=0;sum_permeability1_tmp=0;sum_permeability2_tmp=0;
	sum_Umag_tmp=0;sum_U1mag_tmp=0;sum_U2mag_tmp=0;
	sum_Mu_tmp=0;sum_Mu1_tmp=0;sum_Mu2_tmp=0;
	sum_DeltaPmag_tmp=0;sum_DeltaP1mag_tmp=0;sum_DeltaP2mag_tmp=0;
	sum_Alpha1_tmp=0;sum_Alpha2_tmp=0;


	for (int i=0;i<MarkFluidNode_V075.size();i++){
		alpha=Convert_RhoNToAlpha(RhoNPerm[MarkFluidNode_V075[i]],true);
		alphaM1=1.0-alpha;
		sum_Alpha1_tmp+=alpha;
		sum_Alpha2_tmp+=alphaM1;
		DeltaPmag=std::sqrt(DeltaP[0][MarkFluidNode_V075[i]]*DeltaP[0][MarkFluidNode_V075[i]]+DeltaP[1][MarkFluidNode_V075[i]]*DeltaP[1][MarkFluidNode_V075[i]]);
		sum_DeltaPmag_tmp+=DeltaPmag;
		sum_DeltaP1mag_tmp+=alpha*DeltaPmag;
		sum_DeltaP2mag_tmp+=alphaM1*DeltaPmag;
		Umag=std::sqrt(UPerm[0][MarkFluidNode_V075[i]]*UPerm[0][MarkFluidNode_V075[i]]+UPerm[1][MarkFluidNode_V075[i]]*UPerm[1][MarkFluidNode_V075[i]]);
		sum_Umag_tmp+=Umag;
		visco=PtrViscosityConv->Get_Mu(1,RhoNPerm[MarkFluidNode_V075[i]]);
		sum_Mu_tmp+=visco;
		sum_Mu1_tmp+=alpha*visco;
		sum_Mu2_tmp+=alphaM1*visco;
		if(DeltaPmag>0)
		{
			sum_permeability_tmp+=Permeability(DeltaPmag,Umag,visco);
			sum_permeability1_tmp+=Permeability(DeltaPmag,alpha*Umag,alpha*visco);
			sum_permeability2_tmp+=Permeability(DeltaPmag,alphaM1*Umag,alphaM1*visco);
		}
	}
	sum_permeability+=0.75*sum_permeability_tmp;sum_permeability1+=0.75*sum_permeability1_tmp;sum_permeability2+=0.75*sum_permeability2_tmp;
	sum_Umag+=0.75*sum_Umag_tmp;sum_U1mag+=0.75*sum_U1mag_tmp;sum_U2mag+=0.75*sum_U2mag_tmp;
	sum_Mu+=0.75*sum_Mu_tmp;sum_Mu1+=0.75*sum_Mu1_tmp;sum_Mu2+=0.75*sum_Mu2_tmp;
	sum_DeltaPmag+=0.75*sum_DeltaPmag_tmp;sum_DeltaP1mag+=0.75*sum_DeltaP1mag_tmp;sum_DeltaP2mag+=0.75*sum_DeltaP2mag_tmp;
	sum_Alpha1+=0.75*sum_Alpha1_tmp;sum_Alpha2+=0.75*sum_Alpha2_tmp;

	sum_permeability_tmp=0;sum_permeability1_tmp=0;sum_permeability2_tmp=0;
	sum_Umag_tmp=0;sum_U1mag_tmp=0;sum_U2mag_tmp=0;
	sum_Mu_tmp=0;sum_Mu1_tmp=0;sum_Mu2_tmp=0;
	sum_DeltaPmag_tmp=0;sum_DeltaP1mag_tmp=0;sum_DeltaP2mag_tmp=0;
	sum_Alpha1_tmp=0;sum_Alpha2_tmp=0;

	for (int i=0;i<MarkFluidNode_V05.size();i++){
		alpha=Convert_RhoNToAlpha(RhoNPerm[MarkFluidNode_V05[i]],true);
		alphaM1=1.0-alpha;
		sum_Alpha1_tmp+=alpha;
		sum_Alpha2_tmp+=alphaM1;
		DeltaPmag=std::sqrt(DeltaP[0][MarkFluidNode_V05[i]]*DeltaP[0][MarkFluidNode_V05[i]]+DeltaP[1][MarkFluidNode_V05[i]]*DeltaP[1][MarkFluidNode_V05[i]]);
		sum_DeltaPmag_tmp+=DeltaPmag;
		sum_DeltaP1mag_tmp+=alpha*DeltaPmag;
		sum_DeltaP2mag_tmp+=alphaM1*DeltaPmag;
		Umag=std::sqrt(UPerm[0][MarkFluidNode_V05[i]]*UPerm[0][MarkFluidNode_V05[i]]+UPerm[1][MarkFluidNode_V05[i]]*UPerm[1][MarkFluidNode_V05[i]]);
		sum_Umag_tmp+=Umag;
		visco=PtrViscosityConv->Get_Mu(1,RhoNPerm[MarkFluidNode_V05[i]]);
		sum_Mu_tmp+=visco;
		sum_Mu1_tmp+=alpha*visco;
		sum_Mu2_tmp+=alphaM1*visco;
		if(DeltaPmag>0)
		{
			sum_permeability_tmp+=Permeability(DeltaPmag,Umag,visco);
			sum_permeability1_tmp+=Permeability(DeltaPmag,alpha*Umag,alpha*visco);
			sum_permeability2_tmp+=Permeability(DeltaPmag,alphaM1*Umag,alphaM1*visco);
		}
	}
	sum_permeability+=0.5*sum_permeability_tmp;sum_permeability1+=0.5*sum_permeability1_tmp;sum_permeability2+=0.5*sum_permeability2_tmp;
	sum_Umag+=0.5*sum_Umag_tmp;sum_U1mag+=0.5*sum_U1mag_tmp;sum_U2mag+=0.5*sum_U2mag_tmp;
	sum_Mu+=0.5*sum_Mu_tmp;sum_Mu1+=0.5*sum_Mu1_tmp;sum_Mu2+=0.5*sum_Mu2_tmp;
	sum_DeltaPmag+=0.5*sum_DeltaPmag_tmp;sum_DeltaP1mag+=0.5*sum_DeltaP1mag_tmp;sum_DeltaP2mag+=0.5*sum_DeltaP2mag_tmp;
	sum_Alpha1+=0.5*sum_Alpha1_tmp;sum_Alpha2+=0.5*sum_Alpha2_tmp;

	sum_permeability_tmp=0;sum_permeability1_tmp=0;sum_permeability2_tmp=0;
	sum_Umag_tmp=0;sum_U1mag_tmp=0;sum_U2mag_tmp=0;
	sum_Mu_tmp=0;sum_Mu1_tmp=0;sum_Mu2_tmp=0;
	sum_DeltaPmag_tmp=0;sum_DeltaP1mag_tmp=0;sum_DeltaP2mag_tmp=0;
	sum_Alpha1_tmp=0;sum_Alpha2_tmp=0;

	for (int i=0;i<MarkFluidNode_V025.size();i++){
		alpha=Convert_RhoNToAlpha(RhoNPerm[MarkFluidNode_V025[i]],true);
		alphaM1=1.0-alpha;
		sum_Alpha1_tmp+=alpha;
		sum_Alpha2_tmp+=alphaM1;
		DeltaPmag=std::sqrt(DeltaP[0][MarkFluidNode_V025[i]]*DeltaP[0][MarkFluidNode_V025[i]]+DeltaP[1][MarkFluidNode_V025[i]]*DeltaP[1][MarkFluidNode_V025[i]]);
		sum_DeltaPmag_tmp+=DeltaPmag;
		sum_DeltaP1mag_tmp+=alpha*DeltaPmag;
		sum_DeltaP2mag_tmp+=alphaM1*DeltaPmag;
		Umag=std::sqrt(UPerm[0][MarkFluidNode_V025[i]]*UPerm[0][MarkFluidNode_V025[i]]+UPerm[1][MarkFluidNode_V025[i]]*UPerm[1][MarkFluidNode_V025[i]]);
		sum_Umag_tmp+=Umag;
		visco=PtrViscosityConv->Get_Mu(1,RhoNPerm[MarkFluidNode_V025[i]]);
		sum_Mu_tmp+=visco;
		sum_Mu1_tmp+=alpha*visco;
		sum_Mu2_tmp+=alphaM1*visco;
		if(DeltaPmag>0)
		{
			sum_permeability_tmp+=Permeability(DeltaPmag,Umag,visco);
			sum_permeability1_tmp+=Permeability(DeltaPmag,alpha*Umag,alpha*visco);
			sum_permeability2_tmp+=Permeability(DeltaPmag,alphaM1*Umag,alphaM1*visco);
		}
	}
	sum_permeability+=0.25*sum_permeability_tmp;sum_permeability1+=0.25*sum_permeability1_tmp;sum_permeability2+=0.25*sum_permeability2_tmp;
	sum_Umag+=0.25*sum_Umag_tmp;sum_U1mag+=0.25*sum_U1mag_tmp;sum_U2mag+=0.25*sum_U2mag_tmp;
	sum_Mu+=0.25*sum_Mu_tmp;sum_Mu1+=0.25*sum_Mu1_tmp;sum_Mu2+=0.25*sum_Mu2_tmp;
	sum_DeltaPmag+=0.25*sum_DeltaPmag_tmp;sum_DeltaP1mag+=0.25*sum_DeltaP1mag_tmp;sum_DeltaP2mag+=0.25*sum_DeltaP2mag_tmp;
	sum_Alpha1+=0.25*sum_Alpha1_tmp;sum_Alpha2+=0.25*sum_Alpha2_tmp;

	sum_permeability_tmp=0;sum_permeability1_tmp=0;sum_permeability2_tmp=0;
	sum_Umag_tmp=0;sum_U1mag_tmp=0;sum_U2mag_tmp=0;
	sum_Mu_tmp=0;sum_Mu1_tmp=0;sum_Mu2_tmp=0;
	sum_DeltaPmag_tmp=0;sum_DeltaP1mag_tmp=0;sum_DeltaP2mag_tmp=0;
	sum_Alpha1_tmp=0;sum_Alpha2_tmp=0;

	sum_permeability=PtrMultiBlockConv->SumAllProcessors(&sum_permeability);sum_permeability1=PtrMultiBlockConv->SumAllProcessors(&sum_permeability1);sum_permeability2=PtrMultiBlockConv->SumAllProcessors(&sum_permeability2);
	sum_Umag=PtrMultiBlockConv->SumAllProcessors(&sum_Umag);sum_U1mag=PtrMultiBlockConv->SumAllProcessors(&sum_U1mag);sum_U2mag=PtrMultiBlockConv->SumAllProcessors(&sum_U2mag);
	sum_Mu=PtrMultiBlockConv->SumAllProcessors(&sum_Mu);sum_Mu1=PtrMultiBlockConv->SumAllProcessors(&sum_Mu1);sum_Mu2=PtrMultiBlockConv->SumAllProcessors(&sum_Mu2);
	sum_DeltaPmag=PtrMultiBlockConv->SumAllProcessors(&sum_DeltaPmag);sum_DeltaP1mag=PtrMultiBlockConv->SumAllProcessors(&sum_DeltaP1mag);sum_DeltaP2mag=PtrMultiBlockConv->SumAllProcessors(&sum_DeltaP2mag);
	sum_Alpha1=PtrMultiBlockConv->SumAllProcessors(&sum_Alpha1);sum_Alpha2=PtrMultiBlockConv->SumAllProcessors(&sum_Alpha2);


	double Avg_Umag=sum_Umag/fluidVolumesum;double Avg_U1mag=sum_U1mag/fluidVolumesum;double Avg_U2mag=sum_U2mag/fluidVolumesum;
	double Avg_Mu=sum_Mu/fluidVolumesum;double Avg_Mu1=sum_Mu1/fluidVolumesum;double Avg_Mu2=sum_Mu2/fluidVolumesum;
	double Avg_DeltaPmag=sum_DeltaPmag/fluidVolumesum;double Avg_DeltaP1mag=sum_DeltaP1mag/fluidVolumesum;double Avg_DeltaP2mag=sum_DeltaP2mag/fluidVolumesum;
	double Avg_Alpha1=sum_Alpha1/fluidVolumesum;double Avg_Alpha2=sum_Alpha2/fluidVolumesum;

	//double PorePermeability=sum_permeability/fluidVolumesum;
	double PorePermeability=sum_permeability/fluidVolumesum;double PorePermeability1=sum_permeability1/fluidVolumesum;double PorePermeability2=sum_permeability2/fluidVolumesum;
	double PoreGlobalPermeability=Avg_Umag*Avg_Mu/Avg_DeltaPmag;double PoreGlobalPermeability1=Avg_U1mag*Avg_Mu1/Avg_DeltaP1mag;double PoreGlobalPermeability2=Avg_U2mag*Avg_Mu2/Avg_DeltaP2mag;
	double Permeability=LuToPhy2*porosity*PorePermeability;double Permeability1=LuToPhy2*porosity*PorePermeability1;double Permeability2=LuToPhy2*porosity*PorePermeability2;
	double GlobalPermeability=LuToPhy2*porosity*PoreGlobalPermeability;double GlobalPermeability1=LuToPhy2*porosity*PoreGlobalPermeability1;double GlobalPermeability2=LuToPhy2*porosity*PoreGlobalPermeability2;
	double ScaleAvg_Umag=LuToPhy2*porosity*Avg_Umag;double ScaleAvg_U1mag=LuToPhy2*porosity*Avg_U1mag;double ScaleAvg_U2mag=LuToPhy2*porosity*Avg_U2mag;
	if(PtrMultiBlockConv->IsMainProcessor())
	{
		ofstream Permeabilityfile;
		Permeabilityfile.open("Permeability.txt",ios::out | ios::app);
		Permeabilityfile<<Time<<","<<Permeability<<","<<Permeability1<<","<<Permeability2<<","<<GlobalPermeability<<","<<GlobalPermeability1<<","<<GlobalPermeability2
				<<","<<PorePermeability<<","<<PorePermeability1<<","<<PorePermeability2<<","<<PoreGlobalPermeability<<","<<PoreGlobalPermeability1<<","<<PoreGlobalPermeability2
				<<","<<ScaleAvg_Umag<<","<<ScaleAvg_U1mag<<","<<ScaleAvg_U2mag<<","<<Avg_Umag<<","<<Avg_U1mag<<","<<Avg_U2mag<<","
				<<Avg_Mu<<","<<Avg_Mu1<<","<<Avg_Mu2
				<<","<<Avg_DeltaPmag<<","<<Avg_DeltaP1mag<<","<<Avg_DeltaP2mag<<","<<Avg_Alpha1<<","<<Avg_Alpha2<<std::endl;
		Permeabilityfile.close();
		std::cout<<"Average Permeability [m2]: "<<Permeability<<" Global Permeability [m2]: "<<GlobalPermeability<<" Pore-scale Average Permeability: "<<PorePermeability<<" Pore-scale Global Permeability: "<<PoreGlobalPermeability<<" Average Velocity: "<<Avg_Umag<<" Average DeltaP: "<<Avg_DeltaPmag<<" Average viscosity: "<<Avg_Mu<<" Average Alpha: "<<Avg_Alpha1<<std::endl;
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
	fluidVolumesum=0;solidVolumesum=0;
	fluidVolumesum=MarkFluidNode_V1.size()+0.75*MarkFluidNode_V075.size()+0.5*MarkFluidNode_V05.size()+0.25*MarkFluidNode_V025.size();
	//SumFluidVolume(fluidVolumesum);
	SumSolidVolume(solidVolumesum);
	fluidVolumesum=PtrMultiBlockConv->SumAllProcessors(&fluidVolumesum);
	solidVolumesum=PtrMultiBlockConv->SumAllProcessors(&solidVolumesum);
	if(PtrMultiBlockConv->IsMainProcessor())
		std::cout<<"Porosity of the domain is: "<<fluidVolumesum/(fluidVolumesum+solidVolumesum)
		<<" Fluid volumes: "<<fluidVolumesum<<" Solid Volumes: "<<solidVolumesum<<std::endl;
	return fluidVolumesum/(fluidVolumesum+solidVolumesum);
}
void Convergence::Set_MarkFluidNodes(){
	MarkFluidNode_V1.clear();MarkFluidNode_V075.clear();MarkFluidNode_V05.clear();MarkFluidNode_V025.clear();
//treatment Interior of the domain
	for (int i=0;i<PtrNodeArraysConv->Get_SizeNodeIdInterior();i++)
		if(IsGlobalCornerASpecialWall(PtrNodeArraysConv->Get_NodeIdInterior(i)))
			MarkFluidNode_V025.push_back(PtrNodeArraysConv->Get_NodeIdInterior(i));
		else if(IsInsideDomain(PtrNodeArraysConv->Get_NodeIdInterior(i)))
			MarkFluidNode_V1.push_back(PtrNodeArraysConv->Get_NodeIdInterior(i));
		else if(IsInDomain(PtrNodeArraysConv->Get_NodeIdInterior(i)))
			MarkFluidNode_V05.push_back(PtrNodeArraysConv->Get_NodeIdInterior(i));


//Treatment boundary conditions of the domain
	for (int i=0;i<PtrNodeArraysConv->Get_SizeNodeIdPressure();i++)
			if(IsGlobalCornerASpecialWall(PtrNodeArraysConv->Get_NodeIdPressure(i)))
				MarkFluidNode_V025.push_back(PtrNodeArraysConv->Get_NodeIdPressure(i));
			else if(IsInDomain(PtrNodeArraysConv->Get_NodeIdPressure(i)))
				MarkFluidNode_V05.push_back(PtrNodeArraysConv->Get_NodeIdPressure(i));
	for (int i=0;i<PtrNodeArraysConv->Get_SizeNodeIdVelocity();i++)
			if(IsGlobalCornerASpecialWall(PtrNodeArraysConv->Get_NodeIdVelocity(i)))
				MarkFluidNode_V025.push_back(PtrNodeArraysConv->Get_NodeIdVelocity(i));
			else if(IsInDomain(PtrNodeArraysConv->Get_NodeIdVelocity(i)))
				MarkFluidNode_V05.push_back(PtrNodeArraysConv->Get_NodeIdVelocity(i));
	for (int i=0;i<PtrNodeArraysConv->Get_SizeNodeIdPeriodic();i++)
			if(IsGlobalCornerASpecialWall(PtrNodeArraysConv->Get_NodeIdPeriodic(i)))
				MarkFluidNode_V025.push_back(PtrNodeArraysConv->Get_NodeIdPeriodic(i));
			else if(IsInDomain(PtrNodeArraysConv->Get_NodeIdPeriodic(i)))
				MarkFluidNode_V05.push_back(PtrNodeArraysConv->Get_NodeIdPeriodic(i));
	for (int i=0;i<PtrNodeArraysConv->Get_SizeNodeIdSymmetry();i++)
			if(IsGlobalCornerASpecialWall(PtrNodeArraysConv->Get_NodeIdSymmetry(i)))
				MarkFluidNode_V025.push_back(PtrNodeArraysConv->Get_NodeIdSymmetry(i));
			else if(IsInDomain(PtrNodeArraysConv->Get_NodeIdSymmetry(i)))
				MarkFluidNode_V05.push_back(PtrNodeArraysConv->Get_NodeIdSymmetry(i));

// treatment of walls
	for (int i=0;i<PtrNodeArraysConv->Get_SizeNodeIdWall();i++)
		if(IsInDomain(PtrNodeArraysConv->Get_NodeIdWall(i)))
			if(IsInsideDomain(PtrNodeArraysConv->Get_NodeIdWall(i)))
				MarkFluidNode_V05.push_back(PtrNodeArraysConv->Get_NodeIdWall(i));
	//treatment border of the domain
			else if(IsGlobalCornerASpecialWall(PtrNodeArraysConv->Get_NodeIdWall(i)))
			{
				//need the normal toward the domain
				if(IsNormalLimitDomain(PtrNodeArraysConv->Get_NodeIdWall(i)))
					MarkFluidNode_V025.push_back(PtrNodeArraysConv->Get_NodeIdWall(i));
			}
			else
			{

				if(IsWrongSideDomain(PtrNodeArraysConv->Get_NodeIdWall(i)))
				{
					//treat as a "special wall" of the porous media domain
					if(IsNormalLimitDomain(PtrNodeArraysConv->Get_NodeIdWall(i)))
						MarkFluidNode_V025.push_back(PtrNodeArraysConv->Get_NodeIdWall(i));
				}
				else
					//wall on the boundary of the porous media domain but toward the domain
					MarkFluidNode_V05.push_back(PtrNodeArraysConv->Get_NodeIdWall(i));
			}

// treatment global corner of the computational domain
	for (int i=0;i<PtrNodeArraysConv->Get_SizeNodeIdGlobalCorner();i++)
			if(IsInDomain(PtrNodeArraysConv->Get_NodeIdGlobalCorner(i)))
				MarkFluidNode_V025.push_back(PtrNodeArraysConv->Get_NodeIdGlobalCorner(i));

//treatment concave corner
	for (int i=0;i<PtrNodeArraysConv->Get_SizeNodeIdCornerConcave();i++)
		//Concave corner in the domain
		if(IsInDomain(PtrNodeArraysConv->Get_NodeIdCornerConcave(i)))
			if(IsInsideDomain(PtrNodeArraysConv->Get_NodeIdCornerConcave(i)))
				MarkFluidNode_V025.push_back(PtrNodeArraysConv->Get_NodeIdCornerConcave(i));
			//treatment border of the domain
			//check the orientation toward or not the domain
			else if(IsWrongSideDomain(PtrNodeArraysConv->Get_NodeIdCornerConcave(i)))
			{
				//keep corner in the domain
				if(IsNormalLimitDomain(PtrNodeArraysConv->Get_NodeIdCornerConcave(i)))
					MarkFluidNode_V025.push_back(PtrNodeArraysConv->Get_NodeIdCornerConcave(i));
			}
			else
				MarkFluidNode_V025.push_back(PtrNodeArraysConv->Get_NodeIdCornerConcave(i));


//treatment of the special wall of the computational domain
	for (int i=0;i<PtrNodeArraysConv->Get_SizeNodeIdSpecialWall();i++)
		if(IsGlobalCornerASpecialWall(PtrNodeArraysConv->Get_NodeIdSpecialWall(i)))
		{
			//keep only the special node oriented toward the domain at the border of the porous media domain
			if(IsNormalLimitDomain(PtrNodeArraysConv->Get_NodeIdSpecialWall(i)))
				MarkFluidNode_V025.push_back(PtrNodeArraysConv->Get_NodeIdSpecialWall(i));
		}
		else if(IsInDomain(PtrNodeArraysConv->Get_NodeIdSpecialWall(i)))
			MarkFluidNode_V025.push_back(PtrNodeArraysConv->Get_NodeIdSpecialWall(i));


//	sumfluidvolume+=0.75*PtrNodeArraysConv->Get_SizeNodeIdCornerConvex();
	for (int i=0;i<PtrNodeArraysConv->Get_SizeNodeIdCornerConvex();i++)
			if(IsInDomain(PtrNodeArraysConv->Get_NodeIdCornerConvex(i)))
				if(IsInsideDomain(PtrNodeArraysConv->Get_NodeIdCornerConvex(i)))
					MarkFluidNode_V075.push_back(PtrNodeArraysConv->Get_NodeIdCornerConvex(i));
	//exclude convex corner in the corner of the porous media and with normal toward the outside diagonal
				else if(IsGlobalCornerASpecialWall(PtrNodeArraysConv->Get_NodeIdCornerConvex(i)))
						{
							if(!IsConvexCornerNormalOutside(PtrNodeArraysConv->Get_NodeIdCornerConvex(i)))
								MarkFluidNode_V025.push_back(PtrNodeArraysConv->Get_NodeIdCornerConvex(i));
						}
				else if(IsWrongSideDomain(PtrNodeArraysConv->Get_NodeIdCornerConvex(i)))
					MarkFluidNode_V025.push_back(PtrNodeArraysConv->Get_NodeIdCornerConvex(i));
				else
					MarkFluidNode_V05.push_back(PtrNodeArraysConv->Get_NodeIdCornerConvex(i));
}
/*
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
*/
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
void Convergence::Calcul_localDeltaP(){
	double tmp[2];
	for (int j=0;j<PtrNodeArraysConv->Get_SizeNodeIdInterior();j++)
	{
	// Calculate gradients
		DensityGradient.Grad(&tmp[0],&PressureConv[0],PtrNodeArraysConv->Get_NodeConnectInterior(j),PtrNodeArraysConv->Get_NodeNormalInterior(j));
		DeltaP[0][PtrNodeArraysConv->Get_NodeIdInterior(j)]=tmp[0];DeltaP[1][PtrNodeArraysConv->Get_NodeIdInterior(j)]=tmp[1];
	//	NormDeltaP[PtrNodeArraysConv->Get_NodeIdInterior(j)]=std::sqrt(tmp[0]*tmp[0]+tmp[1]*tmp[1]);
	}

	for (int j=0;j<PtrNodeArraysConv->Get_SizeNodeIdCorner();j++)
	{
		// Calculate gradients
			DensityGradient.GradBc(&tmp[0],&PressureConv[0],PtrNodeArraysConv->Get_NodeConnectCorner(j),PtrNodeArraysConv->Get_NodeNormalCorner(j));
			DeltaP[0][PtrNodeArraysConv->Get_NodeIdCorner(j)]=tmp[0];DeltaP[1][PtrNodeArraysConv->Get_NodeIdCorner(j)]=tmp[1];
//			NormDeltaP[PtrNodeArraysConv->Get_NodeIdCorner(j)]=std::sqrt(tmp[0]*tmp[0]+tmp[1]*tmp[1]);
	}

	for (int j=0;j<PtrNodeArraysConv->Get_SizeNodeIdGlobalCorner();j++)
	{
		// Calculate gradients
			DensityGradient.GradBc(&tmp[0],&PressureConv[0],PtrNodeArraysConv->Get_NodeConnectGlobalCorner(j),PtrNodeArraysConv->Get_NodeNormalGlobalCorner(j));
			DeltaP[0][PtrNodeArraysConv->Get_NodeIdGlobalCorner(j)]=tmp[0];DeltaP[1][PtrNodeArraysConv->Get_NodeIdGlobalCorner(j)]=tmp[1];
	//		NormDeltaP[PtrNodeArraysConv->Get_NodeIdGlobalCorner(j)]=std::sqrt(tmp[0]*tmp[0]+tmp[1]*tmp[1]);
	}

	for (int j=0;j<PtrNodeArraysConv->Get_SizeNodeIdVelocity();j++)
	{
		// Calculate gradients
			DensityGradient.GradBc(&tmp[0],&PressureConv[0],PtrNodeArraysConv->Get_NodeConnectVelocity(j),PtrNodeArraysConv->Get_NodeNormalVelocity(j));
			DeltaP[0][PtrNodeArraysConv->Get_NodeIdVelocity(j)]=tmp[0];DeltaP[1][PtrNodeArraysConv->Get_NodeIdVelocity(j)]=tmp[1];
//			NormDeltaP[PtrNodeArraysConv->Get_NodeIdVelocity(j)]=std::sqrt(tmp[0]*tmp[0]+tmp[1]*tmp[1]);
	}

	for (int j=0;j<PtrNodeArraysConv->Get_SizeNodeIdPressure();j++)
	{
		// Calculate gradients
			DensityGradient.GradBc(&tmp[0],&PressureConv[0],PtrNodeArraysConv->Get_NodeConnectPressure(j),PtrNodeArraysConv->Get_NodeNormalPressure(j));
			DeltaP[0][PtrNodeArraysConv->Get_NodeIdPressure(j)]=tmp[0];DeltaP[1][PtrNodeArraysConv->Get_NodeIdPressure(j)]=tmp[1];
	//		NormDeltaP[PtrNodeArraysConv->Get_NodeIdPressure(j)]=std::sqrt(tmp[0]*tmp[0]+tmp[1]*tmp[1]);
	}

	for (int j=0;j<PtrNodeArraysConv->Get_SizeNodeIdWall();j++)
	{
		// Calculate gradients
			DensityGradient.GradBc(&tmp[0],&PressureConv[0],PtrNodeArraysConv->Get_NodeConnectWall(j),PtrNodeArraysConv->Get_NodeNormalWall(j));
			DeltaP[0][PtrNodeArraysConv->Get_NodeIdWall(j)]=tmp[0];DeltaP[1][PtrNodeArraysConv->Get_NodeIdWall(j)]=tmp[1];
//			NormDeltaP[PtrNodeArraysConv->Get_NodeIdWall(j)]=std::sqrt(tmp[0]*tmp[0]+tmp[1]*tmp[1]);
	}
	for (int j=0;j<PtrNodeArraysConv->Get_SizeNodeIdSpecialWall();j++)
	{
		// Calculate gradients
			DensityGradient.GradBc(&tmp[0],&PressureConv[0],PtrNodeArraysConv->Get_NodeConnectSpecialWall(j),PtrNodeArraysConv->Get_NodeNormalSpecialWall(j));
			DeltaP[0][PtrNodeArraysConv->Get_NodeIdSpecialWall(j)]=tmp[0];DeltaP[1][PtrNodeArraysConv->Get_NodeIdSpecialWall(j)]=tmp[1];
//			NormDeltaP[PtrNodeArraysConv->Get_NodeIdSpecialWall(j)]=std::sqrt(tmp[0]*tmp[0]+tmp[1]*tmp[1]);
	}
	for (int j=0;j<PtrNodeArraysConv->Get_SizeNodeIdSymmetry();j++)
	{
		// Calculate gradients
			DensityGradient.GradBc(&tmp[0],&PressureConv[0],PtrNodeArraysConv->Get_NodeConnectSymmetry(j),PtrNodeArraysConv->Get_NodeNormalSymmetry(j));
			DeltaP[0][PtrNodeArraysConv->Get_NodeIdSymmetry(j)]=tmp[0];DeltaP[1][PtrNodeArraysConv->Get_NodeIdSymmetry(j)]=tmp[1];
//			NormDeltaP[PtrNodeArraysConv->Get_NodeIdSymmetry(j)]=std::sqrt(tmp[0]*tmp[0]+tmp[1]*tmp[1]);
	}
	for (int j=0;j<PtrNodeArraysConv->Get_SizeNodeIdPeriodic();j++)
	{
		// Calculate gradients
			DensityGradient.GradBc(&tmp[0],&PressureConv[0],PtrNodeArraysConv->Get_NodeConnectPeriodic(j),PtrNodeArraysConv->Get_NodeNormalPeriodic(j));
			DeltaP[0][PtrNodeArraysConv->Get_NodeIdPeriodic(j)]=tmp[0];DeltaP[1][PtrNodeArraysConv->Get_NodeIdPeriodic(j)]=tmp[1];
//			NormDeltaP[PtrNodeArraysConv->Get_NodeIdPeriodic(j)]=std::sqrt(tmp[0]*tmp[0]+tmp[1]*tmp[1]);
	}
}
