/*
 * Solver.cpp
 *
 *  Created on: 17 Apr 2015
 *      Author: thomas
 */

#include "Solver.h"

Solver::Solver() {
	PtrParameters=0;
//	InvTau=1; /// Default  value of Tau=1

	nbvelo=1;
	nbnode=1;
}

Solver::~Solver() {

}

SolverSinglePhase::SolverSinglePhase() {
	PtrParameters=0;
	f=0;
	ftmp=0;

}

SolverSinglePhase::~SolverSinglePhase() {
	delete f;
	delete [] ftmp;
}

SolverTwoPhases::SolverTwoPhases() {
	PtrParameters=0;
	f=0;
	ftmp=0;

}

SolverTwoPhases::~SolverTwoPhases() {
	if(f!=0)
	{
	delete f[0];
	delete f[1];
	delete f;
	}
	delete [] ftmp;
}

SolverSinglePhaseLowOrder::SolverSinglePhaseLowOrder() {
	PtrParameters=0;
}

SolverSinglePhaseLowOrder::~SolverSinglePhaseLowOrder() {
	// TODO Auto-generated destructor stub
}
SolverTwoPhasesLowOrder::SolverTwoPhasesLowOrder() {
	PtrParameters=0;
}

SolverTwoPhasesLowOrder::~SolverTwoPhasesLowOrder() {
	// TODO Auto-generated destructor stub
}
SolverSinglePhaseLowOrder2D::SolverSinglePhaseLowOrder2D() {
	Solver::PtrParameters=0;

}

SolverSinglePhaseLowOrder2D::~SolverSinglePhaseLowOrder2D() {
	// TODO Auto-generated destructor stub
}

void SolverSinglePhaseLowOrder2D::get_time() {

}
void SolverSinglePhaseLowOrder2D::Set_Solver(MultiBlock* PtrMultiBlock_,ParallelManager* PtrParallel_,WriterManager* PtrWriter_, Parameters* PtrParameters_){
	MultiBlock_=PtrMultiBlock_;
	PtrMultiBlockConv=PtrMultiBlock_;
	parallel=PtrParallel_;
	Writer=PtrWriter_;
	Solver::PtrParameters=PtrParameters_;
	IniTau(PtrParameters);
	PtrBlockStream=MultiBlock_->Get_Block();
	PtrBlockCollide=PtrBlockStream;
	//Node=MultiBlock_->Get_Block()->Get_PtrNode();
	NodeArrays=MultiBlock_->Get_Block()->Get_NodeArrays2D();
	PatchsBc=MultiBlock_->Get_Block()->Get_PatchBc();
	IdBoundaries=MultiBlock_->Get_Block()->Get_PtrIdBc();
	//InvTau=1.0/PtrParameters->Get_Tau();
	nbvelo=Solver::PtrParameters->Get_NbVelocities();
	nbnode=MultiBlock_->Get_nnodes();
//	Gradients::initGradients(2,nbvelo,PtrParameters->Get_GradientType());
	Extrapolation::initExtrapolation(2,nbvelo,Solver::PtrParameters->Get_ExtrapolationType());
	Dic=new Dictionary(2,nbnode);
	PtrDicConv=Dic;
	PtrParmConv=PtrParameters_;
	PtrPatchBcConv=PatchsBc;
	PtrNodeArraysConv=NodeArrays;
	PtrViscosityConv=Get_viscosity();
	Add_OneDistributionToDictionary();
	Set_Solution(PtrParameters);

}
void SolverSinglePhaseLowOrder2D::Add_OneDistributionToDictionary(){
	//First distribution
	for(int i=0;i<nbvelo;i++)
	{
		std::stringstream sstm;
		sstm<<"f[0]_"<<i;
		Dic->AddDistributionBreakpoint(sstm.str(),f->f[i]);
	}
}
void SolverSinglePhaseLowOrder2D::Updated_OneDistributionToDictionary(){
	//First distribution
	for(int i=0;i<nbvelo;i++)
	{
		Dic->UpdatedDistributionBreakpoint(i,f->f[i]);
	}
}
SolverTwoPhasesLowOrder2D::SolverTwoPhasesLowOrder2D() {
	Solver::PtrParameters=0;


}

SolverTwoPhasesLowOrder2D::~SolverTwoPhasesLowOrder2D() {
	// TODO Auto-generated destructor stub
}

void SolverTwoPhasesLowOrder2D::get_time() {

}
void SolverTwoPhasesLowOrder2D::Set_Solver(MultiBlock* PtrMultiBlock_,ParallelManager* PtrParallel_,WriterManager* PtrWriter_, Parameters* PtrParameters_){
	MultiBlock_=PtrMultiBlock_;
	PtrMultiBlockConv=PtrMultiBlock_;
	parallel=PtrParallel_;
	Writer=PtrWriter_;
	Solver::PtrParameters=PtrParameters_;
	IniTau(PtrParameters);
	PtrBlockStream=MultiBlock_->Get_Block();
	NodeArrays=MultiBlock_->Get_Block()->Get_NodeArrays2D();
	PatchsBc=MultiBlock_->Get_Block()->Get_PatchBc();
	IdBoundaries=MultiBlock_->Get_Block()->Get_PtrIdBc();
	//InvTau=1.0/PtrParameters->Get_Tau();
	nbvelo=Solver::PtrParameters->Get_NbVelocities();
	nbnode=MultiBlock_->Get_nnodes();
	Gradients::initGradients(2,nbvelo,Solver::PtrParameters->Get_GradientType());
	Extrapolation::initExtrapolation(2,nbvelo,Solver::PtrParameters->Get_ExtrapolationType());
	Dic=new Dictionary(2,nbnode);
	PtrDicConv=Dic;
	PtrParmConv=PtrParameters_;
	PtrPatchBcConv=PatchsBc;
	PtrNodeArraysConv=NodeArrays;
	PtrViscosityConv=Get_viscosity();
	Add_TwoDistributionsToDictionary();
	Set_Solution(PtrParameters);
}
void SolverTwoPhasesLowOrder2D::Add_TwoDistributionsToDictionary(){

	//First distribution
	for(int i=0;i<nbvelo;i++)
	{
		std::stringstream sstm;
		sstm<<"f[0]_"<<i;
		Dic->AddDistributionBreakpoint(sstm.str(),f[0]->f[i]);
	}
	//Second distribution
	for(int i=0;i<nbvelo;i++)
	{
		std::stringstream sstm;
		sstm<<"f[1]_"<<i;
		Dic->AddDistributionBreakpoint(sstm.str(),f[1]->f[i]);
	}
}
void SolverTwoPhasesLowOrder2D::Updated_TwoDistributionsToDictionary(){

	//First distribution
	for(int i=0;i<nbvelo;i++)
	{
		Dic->UpdatedDistributionBreakpoint(i,f[0]->f[i]);
	}
	//Second distribution
	for(int i=0;i<nbvelo;i++)
	{
		Dic->UpdatedDistributionBreakpoint(i+nbvelo,f[1]->f[i]);
	}
}
void SolverSinglePhaseLowOrder2D::Write_Breakpoint(Parameters *Param){
	Updated_OneDistributionToDictionary();
	Writer->Set_breakpoint(Dic->Get_PtrExportVarBreakpoint(),Dic->Get_PtrExportVarBreakpointName(),Dic->Get_NbExportVarBreakpoint());
	Writer->Write_breakpoint(*Param);
}
void SolverTwoPhasesLowOrder2D::Write_Breakpoint(Parameters *Param){
	Updated_TwoDistributionsToDictionary();
	Writer->Set_breakpoint(Dic->Get_PtrExportVarBreakpoint(),Dic->Get_PtrExportVarBreakpointName(),Dic->Get_NbExportVarBreakpoint());
	Writer->Write_breakpoint(*Param);
}
