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
	delete [] f[0];
	delete [] f[1];
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
	parallel=PtrParallel_;
	Writer=PtrWriter_;
	Solver::PtrParameters=PtrParameters_;
	PtrBlockStream=MultiBlock_->Get_Block();
	PtrBlockCollide=PtrBlockStream;
	//Node=MultiBlock_->Get_Block()->Get_PtrNode();
	NodeArrays=MultiBlock_->Get_Block()->Get_NodeArrays2D();
	IdBoundaries=MultiBlock_->Get_Block()->Get_PtrIdBc();
	//InvTau=1.0/PtrParameters->Get_Tau();
	nbvelo=Solver::PtrParameters->Get_NbVelocities();
	nbnode=MultiBlock_->Get_nnodes();
//	Gradients::initGradients(2,nbvelo,PtrParameters->Get_GradientType());
	Dic=new Dictionary(2,nbnode);
	Add_OneDistributionToDictionary();
	Set_Solution(PtrParameters);
/*
	//PtrVariablesOutput=new double* [PtrParameters->Get_NbVariablesOutput()];
	PtrVariablesBreakpoint=new double* [nbvelo+3];
	Solution2D::Set_output();
	//Solution2D::Set_output(PtrParameters->Get_PtrVariablesOutput(),PtrParameters->Get_NbVariablesOutput());
	SolverSinglePhase::set_f_name();///save name of variables (distribution function and macroscopic variables)
	SolverSinglePhase::set_f_ini();///save pointers of the distribution function
	///Set the variables names and the variable pointers for breakpoints in solution (the set of solution and writer are separated due to macroscopic variables are not set in the same class)
	Solution2D::Set_breakpoint(get_f_name(),nbvelo+3,get_f_ini()); //add macroscopic variables to save it
///Set the variables names and the variable pointers for breakpoints in writer
	Writer->Set_breakpoint(PtrVariablesBreakpoint, get_f_name(),  nbvelo+3);
*/
/*
///Set the variables names and the variable pointers for output in solution
	Solution2D::Set_output();
///Set the variables names and the variable pointers for breakpoints in solution
	Solution2D::Set_breakpoint();
*/
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
	parallel=PtrParallel_;
	Writer=PtrWriter_;
	Solver::PtrParameters=PtrParameters_;
	PtrBlockStream=MultiBlock_->Get_Block();
	NodeArrays=MultiBlock_->Get_Block()->Get_NodeArrays2D();
	IdBoundaries=MultiBlock_->Get_Block()->Get_PtrIdBc();
	//InvTau=1.0/PtrParameters->Get_Tau();
	nbvelo=Solver::PtrParameters->Get_NbVelocities();
	nbnode=MultiBlock_->Get_nnodes();
	Gradients::initGradients(2,nbvelo,Solver::PtrParameters->Get_GradientType());
	Dic=new Dictionary(2,nbnode);
	Add_TwoDistributionsToDictionary();
	Set_Solution(PtrParameters);
//	PtrVariablesBreakpoint=new double* [2*nbvelo+3];
/*	PtrVariablesOutput=new double* [PtrParameters->Get_NbVariablesOutput()+7];

	std::string* TmpPtrVarOut=new std::string [PtrParameters->Get_NbVariablesOutput()+7];
	for(int i=0;i<PtrParameters->Get_NbVariablesOutput();i++)
		TmpPtrVarOut[i]=PtrParameters->Get_PtrVariablesOutput()[i];
	TmpPtrVarOut[PtrParameters->Get_NbVariablesOutput()]="RhoN";
	TmpPtrVarOut[PtrParameters->Get_NbVariablesOutput()+1]="RhoRed";
	TmpPtrVarOut[PtrParameters->Get_NbVariablesOutput()+2]="RhoBlue";
	TmpPtrVarOut[PtrParameters->Get_NbVariablesOutput()+3]="ColourGradX";
	TmpPtrVarOut[PtrParameters->Get_NbVariablesOutput()+4]="ColourGradY";
	TmpPtrVarOut[PtrParameters->Get_NbVariablesOutput()+5]="InterfaceForceX";
	TmpPtrVarOut[PtrParameters->Get_NbVariablesOutput()+6]="InterfaceForceY";*/
/*
	Rhor=new double [nbnodes_total];
	Rhob=new double [nbnodes_total];
	RhoN=new double [nbnodes_total];
*/
/*	V1=new double* [2];
	V1[0]=new double [nbnodes_total];
	V1[1]=new double [nbnodes_total];

	V2=new double* [2];
	V2[0]=new double [nbnodes_total];
	V2[1]=new double [nbnodes_total];*/
/*
///Set the variables names and the variable pointers for output in solution
	Solution2D::Set_output();
///Set the variables names and the variable pointers for breakpoints in solution
	Solution2D::Set_breakpoint();
*/
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
