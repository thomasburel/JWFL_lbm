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
	f_name=0;
	f_ini=0;
}

SolverSinglePhase::~SolverSinglePhase() {
	delete f;
	delete [] ftmp;
}
double** SolverSinglePhase::get_f_ini(){
	return f_ini;
}
void SolverSinglePhase::set_f_ini(){
	if(f_ini!=0)
		delete f_ini;
	f_ini=new double*[nbvelo];
	for (int i=0;i<nbvelo;i++)
		f_ini[i]=&f->f[i][0];
}
std::string* SolverSinglePhase::get_f_name(){
	return f_name;
}
void SolverSinglePhase::set_f_name(){
	if(f_name!=0)
		delete f_name;
	f_name=new std::string[nbvelo+3];
	f_name[0]="Rho";
	f_name[1]="VelocityX";
	f_name[2]="VelocityY";
	for (int i=3;i<nbvelo+3;i++)
	{
		std::stringstream sstm;
		sstm<<"f_"<<i-3;
		f_name[i]=sstm.str();
	}
}
SolverTwoPhases::SolverTwoPhases() {
	PtrParameters=0;
	f=0;
	ftmp=0;
	f_name=0;
	f_ini=0;
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
double** SolverTwoPhases::get_f_ini(){
	return f_ini;
}
void SolverTwoPhases::set_f_ini(){
	if(f_ini!=0)
		delete f_ini;
	f_ini=new double*[nbvelo*2];
	for (int i=0;i<nbvelo;i++)
		f_ini[i]=&f[0]->f[i][0];
	for (int i=nbvelo;i<nbvelo*2;i++)
		f_ini[i]=&f[1]->f[i-nbvelo][0];
}
std::string* SolverTwoPhases::get_f_name(){
	return f_name;
}
void SolverTwoPhases::set_f_name(){
	if(f_name!=0)
		delete f_name;
	f_name=new std::string[2*nbvelo+3];
	f_name[0]="Rho";
	f_name[1]="VelocityX";
	f_name[2]="VelocityY";
	for (int i=3;i<2*nbvelo+3;i++)
	{
		std::stringstream sstm;
		sstm<<"f_"<<i-3;
		f_name[i]=sstm.str();
	}

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
	PtrParameters=0;

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
	PtrParameters=PtrParameters_;
	PtrBlockStream=MultiBlock_->Get_Block();
	PtrBlockCollide=PtrBlockStream;
	//Node=MultiBlock_->Get_Block()->Get_PtrNode();
	NodeArrays=MultiBlock_->Get_Block()->Get_NodeArrays2D();
	IdBoundaries=MultiBlock_->Get_Block()->Get_PtrIdBc();
	//InvTau=1.0/PtrParameters->Get_Tau();
	nbvelo=PtrParameters->Get_NbVelocities();
	nbnode=MultiBlock_->Get_nnodes();
	Set_Solution();
	PtrVariablesOutput=new double* [PtrParameters->Get_NbVariablesOutput()];
	PtrVariablesBreakpoint=new double* [nbvelo+3];
	Solution2D::Set_output(PtrParameters->Get_PtrVariablesOutput(),PtrParameters->Get_NbVariablesOutput());
	SolverSinglePhase::set_f_name();///save name of variables (distribution function and macroscopic variables)
	SolverSinglePhase::set_f_ini();///save pointers of the distribution function
	///Set the variables names and the variable pointers for breakpoints in solution (the set of solution and writer are separated due to macroscopic variables are not set in the same class)
	Solution2D::Set_breakpoint(get_f_name(),nbvelo+3,get_f_ini()); //add macroscopic variables to save it
///Set the variables names and the variable pointers for breakpoints in writer
	Writer->Set_breakpoint(PtrVariablesBreakpoint, get_f_name(),  nbvelo+3);

}
SolverTwoPhasesLowOrder2D::SolverTwoPhasesLowOrder2D() {
	PtrParameters=0;
	Gradients(2);

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
	PtrParameters=PtrParameters_;
	PtrBlockStream=MultiBlock_->Get_Block();
//	PtrBlockCollide=PtrBlockStream;
	//Node=MultiBlock_->Get_Block()->Get_PtrNode();
	NodeArrays=MultiBlock_->Get_Block()->Get_NodeArrays2D();
	IdBoundaries=MultiBlock_->Get_Block()->Get_PtrIdBc();
	//InvTau=1.0/PtrParameters->Get_Tau();
	nbvelo=PtrParameters->Get_NbVelocities();
	nbnode=MultiBlock_->Get_nnodes();
	Set_Solution();
	PtrVariablesOutput=new double* [PtrParameters->Get_NbVariablesOutput()];
	PtrVariablesBreakpoint=new double* [2*nbvelo+3];
	Solution2D::Set_output(PtrParameters->Get_PtrVariablesOutput(),PtrParameters->Get_NbVariablesOutput());
	SolverTwoPhases::set_f_name();///save name of variables (distribution function and macroscopic variables)
	SolverTwoPhases::set_f_ini();///save pointers of the distribution function
	///Set the variables names and the variable pointers for breakpoints in solution (the set of solution and writer are separated due to macroscopic variables are not set in the same class)
	Solution2D::Set_breakpoint(get_f_name(),2*nbvelo+3,get_f_ini()); //add macroscopic variables to save it
///Set the variables names and the variable pointers for breakpoints in writer
	Writer->Set_breakpoint(PtrVariablesBreakpoint, get_f_name(),  2*nbvelo+3);

}

