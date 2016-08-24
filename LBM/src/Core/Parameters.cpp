/*
 * Parameters.cpp
 *
 *  Created on: 5 May 2015
 *      Author: thomas
 */

#include "Parameters.h"

Parameters::Parameters()
{
	nx=1;
	ny=1;
	nz=1;
	dimension_=D2;
	scheme=Q9;
	model=SinglePhase;
	fluid=Newtonian;
	UserForce=None;
	argc=0;
	argv=0;
	verbous=false;
	parrallel=Serial;
	VariablesOutput=0;
	SolverParameters::set_SolverParameters();
	Set_VariablesOutput(true,true);//export Rho and U by default
	NbGlobalBcType=4;
	GlobalBcType=new NodeType [NbGlobalBcType];// Bottom: Wall, Outlet: Pressure, Top: Wall, Inlet: Velocity
	GlobalBcType[0]=Wall;
	GlobalBcType[1]=Pressure ;
	GlobalBcType[2]=Wall;
	GlobalBcType[3]=Velocity;
	restart=false;
	Format=CGNSFormat;
	OutputFileName="output";
	Rho_1=1;
	Rho_2=1;
	Tau_1=1;
	Tau_2=1;
	tension=0;
	beta=0.7;
	A1=0; A2=0;
	U_ini=0;V_ini=0;W_ini=0;
	Re=1;
	l=1;t=1;dx=1;dt=1;
	ColourGrad=Gunstensen;
	Recolouring=LatvaKokkoRothman;
	ColourOperator=Grunau;

}

Parameters::~Parameters() {
	delete [] GlobalBcType;
	delete [] VariablesOutput;
}
void Parameters::Set_Arguments(int *argc_,char ***argv_,bool verbous_) {
	argc=argc_;
	argv=argv_;
	verbous=verbous_;
}
int* Parameters::Get_Argc() const {
	return argc;
}

char*** Parameters::Get_Argv() const {
	return argv;
}

bool Parameters::Get_Verbous() const {
	return verbous;
}
void Parameters::Change_Dimension(dimension dim_)
{
	if (dimension_!=dim_)
	{
		delete GlobalBcType;
		dimension_=dim_;
		if (dimension_==D2)
		{
			NbGlobalBcType=4;
			GlobalBcType=new NodeType [NbGlobalBcType];// Bottom: Wall, Outlet: Pressure, Top: Wall, Inlet: Velocity
			GlobalBcType[0]=Wall;
			GlobalBcType[1]=Pressure;
			GlobalBcType[2]=Wall;
			GlobalBcType[3]=Velocity;
		}
		else
		{
			NbGlobalBcType=6;
			GlobalBcType=new NodeType [NbGlobalBcType];// Bottom: Wall, Outlet: Pressure, Top: Wall, Inlet: Velocity, Front: Wall, Back: Wall
			GlobalBcType[0]=Wall;
			GlobalBcType[1]=Pressure;
			GlobalBcType[2]=Wall;
			GlobalBcType[3]=Velocity;
			GlobalBcType[4]=Wall;
			GlobalBcType[5]=Wall;
		}
	}

}

void MeshParameters::Set_Dimension(dimension dim_)
{
		dimension_=dim_;

}
dimension MeshParameters::Get_Dimension() const {
	return dimension_;
}

void Parameters::Set_BcType(NodeType Bc0,NodeType Bc1,NodeType Bc2,NodeType Bc3,NodeType Bc4,NodeType Bc5){
	if (MeshParameters::Get_Dimension()==D2)
			{
				NbGlobalBcType=4;
				GlobalBcType[0]=Bc0;
				GlobalBcType[1]=Bc1;
				GlobalBcType[2]=Bc2;
				GlobalBcType[3]=Bc3;
			}
			else
			{
				NbGlobalBcType=6;
				GlobalBcType[0]=Bc0;
				GlobalBcType[1]=Bc1;
				GlobalBcType[2]=Bc2;
				GlobalBcType[3]=Bc3;
				GlobalBcType[4]=Bc4;
				GlobalBcType[5]=Bc5;
			}
}

void MeshParameters::Set_Domain_Size(int nx_,int ny_,int nz_) {
	nx=nx_;
	ny=ny_;
	nz=nz_;
}
int MeshParameters::Get_Nx() const {
	return nx;
}
int MeshParameters::Get_Ny() const {
	return ny;
}
int MeshParameters::Get_Nz() const {
	return nz;
}

void MeshParameters::Set_MeshFile(string MeshFile_) {
	MeshFile=MeshFile_;
}
string MeshParameters::Get_MeshFile() const {
	return MeshFile;
}

void SolverParameters::Set_Scheme(schemetype scheme_)
{
	scheme=scheme_;
}
schemetype SolverParameters::Get_Scheme() const
{
	return scheme;
}
void SolverParameters::Set_Model(modeltype model_)
{
	model=model_;
}
modeltype SolverParameters::Get_Model() const
{
	return model;
}
void CalculationParameters::Set_Parallel(parralleltype parrallel_) {
	parrallel=parrallel_;
}
parralleltype CalculationParameters::Get_Parallel() const {
	return parrallel;
}

void MeshParameters::set_MeshParameters()
{
	if(dimension_==D2)
	{
		NbNodes=(nx+1)*(ny+1);
	}
	else
	{
		NbNodes=(nx+1)*(ny+1)*(nz+1);

	}
}

void SolverParameters::set_SolverParameters()
{
	switch(scheme)
		{
		case Q9:
			NbVelocities=9;
			break;
		default:
			std::cout<< "Scheme not found" << std::endl;
	        break;
		}

}

int SolverParameters::Get_NbVelocities() const {
	return NbVelocities;
}
void CalculationParameters::Set_NbStep(int NbSteps_){
	NbSteps=NbSteps_;
}
int CalculationParameters::Get_NbStep() const{
	return NbSteps;
}
void CalculationParameters::Set_OutPutNSteps(int NbSteps_){
	IntervalOutput=NbSteps_;
}
int CalculationParameters::Get_OutPutNSteps() const{
	return IntervalOutput;
}

void CalculationParameters::Set_VariablesOutput(int nbvar, std::string * strinput)
{
	NbVariablesOutput=nbvar;
	if(VariablesOutput!=0)
		delete [] VariablesOutput;
	VariablesOutput=new string[NbVariablesOutput];
//	if(VariablesOutputSeria==0)
//		delete VariablesOutputSeria;
//	VariablesOutputSeria=new XStringXml[NbVariablesOutput];
	for (int i=0;i<nbvar;i++)
	{
		VariablesOutput[i]=strinput[i];
//		VariablesOutputSeria[i]=VariablesOutput[i];
	}
}
void Parameters::Set_VariablesOutput (bool Rho, bool U) {
	std::string *strtmp=0;
	int nbvar=0;
	if(MeshParameters::Get_Dimension()==D2)
		if (Rho)
			if(U)
			{
				nbvar=3;
				strtmp=new string[nbvar];
				strtmp[0]="Rho"; strtmp[1]="VelocityX"; strtmp[2]="VelocityY";
				CalculationParameters::Set_VariablesOutput(nbvar,strtmp);
				density=true; velocity=true;
			}
			else
			{
				nbvar=1;
				strtmp=new string[nbvar];
				strtmp[0]="Rho";
				CalculationParameters::Set_VariablesOutput(nbvar,strtmp);
				density=true; velocity=false;
			}
		else
			if(U)
			{
				nbvar=2;
				strtmp=new string[nbvar];
				strtmp[0]="VelocityX"; strtmp[1]="VelocityY";
				CalculationParameters::Set_VariablesOutput(nbvar,strtmp);
				density=false; velocity=true;
			}
			else //no variables
			{
				nbvar=0;
				CalculationParameters::Set_VariablesOutput(nbvar,strtmp);
				density=false; velocity=false;
			}
	else
		if (Rho)
			if(U)
			{
				nbvar=4;
				strtmp=new string[nbvar];
				strtmp[0]="Rho"; strtmp[1]="VelocityX"; strtmp[2]="VelocityY"; strtmp[3]="VelocityZ";
				CalculationParameters::Set_VariablesOutput(nbvar,strtmp);
				density=true; velocity=true;
			}
			else
			{
				nbvar=1;
				strtmp=new string[nbvar];
				strtmp[0]="Rho";
				CalculationParameters::Set_VariablesOutput(nbvar,strtmp);
				density=true; velocity=false;
			}
		else
			if(U)
			{
				nbvar=3;
				strtmp=new string[nbvar];
				strtmp[0]="VelocityX"; strtmp[1]="VelocityY"; strtmp[2]="VelocityZ";
				density=false; velocity=true;
				CalculationParameters::Set_VariablesOutput(nbvar,strtmp);
			}
			else//no variables
			{
				nbvar=0;
				CalculationParameters::Set_VariablesOutput(nbvar,strtmp);
				density=false; velocity=false;
			}

}
std::string* CalculationParameters::Get_PtrVariablesOutput() const{
	return VariablesOutput;
}
int CalculationParameters::Get_NbVariablesOutput() const {
	return NbVariablesOutput;

}

void BoundaryParameters::Set_WallType(WallType WallType_){ WallTypeParam= WallType_;}
WallType BoundaryParameters::Get_WallType() const {return WallTypeParam;}
void BoundaryParameters::Set_SymmetryType(SymmetryType SymmetryType_){ SymmetryTypeParam= SymmetryType_;}
SymmetryType BoundaryParameters::Get_SymmetryType() const {return SymmetryTypeParam;}
NodeType BoundaryParameters::Get_GlobalBcType(int side) const {
	return GlobalBcType[side];
}
string CalculationParameters::Get_OutputFileName() const{ return OutputFileName;}
void CalculationParameters::Set_OutputFileName(string OutputFileName_){OutputFileName=OutputFileName_;}

void  CalculationParameters::Set_listing(int IntervalListing_){IntervalListing=IntervalListing_;}
int  CalculationParameters::Get_listing() const {return IntervalListing;}


void PhysicalParameters::Set_deltax(double dx_)
{
	dx=dx_;
}
double PhysicalParameters::Get_deltax() const {
	return dx;
}
