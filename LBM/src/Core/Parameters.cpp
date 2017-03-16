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
	dimension_=SolverEnum::D2;
	scheme=SolverEnum::Q9;
	deltaT=1;deltaX=1;
	model=SolverEnum::SinglePhase;
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
	InitFromFile=false;
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
	ColourGrad=ColourFluidEnum::Gunstensen;
	Recolouring=ColourFluidEnum::LatvaKokkoRothman;
	ColourOperator=ColourFluidEnum::Grunau;
	ColourExtrapolNoramlDensity=true;ColourExtrapolNoramlInterface=true;
	PressureModelParam=HeZouP;
	PressureTypeParam=FixP;
	VelocityModelParam=HeZouV;
	VelocityTypeParam=FixV;
	CornerModelParam=HoChan;
	CornerPressureTypeParam=ExtrapolCP;
	tetaType=ContactAngleEnum::NoTeta;
	tetaModel=ContactAngleEnum::Interpol;
	NormalExtrapol=ContactAngleEnum::WeightDistanceExtrapol;
	NormalInterpol=ContactAngleEnum::LinearLeastSquareInterpol;
	Switchteta=ContactAngleEnum::Binary;
	teta=std::acos(-1)/2;// pi/2
	ErrorMax=1e-10;
	PeriodicTypeParam=Simple;
	PressureDropParam=0;
	RhoLimiter=1e-7;
	GNormLimiter=0.01;
	nNodesInterpolInSolid=1;
	nNodesInterpolInFluid=1;
	refDensity=1;
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
void Parameters::Change_Dimension(SolverEnum::dimension dim_)
{
	if (dimension_!=dim_)
	{
		delete GlobalBcType;
		dimension_=dim_;
		if (dimension_==SolverEnum::D2)
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

void MeshParameters::Set_Dimension(SolverEnum::dimension dim_)
{
		dimension_=dim_;

}
SolverEnum::dimension MeshParameters::Get_Dimension() const {
	return dimension_;
}

void Parameters::Set_BcType(NodeType Bc0,NodeType Bc1,NodeType Bc2,NodeType Bc3,NodeType Bc4,NodeType Bc5){
	//Force periodic if only one side is setup as periodic
	if (Bc0==Periodic||Bc2==Periodic)
		if(Bc0!=Bc2)
		{
			Bc0=Periodic;
			Bc2=Periodic;
		}
	if (Bc1==Periodic||Bc3==Periodic)
		if(Bc1!=Bc3)
		{
			Bc1=Periodic;
			Bc3=Periodic;
		}
	if (MeshParameters::Get_Dimension()==SolverEnum::D2)
		{
			NbGlobalBcType=4;
			GlobalBcType[0]=Bc0;
			GlobalBcType[1]=Bc1;
			GlobalBcType[2]=Bc2;
			GlobalBcType[3]=Bc3;
		}
		else
		{
			//Force periodic if only one side is setup as periodic
			if (Bc4==Periodic||Bc5==Periodic)
				if(Bc4!=Bc5)
				{
					Bc4=Periodic;
					Bc5=Periodic;
				}
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

void SolverParameters::Set_Scheme(SolverEnum::schemetype scheme_)
{
	scheme=scheme_;
}
SolverEnum::schemetype SolverParameters::Get_Scheme() const
{
	return scheme;
}
void SolverParameters::Set_Model(SolverEnum::modeltype model_)
{
	model=model_;
}
SolverEnum::modeltype SolverParameters::Get_Model() const
{
	return model;
}
void SolverParameters::Set_NumberOfInterpolNode(int nNodes){
	if(nNodes<2)
	{
		std::cout<<"Number of node for interpolation is not enough. It will be set to 2."<<std::endl;
		nNodesInterpolInSolid=1;
		nNodesInterpolInFluid=1;
	}
	nNodesInterpolInSolid=round(nNodes/2.0);
	nNodesInterpolInFluid=nNodes-nNodesInterpolInSolid;
}
void SolverParameters::Set_NumberOfInterpolNodeInSolid(int nNodes){
	nNodesInterpolInSolid=nNodes;
}
void SolverParameters::Set_NumberOfInterpolNodeInFluid(int nNodes){
	nNodesInterpolInFluid=nNodes;
}
void CalculationParameters::Set_Parallel(parralleltype parrallel_) {
	parrallel=parrallel_;
}
parralleltype CalculationParameters::Get_Parallel() const {
	return parrallel;
}

void MeshParameters::set_MeshParameters()
{
	if(dimension_==SolverEnum::D2)
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
		case SolverEnum::Q9:
			NbVelocities=9;
			cs=1/std::sqrt(3)*deltaX/deltaT;
			cs2=cs*cs;
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
void Parameters::Add_VariableToInit(std::string filename,SolverEnum::variablesSolver variablesSolvertype){
	switch(variablesSolvertype)
	{
	case SolverEnum::Density:
		InitFromFile=true;
		NameInitFromFile.push_back(filename);
		VariableInitFromFile.push_back("Density");
		NbVariablesInitFromFile=VariableInitFromFile.size();
		break;
	case  SolverEnum::VelocityX:
		InitFromFile=true;
		NameInitFromFile.push_back(filename);
		VariableInitFromFile.push_back("VelocityX");
		NbVariablesInitFromFile=VariableInitFromFile.size();
		break;
	case  SolverEnum::VelocityY:
		InitFromFile=true;
		NameInitFromFile.push_back(filename);
		VariableInitFromFile.push_back("VelocityY");
		NbVariablesInitFromFile=VariableInitFromFile.size();
		break;
	case  SolverEnum::VelocityZ:
		if(dimension_==SolverEnum::D3){
			InitFromFile=true;
			NameInitFromFile.push_back(filename);
			VariableInitFromFile.push_back("VelocityZ");
			NbVariablesInitFromFile=VariableInitFromFile.size();
		}
		break;
	case  SolverEnum::RhoN:
		if(model==SolverEnum::ColourFluid){
			InitFromFile=true;
			NameInitFromFile.push_back(filename);
			VariableInitFromFile.push_back("RhoN");
			NbVariablesInitFromFile=VariableInitFromFile.size();
		}
		break;
	default:
		std::cout<<"Variable to Init is not found"<<std::endl;
	}
}
void Parameters::CheckParameters()
{
	switch(scheme)
		{
		case SolverEnum::Q9:
			NbVelocities=9;
			cs=1/std::sqrt(3)*deltaX/deltaT;
			cs2=cs*cs;
			break;
		default:
			std::cout<< "Scheme not found" << std::endl;
	        break;
		}
}
