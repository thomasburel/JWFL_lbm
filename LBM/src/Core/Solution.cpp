/*
 * Solution.cpp
 *
 *  Created on: 30 Apr 2015
 *      Author: thomas
 */

#include "Solution.h"

Solution::Solution(){
	MultiBlock_=0;
	parallel=0;
	Writer=0;
	U=0;
	Rho=0;
	P=0;CalPressure=false;CalGradP=false;
	Dic=0;
	nbnodes_total=0;
	nbnodes_real=0;
	}
Solution::~Solution(){

	delete Dic;
	}
Solution2D::Solution2D() {

	NodeArrays=0;
	PatchsBc=0;

}

Solution2D::~Solution2D() {

//	delete NodeArrays;

}


void Solution2D::Set_Solution(Parameters *Param) {
	//std::cout<<"number of node: "<<Node->size()<<std::endl;
	nbnodes_total=NodeArrays->TypeOfNode.size();
	nbnodes_real=MultiBlock_->Get_End_Nodes()-MultiBlock_->Get_Start_Nodes()+1;//nbnodes_total;
//	std::cout<<"Processor: "<<parrallel->getRank()<<" Number of total node is: "<<nbnodes_total<<" Number of real node is: "<<nbnodes_real<<std::endl;
	U=new double* [2];
	U[0]=0;U[1]=0;


	Dic->AddVar(Vector,"Velocity",Param->Get_output_velocity(), true,false,U[0],U[1]);//No synchronisation between processor by default
	Dic->AddVar(Scalar,"Density",Param->Get_output_density(), true,false,Rho);//No synchronisation between processor by default
	if(Param->Get_output_pressure()||Param->IsCalculatePermeability())
	{
		Dic->AddVar(Scalar,"Pressure",Param->Get_output_pressure(), false,false,P);//No synchronisation between processor by default
		CalPressure=true;
		if(Param->IsCalculatePermeability())
		{
			Dic->AddSync("Pressure",P);
			CalGradP=true;
		}
	}

	/*	for (int i=0;i<2;i++)
	{
		U[i]=new double [nbnodes_total];
//		for (int j=0;j<nbnodes_total;j++)
//			U[i][j]=0;
	}*/
//	Rho=new double [nbnodes_total];
	for (int i=0;i<nbnodes_total;i++)
		Rho[i]=1;
	for (int j=0;j<nbnodes_total;j++)
	{
		U[0][j]=0;//Node->at(j)->get_x();
		U[1][j]=0;//Node->at(j)->get_y();
	}

}
void Solution2D::Set_output(){
	Writer->Set_solution(Dic->Get_PtrExportVar(),Dic->Get_PtrExportVarName(),Dic->Get_NbExportVar());
}
void Solution2D::Set_breakpoint()
{
	Writer->Set_breakpoint(Dic->Get_PtrExportVarBreakpoint(),Dic->Get_PtrExportVarBreakpointName(),Dic->Get_NbExportVarBreakpoint());
}

void Solution2D::Read_Variable(std::string variablename, std::string filename){
//	double *d=Dic->Get_PtrVar(Dic->Get_Id_Var(variablename));
	bool Var_found;
	Writer->Read_data(Dic->Get_PtrVar(Dic->Get_Id_Var(variablename,Var_found)),variablename,filename);
}
