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
	Error=0;
	Error_tmp=0;
	Error_sum=0;
	Error_avg=0;
	Sum_Current=1;

	Scalar_last=0;
	Vector_last=0;
	NbNodes=0;
	Scalar_CurrentTime=0;
	Vector_CurrentTime=0;


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
		std::cout<<"Error variable type not found. Density will be used."<<std::endl;
		PtrDicConv->Get_PtrVar("Density",Scalar_CurrentTime);
	}
	NbNodes=PtrDicConv->Get_NbNodes();
	PtrDicConv->AddVar(Scalar,"ScalarError_last",false,true,false,Scalar_last);
	for(int i=0;i<NbNodes;i++)
	{
		Scalar_last[i]=1.0;
	}
}
void Convergence::Calcul_Error(){

	Calcul_Error_ScalarField();

}
void Convergence::Calcul_Error_ScalarField(){
	Error=Calcul_Error_ScalarFieldInOneProc();
	Error=PtrMultiBlockConv->SumAllProcessors(&Error);
	Error=Error/PtrMultiBlockConv->NumberOfProcessors();
}
double Convergence::Calcul_Error_ScalarFieldInOneProc(){
	Sum_Current=0;Error_tmp=0;
	for(int i=0;i<NbNodes;i++)
	{
		Error_tmp+=std::abs(Scalar_CurrentTime[i]-Scalar_last[i]);
		Sum_Current+=std::abs(Scalar_CurrentTime[i]);
		Scalar_last[i]=Scalar_CurrentTime[i];//Save the new value
	}
	Error_tmp=Error_tmp/(Sum_Current+1e-15);//std::abs(Scalar_CurrentTime[idx]-Scalar_last[idx])/(std::abs(Scalar_CurrentTime[idx])+1e-15);

	return Error_tmp;
}
void Convergence::Calcul_PorousMediaConvergence(){
	//Calcul EOR and write it in a file
	double EOR=Calcul_EOR();
	EOR=PtrMultiBlockConv->SumAllProcessors(&EOR);
	EOR=EOR/PtrMultiBlockConv->NumberOfProcessors();
	ofstream EORfile;
	EORfile.open("EOR.txt",ios::out | ios::app);
	EORfile<<EOR<<std::endl;
	EORfile.close();
}
double Convergence::Calcul_EOR(){

	double sum=0;
	for (int i=0;i<OutletPatchId.size();i++)
	{
		PtrPatchBcConv->Get_NodeIndex(OutletPatchId[i],*NodeId,*NodeIdSpeWall,*NodeIdGloCorner);
		Avg_ScalarNodeArray(*NodeId,*NodeIdSpeWall,*NodeIdGloCorner,RhoNEOR,avg);
		sum+=avg;
	}
	return sum/(double)OutletPatchId.size();
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
	for(int i=0;i<NodeArray.size();i++)
		sum+=Var1[NodeArray[i]];
	double sum_tmp=0;
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
