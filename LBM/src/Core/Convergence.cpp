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

	Error=Calcul_Error_ScalarField();
	Error=PtrMultiBlockConv->SumAllProcessors(&Error);
	Error=Error/PtrMultiBlockConv->NumberOfProcessors();

}
double Convergence::Calcul_Error_ScalarField(){
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
