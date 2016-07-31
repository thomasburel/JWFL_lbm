/*
 * Solution.cpp
 *
 *  Created on: 30 Apr 2015
 *      Author: thomas
 */

#include "Solution.h"

Solution2D::Solution2D() {
	MultiBlock_=0;
	parallel=0;
	Writer=0;
	U=0;
	Rho=0;
	//Node=0;
	nbnodes_total=0;
	nbnodes_real=0;
	PtrVariablesOutput=0;
	PtrVariablesBreakpoint=0;
	Ptrvariabletest=0;
	NodeArrays=0;
}

Solution2D::~Solution2D() {
	for (int i=0;i<2;i++)
	{
		for (int j=0;j<nbnodes_total;j++)
			delete U[i];
	}
	delete U,Rho,NodeArrays;
}


void Solution2D::Set_Solution() {
	//std::cout<<"number of node: "<<Node->size()<<std::endl;
	nbnodes_total=NodeArrays->TypeOfNode.size();
	nbnodes_real=MultiBlock_->Get_End_Nodes()-MultiBlock_->Get_Start_Nodes()+1;//nbnodes_total;
//	std::cout<<"Processor: "<<parrallel->getRank()<<" Number of total node is: "<<nbnodes_total<<" Number of real node is: "<<nbnodes_real<<std::endl;
	U=new double* [2];
	for (int i=0;i<2;i++)
	{
		U[i]=new double [nbnodes_total];
//		for (int j=0;j<nbnodes_total;j++)
//			U[i][j]=0;
	}
	Rho=new double [nbnodes_total];
	for (int i=0;i<nbnodes_total;i++)
		Rho[i]=1;
	for (int j=0;j<nbnodes_total;j++)
	{
		U[0][j]=0;//Node->at(j)->get_x();
		U[1][j]=0;//Node->at(j)->get_y();
	}

}
void Solution2D::Set_output(std::string *str, int nbvar){
	for (int i=0;i<nbvar;i++)
		if(str[i]=="Rho")
			PtrVariablesOutput[i]=&Rho[0];
		else if (str[i]=="VelocityX")
			PtrVariablesOutput[i]=&U[0][0];//Ptrvariabletest;//&U[0][0];
		else if(str[i]=="VelocityY")
			PtrVariablesOutput[i]=&U[1][0];
		else
			std::cout<< "Problem in output setup" << std::endl;


}
void Solution2D::Set_breakpoint(std::string *str, int nbvar,double **f_ini){
	int i_tmp=0;
	std::stringstream sstm;
	sstm<<"f_"<<i_tmp;
	for (int i=0;i<nbvar;i++)
	{
		std::string str_tmp=sstm.str();
		if(str[i]=="Rho")
			PtrVariablesBreakpoint[i]=&Rho[0];
		else if (str[i]=="VelocityX")
			PtrVariablesBreakpoint[i]=&U[0][0];
		else if (str[i]=="VelocityY")
			PtrVariablesBreakpoint[i]=&U[1][0];
		else if (str[i]==str_tmp)
		{
			PtrVariablesBreakpoint[i]=f_ini[i_tmp];
			i_tmp++;
			//empty stringstream
			sstm.clear();
			sstm.str(std::string());
			// set the next value for sstm
			sstm<<"f_"<<i_tmp;
		}
		else
			std::cout<< "Problem in breakpoint setup" << std::endl;
	}


}
