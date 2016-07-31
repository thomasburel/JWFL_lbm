/*
 * D2Q9TwoPhases.cpp
 *
 *  Created on: 9 Jun 2015
 *      Author: thomas
 */

#include "D2Q9_TwoPhases.h"
#include <iomanip>
D2Q9TwoPhases::D2Q9TwoPhases() {
	// TODO Auto-generated constructor stub
	MultiBlock_=0;
	parallel=0;
	Writer=0;
	PtrParameters=0;
	f=0;
	nbvelo=9;
	nbnode=0;
	ftmp=0;
	tmp=0;
	tmpreturn=0;
	intTmpReturn=0;
	doubleTmpReturn=0;
	buf_recv=0;
	buf_send=0;
	size_buf=0;
	DiagConnect=0;
	Nd_variables_sync=0;
	buf_MacroSend=0;
	buf_MacroRecv=0;
	size_MacroBuf=0;
	Nd_MacroVariables_sync=0;
	Rhor=0;Rhob=0;
}
D2Q9TwoPhases::D2Q9TwoPhases(MultiBlock* MultiBlock__,ParallelManager* parallel__,WriterManager* Writer__, Parameters* Parameters_ ,InitLBM& ini) {

	// TODO Auto-generated constructor stub	Ptrvariabletest=&f[0]->f[1][0];

	f=new DistriFunct*[2];
	f[0]=new DistriFunct(MultiBlock__->Get_nnodes(),Parameters_->Get_NbVelocities());
	f[1]=new DistriFunct(MultiBlock__->Get_nnodes(),Parameters_->Get_NbVelocities());
//	Ptrvariabletest=&f[0]->f[1][0];
	Set_Solver(MultiBlock__,parallel__,Writer__,Parameters_);
	ftmp=new double[nbnode];
	PtrFiStream=f[0];
	PtrFiCollide=PtrFiStream;
	InvTau=1.0/PtrParameters->Get_Tau();
	intTmpReturn=0;
	doubleTmpReturn=0;

	//Distribution velocity
	Ei[0][0]= 0.0;
	Ei[0][1]= 0.0;
	Ei[1][0]= 1.0;
	Ei[1][1]= 0.0;
	Ei[2][0]= 0.0;
	Ei[2][1]= 1.0;
	Ei[3][0]= -1.0;
	Ei[3][1]= 0.0;
	Ei[4][0]= 0.0;
	Ei[4][1]= -1.0;
	Ei[5][0]= 1.0;
	Ei[5][1]= 1.0;
	Ei[6][0]= -1.0;
	Ei[6][1]= 1.0;
	Ei[7][0]= -1.0;
	Ei[7][1]= -1.0;
	Ei[8][0]= 1.0;
	Ei[8][1]= -1.0;
	//Weigh depending of the direction
	omega[0]=4.0/9.0;
	omega[1]=1.0/9.0;
	omega[2]=1.0/9.0;
	omega[3]=1.0/9.0;
	omega[4]=1.0/9.0;
	omega[5]=1.0/36.0;
	omega[6]=1.0/36.0;
	omega[7]=1.0/36.0;
	omega[8]=1.0/36.0;
	// Opposite direction
	Opposite[0]=0;
	Opposite[1]=3;
	Opposite[2]=4;
	Opposite[3]=1;
	Opposite[4]=2;
	Opposite[5]=7;
	Opposite[6]=8;
	Opposite[7]=5;
	Opposite[8]=6;
	//Diffuse variables
	SumWeightS=omega[4]+omega[7]+omega[8];
	SumWeightE=omega[1]+omega[5]+omega[8];
	SumWeightN=omega[2]+omega[5]+omega[6];
	SumWeightW=omega[3]+omega[6]+omega[7];
	SumWeightConcaveSE=omega[1]+omega[4]+omega[8];
	SumWeightConcaveNE=omega[1]+omega[2]+omega[5];
	SumWeightConcaveNW=omega[2]+omega[3]+omega[6];
	SumWeightConcaveSW=omega[3]+omega[4]+omega[7];
	SumWeightConvexSE=omega[2]+omega[3]+omega[6];
	SumWeightConvexNE=omega[3]+omega[4]+omega[7];
	SumWeightConvexNW=omega[1]+omega[4]+omega[8];
	SumWeightConvexSW=omega[1]+omega[2]+omega[5];
	IniComVariables();
	init(ini);
}
D2Q9TwoPhases::~D2Q9TwoPhases() {
	delete [] size_buf;

	for (int i=0;i<Nd_variables_sync;i++)
	{
		delete [] buf_send[i][0];
		delete [] buf_send[i][1];
		delete [] buf_send[i][2];
		delete [] buf_send[i][3];
		delete [] buf_recv[i][0];
		delete [] buf_recv[i][1];
		delete [] buf_recv[i][2];
		delete [] buf_recv[i][3];
	}
	for (int i=0;i<Nd_variables_sync;i++)
	{
		delete [] buf_send[i];
		delete [] buf_recv[i];
	}
	delete [] Rhor;
	delete [] Rhob;
}
void D2Q9TwoPhases::init(InitLBM& ini){
	// Multiphase variables
	Rhor=new double [nbnodes_total];
	Rhob=new double [nbnodes_total];
	double alpha=0;
	double* pos =new double[2];
	double* U_=new double[2];
	for (int j=0;j<NodeArrays->NodeInterior.size();j++)
	{
		pos[0]=NodeArrays->NodeInterior[j].get_x();
		pos[1]=NodeArrays->NodeInterior[j].get_y();
		ini.IniDomain(parallel->getRank(),NodeArrays->NodeInterior[j],0, NodeArrays->NodeInterior[j].Get_index(),pos,Rho[NodeArrays->NodeInterior[j].Get_index()],U_,alpha);
//		Rho[NodeArrays->NodeInterior[j].Get_index()]=Convert_Alpha_To_Rho(alpha);
		Rhor[NodeArrays->NodeInterior[j].Get_index()]=alpha*Rho[NodeArrays->NodeInterior[j].Get_index()];//PtrParameters->Get_Rho_1();
		Rhob[NodeArrays->NodeInterior[j].Get_index()]=(1- alpha) *Rho[NodeArrays->NodeInterior[j].Get_index()];//PtrParameters->Get_Rho_2();
		U[0][NodeArrays->NodeInterior[j].Get_index()]=U_[0];
		U[1][NodeArrays->NodeInterior[j].Get_index()]=U_[1];
	}
	for (int j=0;j<NodeArrays->NodeCorner.size();j++)
	{
		pos[0]=NodeArrays->NodeCorner[j].get_x();
		pos[1]=NodeArrays->NodeCorner[j].get_y();
		ini.IniDomain(parallel->getRank(),NodeArrays->NodeCorner[j],0, NodeArrays->NodeCorner[j].Get_index(),pos,Rho[NodeArrays->NodeCorner[j].Get_index()],U_,alpha);
//		Rho[NodeArrays->NodeCorner[j].Get_index()]=Convert_Alpha_To_Rho(alpha);
		Rhor[NodeArrays->NodeCorner[j].Get_index()]=alpha*Rho[NodeArrays->NodeCorner[j].Get_index()];//PtrParameters->Get_Rho_1();
		Rhob[NodeArrays->NodeCorner[j].Get_index()]=(1- alpha) *Rho[NodeArrays->NodeCorner[j].Get_index()];//PtrParameters->Get_Rho_2();
		U[0][NodeArrays->NodeCorner[j].Get_index()]=U_[0];
		U[1][NodeArrays->NodeCorner[j].Get_index()]=U_[1];
	}
	for (int j=0;j<NodeArrays->NodeGlobalCorner.size();j++)
	{
		pos[0]=NodeArrays->NodeGlobalCorner[j].get_x();
		pos[1]=NodeArrays->NodeGlobalCorner[j].get_y();
		ini.IniDomain(parallel->getRank(),NodeArrays->NodeGlobalCorner[j],0, NodeArrays->NodeGlobalCorner[j].Get_index(),pos,Rho[NodeArrays->NodeGlobalCorner[j].Get_index()],U_,alpha);
//		Rho[NodeArrays->NodeGlobalCorner[j].Get_index()]=Convert_Alpha_To_Rho(alpha);
		Rhor[NodeArrays->NodeGlobalCorner[j].Get_index()]=alpha*Rho[NodeArrays->NodeGlobalCorner[j].Get_index()];//PtrParameters->Get_Rho_1();
		Rhob[NodeArrays->NodeGlobalCorner[j].Get_index()]=(1- alpha) *Rho[NodeArrays->NodeGlobalCorner[j].Get_index()];//PtrParameters->Get_Rho_2();
		U[0][NodeArrays->NodeGlobalCorner[j].Get_index()]=U_[0];
		U[1][NodeArrays->NodeGlobalCorner[j].Get_index()]=U_[1];
		NodeArrays->NodeGlobalCorner[j].Set_UDef(U_[0],U_[1]);
		NodeArrays->NodeGlobalCorner[j].Set_RhoDef(Rho[NodeArrays->NodeGlobalCorner[j].Get_index()]);
		NodeArrays->NodeGlobalCorner[j].Set_AlphaDef(alpha);
	}
	for (int j=0;j<NodeArrays->NodeVelocity.size();j++)
	{
		pos[0]=NodeArrays->NodeVelocity[j].get_x();
		pos[1]=NodeArrays->NodeVelocity[j].get_y();
		ini.IniDomain(parallel->getRank(),NodeArrays->NodeVelocity[j],0, NodeArrays->NodeVelocity[j].Get_index(),pos,Rho[NodeArrays->NodeVelocity[j].Get_index()],U_,alpha);
//		Rho[NodeArrays->NodeVelocity[j].Get_index()]=Convert_Alpha_To_Rho(alpha);
		Rhor[NodeArrays->NodeVelocity[j].Get_index()]=alpha*Rho[NodeArrays->NodeVelocity[j].Get_index()];//PtrParameters->Get_Rho_1();
		Rhob[NodeArrays->NodeVelocity[j].Get_index()]=(1- alpha) *Rho[NodeArrays->NodeVelocity[j].Get_index()];//PtrParameters->Get_Rho_2();
		U[0][NodeArrays->NodeVelocity[j].Get_index()]=U_[0];
		U[1][NodeArrays->NodeVelocity[j].Get_index()]=U_[1];
		NodeArrays->NodeVelocity[j].Set_UDef(U_[0],U_[1]);
	}

	for (int j=0;j<NodeArrays->NodePressure.size();j++)
	{
		pos[0]=NodeArrays->NodePressure[j].get_x();
		pos[1]=NodeArrays->NodePressure[j].get_y();
		ini.IniDomain(parallel->getRank(),NodeArrays->NodePressure[j],0, NodeArrays->NodePressure[j].Get_index(),pos,Rho[NodeArrays->NodePressure[j].Get_index()],U_,alpha);
//		Rho[NodeArrays->NodePressure[j].Get_index()]=Convert_Alpha_To_Rho(alpha);
		Rhor[NodeArrays->NodePressure[j].Get_index()]=alpha*Rho[NodeArrays->NodePressure[j].Get_index()];//PtrParameters->Get_Rho_1();
		Rhob[NodeArrays->NodePressure[j].Get_index()]=(1- alpha) *Rho[NodeArrays->NodePressure[j].Get_index()];//PtrParameters->Get_Rho_2();
		U[0][NodeArrays->NodePressure[j].Get_index()]=U_[0];
		U[1][NodeArrays->NodePressure[j].Get_index()]=U_[1];
		NodeArrays->NodePressure[j].Set_RhoDef(Rho[NodeArrays->NodePressure[j].Get_index()]);
		NodeArrays->NodePressure[j].Set_AlphaDef(alpha);
	}
	for (int j=0;j<NodeArrays->NodeWall.size();j++)
	{
		pos[0]=NodeArrays->NodeWall[j].get_x();
		pos[1]=NodeArrays->NodeWall[j].get_y();
		ini.IniDomain(parallel->getRank(),NodeArrays->NodeWall[j],0, NodeArrays->NodeWall[j].Get_index(),pos,Rho[NodeArrays->NodeWall[j].Get_index()],U_,alpha);
//		Rho[NodeArrays->NodeWall[j].Get_index()]=Convert_Alpha_To_Rho(alpha);
		Rhor[NodeArrays->NodeWall[j].Get_index()]=alpha*Rho[NodeArrays->NodeWall[j].Get_index()];//PtrParameters->Get_Rho_1();
		Rhob[NodeArrays->NodeWall[j].Get_index()]=(1- alpha) *Rho[NodeArrays->NodeWall[j].Get_index()];//PtrParameters->Get_Rho_2();
		U[0][NodeArrays->NodeWall[j].Get_index()]=U_[0];
		U[1][NodeArrays->NodeWall[j].Get_index()]=U_[1];
	}
	for (int j=0;j<NodeArrays->NodeSpecialWall.size();j++)
	{
		pos[0]=NodeArrays->NodeSpecialWall[j].get_x();
		pos[1]=NodeArrays->NodeSpecialWall[j].get_y();
		ini.IniDomain(parallel->getRank(),NodeArrays->NodeSpecialWall[j],0, NodeArrays->NodeSpecialWall[j].Get_index(),pos,Rho[NodeArrays->NodeSpecialWall[j].Get_index()],U_,alpha);
//		Rho[NodeArrays->NodeSpecialWall[j].Get_index()]=Convert_Alpha_To_Rho(alpha);
		Rhor[NodeArrays->NodeSpecialWall[j].Get_index()]=alpha*Rho[NodeArrays->NodeSpecialWall[j].Get_index()];//PtrParameters->Get_Rho_1();
		Rhob[NodeArrays->NodeSpecialWall[j].Get_index()]=(1- alpha) *Rho[NodeArrays->NodeSpecialWall[j].Get_index()];//PtrParameters->Get_Rho_2();
		U[0][NodeArrays->NodeSpecialWall[j].Get_index()]=U_[0];
		U[1][NodeArrays->NodeSpecialWall[j].Get_index()]=U_[1];
	}
	for (int j=0;j<NodeArrays->NodeSymmetry.size();j++)
	{
		pos[0]=NodeArrays->NodeSymmetry[j].get_x();
		pos[1]=NodeArrays->NodeSymmetry[j].get_y();
		ini.IniDomain(parallel->getRank(),NodeArrays->NodeSymmetry[j],0, NodeArrays->NodeSymmetry[j].Get_index(),pos,Rho[NodeArrays->NodeSymmetry[j].Get_index()],U_,alpha);
//		Rho[NodeArrays->NodeSymmetry[j].Get_index()]=Convert_Alpha_To_Rho(alpha);
		Rhor[NodeArrays->NodeSymmetry[j].Get_index()]=alpha*Rho[NodeArrays->NodeSymmetry[j].Get_index()];//PtrParameters->Get_Rho_1();
		Rhob[NodeArrays->NodeSymmetry[j].Get_index()]=(1- alpha) *Rho[NodeArrays->NodeSymmetry[j].Get_index()];//PtrParameters->Get_Rho_2();
		U[0][NodeArrays->NodeSymmetry[j].Get_index()]=U_[0];
		U[1][NodeArrays->NodeSymmetry[j].Get_index()]=U_[1];
	}
	for (int j=0;j<NodeArrays->NodeGhost.size();j++)
	{
		pos[0]=NodeArrays->NodeGhost[j].get_x();
		pos[1]=NodeArrays->NodeGhost[j].get_y();
		ini.IniDomain(parallel->getRank(),NodeArrays->NodeGhost[j],0, NodeArrays->NodeGhost[j].Get_index(),pos,Rho[NodeArrays->NodeGhost[j].Get_index()],U_,alpha);
//		Rho[NodeArrays->NodeGhost[j].Get_index()]=Convert_Alpha_To_Rho(alpha);
		Rhor[NodeArrays->NodeGhost[j].Get_index()]=alpha*Rho[NodeArrays->NodeGhost[j].Get_index()];//PtrParameters->Get_Rho_1();
		Rhob[NodeArrays->NodeGhost[j].Get_index()]=(1- alpha) *Rho[NodeArrays->NodeGhost[j].Get_index()];//PtrParameters->Get_Rho_2();
		U[0][NodeArrays->NodeGhost[j].Get_index()]=U_[0];
		U[1][NodeArrays->NodeGhost[j].Get_index()]=U_[1];
	}
	for (int j=0;j<NodeArrays->NodeSolid.size();j++)
	{
		pos[0]=NodeArrays->NodeSolid[j].get_x();
		pos[1]=NodeArrays->NodeSolid[j].get_y();
		ini.IniDomain(parallel->getRank(),NodeArrays->NodeSolid[j],0, NodeArrays->NodeSolid[j].Get_index(),pos,Rho[NodeArrays->NodeSolid[j].Get_index()],U_,alpha);
//		Rho[NodeArrays->NodeSolid[j].Get_index()]=Convert_Alpha_To_Rho(alpha);
		Rhor[NodeArrays->NodeSolid[j].Get_index()]=alpha*Rho[NodeArrays->NodeSolid[j].Get_index()];//PtrParameters->Get_Rho_1();
		Rhob[NodeArrays->NodeSolid[j].Get_index()]=(1- alpha) *Rho[NodeArrays->NodeSolid[j].Get_index()];//PtrParameters->Get_Rho_2();
		U[0][NodeArrays->NodeSolid[j].Get_index()]=U_[0];
		U[1][NodeArrays->NodeSolid[j].Get_index()]=U_[1];
	}
	for (int i=0;i<nbvelo;i++)
	{
		for (int j=0;j<nbnode;j++)
		{
			f[0]->f[i][j]=CollideLowOrder::EquiDistriFunct2D(Rhor[j], U[0][j], U[1][j], &Ei[i][0], omega[i]);
			f[1]->f[i][j]=CollideLowOrder::EquiDistriFunct2D(Rhob[j], U[0][j], U[1][j], &Ei[i][0], omega[i]);
		}
	}
/*	for (int j=0;j<NodeArrays->NodeSolid.size();j++)
	{
		Rho[NodeArrays->NodeSolid[j].Get_index()]=NULL;;
		U[0][NodeArrays->NodeSolid[j].Get_index()]=NAN;
		U[1][NodeArrays->NodeSolid[j].Get_index()]=NAN;
	}*/

	Set_BcType();
	if(PtrParameters->Get_Verbous())
		for (int j=0;j<NodeArrays->NodeCorner.size();j++)
		{
			std::cout<<"Processor: "<<parallel->getRank()<<" corner number: "<<j<<" Node index: "<<NodeArrays->NodeCorner[j].Get_index() <<" x: "<<NodeArrays->NodeCorner[j].get_x()<<" y: "<<NodeArrays->NodeCorner[j].get_y()<<" orientation: "<<NodeArrays->NodeCorner[j].Get_BcNormal()<<std::endl;
		}
	/*if(parallel->getRank()==0){
	 int testindexcorner;
	testindexcorner=4;
	std::cout<<"Processor: "<<parallel->getRank()<<" Node index: "<<NodeArrays->NodeCorner[testindexcorner].Get_index()<<std::endl;
	for(unsigned int j=0;j<(unsigned int)nbvelo;j++)
	  {
		std::cout<<" streaming: "<<NodeArrays->NodeCorner[testindexcorner].stream()[j]<<" in direction: "<<j<<std::endl;
	  }
	std::cout<<"x Velocity: "<<U[0][NodeArrays->NodeCorner[testindexcorner].Get_index()]<<std::endl;
	testindexcorner=8;
	std::cout<<"Processor: "<<parallel->getRank()<<" Node index: "<<NodeArrays->NodeCorner[testindexcorner].Get_index()<<std::endl;
	for(unsigned int j=0;j<(unsigned int)nbvelo;j++)
	  {
		std::cout<<" streaming: "<<NodeArrays->NodeCorner[testindexcorner].stream()[j]<<" in direction: "<<j<<std::endl;
	  }
	std::cout<<"x Velocity: "<<U[0][NodeArrays->NodeCorner[testindexcorner].Get_index()]<<std::endl;
	testindexcorner=5;
	std::cout<<"Processor: "<<parallel->getRank()<<" Node index: "<<NodeArrays->NodeCorner[testindexcorner].Get_index()<<std::endl;
	for(unsigned int j=0;j<(unsigned int)nbvelo;j++)
	  {
		std::cout<<" streaming: "<<NodeArrays->NodeCorner[testindexcorner].stream()[j]<<" in direction: "<<j<<std::endl;
	  }
	std::cout<<"x Velocity: "<<U[0][NodeArrays->NodeCorner[testindexcorner].Get_index()]<<std::endl;
	testindexcorner=9;
	std::cout<<"Processor: "<<parallel->getRank()<<" Node index: "<<NodeArrays->NodeCorner[testindexcorner].Get_index()<<std::endl;
	for(unsigned int j=0;j<(unsigned int)nbvelo;j++)
	  {
		std::cout<<" streaming: "<<NodeArrays->NodeCorner[testindexcorner].stream()[j]<<" in direction: "<<j<<std::endl;
	  }
	std::cout<<"x Velocity: "<<U[0][NodeArrays->NodeCorner[testindexcorner].Get_index()]<<std::endl;


	  }*/
	//std::cout<<std::endl;
	delete [] pos;
	delete [] U_;
	Writer->Set_solution(PtrVariablesOutput,PtrParameters->Get_PtrVariablesOutput(),PtrParameters->Get_NbVariablesOutput());
	parallel->barrier();

}

void D2Q9TwoPhases::run(){

	int NbStep=PtrParameters->Get_NbStep();
	int OutPutNStep=PtrParameters->Get_OutPutNSteps();
	int listing=PtrParameters->Get_listing();
	double time_inirun=parallel->getTime();
	double time_run=0;
	int it=0;

	if(parallel->getSize()>1)
		SyncToGhost();
	ApplyBc();
	UpdateMacroVariables();
	Writer->Write_Output(it);

	if(parallel->getSize()>1)
	{
		for (int i=1;i<NbStep+1;i++)
		{

			TwoPhase_Collision();
			SyncToGhost();
			StreamD2Q9();;
			ApplyBc();
			UpdateMacroVariables();
			SyncMacroVarToGhost();
			if(i%OutPutNStep==0)
			{
				Writer->Write_Output(i);
			}
			if(i%listing==0 && parallel->isMainProcessor() )
			{
				time_run=parallel->getTime()-time_inirun;

				std::cout<<"Iteration number: "<<i<< " Running time: ";
				if(time_run>60.0)
					std::cout<<trunc(time_run/60)<<"Minutes ";
				else //less than 1 minute
				std::cout<<trunc(time_run)<<"Seconds ";
				std::cout<< " Time per iteration: "<<time_run/i<<"s"<<std::endl;
			}
		}
	}
	else
	{
		for (int i=1;i<NbStep+1;i++)
		{

			TwoPhase_Collision();
			StreamD2Q9();
			ApplyBc();
			UpdateMacroVariables();
			if(i%OutPutNStep==0)
			{
				Writer->Write_Output(i);
			}
			if(i%listing==0)
			{
				time_run=parallel->getTime()-time_inirun;
				std::cout<<"Iteration number: "<<i<< " Running time: ";
				if(time_run>60.0)
					std::cout<<trunc(time_run/60)<<"Minutes ";
				else //less than 1 minute
				std::cout<<trunc(time_run)<<"Seconds ";
				std::cout<< " Time per iteration: "<<time_run/i<<"s"<<std::endl;			}
		}
	}
	Writer->Write_breakpoint(*PtrParameters);
}

void D2Q9TwoPhases::StreamD2Q9() {


	for (int k=0;k<2;k++)
	{
	for (unsigned int i=1;i<(unsigned int)nbvelo;i++) //No need to stream direction 0
	{
		for (int j=0;j<NodeArrays->NodeInterior.size();j++)
		{
			if (NodeArrays->NodeInterior[j].Get_connect()[i]!=NodeArrays->NodeInterior[j].Get_index())
			{
				ftmp[NodeArrays->NodeInterior[j].Get_connect()[i]]=f[k]->f[i][NodeArrays->NodeInterior[j].Get_index()];
			}
		}

		for (int j=0;j<NodeArrays->NodeCorner.size();j++)
		{
			if (NodeArrays->NodeCorner[j].stream()[i])
			{
				ftmp[NodeArrays->NodeCorner[j].Get_connect()[i]]=f[k]->f[i][NodeArrays->NodeCorner[j].Get_index()];
			}		}
		for (int j=0;j<NodeArrays->NodeGlobalCorner.size();j++)
		{
			if (NodeArrays->NodeGlobalCorner[j].stream()[i])
			{
				ftmp[NodeArrays->NodeGlobalCorner[j].Get_connect()[i]]=f[k]->f[i][NodeArrays->NodeGlobalCorner[j].Get_index()];
			}

		}
		for (int j=0;j<NodeArrays->NodeVelocity.size();j++)
		{
			if (NodeArrays->NodeVelocity[j].stream()[i])
					{
						ftmp[NodeArrays->NodeVelocity[j].Get_connect()[i]]=f[k]->f[i][NodeArrays->NodeVelocity[j].Get_index()];
					}
		}

		for (int j=0;j<NodeArrays->NodePressure.size();j++)
		{
			if (NodeArrays->NodePressure[j].stream()[i])
					{
						ftmp[NodeArrays->NodePressure[j].Get_connect()[i]]=f[k]->f[i][NodeArrays->NodePressure[j].Get_index()];
					}
		}
		for (int j=0;j<NodeArrays->NodeWall.size();j++)
		{
			if (NodeArrays->NodeWall[j].stream()[i])
			{
				ftmp[NodeArrays->NodeWall[j].Get_connect()[i]]=f[k]->f[i][NodeArrays->NodeWall[j].Get_index()];
			}
		}
		for (int j=0;j<NodeArrays->NodeSpecialWall.size();j++)
		{
			if (NodeArrays->NodeSpecialWall[j].stream()[i])
			{
				ftmp[NodeArrays->NodeSpecialWall[j].Get_connect()[i]]=f[k]->f[i][NodeArrays->NodeSpecialWall[j].Get_index()];
			}
		}
		for (int j=0;j<NodeArrays->NodeSymmetry.size();j++)
		{
			if (NodeArrays->NodeSymmetry[j].stream()[i])
			{
				ftmp[NodeArrays->NodeSymmetry[j].Get_connect()[i]]=f[k]->f[i][NodeArrays->NodeSymmetry[j].Get_index()];
			}
		}
		for (int j=0;j<NodeArrays->NodeGhost.size();j++)
		{
			if (NodeArrays->NodeGhost[j].stream()[i])
			{
				ftmp[NodeArrays->NodeGhost[j].Get_connect()[i]]=f[k]->f[i][NodeArrays->NodeGhost[j].Get_index()];
			}
		}/*for (int j=0;j<nbnode;j++)
		{
			D2Q9TwoPhases::SelectStream(j,i);
		}*/


		D2Q9TwoPhases::TmptoDistri(i,k);
	}
	}
}
void D2Q9TwoPhases::CollideD2Q9(){
	double wtmp=0;
	double *fi_tmp,*U_tmp, *V_tmp, *F_tmp;
	double Rhoi_tmp;
	fi_tmp=new double [9];
	F_tmp=new double [2];
	int idx_tmp;
	for (int i=0;i<nbvelo;i++)
	{
		for (int j=0;j<NodeArrays->NodeInterior.size();j++)
		{
			Collide_SinglePhase(f[0]->f[i][NodeArrays->NodeInterior[j].Get_index()],Rhor[NodeArrays->NodeInterior[j].Get_index()], U[0][NodeArrays->NodeInterior[j].Get_index()], U[1][NodeArrays->NodeInterior[j].Get_index()], &Ei[i][0], omega[i]);
			Collide_SinglePhase(f[1]->f[i][NodeArrays->NodeInterior[j].Get_index()],Rhob[NodeArrays->NodeInterior[j].Get_index()], U[0][NodeArrays->NodeInterior[j].Get_index()], U[1][NodeArrays->NodeInterior[j].Get_index()], &Ei[i][0], omega[i]);

			//f[0]->f[i][NodeArrays->NodeInterior[j].Get_index()]= f[0]->f[i][NodeArrays->NodeInterior[j].Get_index()]-InvTau*(f[0]->f[i][NodeArrays->NodeInterior[j].Get_index()]-CollideLowOrder::EquiDistriFunct2D(Rho[NodeArrays->NodeInterior[j].Get_index()], U[0][NodeArrays->NodeInterior[j].Get_index()], U[1][NodeArrays->NodeInterior[j].Get_index()], &Ei[i][0], omega[i]))
					;//+LocalForce(i,Rho[NodeArrays->NodeCorner[j].Get_index()],U[0][NodeArrays->NodeInterior[j].Get_index()],U[1][NodeArrays->NodeInterior[j].Get_index()],wtmp);
		}
		for (int j=0;j<NodeArrays->NodeCorner.size();j++)
		{
			Collide_SinglePhase(f[0]->f[i][NodeArrays->NodeCorner[j].Get_index()],Rhor[NodeArrays->NodeCorner[j].Get_index()], U[0][NodeArrays->NodeCorner[j].Get_index()], U[1][NodeArrays->NodeCorner[j].Get_index()], &Ei[i][0], omega[i]);
			Collide_SinglePhase(f[1]->f[i][NodeArrays->NodeCorner[j].Get_index()],Rhob[NodeArrays->NodeCorner[j].Get_index()], U[0][NodeArrays->NodeCorner[j].Get_index()], U[1][NodeArrays->NodeCorner[j].Get_index()], &Ei[i][0], omega[i]);

			//f[0]->f[i][NodeArrays->NodeCorner[j].Get_index()]= f[0]->f[i][NodeArrays->NodeCorner[j].Get_index()]-InvTau*(f[0]->f[i][NodeArrays->NodeCorner[j].Get_index()]-CollideLowOrder::EquiDistriFunct2D(Rho[NodeArrays->NodeCorner[j].Get_index()], U[0][NodeArrays->NodeCorner[j].Get_index()], U[1][NodeArrays->NodeCorner[j].Get_index()], &Ei[i][0], omega[i]))
					;//+LocalForce(i,Rho[NodeArrays->NodeCorner[j].Get_index()],U[0][NodeArrays->NodeCorner[j].Get_index()],U[1][NodeArrays->NodeCorner[j].Get_index()],wtmp);
		}
		for (int j=0;j<NodeArrays->NodeGlobalCorner.size();j++)
		{
			Collide_SinglePhase(f[0]->f[i][NodeArrays->NodeGlobalCorner[j].Get_index()],Rhor[NodeArrays->NodeGlobalCorner[j].Get_index()], U[0][NodeArrays->NodeGlobalCorner[j].Get_index()], U[1][NodeArrays->NodeGlobalCorner[j].Get_index()], &Ei[i][0], omega[i]);
			Collide_SinglePhase(f[1]->f[i][NodeArrays->NodeGlobalCorner[j].Get_index()],Rhob[NodeArrays->NodeGlobalCorner[j].Get_index()], U[0][NodeArrays->NodeGlobalCorner[j].Get_index()], U[1][NodeArrays->NodeGlobalCorner[j].Get_index()], &Ei[i][0], omega[i]);

			//f[0]->f[i][NodeArrays->NodeGlobalCorner[j].Get_index()]= f[0]->f[i][NodeArrays->NodeGlobalCorner[j].Get_index()]-InvTau*(f[0]->f[i][NodeArrays->NodeGlobalCorner[j].Get_index()]-CollideLowOrder::EquiDistriFunct2D(Rho[NodeArrays->NodeGlobalCorner[j].Get_index()], U[0][NodeArrays->NodeGlobalCorner[j].Get_index()], U[1][NodeArrays->NodeGlobalCorner[j].Get_index()], &Ei[i][0], omega[i]))
					;//+LocalForce(i,Rho[NodeArrays->NodeGlobalCorner[j].Get_index()],U[0][NodeArrays->NodeGlobalCorner[j].Get_index()],U[1][NodeArrays->NodeGlobalCorner[j].Get_index()],wtmp);
		}
		for (int j=0;j<NodeArrays->NodeVelocity.size();j++)
		{
			Collide_SinglePhase(f[0]->f[i][NodeArrays->NodeVelocity[j].Get_index()],Rhor[NodeArrays->NodeVelocity[j].Get_index()], U[0][NodeArrays->NodeVelocity[j].Get_index()], U[1][NodeArrays->NodeVelocity[j].Get_index()], &Ei[i][0], omega[i]);
			Collide_SinglePhase(f[1]->f[i][NodeArrays->NodeVelocity[j].Get_index()],Rhob[NodeArrays->NodeVelocity[j].Get_index()], U[0][NodeArrays->NodeVelocity[j].Get_index()], U[1][NodeArrays->NodeVelocity[j].Get_index()], &Ei[i][0], omega[i]);

			//f[0]->f[i][NodeArrays->NodeVelocity[j].Get_index()]= f[0]->f[i][NodeArrays->NodeVelocity[j].Get_index()]-InvTau*(f[0]->f[i][NodeArrays->NodeVelocity[j].Get_index()]-CollideLowOrder::EquiDistriFunct2D(Rho[NodeArrays->NodeVelocity[j].Get_index()], U[0][NodeArrays->NodeVelocity[j].Get_index()], U[1][NodeArrays->NodeVelocity[j].Get_index()], &Ei[i][0], omega[i]))
					;//+LocalForce(i,Rho[NodeArrays->NodeVelocity[j].Get_index()],U[0][NodeArrays->NodeVelocity[j].Get_index()],U[1][NodeArrays->NodeVelocity[j].Get_index()],wtmp);
		}

		for (int j=0;j<NodeArrays->NodePressure.size();j++)
		{
			Collide_SinglePhase(f[0]->f[i][NodeArrays->NodePressure[j].Get_index()],Rhor[NodeArrays->NodePressure[j].Get_index()], U[0][NodeArrays->NodePressure[j].Get_index()], U[1][NodeArrays->NodePressure[j].Get_index()], &Ei[i][0], omega[i]);
			Collide_SinglePhase(f[1]->f[i][NodeArrays->NodePressure[j].Get_index()],Rhob[NodeArrays->NodePressure[j].Get_index()], U[0][NodeArrays->NodePressure[j].Get_index()], U[1][NodeArrays->NodePressure[j].Get_index()], &Ei[i][0], omega[i]);

			//f[0]->f[i][NodeArrays->NodePressure[j].Get_index()]= f[0]->f[i][NodeArrays->NodePressure[j].Get_index()]-InvTau*(f[0]->f[i][NodeArrays->NodePressure[j].Get_index()]-CollideLowOrder::EquiDistriFunct2D(Rho[NodeArrays->NodePressure[j].Get_index()], U[0][NodeArrays->NodePressure[j].Get_index()], U[1][NodeArrays->NodePressure[j].Get_index()], &Ei[i][0], omega[i]))
					;//+LocalForce(i,Rho[NodeArrays->NodePressure[j].Get_index()],U[0][NodeArrays->NodePressure[j].Get_index()],U[1][NodeArrays->NodePressure[j].Get_index()],wtmp);
		}
		for (int j=0;j<NodeArrays->NodeWall.size();j++)
		{
			Collide_SinglePhase(f[0]->f[i][NodeArrays->NodeWall[j].Get_index()],Rhor[NodeArrays->NodeWall[j].Get_index()], U[0][NodeArrays->NodeWall[j].Get_index()], U[1][NodeArrays->NodeWall[j].Get_index()], &Ei[i][0], omega[i]);
			Collide_SinglePhase(f[1]->f[i][NodeArrays->NodeWall[j].Get_index()],Rhob[NodeArrays->NodeWall[j].Get_index()], U[0][NodeArrays->NodeWall[j].Get_index()], U[1][NodeArrays->NodeWall[j].Get_index()], &Ei[i][0], omega[i]);

			//f[0]->f[i][NodeArrays->NodeWall[j].Get_index()] = f[0]->f[i][NodeArrays->NodeWall[j].Get_index()]-InvTau*(f[0]->f[i][NodeArrays->NodeWall[j].Get_index()]-CollideLowOrder::EquiDistriFunct2D(Rho[NodeArrays->NodeWall[j].Get_index()], U[0][NodeArrays->NodeWall[j].Get_index()], U[1][NodeArrays->NodeWall[j].Get_index()], &Ei[i][0], omega[i]))
					;//+LocalForce(i,Rho[NodeArrays->NodeWall[j].Get_index()],U[0][NodeArrays->NodeWall[j].Get_index()],U[1][NodeArrays->NodeWall[j].Get_index()],wtmp);
		}
		for (int j=0;j<NodeArrays->NodeSymmetry.size();j++)
		{
			Collide_SinglePhase(f[0]->f[i][NodeArrays->NodeSymmetry[j].Get_index()],Rhob[NodeArrays->NodeSymmetry[j].Get_index()], U[0][NodeArrays->NodeSymmetry[j].Get_index()], U[1][NodeArrays->NodeSymmetry[j].Get_index()], &Ei[i][0], omega[i]);
			Collide_SinglePhase(f[1]->f[i][NodeArrays->NodeSymmetry[j].Get_index()],Rhob[NodeArrays->NodeSymmetry[j].Get_index()], U[0][NodeArrays->NodeSymmetry[j].Get_index()], U[1][NodeArrays->NodeSymmetry[j].Get_index()], &Ei[i][0], omega[i]);

			//f[0]->f[i][NodeArrays->NodeSymmetry[j].Get_index()] = f[0]->f[i][NodeArrays->NodeSymmetry[j].Get_index()]-InvTau*(f[0]->f[i][NodeArrays->NodeSymmetry[j].Get_index()]-CollideLowOrder::EquiDistriFunct2D(Rho[NodeArrays->NodeSymmetry[j].Get_index()], U[0][NodeArrays->NodeSymmetry[j].Get_index()], U[1][NodeArrays->NodeSymmetry[j].Get_index()], &Ei[i][0], omega[i]))
					;//+LocalForce(i,Rho[NodeArrays->NodeSymmetry[j].Get_index()],U[0][NodeArrays->NodeSymmetry[j].Get_index()],U[1][NodeArrays->NodeSymmetry[j].Get_index()],wtmp);
		}
	//	for (int j=0;j<nbnodes_real;j++)
	//		f[0]->f[i][j]=f[0]->f[i][j]-InvTau*(f[0]->f[i][j]-CollideLowOrder::EquiDistriFunct2D(Rho[j], U[0][j], U[1][j], &Ei[i][0], omega[i]));
	}
	delete [] fi_tmp;
}
void D2Q9TwoPhases::UpdateMacroVariables(){
	//for (int i=0;i<nbvelo;i++)
	//{
		for (int j=0;j<NodeArrays->NodeInterior.size();j++)
		{
			MacroVariables(NodeArrays->NodeInterior[j].Get_index());
		}
		for (int j=0;j<NodeArrays->NodeCorner.size();j++)
		{
			MacroVariables(NodeArrays->NodeCorner[j].Get_index());
		}
		for (int j=0;j<NodeArrays->NodeGlobalCorner.size();j++)
		{
			MacroVariables(NodeArrays->NodeGlobalCorner[j].Get_index());
		}
		for (int j=0;j<NodeArrays->NodeVelocity.size();j++)
		{
			MacroVariables(NodeArrays->NodeVelocity[j].Get_index());
		}
		for (int j=0;j<NodeArrays->NodePressure.size();j++)
		{
			MacroVariables(NodeArrays->NodePressure[j].Get_index());
		}
		for (int j=0;j<NodeArrays->NodeWall.size();j++)
		{
			MacroVariables(NodeArrays->NodeWall[j].Get_index());
		}
		for (int j=0;j<NodeArrays->NodeSymmetry.size();j++)
		{
			MacroVariables(NodeArrays->NodeSymmetry[j].Get_index());
		}
	//	for (int j=0;j<nbnodes_real;j++)
	//		f[0]->f[i][j]=f[0]->f[i][j]-InvTau*(f[0]->f[i][j]-CollideLowOrder::EquiDistriFunct2D(Rho[j], U[0][j], U[1][j], &Ei[i][0], omega[i]));
	//}
}
void D2Q9TwoPhases::MacroVariables(int& idx){

		U[0][idx]=0;
		U[1][idx]=0;
		Rhor[idx]=0,Rhob[idx]=0;
		//Rho[idx]=0;
		for (int k=0; k<nbvelo;k++)
		{
			Rhor[idx]+=f[0]->f[k][idx];
			Rhob[idx]+=f[1]->f[k][idx];
			for (int j=0;j<2;j++)
			{
				U[j][idx]+=(f[0]->f[k][idx]+f[1]->f[k][idx])*Ei[k][j];
			}
		}
		Rho[idx]=Rhor[idx]+Rhob[idx];
		U[0][idx]=U[0][idx]/Rho[idx];
		U[1][idx]=U[1][idx]/Rho[idx];

}
void D2Q9TwoPhases::TmptoDistri(unsigned int& direction, int& IdDistri){
	tmp= f[IdDistri]->f[direction];
	f[IdDistri]->f[direction]=ftmp;
	ftmp=tmp;
}
void D2Q9TwoPhases::TwoPhase_Collision()
{
	double Ak=0.65;
	double F[2];
	double F_Norm=0;
	int idx_tmp;
	double wtmp=0;
//	double *U_tmp, *V_tmp, *F_tmp;
	double  InvTau_,fi_tmp;
	InvTau_=InvTau; //Tmp

	for (int j=0;j<NodeArrays->NodeInterior.size();j++)
	{
	// Common variables
		idx_tmp=NodeArrays->NodeInterior[j].Get_index();
	// Calculate gradients
		Colour_gradient(idx_tmp,&F[0]);
	//Calculate Norms
		F_Norm=std::sqrt(F[0]*F[0]+F[1]*F[1]);
	//Model
		for (int i=0;i<9;i++)
		{
			fi_tmp=f[0]->f[i][idx_tmp]+f[1]->f[i][idx_tmp];
			Collide_ColorFluid(i,fi_tmp,Rho[idx_tmp],&F[0],F_Norm, InvTau_, U[0][idx_tmp], U[1][idx_tmp]);
			Recoloring(fi_tmp, f[0]->f[i][idx_tmp], f[1]->f[i][idx_tmp], Rho[idx_tmp], Rhor[idx_tmp], Rhob[idx_tmp]);
		}
	}
	for (int j=0;j<NodeArrays->NodeCorner.size();j++)
	{
	// Common variables
		idx_tmp=NodeArrays->NodeCorner[j].Get_index();
	// Calculate gradients
		Colour_gradient(idx_tmp,&F[0]);
	//Calculate Norms
		F_Norm=std::sqrt(F[0]*F[0]+F[1]*F[1]);
	//Model
		for (int i=0;i<9;i++)
		{
			fi_tmp=f[0]->f[i][idx_tmp]+f[1]->f[i][idx_tmp];
			Collide_ColorFluid(i,fi_tmp,Rho[idx_tmp],&F[0],F_Norm, InvTau_, U[0][idx_tmp], U[1][idx_tmp]);
			Recoloring(fi_tmp, f[0]->f[i][idx_tmp], f[1]->f[i][idx_tmp], Rho[idx_tmp], Rhor[idx_tmp], Rhob[idx_tmp]);
		}
	}
	for (int j=0;j<NodeArrays->NodeGlobalCorner.size();j++)
	{
	// Common variables
		idx_tmp=NodeArrays->NodeGlobalCorner[j].Get_index();
	// Calculate gradients
		Colour_gradient(idx_tmp,&F[0]);
	//Calculate Norms
		F_Norm=std::sqrt(F[0]*F[0]+F[1]*F[1]);
	//Model
		for (int i=0;i<9;i++)
		{
			fi_tmp=f[0]->f[i][idx_tmp]+f[1]->f[i][idx_tmp];
			Collide_ColorFluid(i,fi_tmp,Rho[idx_tmp],&F[0],F_Norm, InvTau_, U[0][idx_tmp], U[1][idx_tmp]);
			Recoloring(fi_tmp, f[0]->f[i][idx_tmp], f[1]->f[i][idx_tmp], Rho[idx_tmp], Rhor[idx_tmp], Rhob[idx_tmp]);
		}
	}
	for (int j=0;j<NodeArrays->NodeVelocity.size();j++)
	{
	// Common variables
		idx_tmp=NodeArrays->NodeVelocity[j].Get_index();
	// Calculate gradients
		Colour_gradient(idx_tmp,&F[0]);
	//Calculate Norms
		F_Norm=std::sqrt(F[0]*F[0]+F[1]*F[1]);
	//Model
		for (int i=0;i<9;i++)
		{
			fi_tmp=f[0]->f[i][idx_tmp]+f[1]->f[i][idx_tmp];
			Collide_ColorFluid(i,fi_tmp,Rho[idx_tmp],&F[0],F_Norm, InvTau_, U[0][idx_tmp], U[1][idx_tmp]);
			Recoloring(fi_tmp, f[0]->f[i][idx_tmp], f[1]->f[i][idx_tmp], Rho[idx_tmp], Rhor[idx_tmp], Rhob[idx_tmp]);
		}
	}
	for (int j=0;j<NodeArrays->NodePressure.size();j++)
	{
	// Common variables
		idx_tmp=NodeArrays->NodePressure[j].Get_index();
	// Calculate gradients
		Colour_gradient(idx_tmp,&F[0]);
	//Calculate Norms
		F_Norm=std::sqrt(F[0]*F[0]+F[1]*F[1]);
	//Model
		for (int i=0;i<9;i++)
		{
			fi_tmp=f[0]->f[i][idx_tmp]+f[1]->f[i][idx_tmp];
			Collide_ColorFluid(i,fi_tmp,Rho[idx_tmp],&F[0],F_Norm, InvTau_, U[0][idx_tmp], U[1][idx_tmp]);
			Recoloring(fi_tmp, f[0]->f[i][idx_tmp], f[1]->f[i][idx_tmp], Rho[idx_tmp], Rhor[idx_tmp], Rhob[idx_tmp]);
		}
	}
	for (int j=0;j<NodeArrays->NodeWall.size();j++)
	{
	// Common variables
		idx_tmp=NodeArrays->NodeWall[j].Get_index();
	// Calculate gradients
		Colour_gradient(idx_tmp,&F[0]);
	//Calculate Norms
		F_Norm=std::sqrt(F[0]*F[0]+F[1]*F[1]);
	//Model
		for (int i=0;i<9;i++)
		{
			fi_tmp=f[0]->f[i][idx_tmp]+f[1]->f[i][idx_tmp];
			Collide_ColorFluid(i,fi_tmp,Rho[idx_tmp],&F[0],F_Norm, InvTau_, U[0][idx_tmp], U[1][idx_tmp]);
			Recoloring(fi_tmp, f[0]->f[i][idx_tmp], f[1]->f[i][idx_tmp], Rho[idx_tmp], Rhor[idx_tmp], Rhob[idx_tmp]);
		}
	}
	for (int j=0;j<NodeArrays->NodeSymmetry.size();j++)
	{
	// Common variables
		idx_tmp=NodeArrays->NodeSymmetry[j].Get_index();
	// Calculate gradients
		Colour_gradient(idx_tmp,&F[0]);
	//Calculate Norms
		F_Norm=std::sqrt(F[0]*F[0]+F[1]*F[1]);
	//Model
		for (int i=0;i<9;i++)
		{
			fi_tmp=f[0]->f[i][idx_tmp]+f[1]->f[i][idx_tmp];
			Collide_ColorFluid(i,fi_tmp,Rho[idx_tmp],&F[0],F_Norm, InvTau_, U[0][idx_tmp], U[1][idx_tmp]);
			Recoloring(fi_tmp, f[0]->f[i][idx_tmp], f[1]->f[i][idx_tmp], Rho[idx_tmp], Rhor[idx_tmp], Rhob[idx_tmp]);
		}
	}
}
void D2Q9TwoPhases::Colour_gradient(int & nodenumber, double* F){
	for (int k=0; k<nbvelo;k++)
	{
		for (int j=0;j<2;j++)
		{
			F[j]=(Rhor[nodenumber]-Rhob[nodenumber])*Ei[k][j];
		}
	}
}
double D2Q9TwoPhases::TwoPhase_Collision_operator(int & nodenumber, int & i, double & Ak, double* F, double & F_Norm){
 return Ak*0.5*F_Norm*(((F[0]*Ei[i][0]+F[1]*Ei[i][1])/F_Norm)*((F[0]*Ei[i][0]+F[1]*Ei[i][1])/F_Norm)-3/4);
}

void D2Q9TwoPhases::Set_BcType(){
	for (int j=0;j<NodeArrays->NodeCorner.size();j++)
	{
		Set_CornerType(NodeArrays->NodeCorner[j]);
	}
	for (int j=0;j<NodeArrays->NodeGlobalCorner.size();j++)
	{
		Set_CornerType(NodeArrays->NodeGlobalCorner[j]);
	}
	for (int j=0;j<NodeArrays->NodeVelocity.size();j++)
	{
		Set_VelocityType(NodeArrays->NodeVelocity[j]);
	}

	for (int j=0;j<NodeArrays->NodePressure.size();j++)
	{
		Set_PressureType(NodeArrays->NodePressure[j]);
	}
	for (int j=0;j<NodeArrays->NodeWall.size();j++)
	{
		Set_WallType(NodeArrays->NodeWall[j]);
	}
	for (int j=0;j<NodeArrays->NodeSpecialWall.size();j++)
	{
		Set_WallType(NodeArrays->NodeSpecialWall[j]);
	}
	for (int j=0;j<NodeArrays->NodeSymmetry.size();j++)
	{
		Set_SymmetryType(NodeArrays->NodeSymmetry[j]);
	}
	for (int j=0;j<NodeArrays->NodeGhost.size();j++)
	{
		Set_GhostType(NodeArrays->NodeGhost[j]);
	}
}
void D2Q9TwoPhases::Set_GhostType(NodeGhost2D& NodeIn){
	bool GhostStreaming[9];
	StreamingOrientation(NodeIn,GhostStreaming);
	NodeIn.Set_stream(GhostStreaming,nbvelo);

}
void D2Q9TwoPhases::Set_WallType(NodeWall2D& NodeIn){
	bool WallStreaming[9];
	StreamingOrientation(NodeIn,WallStreaming);
	NodeIn.Set_stream(WallStreaming,nbvelo);
}
void D2Q9TwoPhases::Set_SymmetryType(NodeSymmetry2D& NodeIn){
	bool SymmetryStreaming[9];
	StreamingOrientation(NodeIn,SymmetryStreaming);
	NodeIn.Set_stream(SymmetryStreaming,nbvelo);

}
void D2Q9TwoPhases::Set_CornerType(NodeCorner2D& NodeIn){
	bool CornerStreaming[9];
	StreamingOrientation(NodeIn,CornerStreaming);
	NodeIn.Set_stream(CornerStreaming,nbvelo);
	NodeIn.Set_UDef(U[0][NodeIn.Get_index()],U[1][NodeIn.Get_index()]);
	NodeIn.Set_RhoDef(Rho[NodeIn.Get_index()]);
}
void D2Q9TwoPhases::Set_VelocityType(NodeVelocity2D& NodeIn){
	bool VelocityStreaming[9];
	StreamingOrientation(NodeIn,VelocityStreaming);
	NodeIn.Set_stream(VelocityStreaming,nbvelo);
	NodeIn.Set_UDef(U[0][NodeIn.Get_index()],U[1][NodeIn.Get_index()]);
}
void D2Q9TwoPhases::Set_PressureType(NodePressure2D& NodeIn){
	bool PressureStreaming[9];
	StreamingOrientation(NodeIn,PressureStreaming);
	NodeIn.Set_stream(PressureStreaming,nbvelo);
	NodeIn.Set_RhoDef(Rho[NodeIn.Get_index()]);
}


///Select and apply boundary conditions
void D2Q9TwoPhases::ApplyBc(){
	for (int j=0;j<NodeArrays->NodeVelocity.size();j++)
	{
		ApplyHeZou_U(NodeArrays->NodeVelocity[j]);
	}

	for (int j=0;j<NodeArrays->NodePressure.size();j++)
	{
		ApplyHeZou_P(NodeArrays->NodePressure[j]);
	}
	switch(PtrParameters->Get_WallType())
	{
		case BounceBack:
		for (int j=0;j<NodeArrays->NodeWall.size();j++)
		{
			ApplyBounceBack(NodeArrays->NodeWall[j]);
		}
		for (int j=0;j<NodeArrays->NodeSpecialWall.size();j++)
		{
			ApplyBounceBackSymmetry(NodeArrays->NodeSpecialWall[j]);
		}
		for (int j=0;j<NodeArrays->NodeCorner.size();j++)
		{
	//		ApplyCorner(NodeArrays->NodeCorner[j]);
			ApplyBounceBack(NodeArrays->NodeCorner[j]);
		}
		break;
		case Diffuse:
		for (int j=0;j<NodeArrays->NodeWall.size();j++)
		{
			ApplyDiffuseWall(NodeArrays->NodeWall[j]);
		}
		for (int j=0;j<NodeArrays->NodeSpecialWall.size();j++)
		{
			ApplyDiffuseWallSymmetry(NodeArrays->NodeSpecialWall[j]);
		}
		for (int j=0;j<NodeArrays->NodeCorner.size();j++)
		{
			ApplyDiffuseWall(NodeArrays->NodeCorner[j]);
		}
		break;
	}
	for (int j=0;j<NodeArrays->NodeGlobalCorner.size();j++)
	{
//		ApplyCorner(NodeArrays->NodeCorner[j]);
		ApplyGlobalCorner(NodeArrays->NodeGlobalCorner[j]);
	}
	for (int j=0;j<NodeArrays->NodeSymmetry.size();j++)
	{
		if(NodeArrays->NodeSymmetry[j].Get_SymmetryType()==OnNode)
			ApplySymmetryOnNode(NodeArrays->NodeSymmetry[j]);
		if(NodeArrays->NodeSymmetry[j].Get_SymmetryType()==SymmetryPressureOnNode)
			ApplySymmetryPressureOnNode(NodeArrays->NodeSymmetry[j]);
	}
   parallel->barrier();

}
/// Impose velocity by HeZou
void D2Q9TwoPhases::ApplyHeZou_U(NodeVelocity2D& NodeIn){
	double* fi=new double [9];
	for(int j=0;j<2;j++)
	{
		for (int i=0;i<9;i++)
		fi[i]=f[j]->f[i][NodeIn.Get_index()];

		BC_HeZou_U(NodeIn.Get_BcNormal(),fi, NodeIn.Get_UDef()[0],NodeIn.Get_UDef()[1]);

		for (int i=0;i<9;i++)
		f[j]->f[i][NodeIn.Get_index()]=fi[i];
	}
	delete [] fi;
}
void D2Q9TwoPhases::ApplyHeZou_U(NodeCorner2D& NodeIn, int normal){
	double* fi=new double [9];
	for(int j=0;j<2;j++)
	{
		for (int i=0;i<9;i++)
		fi[i]=f[j]->f[i][NodeIn.Get_index()];

		BC_HeZou_U(normal,fi, NodeIn.Get_UDef()[0],NodeIn.Get_UDef()[1]);

		for (int i=0;i<9;i++)
		f[j]->f[i][NodeIn.Get_index()]=fi[i];
	}
	delete [] fi;
}
/// Impose pressure by HeZou
void D2Q9TwoPhases::ApplyHeZou_P(NodePressure2D& NodeIn){
	double* fi=new double [9];

	if(NodeIn.Get_AlphaDef()==0)
	{
		for (int i=0;i<9;i++)
					f[0]->f[i][NodeIn.Get_index()]=0;
	}
	else
	{
		for (int i=0;i<9;i++)
			fi[i]=f[0]->f[i][NodeIn.Get_index()];
		BC_HeZou_P(NodeIn.Get_BcNormal(),fi,NodeIn.Get_AlphaDef()*NodeIn.Get_RhoDef(), U[0][NodeIn.Get_index()],U[1][NodeIn.Get_index()]);
		for (int i=0;i<9;i++)
			f[0]->f[i][NodeIn.Get_index()]=fi[i];
	}
	if(NodeIn.Get_AlphaDef()==1)
	{
		for (int i=0;i<9;i++)
					f[1]->f[i][NodeIn.Get_index()]=0;
	}
	else
	{
		for (int i=0;i<9;i++)
			fi[i]=f[1]->f[i][NodeIn.Get_index()];
		BC_HeZou_P(NodeIn.Get_BcNormal(),fi,(1-NodeIn.Get_AlphaDef())*NodeIn.Get_RhoDef(), U[0][NodeIn.Get_index()],U[1][NodeIn.Get_index()]);
		for (int i=0;i<9;i++)
			f[1]->f[i][NodeIn.Get_index()]=fi[i];
	}

	delete [] fi;
}
/// Impose pressure by HeZou
void D2Q9TwoPhases::ApplyHeZou_P(NodeCorner2D& NodeIn, int normal){
	double* fi=new double [9];
	if(NodeIn.Get_AlphaDef()==0)
	{
		for (int i=0;i<9;i++)
					f[0]->f[i][NodeIn.Get_index()]=0;
	}
	else
	{
		for (int i=0;i<9;i++)
			fi[i]=f[0]->f[i][NodeIn.Get_index()];
		BC_HeZou_P(normal,fi,NodeIn.Get_AlphaDef()*NodeIn.Get_RhoDef(), U[0][NodeIn.Get_index()],U[1][NodeIn.Get_index()]);
		for (int i=0;i<9;i++)
			f[0]->f[i][NodeIn.Get_index()]=fi[i];
	}
	if(NodeIn.Get_AlphaDef()==1)
	{
		for (int i=0;i<9;i++)
					f[1]->f[i][NodeIn.Get_index()]=0;
	}
	else
	{
		for (int i=0;i<9;i++)
			fi[i]=f[1]->f[i][NodeIn.Get_index()];
		BC_HeZou_P(normal,fi,(1-NodeIn.Get_AlphaDef())*NodeIn.Get_RhoDef(), U[0][NodeIn.Get_index()],U[1][NodeIn.Get_index()]);
		for (int i=0;i<9;i++)
			f[1]->f[i][NodeIn.Get_index()]=fi[i];
	}
	delete [] fi;
}
//Symmetry treatment on node
void D2Q9TwoPhases::ApplySymmetryOnNode(NodeSymmetry2D& NodeIn){
	switch(NodeIn.Get_BcNormal())
			{
			case 2:
				f[0]->f[2][NodeIn.Get_index()]=f[0]->f[4][NodeIn.Get_index()];
				f[0]->f[5][NodeIn.Get_index()]=f[0]->f[8][NodeIn.Get_index()];
				f[0]->f[6][NodeIn.Get_index()]=f[0]->f[7][NodeIn.Get_index()];

				f[1]->f[2][NodeIn.Get_index()]=f[1]->f[4][NodeIn.Get_index()];
				f[1]->f[5][NodeIn.Get_index()]=f[1]->f[8][NodeIn.Get_index()];
				f[1]->f[6][NodeIn.Get_index()]=f[1]->f[7][NodeIn.Get_index()];
				break;
			case 4:
				f[0]->f[4][NodeIn.Get_index()]=f[0]->f[2][NodeIn.Get_index()];
				f[0]->f[7][NodeIn.Get_index()]=f[0]->f[6][NodeIn.Get_index()];
				f[0]->f[8][NodeIn.Get_index()]=f[0]->f[5][NodeIn.Get_index()];

				f[1]->f[4][NodeIn.Get_index()]=f[1]->f[2][NodeIn.Get_index()];
				f[1]->f[7][NodeIn.Get_index()]=f[1]->f[6][NodeIn.Get_index()];
				f[1]->f[8][NodeIn.Get_index()]=f[1]->f[5][NodeIn.Get_index()];
				break;
			case 1:
				f[0]->f[1][NodeIn.Get_index()]=f[0]->f[3][NodeIn.Get_index()];
				f[0]->f[5][NodeIn.Get_index()]=f[0]->f[6][NodeIn.Get_index()];
				f[0]->f[8][NodeIn.Get_index()]=f[0]->f[7][NodeIn.Get_index()];

				f[1]->f[1][NodeIn.Get_index()]=f[1]->f[3][NodeIn.Get_index()];
				f[1]->f[5][NodeIn.Get_index()]=f[1]->f[6][NodeIn.Get_index()];
				f[1]->f[8][NodeIn.Get_index()]=f[1]->f[7][NodeIn.Get_index()];
				break;
			case 3:
				f[0]->f[3][NodeIn.Get_index()]=f[0]->f[1][NodeIn.Get_index()];
				f[0]->f[6][NodeIn.Get_index()]=f[0]->f[5][NodeIn.Get_index()];
				f[0]->f[7][NodeIn.Get_index()]=f[0]->f[8][NodeIn.Get_index()];

				f[1]->f[3][NodeIn.Get_index()]=f[1]->f[1][NodeIn.Get_index()];
				f[1]->f[6][NodeIn.Get_index()]=f[1]->f[5][NodeIn.Get_index()];
				f[1]->f[7][NodeIn.Get_index()]=f[1]->f[8][NodeIn.Get_index()];
				break;
			default:
				std::cerr<<"Direction symmetry not found"<<std::endl;
				break;
			}
}

void D2Q9TwoPhases::ApplySymmetryPressureOnNode(NodeSymmetry2D& NodeIn){
	double* fi=new double [9];
	double rhosym=0;
	double rhoadd=0;
	switch(NodeIn.Get_BcNormal())
			{
			case 2:
				for(int j=0;j<2;j++)
				{
					//Applying symmetry
					fi[0]=f[j]->f[0][NodeIn.Get_index()];
					fi[1]=f[j]->f[1][NodeIn.Get_index()];
					fi[2]=f[j]->f[4][NodeIn.Get_index()];
					fi[3]=f[j]->f[3][NodeIn.Get_index()];
					fi[4]=f[j]->f[4][NodeIn.Get_index()];
					fi[5]=f[j]->f[8][NodeIn.Get_index()];
					fi[6]=f[j]->f[7][NodeIn.Get_index()];
					fi[7]=f[j]->f[7][NodeIn.Get_index()];
					fi[8]=f[j]->f[8][NodeIn.Get_index()];
					//Applying He Zou Pressure
					BC_HeZou_P(NodeIn.Get_BcNormal(),fi,NodeIn.Get_RhoDef(), U[0][NodeIn.Get_index()],U[1][NodeIn.Get_index()]);
					for (int i=0;i<9;i++)
						f[j]->f[i][NodeIn.Get_index()]=fi[i];
				}
				break;
			case 4:
				for(int j=0;j<2;j++)
				{
					//Applying symmetry
					fi[0]=f[j]->f[0][NodeIn.Get_index()];
					fi[1]=f[j]->f[1][NodeIn.Get_index()];
					fi[2]=f[j]->f[2][NodeIn.Get_index()];
					fi[3]=f[j]->f[3][NodeIn.Get_index()];
					fi[4]=f[j]->f[2][NodeIn.Get_index()];
					fi[5]=f[j]->f[5][NodeIn.Get_index()];
					fi[6]=f[j]->f[6][NodeIn.Get_index()];
					fi[7]=f[j]->f[6][NodeIn.Get_index()];
					fi[8]=f[j]->f[5][NodeIn.Get_index()];
					//Applying He Zou Pressure
					BC_HeZou_P(NodeIn.Get_BcNormal(),fi,NodeIn.Get_RhoDef(), U[0][NodeIn.Get_index()],U[1][NodeIn.Get_index()]);
					for (int i=0;i<9;i++)
						f[j]->f[i][NodeIn.Get_index()]=fi[i];
				}
				break;
			case 1:
				for(int j=0;j<2;j++)
				{
					//Applying symmetry
					fi[0]=f[j]->f[0][NodeIn.Get_index()];
					fi[1]=f[j]->f[3][NodeIn.Get_index()];
					fi[2]=f[j]->f[2][NodeIn.Get_index()];
					fi[3]=f[j]->f[3][NodeIn.Get_index()];
					fi[4]=f[j]->f[4][NodeIn.Get_index()];
					fi[5]=f[j]->f[6][NodeIn.Get_index()];
					fi[6]=f[j]->f[6][NodeIn.Get_index()];
					fi[7]=f[j]->f[7][NodeIn.Get_index()];
					fi[8]=f[j]->f[7][NodeIn.Get_index()];
					//Applying He Zou Pressure
					BC_HeZou_P(NodeIn.Get_BcNormal(),fi,NodeIn.Get_RhoDef(), U[0][NodeIn.Get_index()],U[1][NodeIn.Get_index()]);
					for (int i=0;i<9;i++)
						f[j]->f[i][NodeIn.Get_index()]=fi[i];
				}
				break;
			case 3:
				for(int j=0;j<2;j++)
				{
					//Applying symmetry
					fi[0]=f[j]->f[0][NodeIn.Get_index()];
					fi[1]=f[j]->f[1][NodeIn.Get_index()];
					fi[2]=f[j]->f[2][NodeIn.Get_index()];
					fi[3]=f[j]->f[1][NodeIn.Get_index()];
					fi[4]=f[j]->f[4][NodeIn.Get_index()];
					fi[5]=f[j]->f[5][NodeIn.Get_index()];
					fi[6]=f[j]->f[5][NodeIn.Get_index()];
					fi[7]=f[j]->f[8][NodeIn.Get_index()];
					fi[8]=f[j]->f[8][NodeIn.Get_index()];
					//Applying He Zou Pressure
					BC_HeZou_P(NodeIn.Get_BcNormal(),fi,NodeIn.Get_RhoDef(), U[0][NodeIn.Get_index()],U[1][NodeIn.Get_index()]);
					for (int i=0;i<9;i++)
						f[j]->f[i][NodeIn.Get_index()]=fi[i];
				}
				break;
			default:
				std::cerr<<"Direction symmetry not found"<<std::endl;
				break;
			}
	delete [] fi;
}
///Bounceback Wall treatment
void D2Q9TwoPhases::ApplyBounceBack(NodeWall2D& NodeIn){

	switch(NodeIn.Get_BcNormal())
			{
			case 2:
				f[0]->f[2][NodeIn.Get_index()]=f[0]->f[Opposite[2]][NodeIn.Get_index()];
				f[0]->f[5][NodeIn.Get_index()]=f[0]->f[Opposite[5]][NodeIn.Get_index()];
				f[0]->f[6][NodeIn.Get_index()]=f[0]->f[Opposite[6]][NodeIn.Get_index()];

				f[1]->f[2][NodeIn.Get_index()]=f[1]->f[Opposite[2]][NodeIn.Get_index()];
				f[1]->f[5][NodeIn.Get_index()]=f[1]->f[Opposite[5]][NodeIn.Get_index()];
				f[1]->f[6][NodeIn.Get_index()]=f[1]->f[Opposite[6]][NodeIn.Get_index()];
				break;
			case 4:
				f[0]->f[4][NodeIn.Get_index()]=f[0]->f[Opposite[4]][NodeIn.Get_index()];
				f[0]->f[7][NodeIn.Get_index()]=f[0]->f[Opposite[7]][NodeIn.Get_index()];
				f[0]->f[8][NodeIn.Get_index()]=f[0]->f[Opposite[8]][NodeIn.Get_index()];

				f[1]->f[4][NodeIn.Get_index()]=f[1]->f[Opposite[4]][NodeIn.Get_index()];
				f[1]->f[7][NodeIn.Get_index()]=f[1]->f[Opposite[7]][NodeIn.Get_index()];
				f[1]->f[8][NodeIn.Get_index()]=f[1]->f[Opposite[8]][NodeIn.Get_index()];
				break;
			case 1:
				f[0]->f[1][NodeIn.Get_index()]=f[0]->f[Opposite[1]][NodeIn.Get_index()];
				f[0]->f[5][NodeIn.Get_index()]=f[0]->f[Opposite[5]][NodeIn.Get_index()];
				f[0]->f[8][NodeIn.Get_index()]=f[0]->f[Opposite[8]][NodeIn.Get_index()];

				f[1]->f[1][NodeIn.Get_index()]=f[1]->f[Opposite[1]][NodeIn.Get_index()];
				f[1]->f[5][NodeIn.Get_index()]=f[1]->f[Opposite[5]][NodeIn.Get_index()];
				f[1]->f[8][NodeIn.Get_index()]=f[1]->f[Opposite[8]][NodeIn.Get_index()];
				break;
			case 3:
				f[0]->f[3][NodeIn.Get_index()]=f[0]->f[Opposite[3]][NodeIn.Get_index()];
				f[0]->f[6][NodeIn.Get_index()]=f[0]->f[Opposite[6]][NodeIn.Get_index()];
				f[0]->f[7][NodeIn.Get_index()]=f[0]->f[Opposite[7]][NodeIn.Get_index()];


				f[1]->f[3][NodeIn.Get_index()]=f[1]->f[Opposite[3]][NodeIn.Get_index()];
				f[1]->f[6][NodeIn.Get_index()]=f[1]->f[Opposite[6]][NodeIn.Get_index()];
				f[1]->f[7][NodeIn.Get_index()]=f[1]->f[Opposite[7]][NodeIn.Get_index()];				break;
			default:
				std::cerr<<"Direction wall bounce back not found. Index: "<<NodeIn.Get_index()<<" x: "<<NodeIn.get_x()<<" y: "<<NodeIn.get_y()<<" direction not found: "<<NodeIn.Get_BcNormal()<<std::endl;
				break;
			}
}
void D2Q9TwoPhases::ApplyBounceBackSymmetry(NodeWall2D& NodeIn){

	switch(NodeIn.Get_BcNormal())
			{
			case 2:
				f[0]->f[2][NodeIn.Get_index()]=f[0]->f[Opposite[2]][NodeIn.Get_index()];
				f[0]->f[5][NodeIn.Get_index()]=f[0]->f[Opposite[5]][NodeIn.Get_index()];
				f[0]->f[6][NodeIn.Get_index()]=f[0]->f[Opposite[6]][NodeIn.Get_index()];
				break;
			case 4:
				f[0]->f[4][NodeIn.Get_index()]=f[0]->f[Opposite[4]][NodeIn.Get_index()];
				f[0]->f[7][NodeIn.Get_index()]=f[0]->f[Opposite[7]][NodeIn.Get_index()];
				f[0]->f[8][NodeIn.Get_index()]=f[0]->f[Opposite[8]][NodeIn.Get_index()];
				break;
			case 1:
				//apply Symmetry
				if(NodeArrays->TypeOfNode[NodeIn.Get_connect(0)]==Solid)//Bottom
				{
					f[0]->f[2][NodeIn.Get_index()]=f[0]->f[Opposite[2]][NodeIn.Get_index()];
					f[0]->f[6][NodeIn.Get_index()]=f[0]->f[7][NodeIn.Get_index()];

					f[1]->f[2][NodeIn.Get_index()]=f[1]->f[Opposite[2]][NodeIn.Get_index()];
					f[1]->f[6][NodeIn.Get_index()]=f[1]->f[7][NodeIn.Get_index()];
				}
				else //Top
				{
					f[0]->f[4][NodeIn.Get_index()]=f[0]->f[Opposite[4]][NodeIn.Get_index()];
					f[0]->f[7][NodeIn.Get_index()]=f[0]->f[6][NodeIn.Get_index()];

					f[1]->f[4][NodeIn.Get_index()]=f[1]->f[Opposite[4]][NodeIn.Get_index()];
					f[1]->f[7][NodeIn.Get_index()]=f[1]->f[6][NodeIn.Get_index()];
				}
				f[0]->f[1][NodeIn.Get_index()]=f[0]->f[Opposite[1]][NodeIn.Get_index()];
				f[0]->f[5][NodeIn.Get_index()]=f[0]->f[Opposite[5]][NodeIn.Get_index()];
				f[0]->f[8][NodeIn.Get_index()]=f[0]->f[Opposite[8]][NodeIn.Get_index()];

				f[1]->f[1][NodeIn.Get_index()]=f[1]->f[Opposite[1]][NodeIn.Get_index()];
				f[1]->f[5][NodeIn.Get_index()]=f[1]->f[Opposite[5]][NodeIn.Get_index()];
				f[1]->f[8][NodeIn.Get_index()]=f[1]->f[Opposite[8]][NodeIn.Get_index()];
				break;
			case 3:
				//apply Symmetry
				if(NodeArrays->TypeOfNode[NodeIn.Get_connect(0)]==Solid)//Bottom
				{
					f[0]->f[2][NodeIn.Get_index()]=f[0]->f[Opposite[2]][NodeIn.Get_index()];
					f[0]->f[5][NodeIn.Get_index()]=f[0]->f[8][NodeIn.Get_index()];

					f[1]->f[2][NodeIn.Get_index()]=f[1]->f[Opposite[2]][NodeIn.Get_index()];
					f[1]->f[5][NodeIn.Get_index()]=f[1]->f[8][NodeIn.Get_index()];
				}
				else //Top
				{
					f[0]->f[4][NodeIn.Get_index()]=f[0]->f[Opposite[4]][NodeIn.Get_index()];
					f[0]->f[8][NodeIn.Get_index()]=f[0]->f[5][NodeIn.Get_index()];

					f[1]->f[4][NodeIn.Get_index()]=f[1]->f[Opposite[4]][NodeIn.Get_index()];
					f[1]->f[8][NodeIn.Get_index()]=f[1]->f[5][NodeIn.Get_index()];
				}
				f[0]->f[3][NodeIn.Get_index()]=f[0]->f[Opposite[3]][NodeIn.Get_index()];
				f[0]->f[6][NodeIn.Get_index()]=f[0]->f[Opposite[6]][NodeIn.Get_index()];
				f[0]->f[7][NodeIn.Get_index()]=f[0]->f[Opposite[7]][NodeIn.Get_index()];

				f[1]->f[3][NodeIn.Get_index()]=f[1]->f[Opposite[3]][NodeIn.Get_index()];
				f[1]->f[6][NodeIn.Get_index()]=f[1]->f[Opposite[6]][NodeIn.Get_index()];
				f[1]->f[7][NodeIn.Get_index()]=f[1]->f[Opposite[7]][NodeIn.Get_index()];
				break;
			default:
				std::cerr<<"Direction wall bounce back not found. Index: "<<NodeIn.Get_index()<<std::endl;
				break;
			}
}
void D2Q9TwoPhases::ApplyBounceBack(NodeCorner2D& NodeIn){
	double tmp=0;
	switch(NodeIn.Get_BcNormal())
			{
			case 5:
				f[0]->f[5][NodeIn.Get_index()]=f[0]->f[Opposite[5]][NodeIn.Get_index()];
				f[0]->f[1][NodeIn.Get_index()]=f[0]->f[Opposite[1]][NodeIn.Get_index()];
				f[0]->f[2][NodeIn.Get_index()]=f[0]->f[Opposite[2]][NodeIn.Get_index()];

				f[1]->f[5][NodeIn.Get_index()]=f[1]->f[Opposite[5]][NodeIn.Get_index()];
				f[1]->f[1][NodeIn.Get_index()]=f[1]->f[Opposite[1]][NodeIn.Get_index()];
				f[1]->f[2][NodeIn.Get_index()]=f[1]->f[Opposite[2]][NodeIn.Get_index()];
				if (NodeIn.stream()[3]==false)
				{

					tmp=f[0]->f[0][NodeIn.Get_index()]+f[0]->f[1][NodeIn.Get_index()]+f[0]->f[2][NodeIn.Get_index()]+f[0]->f[3][NodeIn.Get_index()]+f[0]->f[4][NodeIn.Get_index()]+f[0]->f[5][NodeIn.Get_index()]+f[0]->f[7][NodeIn.Get_index()];
					f[0]->f[6][NodeIn.Get_index()]=(Cal_RhoR_Corner(NodeIn)-tmp)/2;
					f[0]->f[8][NodeIn.Get_index()]=(Cal_RhoR_Corner(NodeIn)-tmp)/2;

					tmp=f[1]->f[0][NodeIn.Get_index()]+f[1]->f[1][NodeIn.Get_index()]+f[1]->f[2][NodeIn.Get_index()]+f[1]->f[3][NodeIn.Get_index()]+f[1]->f[4][NodeIn.Get_index()]+f[1]->f[5][NodeIn.Get_index()]+f[1]->f[7][NodeIn.Get_index()];
					f[1]->f[6][NodeIn.Get_index()]=(Cal_RhoB_Corner(NodeIn)-tmp)/2;
					f[1]->f[8][NodeIn.Get_index()]=(Cal_RhoB_Corner(NodeIn)-tmp)/2;
				}
				else
				{
					tmp=f[0]->f[6][NodeIn.Get_index()]+f[0]->f[8][NodeIn.Get_index()];
					f[0]->f[6][NodeIn.Get_index()]=tmp*0.5;
					f[0]->f[8][NodeIn.Get_index()]=tmp*0.5;

					tmp=f[1]->f[6][NodeIn.Get_index()]+f[1]->f[8][NodeIn.Get_index()];
					f[1]->f[6][NodeIn.Get_index()]=tmp*0.5;
					f[1]->f[8][NodeIn.Get_index()]=tmp*0.5;
				}
				break;
			case 6:
				f[0]->f[6][NodeIn.Get_index()]=f[0]->f[Opposite[6]][NodeIn.Get_index()];
				f[0]->f[2][NodeIn.Get_index()]=f[0]->f[Opposite[2]][NodeIn.Get_index()];
				f[0]->f[3][NodeIn.Get_index()]=f[0]->f[Opposite[3]][NodeIn.Get_index()];

				f[1]->f[6][NodeIn.Get_index()]=f[1]->f[Opposite[6]][NodeIn.Get_index()];
				f[1]->f[2][NodeIn.Get_index()]=f[1]->f[Opposite[2]][NodeIn.Get_index()];
				f[1]->f[3][NodeIn.Get_index()]=f[1]->f[Opposite[3]][NodeIn.Get_index()];
				if (NodeIn.stream()[1]==false)
				{
					tmp=f[0]->f[0][NodeIn.Get_index()]+f[0]->f[1][NodeIn.Get_index()]+f[0]->f[2][NodeIn.Get_index()]+f[0]->f[3][NodeIn.Get_index()]+f[0]->f[4][NodeIn.Get_index()]+f[0]->f[6][NodeIn.Get_index()]+f[0]->f[8][NodeIn.Get_index()];
					f[0]->f[5][NodeIn.Get_index()]=(Cal_RhoR_Corner(NodeIn)-tmp)/2;
					f[0]->f[7][NodeIn.Get_index()]=(Cal_RhoR_Corner(NodeIn)-tmp)/2;

					tmp=f[1]->f[0][NodeIn.Get_index()]+f[1]->f[1][NodeIn.Get_index()]+f[1]->f[2][NodeIn.Get_index()]+f[1]->f[3][NodeIn.Get_index()]+f[1]->f[4][NodeIn.Get_index()]+f[1]->f[6][NodeIn.Get_index()]+f[1]->f[8][NodeIn.Get_index()];
					f[1]->f[5][NodeIn.Get_index()]=(Cal_RhoB_Corner(NodeIn)-tmp)/2;
					f[1]->f[7][NodeIn.Get_index()]=(Cal_RhoB_Corner(NodeIn)-tmp)/2;
				}
				else
				{
					tmp=f[0]->f[7][NodeIn.Get_index()]+f[0]->f[5][NodeIn.Get_index()];
					f[0]->f[5][NodeIn.Get_index()]=tmp*0.5;
					f[0]->f[7][NodeIn.Get_index()]=tmp*0.5;

					tmp=f[1]->f[7][NodeIn.Get_index()]+f[1]->f[5][NodeIn.Get_index()];
					f[1]->f[5][NodeIn.Get_index()]=tmp*0.5;
					f[1]->f[7][NodeIn.Get_index()]=tmp*0.5;
				}
				break;
			case 7:
				f[0]->f[7][NodeIn.Get_index()]=f[0]->f[Opposite[7]][NodeIn.Get_index()];
				f[0]->f[3][NodeIn.Get_index()]=f[0]->f[Opposite[3]][NodeIn.Get_index()];
				f[0]->f[4][NodeIn.Get_index()]=f[0]->f[Opposite[4]][NodeIn.Get_index()];

				f[1]->f[7][NodeIn.Get_index()]=f[1]->f[Opposite[7]][NodeIn.Get_index()];
				f[1]->f[3][NodeIn.Get_index()]=f[1]->f[Opposite[3]][NodeIn.Get_index()];
				f[1]->f[4][NodeIn.Get_index()]=f[1]->f[Opposite[4]][NodeIn.Get_index()];
				if (NodeIn.stream()[1]==false)
				{
					tmp=f[0]->f[0][NodeIn.Get_index()]+f[0]->f[1][NodeIn.Get_index()]+f[0]->f[2][NodeIn.Get_index()]+f[0]->f[3][NodeIn.Get_index()]+f[0]->f[4][NodeIn.Get_index()]+f[0]->f[5][NodeIn.Get_index()]+f[0]->f[7][NodeIn.Get_index()];
					f[0]->f[6][NodeIn.Get_index()]=(Cal_RhoR_Corner(NodeIn)-tmp)/2;
					f[0]->f[8][NodeIn.Get_index()]=(Cal_RhoR_Corner(NodeIn)-tmp)/2;

					tmp=f[1]->f[0][NodeIn.Get_index()]+f[1]->f[1][NodeIn.Get_index()]+f[1]->f[2][NodeIn.Get_index()]+f[1]->f[3][NodeIn.Get_index()]+f[1]->f[4][NodeIn.Get_index()]+f[1]->f[5][NodeIn.Get_index()]+f[1]->f[7][NodeIn.Get_index()];
					f[1]->f[6][NodeIn.Get_index()]=(Cal_RhoB_Corner(NodeIn)-tmp)/2;
					f[1]->f[8][NodeIn.Get_index()]=(Cal_RhoB_Corner(NodeIn)-tmp)/2;
				}
				else
				{
					tmp=f[0]->f[6][NodeIn.Get_index()]+f[0]->f[8][NodeIn.Get_index()];
					f[0]->f[6][NodeIn.Get_index()]=tmp*0.5;
					f[0]->f[8][NodeIn.Get_index()]=tmp*0.5;

					tmp=f[1]->f[6][NodeIn.Get_index()]+f[1]->f[8][NodeIn.Get_index()];
					f[1]->f[6][NodeIn.Get_index()]=tmp*0.5;
					f[1]->f[8][NodeIn.Get_index()]=tmp*0.5;
				}
				break;
			case 8:
				f[0]->f[8][NodeIn.Get_index()]=f[0]->f[Opposite[8]][NodeIn.Get_index()];
				f[0]->f[1][NodeIn.Get_index()]=f[0]->f[Opposite[1]][NodeIn.Get_index()];
				f[0]->f[4][NodeIn.Get_index()]=f[0]->f[Opposite[4]][NodeIn.Get_index()];

				f[1]->f[8][NodeIn.Get_index()]=f[1]->f[Opposite[8]][NodeIn.Get_index()];
				f[1]->f[1][NodeIn.Get_index()]=f[1]->f[Opposite[1]][NodeIn.Get_index()];
				f[1]->f[4][NodeIn.Get_index()]=f[1]->f[Opposite[4]][NodeIn.Get_index()];
				if (NodeIn.stream()[2]==false)
				{
					tmp=f[0]->f[0][NodeIn.Get_index()]+f[0]->f[1][NodeIn.Get_index()]+f[0]->f[2][NodeIn.Get_index()]+f[0]->f[3][NodeIn.Get_index()]+f[0]->f[4][NodeIn.Get_index()]+f[0]->f[6][NodeIn.Get_index()]+f[0]->f[8][NodeIn.Get_index()];
					f[0]->f[5][NodeIn.Get_index()]=(Cal_RhoR_Corner(NodeIn)-tmp)/2;
					f[0]->f[7][NodeIn.Get_index()]=(Cal_RhoR_Corner(NodeIn)-tmp)/2;

					tmp=f[1]->f[0][NodeIn.Get_index()]+f[1]->f[1][NodeIn.Get_index()]+f[1]->f[2][NodeIn.Get_index()]+f[1]->f[3][NodeIn.Get_index()]+f[1]->f[4][NodeIn.Get_index()]+f[1]->f[6][NodeIn.Get_index()]+f[1]->f[8][NodeIn.Get_index()];
					f[1]->f[5][NodeIn.Get_index()]=(Cal_RhoB_Corner(NodeIn)-tmp)/2;
					f[1]->f[7][NodeIn.Get_index()]=(Cal_RhoB_Corner(NodeIn)-tmp)/2;
				}
				else
				{
					tmp=f[0]->f[5][NodeIn.Get_index()]+f[0]->f[7][NodeIn.Get_index()];
					f[0]->f[5][NodeIn.Get_index()]=tmp*0.5;
					f[0]->f[7][NodeIn.Get_index()]=tmp*0.5;

					tmp=f[1]->f[5][NodeIn.Get_index()]+f[1]->f[7][NodeIn.Get_index()];
					f[1]->f[5][NodeIn.Get_index()]=tmp*0.5;
					f[1]->f[7][NodeIn.Get_index()]=tmp*0.5;
				}
				break;
			default:
				std::cerr<<"Direction corner bounce back not found. x:"<<NodeIn.get_x()<<" y: "<<NodeIn.get_y()<<std::endl;
				break;
			}
}

///Diffuse Wall treatment
void D2Q9TwoPhases::ApplyDiffuseWall(NodeWall2D& NodeIn){

	switch(NodeIn.Get_BcNormal())
			{
			case 2:
				rhodiff=(f[0]->f[4][NodeIn.Get_index()]+f[0]->f[7][NodeIn.Get_index()]+f[0]->f[8][NodeIn.Get_index()])/SumWeightS;
				f[0]->f[2][NodeIn.Get_index()]=omega[2]*rhodiff;
				f[0]->f[5][NodeIn.Get_index()]=omega[5]*rhodiff;
				f[0]->f[6][NodeIn.Get_index()]=omega[6]*rhodiff;

				rhodiff=(f[1]->f[4][NodeIn.Get_index()]+f[1]->f[7][NodeIn.Get_index()]+f[1]->f[8][NodeIn.Get_index()])/SumWeightS;
				f[1]->f[2][NodeIn.Get_index()]=omega[2]*rhodiff;
				f[1]->f[5][NodeIn.Get_index()]=omega[5]*rhodiff;
				f[1]->f[6][NodeIn.Get_index()]=omega[6]*rhodiff;
				break;
			case 4:
				rhodiff=(f[0]->f[2][NodeIn.Get_index()]+f[0]->f[5][NodeIn.Get_index()]+f[0]->f[6][NodeIn.Get_index()])/SumWeightN;
				f[0]->f[4][NodeIn.Get_index()]=omega[4]*rhodiff;
				f[0]->f[7][NodeIn.Get_index()]=omega[7]*rhodiff;
				f[0]->f[8][NodeIn.Get_index()]=omega[8]*rhodiff;

				rhodiff=(f[1]->f[2][NodeIn.Get_index()]+f[1]->f[5][NodeIn.Get_index()]+f[1]->f[6][NodeIn.Get_index()])/SumWeightN;
				f[1]->f[4][NodeIn.Get_index()]=omega[4]*rhodiff;
				f[1]->f[7][NodeIn.Get_index()]=omega[7]*rhodiff;
				f[1]->f[8][NodeIn.Get_index()]=omega[8]*rhodiff;
				break;
			case 1:
				rhodiff=(f[0]->f[3][NodeIn.Get_index()]+f[0]->f[6][NodeIn.Get_index()]+f[0]->f[7][NodeIn.Get_index()])/SumWeightW;
				f[0]->f[1][NodeIn.Get_index()]=omega[1]*rhodiff;
				f[0]->f[5][NodeIn.Get_index()]=omega[5]*rhodiff;
				f[0]->f[8][NodeIn.Get_index()]=omega[8]*rhodiff;

				rhodiff=(f[1]->f[3][NodeIn.Get_index()]+f[1]->f[6][NodeIn.Get_index()]+f[1]->f[7][NodeIn.Get_index()])/SumWeightW;
				f[1]->f[1][NodeIn.Get_index()]=omega[1]*rhodiff;
				f[1]->f[5][NodeIn.Get_index()]=omega[5]*rhodiff;
				f[1]->f[8][NodeIn.Get_index()]=omega[8]*rhodiff;
				break;
			case 3:
				rhodiff=(f[0]->f[1][NodeIn.Get_index()]+f[0]->f[5][NodeIn.Get_index()]+f[0]->f[8][NodeIn.Get_index()])/SumWeightE;
				f[0]->f[3][NodeIn.Get_index()]=omega[3]*rhodiff;
				f[0]->f[6][NodeIn.Get_index()]=omega[6]*rhodiff;
				f[0]->f[7][NodeIn.Get_index()]=omega[7]*rhodiff;

				rhodiff=(f[1]->f[1][NodeIn.Get_index()]+f[1]->f[5][NodeIn.Get_index()]+f[1]->f[8][NodeIn.Get_index()])/SumWeightE;
				f[1]->f[3][NodeIn.Get_index()]=omega[3]*rhodiff;
				f[1]->f[6][NodeIn.Get_index()]=omega[6]*rhodiff;
				f[1]->f[7][NodeIn.Get_index()]=omega[7]*rhodiff;
				break;
			default:
				std::cerr<<"Direction: "<< NodeIn.Get_BcNormal()<<" (Wall diffuse boundary condition) not found"<<" x: "<<NodeIn.get_x()<<" y: "<<NodeIn.get_y()<<" processor: "<<parallel->getRank()<<std::endl;
				break;
			}
}
void D2Q9TwoPhases::ApplyDiffuseWallSymmetry(NodeWall2D& NodeIn){

	switch(NodeIn.Get_BcNormal())
			{
			case 2:
				rhodiff=(f[0]->f[4][NodeIn.Get_index()]+f[0]->f[7][NodeIn.Get_index()]+f[0]->f[8][NodeIn.Get_index()])/SumWeightS;
				f[0]->f[2][NodeIn.Get_index()]=omega[2]*rhodiff;
				f[0]->f[5][NodeIn.Get_index()]=omega[5]*rhodiff;
				f[0]->f[6][NodeIn.Get_index()]=omega[6]*rhodiff;

				rhodiff=(f[1]->f[4][NodeIn.Get_index()]+f[1]->f[7][NodeIn.Get_index()]+f[1]->f[8][NodeIn.Get_index()])/SumWeightS;
				f[1]->f[2][NodeIn.Get_index()]=omega[2]*rhodiff;
				f[1]->f[5][NodeIn.Get_index()]=omega[5]*rhodiff;
				f[1]->f[6][NodeIn.Get_index()]=omega[6]*rhodiff;
				break;
			case 4:
				rhodiff=(f[0]->f[2][NodeIn.Get_index()]+f[0]->f[5][NodeIn.Get_index()]+f[0]->f[6][NodeIn.Get_index()])/SumWeightN;
				f[0]->f[4][NodeIn.Get_index()]=omega[4]*rhodiff;
				f[0]->f[7][NodeIn.Get_index()]=omega[7]*rhodiff;
				f[0]->f[8][NodeIn.Get_index()]=omega[8]*rhodiff;

				rhodiff=(f[1]->f[2][NodeIn.Get_index()]+f[1]->f[5][NodeIn.Get_index()]+f[1]->f[6][NodeIn.Get_index()])/SumWeightN;
				f[1]->f[4][NodeIn.Get_index()]=omega[4]*rhodiff;
				f[1]->f[7][NodeIn.Get_index()]=omega[7]*rhodiff;
				f[1]->f[8][NodeIn.Get_index()]=omega[8]*rhodiff;
				break;
			case 1:
				//apply Symmetry
				if(NodeArrays->TypeOfNode[NodeIn.Get_connect(0)]==Solid)//Bottom
				{
					f[0]->f[2][NodeIn.Get_index()]=f[0]->f[Opposite[2]][NodeIn.Get_index()];
					f[0]->f[6][NodeIn.Get_index()]=f[0]->f[7][NodeIn.Get_index()];

					f[1]->f[2][NodeIn.Get_index()]=f[1]->f[Opposite[2]][NodeIn.Get_index()];
					f[1]->f[6][NodeIn.Get_index()]=f[1]->f[7][NodeIn.Get_index()];
				}
				else //Top
				{
					f[0]->f[4][NodeIn.Get_index()]=f[0]->f[Opposite[4]][NodeIn.Get_index()];
					f[0]->f[7][NodeIn.Get_index()]=f[0]->f[6][NodeIn.Get_index()];

					f[1]->f[4][NodeIn.Get_index()]=f[1]->f[Opposite[4]][NodeIn.Get_index()];
					f[1]->f[7][NodeIn.Get_index()]=f[1]->f[6][NodeIn.Get_index()];
				}
				rhodiff=(f[0]->f[3][NodeIn.Get_index()]+f[0]->f[6][NodeIn.Get_index()]+f[0]->f[7][NodeIn.Get_index()])/SumWeightW;
				f[0]->f[1][NodeIn.Get_index()]=omega[1]*rhodiff;
				f[0]->f[5][NodeIn.Get_index()]=omega[5]*rhodiff;
				f[0]->f[8][NodeIn.Get_index()]=omega[8]*rhodiff;

				rhodiff=(f[1]->f[3][NodeIn.Get_index()]+f[1]->f[6][NodeIn.Get_index()]+f[1]->f[7][NodeIn.Get_index()])/SumWeightW;
				f[1]->f[1][NodeIn.Get_index()]=omega[1]*rhodiff;
				f[1]->f[5][NodeIn.Get_index()]=omega[5]*rhodiff;
				f[1]->f[8][NodeIn.Get_index()]=omega[8]*rhodiff;
				break;
			case 3:
				//apply Symmetry
				if(NodeArrays->TypeOfNode[NodeIn.Get_connect(0)]==Solid)//Bottom
				{
					f[0]->f[2][NodeIn.Get_index()]=f[0]->f[Opposite[2]][NodeIn.Get_index()];
					f[0]->f[5][NodeIn.Get_index()]=f[0]->f[8][NodeIn.Get_index()];

					f[1]->f[2][NodeIn.Get_index()]=f[1]->f[Opposite[2]][NodeIn.Get_index()];
					f[1]->f[5][NodeIn.Get_index()]=f[1]->f[8][NodeIn.Get_index()];
				}
				else //Top
				{
					f[0]->f[4][NodeIn.Get_index()]=f[0]->f[Opposite[4]][NodeIn.Get_index()];
					f[0]->f[8][NodeIn.Get_index()]=f[0]->f[5][NodeIn.Get_index()];

					f[1]->f[4][NodeIn.Get_index()]=f[1]->f[Opposite[4]][NodeIn.Get_index()];
					f[1]->f[8][NodeIn.Get_index()]=f[1]->f[5][NodeIn.Get_index()];
				}
				rhodiff=(f[0]->f[1][NodeIn.Get_index()]+f[0]->f[5][NodeIn.Get_index()]+f[0]->f[8][NodeIn.Get_index()])/SumWeightE;
				f[0]->f[3][NodeIn.Get_index()]=omega[3]*rhodiff;
				f[0]->f[6][NodeIn.Get_index()]=omega[6]*rhodiff;
				f[0]->f[7][NodeIn.Get_index()]=omega[7]*rhodiff;

				rhodiff=(f[1]->f[1][NodeIn.Get_index()]+f[1]->f[5][NodeIn.Get_index()]+f[1]->f[8][NodeIn.Get_index()])/SumWeightE;
				f[1]->f[3][NodeIn.Get_index()]=omega[3]*rhodiff;
				f[1]->f[6][NodeIn.Get_index()]=omega[6]*rhodiff;
				f[1]->f[7][NodeIn.Get_index()]=omega[7]*rhodiff;
				break;
			default:
				std::cerr<<"Direction: "<< NodeIn.Get_BcNormal()<<" (Wall diffuse boundary condition) not found"<<std::endl;
				break;
			}
}
///Diffuse boundary condition for Concave and Convex corners
///For Concave corners, the incoming distributions go to outcoming so divide by 0.5 and the diagonal distribution doesn't participate for the streaming is approximate by the density at the neibourgh nodes
///For Convex corners, the diffusion goes in all directions.
void D2Q9TwoPhases::ApplyDiffuseWall(NodeCorner2D& NodeIn){
	double tmp=0;
	switch(NodeIn.Get_BcNormal())
			{
			case 5:
			/*	if (NodeIn.stream()[3]==false)
				{
					rhodiff=(f[0]->f[3][NodeIn.Get_index()]+f[0]->f[4][NodeIn.Get_index()]+f[0]->f[7][NodeIn.Get_index()])*2;
					f[0]->f[5][NodeIn.Get_index()]=omega[5]*rhodiff;
					f[0]->f[1][NodeIn.Get_index()]=omega[1]*rhodiff;
					f[0]->f[2][NodeIn.Get_index()]=omega[2]*rhodiff;
					double sumfi=f[0]->f[0][NodeIn.Get_index()]+f[0]->f[1][NodeIn.Get_index()]+f[0]->f[2][NodeIn.Get_index()]+f[0]->f[3][NodeIn.Get_index()]+f[0]->f[4][NodeIn.Get_index()]+f[0]->f[5][NodeIn.Get_index()]+f[0]->f[7][NodeIn.Get_index()];
					f[0]->f[6][NodeIn.Get_index()]=(Cal_RhoCorner(NodeIn)-sumfi)/2;
					f[0]->f[8][NodeIn.Get_index()]=(Cal_RhoCorner(NodeIn)-sumfi)/2;
				}
				else
				{*/
					rhodiff=(f[0]->f[3][NodeIn.Get_index()]+f[0]->f[4][NodeIn.Get_index()]+f[0]->f[7][NodeIn.Get_index()])/SumWeightConvexNE;
					f[0]->f[5][NodeIn.Get_index()]=omega[5]*rhodiff;
					f[0]->f[1][NodeIn.Get_index()]=omega[1]*rhodiff;
					f[0]->f[2][NodeIn.Get_index()]=omega[2]*rhodiff;
					f[0]->f[6][NodeIn.Get_index()]=omega[6]*rhodiff;
					f[0]->f[8][NodeIn.Get_index()]=omega[8]*rhodiff;

					rhodiff=(f[1]->f[3][NodeIn.Get_index()]+f[1]->f[4][NodeIn.Get_index()]+f[1]->f[7][NodeIn.Get_index()])/SumWeightConvexNE;
					f[1]->f[5][NodeIn.Get_index()]=omega[5]*rhodiff;
					f[1]->f[1][NodeIn.Get_index()]=omega[1]*rhodiff;
					f[1]->f[2][NodeIn.Get_index()]=omega[2]*rhodiff;
					f[1]->f[6][NodeIn.Get_index()]=omega[6]*rhodiff;
					f[1]->f[8][NodeIn.Get_index()]=omega[8]*rhodiff;
				//}
				break;
			case 6:
			/*	if (NodeIn.stream()[1]==false)
				{
					rhodiff=(f[0]->f[1][NodeIn.Get_index()]+f[0]->f[4][NodeIn.Get_index()]+f[0]->f[8][NodeIn.Get_index()])*2;
					f[0]->f[6][NodeIn.Get_index()]=omega[6]*rhodiff;
					f[0]->f[2][NodeIn.Get_index()]=omega[2]*rhodiff;
					f[0]->f[3][NodeIn.Get_index()]=omega[3]*rhodiff;
					double sumfi=f[0]->f[0][NodeIn.Get_index()]+f[0]->f[1][NodeIn.Get_index()]+f[0]->f[2][NodeIn.Get_index()]+f[0]->f[3][NodeIn.Get_index()]+f[0]->f[4][NodeIn.Get_index()]+f[0]->f[6][NodeIn.Get_index()]+f[0]->f[8][NodeIn.Get_index()];
					f[0]->f[5][NodeIn.Get_index()]=(Cal_RhoCorner(NodeIn)-sumfi)/2;
					f[0]->f[7][NodeIn.Get_index()]=(Cal_RhoCorner(NodeIn)-sumfi)/2;
				}
				else
				{*/
					rhodiff=(f[0]->f[1][NodeIn.Get_index()]+f[0]->f[4][NodeIn.Get_index()]+f[0]->f[8][NodeIn.Get_index()])/SumWeightConvexNW;
					f[0]->f[6][NodeIn.Get_index()]=omega[6]*rhodiff;
					f[0]->f[2][NodeIn.Get_index()]=omega[2]*rhodiff;
					f[0]->f[3][NodeIn.Get_index()]=omega[3]*rhodiff;
					f[0]->f[5][NodeIn.Get_index()]=omega[5]*rhodiff;
					f[0]->f[7][NodeIn.Get_index()]=omega[7]*rhodiff;

					rhodiff=(f[1]->f[1][NodeIn.Get_index()]+f[1]->f[4][NodeIn.Get_index()]+f[1]->f[8][NodeIn.Get_index()])/SumWeightConvexNW;
					f[1]->f[6][NodeIn.Get_index()]=omega[6]*rhodiff;
					f[1]->f[2][NodeIn.Get_index()]=omega[2]*rhodiff;
					f[1]->f[3][NodeIn.Get_index()]=omega[3]*rhodiff;
					f[1]->f[5][NodeIn.Get_index()]=omega[5]*rhodiff;
					f[1]->f[7][NodeIn.Get_index()]=omega[7]*rhodiff;
				//}
				break;
			case 7:
				/*if (NodeIn.stream()[1]==false)
				{
					rhodiff=(f[0]->f[1][NodeIn.Get_index()]+f[0]->f[2][NodeIn.Get_index()]+f[0]->f[5][NodeIn.Get_index()])*2;
					f[0]->f[7][NodeIn.Get_index()]=omega[7]*rhodiff;
					f[0]->f[3][NodeIn.Get_index()]=omega[3]*rhodiff;
					f[0]->f[4][NodeIn.Get_index()]=omega[4]*rhodiff;
					double sumfi=f[0]->f[0][NodeIn.Get_index()]+f[0]->f[1][NodeIn.Get_index()]+f[0]->f[2][NodeIn.Get_index()]+f[0]->f[3][NodeIn.Get_index()]+f[0]->f[4][NodeIn.Get_index()]+f[0]->f[5][NodeIn.Get_index()]+f[0]->f[7][NodeIn.Get_index()];
					f[0]->f[6][NodeIn.Get_index()]=(Cal_RhoCorner(NodeIn)-sumfi)/2;
					f[0]->f[8][NodeIn.Get_index()]=(Cal_RhoCorner(NodeIn)-sumfi)/2;
				}
				else
				{*/
					rhodiff=(f[0]->f[1][NodeIn.Get_index()]+f[0]->f[2][NodeIn.Get_index()]+f[0]->f[5][NodeIn.Get_index()])/SumWeightConvexSW;
					f[0]->f[7][NodeIn.Get_index()]=omega[7]*rhodiff;
					f[0]->f[3][NodeIn.Get_index()]=omega[3]*rhodiff;
					f[0]->f[4][NodeIn.Get_index()]=omega[4]*rhodiff;
					f[0]->f[6][NodeIn.Get_index()]=omega[6]*rhodiff;
					f[0]->f[8][NodeIn.Get_index()]=omega[8]*rhodiff;

					rhodiff=(f[1]->f[1][NodeIn.Get_index()]+f[1]->f[2][NodeIn.Get_index()]+f[1]->f[5][NodeIn.Get_index()])/SumWeightConvexSW;
					f[1]->f[7][NodeIn.Get_index()]=omega[7]*rhodiff;
					f[1]->f[3][NodeIn.Get_index()]=omega[3]*rhodiff;
					f[1]->f[4][NodeIn.Get_index()]=omega[4]*rhodiff;
					f[1]->f[6][NodeIn.Get_index()]=omega[6]*rhodiff;
					f[1]->f[8][NodeIn.Get_index()]=omega[8]*rhodiff;
				//}
				break;
			case 8:
				/*if (NodeIn.stream()[2]==false)
				{
					rhodiff=(f[0]->f[2][NodeIn.Get_index()]+f[0]->f[3][NodeIn.Get_index()]+f[0]->f[6][NodeIn.Get_index()])*2;
					f[0]->f[8][NodeIn.Get_index()]=omega[8]*rhodiff;
					f[0]->f[1][NodeIn.Get_index()]=omega[1]*rhodiff;
					f[0]->f[4][NodeIn.Get_index()]=omega[4]*rhodiff;
					double sumfi=f[0]->f[0][NodeIn.Get_index()]+f[0]->f[1][NodeIn.Get_index()]+f[0]->f[2][NodeIn.Get_index()]+f[0]->f[3][NodeIn.Get_index()]+f[0]->f[4][NodeIn.Get_index()]+f[0]->f[6][NodeIn.Get_index()]+f[0]->f[8][NodeIn.Get_index()];
					f[0]->f[5][NodeIn.Get_index()]=(Cal_RhoCorner(NodeIn)-sumfi)/2;
					f[0]->f[7][NodeIn.Get_index()]=(Cal_RhoCorner(NodeIn)-sumfi)/2;
				}
				else
				{*/
					rhodiff=(f[0]->f[2][NodeIn.Get_index()]+f[0]->f[3][NodeIn.Get_index()]+f[0]->f[6][NodeIn.Get_index()])/SumWeightConvexSE;
					f[0]->f[8][NodeIn.Get_index()]=omega[8]*rhodiff;
					f[0]->f[1][NodeIn.Get_index()]=omega[1]*rhodiff;
					f[0]->f[4][NodeIn.Get_index()]=omega[4]*rhodiff;
					f[0]->f[5][NodeIn.Get_index()]=omega[5]*rhodiff;
					f[0]->f[7][NodeIn.Get_index()]=omega[7]*rhodiff;

					rhodiff=(f[1]->f[2][NodeIn.Get_index()]+f[1]->f[3][NodeIn.Get_index()]+f[1]->f[6][NodeIn.Get_index()])/SumWeightConvexSE;
					f[1]->f[8][NodeIn.Get_index()]=omega[8]*rhodiff;
					f[1]->f[1][NodeIn.Get_index()]=omega[1]*rhodiff;
					f[1]->f[4][NodeIn.Get_index()]=omega[4]*rhodiff;
					f[1]->f[5][NodeIn.Get_index()]=omega[5]*rhodiff;
					f[1]->f[7][NodeIn.Get_index()]=omega[7]*rhodiff;
				//}
				break;
			default:
				std::cerr<<"Direction: "<< NodeIn.Get_BcNormal()<<" (Corner diffuse boundary conditions) not found"<<std::endl;
				break;
			}
}


/// Corner treat by Chih-Fung Ho, Cheng Chang, Kuen-Hau Lin and Chao-An Lin
/// Consistent Boundary Conditions for 2D and 3D Lattice Boltzmann Simulations
void D2Q9TwoPhases::ApplyCorner(NodeCorner2D& NodeIn){
	double* fi=new double [9];
	for (int i=0;i<9;i++)
		fi[i]=f[0]->f[i][NodeIn.Get_index()];
	BC_corner(NodeIn.Get_BcNormal(),fi,Cal_RhoR_Corner(NodeIn) , NodeIn.Get_UDef()[0],NodeIn.Get_UDef()[1]);
	for (int i=0;i<9;i++)
	f[0]->f[i][NodeIn.Get_index()]=fi[i];

	for (int i=0;i<9;i++)
		fi[i]=f[1]->f[i][NodeIn.Get_index()];
	BC_corner(NodeIn.Get_BcNormal(),fi,Cal_RhoB_Corner(NodeIn) , NodeIn.Get_UDef()[0],NodeIn.Get_UDef()[1]);
	for (int i=0;i<9;i++)
		f[1]->f[i][NodeIn.Get_index()]=fi[i];
	delete [] fi;
}
void D2Q9TwoPhases::ApplyGlobalCorner(NodeCorner2D& NodeIn){
	NodeType NodeTypeTmp1;
	NodeType NodeTypeTmp2;
	int normal=0;
	for(int j=0;j<2;j++)//Not an efficient way to treat the two distribution function but that affects only 4 nodes maximum
	{
		switch(NodeIn.Get_BcNormal())
		{
		case 5:
			NodeTypeTmp1=NodeArrays->TypeOfNode[NodeIn.Get_connect()[2]];
			NodeTypeTmp2=NodeArrays->TypeOfNode[NodeIn.Get_connect()[1]];
			if(NodeTypeTmp1==Symmetry||NodeTypeTmp2==Symmetry)
			{
				if(NodeTypeTmp1==NodeTypeTmp2)
				{
					f[j]->f[1][NodeIn.Get_index()]=f[j]->f[3][NodeIn.Get_index()];//+0.00007/3;
					f[j]->f[2][NodeIn.Get_index()]=f[j]->f[4][NodeIn.Get_index()];
					f[j]->f[6][NodeIn.Get_index()]=f[j]->f[7][NodeIn.Get_index()];
					f[j]->f[8][NodeIn.Get_index()]=f[j]->f[7][NodeIn.Get_index()];//+0.00007/3;
					f[j]->f[5][NodeIn.Get_index()]=f[j]->f[7][NodeIn.Get_index()];//+0.00007/3;
					for (int i=0;i<9;i++)
						f[j]->f[i][NodeIn.Get_index()]=f[j]->f[i][NodeIn.Get_index()]+omega[i]*Ei[i][0]*0.00007;
					//if(NodeTypeTmp1!=Symmetry)
					//	{
					/*normal=1;
					ApplyHeZou_P(NodeIn,normal);
					+0.00007/3*/
				}
	//				std::cerr<<"Symmetry-Symmetry BC not yet implemented"<<std::endl;
				else
				{
				if(NodeTypeTmp1!=Symmetry)
					{
						normal=1;
						f[j]->f[2][NodeIn.Get_index()]=f[j]->f[4][NodeIn.Get_index()];
						f[j]->f[6][NodeIn.Get_index()]=f[j]->f[7][NodeIn.Get_index()];
						if(NodeTypeTmp1==Pressure)
							ApplyHeZou_P(NodeIn,normal);
						else
							ApplyHeZou_U(NodeIn,normal);

					}
					else
					{
						normal=2;
						f[j]->f[8][NodeIn.Get_index()]=f[j]->f[7][NodeIn.Get_index()];
						f[j]->f[1][NodeIn.Get_index()]=f[j]->f[3][NodeIn.Get_index()];
						if(NodeTypeTmp2==Pressure)
							ApplyHeZou_P(NodeIn,normal);
						else
							ApplyHeZou_U(NodeIn,normal);
					}
				}
			}
			else if(NodeTypeTmp1==Wall||NodeTypeTmp2==Wall)
			{
				switch(PtrParameters->Get_WallType())
				{
					case BounceBack:
						ApplyBounceBack(NodeIn);
						break;
					case Diffuse:
						ApplyDiffuseWall(NodeIn);
						break;
				}
				if(NodeTypeTmp1==Pressure)
				{
					normal=1;
					ApplyHeZou_P(NodeIn,normal);
				}
				if(NodeTypeTmp2==Pressure)
				{
					normal=2;
					ApplyHeZou_P(NodeIn,normal);
				}
				if(NodeTypeTmp1==Velocity)
				{
					normal=1;
					ApplyHeZou_U(NodeIn,normal);
				}
				if(NodeTypeTmp2==Velocity)
				{
					normal=2;
					ApplyHeZou_U(NodeIn,normal);
				}
			}
			break;
		case 6:
			NodeTypeTmp1=NodeArrays->TypeOfNode[NodeIn.Get_connect()[2]];
			NodeTypeTmp2=NodeArrays->TypeOfNode[NodeIn.Get_connect()[3]];
			if(NodeTypeTmp1==Symmetry||NodeTypeTmp2==Symmetry)
			{
				if(NodeTypeTmp1==NodeTypeTmp2)
				{
					f[j]->f[2][NodeIn.Get_index()]=f[j]->f[4][NodeIn.Get_index()];
					f[j]->f[3][NodeIn.Get_index()]=f[j]->f[1][NodeIn.Get_index()];//-0.00007/3;
					f[j]->f[6][NodeIn.Get_index()]=f[j]->f[8][NodeIn.Get_index()];//-0.00007/3;
					f[j]->f[5][NodeIn.Get_index()]=f[j]->f[8][NodeIn.Get_index()];
					f[j]->f[7][NodeIn.Get_index()]=f[j]->f[8][NodeIn.Get_index()];//-0.00007/3;
					for (int i=0;i<9;i++)
						f[j]->f[i][NodeIn.Get_index()]=f[j]->f[i][NodeIn.Get_index()]-omega[i]*Ei[i][0]*0.00007;
					/*normal=3;
					ApplyHeZou_P(NodeIn,normal);*/
				}
					//std::cerr<<"Symmetry-Symmetry BC not yet implemented"<<std::endl;
				else
				{
				if(NodeTypeTmp1!=Symmetry)
					{
						normal=3;
						f[j]->f[2][NodeIn.Get_index()]=f[j]->f[4][NodeIn.Get_index()];
						f[j]->f[5][NodeIn.Get_index()]=f[j]->f[8][NodeIn.Get_index()];
						if(NodeTypeTmp1==Pressure)
							ApplyHeZou_P(NodeIn,normal);
						else
							ApplyHeZou_U(NodeIn,normal);

					}
					else
					{
						normal=2;
						f[j]->f[7][NodeIn.Get_index()]=f[j]->f[8][NodeIn.Get_index()];
						f[j]->f[3][NodeIn.Get_index()]=f[j]->f[1][NodeIn.Get_index()];
						if(NodeTypeTmp2==Pressure)
							ApplyHeZou_P(NodeIn,normal);
						else
							ApplyHeZou_U(NodeIn,normal);
					}
				}
			}
			else if(NodeTypeTmp1==Wall||NodeTypeTmp2==Wall)
			{
				switch(PtrParameters->Get_WallType())
				{
					case BounceBack:
						ApplyBounceBack(NodeIn);
						break;
					case Diffuse:
						ApplyDiffuseWall(NodeIn);
						break;
				}
				if(NodeTypeTmp1==Pressure)
				{
					normal=3;
					ApplyHeZou_P(NodeIn,normal);
				}
				if(NodeTypeTmp2==Pressure)
				{
					normal=2;
					ApplyHeZou_P(NodeIn,normal);
				}
				if(NodeTypeTmp1==Velocity)
				{
					normal=3;
					ApplyHeZou_U(NodeIn,normal);
				}
				if(NodeTypeTmp2==Velocity)
				{
					normal=2;
					ApplyHeZou_U(NodeIn,normal);
				}
			}
			break;
		case 7:
			NodeTypeTmp1=NodeArrays->TypeOfNode[NodeIn.Get_connect()[4]];
			NodeTypeTmp2=NodeArrays->TypeOfNode[NodeIn.Get_connect()[3]];
			if(NodeTypeTmp1==Symmetry||NodeTypeTmp2==Symmetry)
			{
				if(NodeTypeTmp1==NodeTypeTmp2)
				{
					f[j]->f[3][NodeIn.Get_index()]=f[j]->f[1][NodeIn.Get_index()];//-0.00007/3;
					f[j]->f[4][NodeIn.Get_index()]=f[j]->f[2][NodeIn.Get_index()];
					f[j]->f[7][NodeIn.Get_index()]=f[j]->f[5][NodeIn.Get_index()];//-0.00007/3;
					f[j]->f[6][NodeIn.Get_index()]=f[j]->f[5][NodeIn.Get_index()];//-0.00007/3;
					f[j]->f[8][NodeIn.Get_index()]=f[j]->f[5][NodeIn.Get_index()];
					for (int i=0;i<9;i++)
						f[j]->f[i][NodeIn.Get_index()]=f[j]->f[i][NodeIn.Get_index()]-omega[i]*Ei[i][0]*0.00007;
					/*normal=3;
					ApplyHeZou_P(NodeIn,normal);*/
				}
				//	std::cerr<<"Symmetry-Symmetry BC not yet implemented"<<std::endl;
				else
				{
				if(NodeTypeTmp1!=Symmetry)
					{
						normal=3;
						f[j]->f[4][NodeIn.Get_index()]=f[j]->f[2][NodeIn.Get_index()];
						f[j]->f[8][NodeIn.Get_index()]=f[j]->f[5][NodeIn.Get_index()];
						if(NodeTypeTmp1==Pressure)
							ApplyHeZou_P(NodeIn,normal);
						else
							ApplyHeZou_U(NodeIn,normal);
					}
					else
					{
						normal=4;
						f[j]->f[3][NodeIn.Get_index()]=f[j]->f[1][NodeIn.Get_index()];
						f[j]->f[6][NodeIn.Get_index()]=f[j]->f[5][NodeIn.Get_index()];
						if(NodeTypeTmp2==Pressure)
							ApplyHeZou_P(NodeIn,normal);
						else
							ApplyHeZou_U(NodeIn,normal);
					}
				}
			}
			else if(NodeTypeTmp1==Wall||NodeTypeTmp2==Wall)
			{
				switch(PtrParameters->Get_WallType())
				{
					case BounceBack:
						ApplyBounceBack(NodeIn);
						break;
					case Diffuse:
						ApplyDiffuseWall(NodeIn);
						break;
				}
				if(NodeTypeTmp1==Pressure)
				{
					normal=3;
					ApplyHeZou_P(NodeIn,normal);
				}
				if(NodeTypeTmp2==Pressure)
				{
					normal=4;
					ApplyHeZou_P(NodeIn,normal);
				}
				if(NodeTypeTmp1==Velocity)
				{
					normal=3;
					ApplyHeZou_U(NodeIn,normal);
				}
				if(NodeTypeTmp2==Velocity)
				{
					normal=4;
					ApplyHeZou_U(NodeIn,normal);
				}
			}
			break;
		case 8:
			NodeTypeTmp1=NodeArrays->TypeOfNode[NodeIn.Get_connect()[4]];
			NodeTypeTmp2=NodeArrays->TypeOfNode[NodeIn.Get_connect()[1]];
			if(NodeTypeTmp1==Symmetry||NodeTypeTmp2==Symmetry)
			{
				if(NodeTypeTmp1==NodeTypeTmp2)
				{
					f[j]->f[1][NodeIn.Get_index()]=f[j]->f[3][NodeIn.Get_index()];//+0.00007/3;
					f[j]->f[4][NodeIn.Get_index()]=f[j]->f[2][NodeIn.Get_index()];
					f[j]->f[8][NodeIn.Get_index()]=f[j]->f[6][NodeIn.Get_index()];//+0.00007/3;
					f[j]->f[5][NodeIn.Get_index()]=f[j]->f[6][NodeIn.Get_index()];//+0.00007/3;
					f[j]->f[7][NodeIn.Get_index()]=f[j]->f[6][NodeIn.Get_index()];
					for (int i=0;i<9;i++)
						f[j]->f[i][NodeIn.Get_index()]=f[j]->f[i][NodeIn.Get_index()]+omega[i]*Ei[i][0]*0.00007;
					/*normal=1;
					ApplyHeZou_P(NodeIn,normal);*/
				}
				//	std::cerr<<"Symmetry-Symmetry BC not yet implemented"<<std::endl;
				else
				{
				if(NodeTypeTmp1!=Symmetry)
					{
						normal=1;
						f[j]->f[4][NodeIn.Get_index()]=f[j]->f[2][NodeIn.Get_index()];
						f[j]->f[7][NodeIn.Get_index()]=f[j]->f[6][NodeIn.Get_index()];
						if(NodeTypeTmp1==Pressure)
							ApplyHeZou_P(NodeIn,normal);
						else
							ApplyHeZou_U(NodeIn,normal);

					}
					else
					{
						normal=4;
						f[j]->f[5][NodeIn.Get_index()]=f[j]->f[6][NodeIn.Get_index()];
						f[j]->f[1][NodeIn.Get_index()]=f[j]->f[3][NodeIn.Get_index()];
						if(NodeTypeTmp2==Pressure)
							ApplyHeZou_P(NodeIn,normal);
						else
							ApplyHeZou_U(NodeIn,normal);
					}
				}
			}
			else if(NodeTypeTmp1==Wall||NodeTypeTmp2==Wall)
			{
				switch(PtrParameters->Get_WallType())
				{
					case BounceBack:
						ApplyBounceBack(NodeIn);
						break;
					case Diffuse:
						ApplyDiffuseWall(NodeIn);
						break;
				}
				if(NodeTypeTmp1==Pressure)
				{
					normal=1;
					ApplyHeZou_P(NodeIn,normal);
				}
				if(NodeTypeTmp2==Pressure)
				{
					normal=4;
					ApplyHeZou_P(NodeIn,normal);
				}
				if(NodeTypeTmp1==Velocity)
				{
					normal=1;
					ApplyHeZou_U(NodeIn,normal);
				}
				if(NodeTypeTmp2==Velocity)
				{
					normal=4;
					ApplyHeZou_U(NodeIn,normal);
				}
			}
			break;
		default:
			std::cerr<<"Direction: "<< NodeIn.Get_BcNormal()<<" not found for Global Corner"<<std::endl;
		}
	}
}
double D2Q9TwoPhases::Cal_RhoCorner(NodeCorner2D& nodeIn){

	unsigned int direction1,direction2;
	switch(nodeIn.Get_BcNormal())
	{
	case 5:
		direction1=1;
		direction2=2;
		doubleTmpReturn=(Rho[nodeIn.Get_connect()[direction1]]+Rho[nodeIn.Get_connect()[direction2]])*0.5;//Connect(nodenumber,direction1)]+Rho[Connect(nodenumber,direction2)])*0.5;
		break;
	case 6:
		direction1=2;
		direction2=3;
		doubleTmpReturn=(Rho[nodeIn.Get_connect()[direction1]]+Rho[nodeIn.Get_connect()[direction2]])*0.5;//(Rho[Connect(nodenumber,direction1)]+Rho[Connect(nodenumber,direction2)])*0.5;
		break;
	case 7:
		direction1=3;
		direction2=4;
		doubleTmpReturn=(Rho[nodeIn.Get_connect()[direction1]]+Rho[nodeIn.Get_connect()[direction2]])*0.5;//(Rho[Connect(nodenumber,direction1)]+Rho[Connect(nodenumber,direction2)])*0.5;
		break;
	case 8:
		direction1=1;
		direction2=4;
		doubleTmpReturn=(Rho[nodeIn.Get_connect()[direction1]]+Rho[nodeIn.Get_connect()[direction2]])*0.5;//(Rho[Connect(nodenumber,direction1)]+Rho[Connect(nodenumber,direction2)])*0.5;
		break;
	}

	return doubleTmpReturn;
}
double D2Q9TwoPhases::Cal_RhoR_Corner(NodeCorner2D& nodeIn){

	unsigned int direction1,direction2;
	switch(nodeIn.Get_BcNormal())
	{
	case 5:
		direction1=1;
		direction2=2;
		doubleTmpReturn=(Rhor[nodeIn.Get_connect()[direction1]]+Rhor[nodeIn.Get_connect()[direction2]])*0.5;//Connect(nodenumber,direction1)]+Rho[Connect(nodenumber,direction2)])*0.5;
		break;
	case 6:
		direction1=2;
		direction2=3;
		doubleTmpReturn=(Rhor[nodeIn.Get_connect()[direction1]]+Rhor[nodeIn.Get_connect()[direction2]])*0.5;//(Rho[Connect(nodenumber,direction1)]+Rho[Connect(nodenumber,direction2)])*0.5;
		break;
	case 7:
		direction1=3;
		direction2=4;
		doubleTmpReturn=(Rhor[nodeIn.Get_connect()[direction1]]+Rhor[nodeIn.Get_connect()[direction2]])*0.5;//(Rho[Connect(nodenumber,direction1)]+Rho[Connect(nodenumber,direction2)])*0.5;
		break;
	case 8:
		direction1=1;
		direction2=4;
		doubleTmpReturn=(Rhor[nodeIn.Get_connect()[direction1]]+Rhor[nodeIn.Get_connect()[direction2]])*0.5;//(Rho[Connect(nodenumber,direction1)]+Rho[Connect(nodenumber,direction2)])*0.5;
		break;
	}

	return doubleTmpReturn;
}
double D2Q9TwoPhases::Cal_RhoB_Corner(NodeCorner2D& nodeIn){

	unsigned int direction1,direction2;
	switch(nodeIn.Get_BcNormal())
	{
	case 5:
		direction1=1;
		direction2=2;
		doubleTmpReturn=(Rhob[nodeIn.Get_connect()[direction1]]+Rhob[nodeIn.Get_connect()[direction2]])*0.5;//Connect(nodenumber,direction1)]+Rho[Connect(nodenumber,direction2)])*0.5;
		break;
	case 6:
		direction1=2;
		direction2=3;
		doubleTmpReturn=(Rhob[nodeIn.Get_connect()[direction1]]+Rhob[nodeIn.Get_connect()[direction2]])*0.5;//(Rho[Connect(nodenumber,direction1)]+Rho[Connect(nodenumber,direction2)])*0.5;
		break;
	case 7:
		direction1=3;
		direction2=4;
		doubleTmpReturn=(Rhob[nodeIn.Get_connect()[direction1]]+Rhob[nodeIn.Get_connect()[direction2]])*0.5;//(Rho[Connect(nodenumber,direction1)]+Rho[Connect(nodenumber,direction2)])*0.5;
		break;
	case 8:
		direction1=1;
		direction2=4;
		doubleTmpReturn=(Rhob[nodeIn.Get_connect()[direction1]]+Rhob[nodeIn.Get_connect()[direction2]])*0.5;//(Rho[Connect(nodenumber,direction1)]+Rho[Connect(nodenumber,direction2)])*0.5;
		break;
	}

	return doubleTmpReturn;
}
void D2Q9TwoPhases::StreamingOrientation(NodeGhost2D& nodeIn, bool Streaming[9]){
	Streaming[0]=false;
	for (unsigned int i=1;i<(unsigned int)nbvelo;i++)
		if((nodeIn.Get_connect()[i]==nodeIn.Get_index())||(nodeIn.get_NodeType()==Solid)||(NodeArrays->TypeOfNode[nodeIn.Get_connect()[i]]==Solid)||(NodeArrays->TypeOfNode[nodeIn.Get_connect()[i]]==Ghost))
			Streaming[i]=false;
		else
			Streaming[i]=true;

	if(nodeIn.get_x()==20 && nodeIn.get_y()==24)
	{
		std::cout<<"streaming ghost: "<<std::endl;
		for (unsigned int i=1;i<(unsigned int)nbvelo;i++)
			std::cout<<Streaming[i];
		std::cout<<std::endl;
	}
}
void D2Q9TwoPhases::StreamingOrientation(NodeCorner2D& nodeIn, bool Streaming[9]){

	if(nodeIn.Get_CornerType()==Concave)
	{
		switch(nodeIn.Get_BcNormal())
			{
			case 5:
				Streaming[0]=false;
				Streaming[1]=true;
				Streaming[2]=true;
				Streaming[3]=false;
				Streaming[4]=false;
				Streaming[5]=true;
				Streaming[6]=false;
				Streaming[7]=false;
				Streaming[8]=false;
				break;
			case 6:
				Streaming[0]=false;
				Streaming[1]=false;
				Streaming[2]=true;
				Streaming[3]=true;
				Streaming[4]=false;
				Streaming[5]=false;
				Streaming[6]=true;
				Streaming[7]=false;
				Streaming[8]=false;
				break;
			case 7:
				Streaming[0]=false;
				Streaming[1]=false;
				Streaming[2]=false;
				Streaming[3]=true;
				Streaming[4]=true;
				Streaming[5]=false;
				Streaming[6]=false;
				Streaming[7]=true;
				Streaming[8]=false;
				break;
			case 8:
				Streaming[0]=false;
				Streaming[1]=true;
				Streaming[2]=false;
				Streaming[3]=false;
				Streaming[4]=true;
				Streaming[5]=false;
				Streaming[6]=false;
				Streaming[7]=false;
				Streaming[8]=true;
				break;
			}
	}
	else
	{
		switch(nodeIn.Get_BcNormal())
			{
			case 5:
				Streaming[0]=false;
				Streaming[1]=true;
				Streaming[2]=true;
				Streaming[3]=true;
				Streaming[4]=true;
				Streaming[5]=true;
				Streaming[6]=true;
				Streaming[7]=false;
				Streaming[8]=true;
				break;
			case 6:
				Streaming[0]=false;
				Streaming[1]=true;
				Streaming[2]=true;
				Streaming[3]=true;
				Streaming[4]=true;
				Streaming[5]=true;
				Streaming[6]=true;
				Streaming[7]=true;
				Streaming[8]=false;
				break;
			case 7:
				Streaming[0]=false;
				Streaming[1]=true;
				Streaming[2]=true;
				Streaming[3]=true;
				Streaming[4]=true;
				Streaming[5]=false;
				Streaming[6]=true;
				Streaming[7]=true;
				Streaming[8]=true;
				break;
			case 8:
				Streaming[0]=false;
				Streaming[1]=true;
				Streaming[2]=true;
				Streaming[3]=true;
				Streaming[4]=true;
				Streaming[5]=true;
				Streaming[6]=false;
				Streaming[7]=true;
				Streaming[8]=true;
				break;
			}
		}

/*	for (unsigned int i=1;i<(unsigned int)nbvelo;i++)
		if((nodeIn.Get_connect()[i]==nodeIn.Get_index())||(nodeIn.get_NodeType()==Solid)||(NodeArrays->TypeOfNode[nodeIn.Get_connect()[i]]==Solid))
			Streaming[i]=false;
		else
			Streaming[i]=true;*/
	if(nodeIn.get_x()==19 && nodeIn.get_y()==24)
	{
		std::cout<<"streaming corner: "<<std::endl;
		for (unsigned int i=1;i<(unsigned int)nbvelo;i++)
			std::cout<<Streaming[i];
		std::cout<<std::endl;
	}
}
void D2Q9TwoPhases::StreamingOrientation(NodeWall2D& nodeIn, bool Streaming[9]){
	switch(nodeIn.Get_BcNormal())
		{
		case 1:
			Streaming[0]=false;
			Streaming[1]=true;
			Streaming[2]=true;
			Streaming[3]=false;
			Streaming[4]=true;
			Streaming[5]=true;
			Streaming[6]=false;
			Streaming[7]=false;
			Streaming[8]=true;
			break;
		case 2:
			Streaming[0]=false;
			Streaming[1]=true;
			Streaming[2]=true;
			Streaming[3]=true;
			Streaming[4]=false;
			Streaming[5]=true;
			Streaming[6]=true;
			Streaming[7]=false;
			Streaming[8]=false;
			break;
		case 3:
			Streaming[0]=false;
			Streaming[1]=false;
			Streaming[2]=true;
			Streaming[3]=true;
			Streaming[4]=true;
			Streaming[5]=false;
			Streaming[6]=true;
			Streaming[7]=true;
			Streaming[8]=false;
			break;
		case 4:
			Streaming[0]=false;
			Streaming[1]=true;
			Streaming[2]=false;
			Streaming[3]=true;
			Streaming[4]=true;
			Streaming[5]=false;
			Streaming[6]=false;
			Streaming[7]=true;
			Streaming[8]=true;
			break;
		}
}
void D2Q9TwoPhases::StreamingOrientation(NodeSymmetry2D& nodeIn, bool Streaming[9]){
	switch(nodeIn.Get_BcNormal())
		{
		case 1:
			Streaming[0]=false;
			Streaming[1]=true;
			Streaming[2]=true;
			Streaming[3]=false;
			Streaming[4]=true;
			Streaming[5]=true;
			Streaming[6]=false;
			Streaming[7]=false;
			Streaming[8]=true;
			break;
		case 2:
			Streaming[0]=false;
			Streaming[1]=true;
			Streaming[2]=true;
			Streaming[3]=true;
			Streaming[4]=false;
			Streaming[5]=true;
			Streaming[6]=true;
			Streaming[7]=false;
			Streaming[8]=false;
			break;
		case 3:
			Streaming[0]=false;
			Streaming[1]=false;
			Streaming[2]=true;
			Streaming[3]=true;
			Streaming[4]=true;
			Streaming[5]=false;
			Streaming[6]=true;
			Streaming[7]=true;
			Streaming[8]=false;
			break;
		case 4:
			Streaming[0]=false;
			Streaming[1]=true;
			Streaming[2]=false;
			Streaming[3]=true;
			Streaming[4]=true;
			Streaming[5]=false;
			Streaming[6]=false;
			Streaming[7]=true;
			Streaming[8]=true;
			break;
		}
}
void D2Q9TwoPhases::StreamingOrientation(NodeVelocity2D& nodeIn, bool Streaming[9]){
	Streaming[0]=false;
	for (unsigned int i=1;i<(unsigned int)nbvelo;i++)
		if((nodeIn.Get_connect()[i]==nodeIn.Get_index())||(nodeIn.get_NodeType()==Solid)||(NodeArrays->TypeOfNode[nodeIn.Get_connect()[i]]==Solid))
			Streaming[i]=false;
		else
			Streaming[i]=true;
}
void D2Q9TwoPhases::StreamingOrientation(NodePressure2D& nodeIn, bool Streaming[9]){
	Streaming[0]=false;
	for (unsigned int i=1;i<(unsigned int)nbvelo;i++)
		if((nodeIn.Get_connect()[i]==nodeIn.Get_index())||(nodeIn.get_NodeType()==Solid)||(NodeArrays->TypeOfNode[nodeIn.Get_connect()[i]]==Solid))
			Streaming[i]=false;
		else
			Streaming[i]=true;
}




//Communication Functions
void D2Q9TwoPhases::IniComVariables(){
	MultiBlock_->Get_Connect_Node(IdNodeN,IdNodeE,IdNodeS,IdNodeW,IdNodeSW,IdNodeSE,IdNodeNW,IdNodeNE);
	MultiBlock_->Get_Connect_Node(IdRNodeN,IdRNodeE,IdRNodeS,IdRNodeW,IdGNodeN,IdGNodeE,IdGNodeS,IdGNodeW,
			IdRNodeSW,IdRNodeSE,IdRNodeNW,IdRNodeNE,IdGNodeSW,IdGNodeSE,IdGNodeNW,IdGNodeNE);
	Nd_variables_sync=3*2;
	buf_send=new double** [Nd_variables_sync];
	buf_recv=new double** [Nd_variables_sync];
	for (int i=0;i<Nd_variables_sync;i++)
	{
		buf_send[i]=new double* [4];
		buf_recv[i]=new double* [4];
	}

	size_buf=new int [4];

	for (int i=0;i<Nd_variables_sync;i++)
	{
		buf_send[i][0]=new double[IdNodeE.size()];
		buf_send[i][1]=new double[IdNodeW.size()];
		buf_send[i][2]=new double[IdNodeS.size()];
		buf_send[i][3]=new double[IdNodeN.size()];
		buf_recv[i][0]=new double[IdNodeW.size()];
		buf_recv[i][1]=new double[IdNodeE.size()];
		buf_recv[i][2]=new double[IdNodeN.size()];
		buf_recv[i][3]=new double[IdNodeS.size()];
	}
	size_buf[0]=IdNodeE.size();
	size_buf[1]=IdNodeW.size();
	size_buf[2]=IdNodeS.size();
	size_buf[3]=IdNodeN.size();
// Macro sync
	Nd_MacroVariables_sync=3;
	buf_MacroSend=new double** [Nd_MacroVariables_sync];
	buf_MacroRecv=new double** [Nd_MacroVariables_sync];
	for (int i=0;i<Nd_MacroVariables_sync;i++)
	{
		buf_MacroSend[i]=new double* [4];
		buf_MacroRecv[i]=new double* [4];
	}

	size_MacroBuf=new int [8];

	for (int i=0;i<Nd_MacroVariables_sync;i++)
	{
		buf_MacroSend[i][0]=new double[IdRNodeE.size()];
		buf_MacroSend[i][1]=new double[IdRNodeW.size()];
		buf_MacroSend[i][2]=new double[IdRNodeS.size()];
		buf_MacroSend[i][3]=new double[IdRNodeN.size()];
		buf_MacroRecv[i][0]=new double[IdGNodeW.size()];
		buf_MacroRecv[i][1]=new double[IdGNodeE.size()];
		buf_MacroRecv[i][2]=new double[IdGNodeN.size()];
		buf_MacroRecv[i][3]=new double[IdGNodeS.size()];
	}
	size_MacroBuf[0]=IdRNodeE.size();
	size_MacroBuf[1]=IdRNodeW.size();
	size_MacroBuf[2]=IdRNodeS.size();
	size_MacroBuf[3]=IdRNodeN.size();
	size_MacroBuf[4]=IdGNodeW.size();
	size_MacroBuf[5]=IdGNodeE.size();
	size_MacroBuf[6]=IdGNodeN.size();
	size_MacroBuf[7]=IdGNodeS.size();
}
void D2Q9TwoPhases::GhostNodesSyncFromGhost(){

	for (unsigned int i=0;i<IdNodeE.size();i++)
	{
		buf_send[0][0][i]=f[0]->f[1][IdNodeE[i]];
		buf_send[1][0][i]=f[0]->f[5][IdNodeE[i]];
		buf_send[2][0][i]=f[0]->f[8][IdNodeE[i]];

		buf_send[3][0][i]=f[1]->f[1][IdNodeE[i]];
		buf_send[4][0][i]=f[1]->f[5][IdNodeE[i]];
		buf_send[5][0][i]=f[1]->f[8][IdNodeE[i]];
	}

	for (unsigned int i=0;i<IdNodeN.size();i++)
	{
		buf_send[0][3][i]=f[0]->f[2][IdNodeN[i]];
		buf_send[1][3][i]=f[0]->f[5][IdNodeN[i]];
		buf_send[2][3][i]=f[0]->f[6][IdNodeN[i]];

		buf_send[3][3][i]=f[1]->f[2][IdNodeN[i]];
		buf_send[4][3][i]=f[1]->f[5][IdNodeN[i]];
		buf_send[5][3][i]=f[1]->f[6][IdNodeN[i]];
	}

	int itmp=0;
	for (int i=0;i<Nd_variables_sync;i++)
	{
		MultiBlock_->CommunicationFromGhost(buf_send[i],buf_recv[i],size_buf);
		//MultiBlock_->CommunicationFromGhost(buf_send[i],buf_recv[i],size_buf,&status[itmp],&request[itmp]);
		//itmp+=2;
	}
	for (unsigned int i=0;i<IdNodeW.size();i++)
	{
		f[0]->f[1][IdNodeW[i]]=buf_recv[0][0][i];
		f[0]->f[5][IdNodeW[i]]=buf_recv[1][0][i];
		f[0]->f[8][IdNodeW[i]]=buf_recv[2][0][i];
	}

	for (unsigned int i=0;i<IdNodeS.size();i++)
	{
		f[0]->f[2][IdNodeS[i]]=buf_recv[0][3][i];
		f[0]->f[5][IdNodeS[i]]=buf_recv[1][3][i];
		f[0]->f[6][IdNodeS[i]]=buf_recv[2][3][i];
	}

}
void D2Q9TwoPhases::GhostNodesSyncToGhost(){
/*	for (unsigned int i=0;i<IdNodeW.size();i++)
	{
		buf_send[0][1][i]=f[0]->f[3][IdNodeW[i]];
		buf_send[1][1][i]=f[0]->f[6][IdNodeW[i]];
		buf_send[2][1][i]=f[0]->f[7][IdNodeW[i]];
	}
	for (unsigned int i=0;i<IdNodeS.size();i++)
	{
		buf_send[0][2][i]=f[0]->f[4][IdNodeS[i]];
		buf_send[1][2][i]=f[0]->f[7][IdNodeS[i]];
		buf_send[2][2][i]=f[0]->f[8][IdNodeS[i]];
	}
	//MPI_Barrier(MPI_COMM_WORLD);
	for (int i=0;i<Nd_variables_sync;i++)
	{
		//MultiBlock_->Communication(buf_send[i],buf_recv[i],size_buf);
		MultiBlock_->CommunicationToGhost(buf_send[i],buf_recv[i],size_buf);

	}
	for (unsigned int i=0;i<IdNodeE.size();i++)
	{
		f[0]->f[3][IdNodeE[i]]=buf_recv[0][1][i];
		f[0]->f[6][IdNodeE[i]]=buf_recv[1][1][i];
		f[0]->f[7][IdNodeE[i]]=buf_recv[2][1][i];
	}
	for (unsigned int i=0;i<IdNodeN.size();i++)
	{
		f[0]->f[4][IdNodeN[i]]=buf_recv[0][2][i];
		f[0]->f[7][IdNodeN[i]]=buf_recv[1][2][i];
		f[0]->f[8][IdNodeN[i]]=buf_recv[2][2][i];
	}*/

	for (unsigned int i=0;i<IdRNodeW.size();i++)
	{
		buf_send[0][1][i]=f[0]->f[3][IdRNodeW[i]];
		buf_send[1][1][i]=f[0]->f[6][IdRNodeW[i]];
		buf_send[2][1][i]=f[0]->f[7][IdRNodeW[i]];

		buf_send[3][1][i]=f[1]->f[3][IdRNodeW[i]];
		buf_send[4][1][i]=f[1]->f[6][IdRNodeW[i]];
		buf_send[5][1][i]=f[1]->f[7][IdRNodeW[i]];

	}
	for (unsigned int i=0;i<IdRNodeS.size();i++)
	{
		buf_send[0][2][i]=f[0]->f[4][IdRNodeS[i]];
		buf_send[1][2][i]=f[0]->f[7][IdRNodeS[i]];
		buf_send[2][2][i]=f[0]->f[8][IdRNodeS[i]];

		buf_send[3][2][i]=f[1]->f[4][IdRNodeS[i]];
		buf_send[4][2][i]=f[1]->f[7][IdRNodeS[i]];
		buf_send[5][2][i]=f[1]->f[8][IdRNodeS[i]];
	}
	for (unsigned int i=0;i<IdRNodeE.size();i++)
	{
		buf_send[0][0][i]=f[0]->f[1][IdRNodeE[i]];
		buf_send[1][0][i]=f[0]->f[5][IdRNodeE[i]];
		buf_send[2][0][i]=f[0]->f[8][IdRNodeE[i]];

		buf_send[3][0][i]=f[1]->f[1][IdRNodeE[i]];
		buf_send[4][0][i]=f[1]->f[5][IdRNodeE[i]];
		buf_send[5][0][i]=f[1]->f[8][IdRNodeE[i]];
	}
	for (unsigned int i=0;i<IdRNodeN.size();i++)
	{
		buf_send[0][3][i]=f[0]->f[2][IdRNodeN[i]];
		buf_send[1][3][i]=f[0]->f[5][IdRNodeN[i]];
		buf_send[2][3][i]=f[0]->f[6][IdRNodeN[i]];

		buf_send[3][3][i]=f[1]->f[2][IdRNodeN[i]];
		buf_send[4][3][i]=f[1]->f[5][IdRNodeN[i]];
		buf_send[5][3][i]=f[1]->f[6][IdRNodeN[i]];
	}
/*	for (int i=0;i<Nd_variables_sync;i++)
	{
		MultiBlock_->Communication(buf_MacroSend[i],buf_MacroRecv[i],size_buf);
		//MultiBlock_->CommunicationToGhost(buf_MacroSend[i],buf_MacroRecv[i],size_MacroBuf);

	}*/
	int tag_x_r=1;
	int tag_x_l=2;
	int tag_y_t=3;
	int tag_y_b=4;
	MPI_Status status;
	if(IdRNodeE.size()>=1)
	{
		MultiBlock_->Send(&buf_send[0][0][0],IdRNodeE.size(),1,tag_x_r);
		MultiBlock_->Send(&buf_send[1][0][0],IdRNodeE.size(),1,tag_x_r);
		MultiBlock_->Send(&buf_send[2][0][0],IdRNodeE.size(),1,tag_x_r);

		MultiBlock_->Send(&buf_send[3][0][0],IdRNodeE.size(),1,tag_x_r);
		MultiBlock_->Send(&buf_send[4][0][0],IdRNodeE.size(),1,tag_x_r);
		MultiBlock_->Send(&buf_send[5][0][0],IdRNodeE.size(),1,tag_x_r);
	}
	if(IdGNodeW.size()>=1)
	{
		MultiBlock_->Recv(&buf_recv[0][0][0],IdGNodeW.size(),3,tag_x_r,status);
		MultiBlock_->Recv(&buf_recv[1][0][0],IdGNodeW.size(),3,tag_x_r,status);
		MultiBlock_->Recv(&buf_recv[2][0][0],IdGNodeW.size(),3,tag_x_r,status);

		MultiBlock_->Recv(&buf_recv[3][0][0],IdGNodeW.size(),3,tag_x_r,status);
		MultiBlock_->Recv(&buf_recv[4][0][0],IdGNodeW.size(),3,tag_x_r,status);
		MultiBlock_->Recv(&buf_recv[5][0][0],IdGNodeW.size(),3,tag_x_r,status);
	}

	if(IdRNodeW.size()>=1)
	{
		MultiBlock_->Send(&buf_send[0][1][0],IdRNodeW.size(),3,tag_x_l);
		MultiBlock_->Send(&buf_send[1][1][0],IdRNodeW.size(),3,tag_x_l);
		MultiBlock_->Send(&buf_send[2][1][0],IdRNodeW.size(),3,tag_x_l);

		MultiBlock_->Send(&buf_send[3][1][0],IdRNodeW.size(),3,tag_x_l);
		MultiBlock_->Send(&buf_send[4][1][0],IdRNodeW.size(),3,tag_x_l);
		MultiBlock_->Send(&buf_send[5][1][0],IdRNodeW.size(),3,tag_x_l);
	}
	if(IdGNodeE.size()>=1)
	{
		MultiBlock_->Recv(&buf_recv[0][1][0],IdGNodeE.size(),1,tag_x_l,status);
		MultiBlock_->Recv(&buf_recv[1][1][0],IdGNodeE.size(),1,tag_x_l,status);
		MultiBlock_->Recv(&buf_recv[2][1][0],IdGNodeE.size(),1,tag_x_l,status);

		MultiBlock_->Recv(&buf_recv[3][1][0],IdGNodeE.size(),1,tag_x_l,status);
		MultiBlock_->Recv(&buf_recv[4][1][0],IdGNodeE.size(),1,tag_x_l,status);
		MultiBlock_->Recv(&buf_recv[5][1][0],IdGNodeE.size(),1,tag_x_l,status);
	}
		if(IdRNodeN.size()>=1)
	{
		MultiBlock_->Send(&buf_send[0][3][0],IdRNodeN.size(),0,tag_y_t);
		MultiBlock_->Send(&buf_send[1][3][0],IdRNodeN.size(),0,tag_y_t);
		MultiBlock_->Send(&buf_send[2][3][0],IdRNodeN.size(),0,tag_y_t);

		MultiBlock_->Send(&buf_send[3][3][0],IdRNodeN.size(),0,tag_y_t);
		MultiBlock_->Send(&buf_send[4][3][0],IdRNodeN.size(),0,tag_y_t);
		MultiBlock_->Send(&buf_send[5][3][0],IdRNodeN.size(),0,tag_y_t);
	}
	if(IdGNodeS.size()>=1)
	{
		MultiBlock_->Recv(&buf_recv[0][3][0],IdGNodeS.size(),2,tag_y_t,status);
		MultiBlock_->Recv(&buf_recv[1][3][0],IdGNodeS.size(),2,tag_y_t,status);
		MultiBlock_->Recv(&buf_recv[2][3][0],IdGNodeS.size(),2,tag_y_t,status);

		MultiBlock_->Recv(&buf_recv[3][3][0],IdGNodeS.size(),2,tag_y_t,status);
		MultiBlock_->Recv(&buf_recv[4][3][0],IdGNodeS.size(),2,tag_y_t,status);
		MultiBlock_->Recv(&buf_recv[5][3][0],IdGNodeS.size(),2,tag_y_t,status);
	}
	if(IdRNodeS.size()>=1)
	{
		MultiBlock_->Send(&buf_send[0][2][0],IdRNodeS.size(),2,tag_y_b);
		MultiBlock_->Send(&buf_send[1][2][0],IdRNodeS.size(),2,tag_y_b);
		MultiBlock_->Send(&buf_send[2][2][0],IdRNodeS.size(),2,tag_y_b);

		MultiBlock_->Send(&buf_send[3][2][0],IdRNodeS.size(),2,tag_y_b);
		MultiBlock_->Send(&buf_send[4][2][0],IdRNodeS.size(),2,tag_y_b);
		MultiBlock_->Send(&buf_send[5][2][0],IdRNodeS.size(),2,tag_y_b);
	}
	if(IdGNodeN.size()>=1)
	{
		MultiBlock_->Recv(&buf_recv[0][2][0],IdGNodeN.size(),0,tag_y_b,status);
		MultiBlock_->Recv(&buf_recv[1][2][0],IdGNodeN.size(),0,tag_y_b,status);
		MultiBlock_->Recv(&buf_recv[2][2][0],IdGNodeN.size(),0,tag_y_b,status);

		MultiBlock_->Recv(&buf_recv[3][2][0],IdGNodeN.size(),0,tag_y_b,status);
		MultiBlock_->Recv(&buf_recv[4][2][0],IdGNodeN.size(),0,tag_y_b,status);
		MultiBlock_->Recv(&buf_recv[5][2][0],IdGNodeN.size(),0,tag_y_b,status);
	}
	//MultiBlock_->Communication(buf_MacroSend[0],buf_MacroRecv[0],size_buf);
	//int* block=MultiBlock_->get_Block_Connect();
	//if(block[1]>=0)
	for (unsigned int i=0;i<IdGNodeE.size();i++)
	{
		f[0]->f[3][IdGNodeE[i]]=buf_recv[0][1][i];
		f[0]->f[6][IdGNodeE[i]]=buf_recv[1][1][i];
		f[0]->f[7][IdGNodeE[i]]=buf_recv[2][1][i];

		f[1]->f[3][IdGNodeE[i]]=buf_recv[3][1][i];
		f[1]->f[6][IdGNodeE[i]]=buf_recv[4][1][i];
		f[1]->f[7][IdGNodeE[i]]=buf_recv[5][1][i];
	}
	//if(block[0]>=0)
	for (unsigned int i=0;i<IdGNodeN.size();i++)
	{
		f[0]->f[4][IdGNodeN[i]]=buf_recv[0][2][i];
		f[0]->f[7][IdGNodeN[i]]=buf_recv[1][2][i];
		f[0]->f[8][IdGNodeN[i]]=buf_recv[2][2][i];

		f[1]->f[4][IdGNodeN[i]]=buf_recv[3][2][i];
		f[1]->f[7][IdGNodeN[i]]=buf_recv[4][2][i];
		f[1]->f[8][IdGNodeN[i]]=buf_recv[5][2][i];
	}
	//if(block[3]>=0)
	for (unsigned int i=0;i<IdGNodeW.size();i++)
	{
		f[0]->f[1][IdGNodeW[i]]=buf_recv[0][0][i];
		f[0]->f[5][IdGNodeW[i]]=buf_recv[1][0][i];
		f[0]->f[8][IdGNodeW[i]]=buf_recv[2][0][i];

		f[1]->f[1][IdGNodeW[i]]=buf_recv[3][0][i];
		f[1]->f[5][IdGNodeW[i]]=buf_recv[4][0][i];
		f[1]->f[8][IdGNodeW[i]]=buf_recv[5][0][i];
	}
	//if(block[2]>=0)
	for (unsigned int i=0;i<IdGNodeS.size();i++)
	{
		f[0]->f[2][IdGNodeS[i]]=buf_recv[0][3][i];
		f[0]->f[5][IdGNodeS[i]]=buf_recv[1][3][i];
		f[0]->f[6][IdGNodeS[i]]=buf_recv[2][3][i];

		f[1]->f[2][IdGNodeS[i]]=buf_recv[3][3][i];
		f[1]->f[5][IdGNodeS[i]]=buf_recv[4][3][i];
		f[1]->f[6][IdGNodeS[i]]=buf_recv[5][3][i];
	}
}
void D2Q9TwoPhases::CornerNodesSyncFromGhost(){
	 int tag_x_r=1;
	 int tag_y_t=2;
	 int tag_d_tr=3;

	//N=0,E=1,S=2,W=3,NE=4, SE=5, NW=6, SW=7;
	MPI_Status status;
	if(IdNodeNE.size()>=1)
	{
		MultiBlock_->Send(&f[0]->f[5][IdNodeNE[0]],1,4,tag_d_tr);
		MultiBlock_->Send(&f[1]->f[5][IdNodeNE[0]],1,4,tag_d_tr);
	}
	if(IdNodeSW.size()>=1)
	{
		MultiBlock_->Recv(&f[0]->f[5][IdNodeSW[0]],1,7,tag_d_tr,status);
		MultiBlock_->Recv(&f[1]->f[5][IdNodeSW[0]],1,7,tag_d_tr,status);
	}
	if(IdNodeNW.size()>=1)
	{
		MultiBlock_->Send(&f[0]->f[2][IdNodeNW[0]],1,0,tag_y_t);
		MultiBlock_->Send(&f[1]->f[2][IdNodeNW[0]],1,0,tag_y_t);
	}
	if(IdNodeSW.size()>=1)
	{
		MultiBlock_->Recv(&f[0]->f[2][IdNodeSW[0]],1,2,tag_y_t,status);
		MultiBlock_->Recv(&f[1]->f[2][IdNodeSW[0]],1,2,tag_y_t,status);
	}
	if(IdNodeNW.size()>=1)
	{
		MultiBlock_->Send(&f[0]->f[6][IdNodeNW[0]],1,0,tag_y_t);
		MultiBlock_->Send(&f[1]->f[6][IdNodeNW[0]],1,0,tag_y_t);
	}
	if(IdNodeSW.size()>=1)
	{
		MultiBlock_->Recv(&f[0]->f[6][IdNodeSW[0]],1,2,tag_y_t,status);
		MultiBlock_->Recv(&f[1]->f[6][IdNodeSW[0]],1,2,tag_y_t,status);
	}
	if(IdNodeSE.size()>=1)
	{
		MultiBlock_->Send(&f[0]->f[1][IdNodeSE[0]],1,1,tag_x_r);
		MultiBlock_->Send(&f[1]->f[1][IdNodeSE[0]],1,1,tag_x_r);
	}
	if(IdNodeSW.size()>=1)
	{
		MultiBlock_->Recv(&f[0]->f[1][IdNodeSW[0]],1,3,tag_x_r,status);
		MultiBlock_->Recv(&f[1]->f[1][IdNodeSW[0]],1,3,tag_x_r,status);
	}
	if(IdNodeSE.size()>=1)
	{
		MultiBlock_->Send(&f[0]->f[8][IdNodeSE[0]],1,1,tag_x_r);
		MultiBlock_->Send(&f[1]->f[8][IdNodeSE[0]],1,1,tag_x_r);
	}
	if(IdNodeSW.size()>=1)
	{
		MultiBlock_->Recv(&f[0]->f[8][IdNodeSW[0]],1,3,tag_x_r,status);
		MultiBlock_->Recv(&f[1]->f[8][IdNodeSW[0]],1,3,tag_x_r,status);
	}
/*	if(IdGNodeNE.size()>=1)
		MultiBlock_->Send(&f[0]->f[5][IdGNodeNE[0]],1,4,tag_d_tr);
	if(IdRNodeSW.size()>=1)
		MultiBlock_->Recv(&f[0]->f[5][IdRNodeSW[0]],1,7,tag_d_tr,status);
	if(IdGNodeNW.size()>=1)
		MultiBlock_->Send(&f[0]->f[2][IdGNodeNW[0]],1,0,tag_y_t);
	if(IdRNodeSW.size()>=1)
		MultiBlock_->Recv(&f[0]->f[2][IdRNodeSW[0]],1,2,tag_y_t,status);
	if(IdGNodeNW.size()>=1)
		MultiBlock_->Send(&f[0]->f[6][IdGNodeNW[0]],1,0,tag_y_t);
	if(IdRNodeSW.size()>=1)
		MultiBlock_->Recv(&f[0]->f[6][IdRNodeSW[0]],1,2,tag_y_t,status);
	if(IdGNodeSE.size()>=1)
		MultiBlock_->Send(&f[0]->f[1][IdGNodeSE[0]],1,1,tag_x_r);
	if(IdRNodeSW.size()>=1)
		MultiBlock_->Recv(&f[0]->f[1][IdRNodeSW[0]],1,3,tag_x_r,status);
	if(IdGNodeSE.size()>=1)
		MultiBlock_->Send(&f[0]->f[8][IdGNodeSE[0]],1,1,tag_x_r);
	if(IdRNodeSW.size()>=1)
		MultiBlock_->Recv(&f[0]->f[8][IdRNodeSW[0]],1,3,tag_x_r,status);*/
}
void D2Q9TwoPhases::CornerNodesSyncToGhost(){
	 int tag_x_l=1;
	 int tag_y_b=2;
	 int tag_d_bl=3;
	//N=0,E=1,S=2,W=3,NE=4, SE=5, NW=6, SW=7;
	MPI_Status status;

	if(IdRNodeSE.size()>=1)
	{
		MultiBlock_->Send(&f[0]->f[8][IdRNodeSE[0]],1,5,tag_x_l);
		MultiBlock_->Send(&f[1]->f[8][IdRNodeSE[0]],1,5,tag_x_l);
	}
	if(IdGNodeNW.size()>=1)
	{
		MultiBlock_->Recv(&f[0]->f[8][IdGNodeNW[0]],1,6,tag_x_l,status);
		MultiBlock_->Recv(&f[1]->f[8][IdGNodeNW[0]],1,6,tag_x_l,status);
	}
	if(IdRNodeSW.size()>=1)
	{
		MultiBlock_->Send(&f[0]->f[7][IdRNodeSW[0]],1,7,tag_x_l);
		MultiBlock_->Send(&f[1]->f[7][IdRNodeSW[0]],1,7,tag_x_l);
	}
	if(IdGNodeNE.size()>=1)
	{
		MultiBlock_->Recv(&f[0]->f[7][IdGNodeNE[0]],1,4,tag_x_l,status);
		MultiBlock_->Recv(&f[1]->f[7][IdGNodeNE[0]],1,4,tag_x_l,status);
	}
	if(IdRNodeNE.size()>=1)
	{
		MultiBlock_->Send(&f[0]->f[5][IdRNodeNE[0]],1,4,tag_x_l);
		MultiBlock_->Send(&f[1]->f[5][IdRNodeNE[0]],1,4,tag_x_l);
	}
	if(IdGNodeSW.size()>=1)
	{
		MultiBlock_->Recv(&f[0]->f[5][IdGNodeSW[0]],1,7,tag_x_l,status);
		MultiBlock_->Recv(&f[1]->f[5][IdGNodeSW[0]],1,7,tag_x_l,status);
	}
	if(IdRNodeNW.size()>=1)
	{
		MultiBlock_->Send(&f[0]->f[6][IdRNodeNW[0]],1,6,tag_d_bl);
		MultiBlock_->Send(&f[1]->f[6][IdRNodeNW[0]],1,6,tag_d_bl);
	}
	if(IdGNodeSE.size()>=1)
	{
		MultiBlock_->Recv(&f[0]->f[6][IdGNodeSE[0]],1,5,tag_d_bl,status);
		MultiBlock_->Recv(&f[1]->f[6][IdGNodeSE[0]],1,5,tag_d_bl,status);
	}
}
void D2Q9TwoPhases::SyncFromGhost(){
	GhostNodesSyncFromGhost();
	CornerNodesSyncFromGhost();
}
void D2Q9TwoPhases::SyncToGhost(){
	GhostNodesSyncToGhost();
	CornerNodesSyncToGhost();
}

void D2Q9TwoPhases::SyncMacroVarToGhost(){



	for (unsigned int i=0;i<IdRNodeW.size();i++)
	{
		buf_MacroSend[0][1][i]=Rho[IdRNodeW[i]];
//		buf_MacroSend[1][1][i]=U[0][IdRNodeW[i]];
//		buf_MacroSend[2][1][i]=U[0][IdRNodeW[i]];
	}
	for (unsigned int i=0;i<IdRNodeS.size();i++)
	{
		buf_MacroSend[0][2][i]=Rho[IdRNodeS[i]];
//		buf_MacroSend[1][2][i]=U[0][IdRNodeS[i]];
//		buf_MacroSend[2][2][i]=U[1][IdRNodeS[i]];
	}
	for (unsigned int i=0;i<IdRNodeE.size();i++)
	{
		buf_MacroSend[0][0][i]=Rho[IdRNodeE[i]];
	}
	for (unsigned int i=0;i<IdRNodeN.size();i++)
	{
		buf_MacroSend[0][3][i]=Rho[IdRNodeN[i]];
	}
//		MultiBlock_->Communication(buf_MacroSend[0],buf_MacroRecv[0],size_MacroBuf);
	int tag_x_r=1;
	int tag_x_l=2;
	int tag_y_t=3;
	int tag_y_b=4;
	MPI_Status status;
	if(IdRNodeE.size()>=1)
	{
		MultiBlock_->Send(&buf_MacroSend[0][0][0],IdRNodeE.size(),1,tag_x_r);
		//MultiBlock_->Send(&buf_MacroSend[1][0][0],IdRNodeE.size(),1,tag_x_r);
		//MultiBlock_->Send(&buf_MacroSend[2][0][0],IdRNodeE.size(),1,tag_x_r);
	}
	if(IdGNodeW.size()>=1)
	{
		MultiBlock_->Recv(&buf_MacroRecv[0][0][0],IdGNodeW.size(),3,tag_x_r,status);
	//	MultiBlock_->Recv(&buf_MacroRecv[1][0][0],IdGNodeW.size(),3,tag_x_r,status);
	//	MultiBlock_->Recv(&buf_MacroRecv[2][0][0],IdGNodeW.size(),3,tag_x_r,status);
	}

	if(IdRNodeW.size()>=1)
	{
		MultiBlock_->Send(&buf_MacroSend[0][1][0],IdRNodeW.size(),3,tag_x_l);
	//	MultiBlock_->Send(&buf_MacroSend[1][1][0],IdRNodeW.size(),3,tag_x_l);
	//	MultiBlock_->Send(&buf_MacroSend[2][1][0],IdRNodeW.size(),3,tag_x_l);
	}
	if(IdGNodeE.size()>=1)
	{
		MultiBlock_->Recv(&buf_MacroRecv[0][1][0],IdGNodeE.size(),1,tag_x_l,status);
	//	MultiBlock_->Recv(&buf_MacroRecv[1][1][0],IdGNodeE.size(),1,tag_x_l,status);
	//	MultiBlock_->Recv(&buf_MacroRecv[2][1][0],IdGNodeE.size(),1,tag_x_l,status);
	}
		if(IdRNodeN.size()>=1)
	{
		MultiBlock_->Send(&buf_MacroSend[0][3][0],IdRNodeN.size(),0,tag_y_t);
	//	MultiBlock_->Send(&buf_MacroSend[1][3][0],IdRNodeN.size(),0,tag_y_t);
	//	MultiBlock_->Send(&buf_MacroSend[2][3][0],IdRNodeN.size(),0,tag_y_t);
	}
	if(IdGNodeS.size()>=1)
	{
		MultiBlock_->Recv(&buf_MacroRecv[0][3][0],IdGNodeS.size(),2,tag_y_t,status);
	//	MultiBlock_->Recv(&buf_MacroRecv[1][3][0],IdGNodeS.size(),2,tag_y_t,status);
	//	MultiBlock_->Recv(&buf_MacroRecv[2][3][0],IdGNodeS.size(),2,tag_y_t,status);
	}
	if(IdRNodeS.size()>=1)
	{
		MultiBlock_->Send(&buf_MacroSend[0][2][0],IdRNodeS.size(),2,tag_y_b);
	//	MultiBlock_->Send(&buf_MacroSend[1][2][0],IdRNodeS.size(),2,tag_y_b);
	//	MultiBlock_->Send(&buf_MacroSend[2][2][0],IdRNodeS.size(),2,tag_y_b);
	}
	if(IdGNodeN.size()>=1)
	{
		MultiBlock_->Recv(&buf_MacroRecv[0][2][0],IdGNodeN.size(),0,tag_y_b,status);
	//	MultiBlock_->Recv(&buf_MacroRecv[1][2][0],IdGNodeN.size(),0,tag_y_b,status);
	//	MultiBlock_->Recv(&buf_MacroRecv[2][2][0],IdGNodeN.size(),0,tag_y_b,status);
	}
	for (unsigned int i=0;i<IdGNodeE.size();i++)
	{
		Rho[IdGNodeE[i]]=buf_MacroRecv[0][1][i];
	}
	for (unsigned int i=0;i<IdGNodeN.size();i++)
	{
		Rho[IdGNodeN[i]]=buf_MacroRecv[0][2][i];
	}
	for (unsigned int i=0;i<IdGNodeW.size();i++)
	{
		Rho[IdGNodeW[i]]=buf_MacroRecv[0][0][i];
	}

	for (unsigned int i=0;i<IdGNodeS.size();i++)
	{
		Rho[IdGNodeS[i]]=buf_MacroRecv[0][3][i];
	}
/*	 int tag_x_l=1;
	 int tag_y_b=2;
	 int tag_d_bl=3;
	//N=0,E=1,S=2,W=3,NE=4, SE=5, NW=6, SW=7;
	MPI_Status status;
	if(IdRNodeSW.size()>=1)
	{
		MultiBlock_->Send(&Rho[IdRNodeSW[0]],1,7,tag_d_bl);
	}
	if(IdGNodeNE.size()>=1)
	{
		MultiBlock_->Recv(&Rho[IdGNodeNE[0]],1,4,tag_d_bl,status);
		MultiBlock_->Send(&Rho[IdRNodeNW[0]],1,4,tag_d_bl);
	}
	if(IdRNodeSW.size()>=1)
	{
		MultiBlock_->Recv(&Rho[IdGNodeSW[0]],1,7,tag_d_bl,status);
		MultiBlock_->Send(&Rho[IdRNodeSW[0]],1,2,tag_y_b);
	}
	if(IdGNodeNW.size()>=1)
	{
		MultiBlock_->Recv(&Rho[IdGNodeNW[0]],1,0,tag_y_b,status);
	}
	if(IdRNodeSW.size()>=1)
	{
		MultiBlock_->Send(&Rho[IdRNodeSW[0]],1,3,tag_x_l);
	}
	if(IdGNodeSE.size()>=1)
	{
		MultiBlock_->Recv(&Rho[IdGNodeSE[0]],1,1,tag_x_l,status);
	}*/
}

double D2Q9TwoPhases::Convert_Alpha_To_Rho(double alpha)
{

	return PtrParameters->Get_Rho_1()*alpha + PtrParameters->Get_Rho_2()*(1-alpha);
}
double D2Q9TwoPhases::Convert_Rho_To_Alpha(double Rho)
{

	return (Rho-PtrParameters->Get_Rho_2()) /(PtrParameters->Get_Rho_1() - PtrParameters->Get_Rho_2());
}
