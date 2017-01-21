/*
 * InitLBM.cpp
 *
 *  Created on: 5 May 2015
 *      Author: thomas
 */

#include "InitLBM.h"


InitLBM::InitLBM() {
	PtrParameters=0;

}
void InitLBM::IniMPI(ParallelManager* parrallel_,int *argc, char ***argv, bool verbous_) {

	/// Initialise MPI
	parrallel_->init(argc,argv,verbous_);

}


InitLBM::~InitLBM() {
	// TODO Auto-generated destructor stub
}

void InitLBM::Set_Parameters(Parameters *Parameters) {
	PtrParameters=Parameters;
}
void InitLBM::IniDomainSinglePhase(int rank,Node2D & Node,int elem, int nodenumber, double* pos,double& Rho, double* U){
	double alpha=0;
	IniDomainTwoPhases(rank,Node,elem,nodenumber, pos,Rho,U,alpha);
}
void InitLBM::IniDomainTwoPhases(int rank,Node2D & Node,int elem, int nodenumber, double* pos,double& Rho, double* U,double & alpha){
	switch(Node.get_NodeType())
	{
	case Interior:
		UserIc(*PtrParameters,elem, nodenumber, pos,Rho, U,alpha);
		break;
	case Wall:
		UserBc(*PtrParameters,elem, nodenumber, pos,Rho, U,alpha);
		U[0]=0;
		U[1]=0;
		break;
	case SpecialWall:
		UserBc(*PtrParameters,elem, nodenumber, pos,Rho, U,alpha);
		U[0]=0;
		U[1]=0;
		break;
	case Corner:
		UserBc(*PtrParameters,elem, nodenumber, pos,Rho, U,alpha);
		Node.Set_UDef(U[0],U[1]);
		Node.Set_RhoDef(Rho);
		break;
	case GlobalCorner:
		UserBc(*PtrParameters,elem, nodenumber, pos,Rho, U,alpha);
		break;
	case Periodic:
		UserBc(*PtrParameters,elem, nodenumber, pos,Rho, U,alpha);
		Node.Set_UDef(U[0],U[1]);
		Node.Set_RhoDef(Rho);
		break;
	case Velocity:
		UserBc(*PtrParameters,elem, nodenumber, pos,Rho, U,alpha);
		Node.Set_UDef(U[0],U[1]);
		break;
	case Pressure:
		UserBc(*PtrParameters,elem, nodenumber, pos,Rho, U,alpha);
		Node.Set_UDef(U[0],U[1]);
		Node.Set_RhoDef(Rho);
		break;
	case Solid:
		U[0]=0;
		U[1]=0;
		Rho=0;
		alpha=1;
		break;
	case Symmetry:
		UserBc(*PtrParameters,elem, nodenumber, pos,Rho, U,alpha);
		Node.Set_UDef(U[0],U[1]);
		Node.Set_RhoDef(Rho);
		break;
	case Ghost:
		UserIc(*PtrParameters,elem, nodenumber, pos,Rho, U,alpha);
		break;
	default:
		std::cout<<"Processor ID: "<<rank<<" Node number: "<<nodenumber<<" Node type not find in initialisation "<<std::endl;
		std::cout<<"Processor ID: "<<rank<<" Node number: "<<nodenumber<<" Type of node is : "<<Node.get_NodeType()<<std::endl;
		break;
	}
}
void InitLBM::IniContactAngle(int rank,Node2D & Node,int elem, int nodenumber, double* pos,double & teta){
	if(PtrParameters->Get_ContactAngleType()!=ContactAngleEnum::NoTeta)
	{
		switch(Node.get_NodeType())
		{
		case Wall:
			if(PtrParameters->Get_ContactAngleType()==ContactAngleEnum::FixTeta)
				teta=PtrParameters->Get_ContactAngle();
			else
				Set_ContactAngle(*PtrParameters,elem, nodenumber, pos,teta);
			break;
		case SpecialWall:
			if(PtrParameters->Get_ContactAngleType()==ContactAngleEnum::FixTeta)
				teta=PtrParameters->Get_ContactAngle();
			else
				Set_ContactAngle(*PtrParameters,elem, nodenumber, pos,teta);
			break;
		case Corner:
			if(PtrParameters->Get_ContactAngleType()==ContactAngleEnum::FixTeta)
				teta=PtrParameters->Get_ContactAngle();
			else
				Set_ContactAngle(*PtrParameters,elem, nodenumber, pos,teta);
			break;
		default:
			std::cout<<"Processor ID: "<<rank<<" Node number: "<<nodenumber<<" Node type not find in initialisation "<<std::endl;
			std::cout<<"Processor ID: "<<rank<<" Node number: "<<nodenumber<<" Type of node is : "<<Node.get_NodeType()<<std::endl;
			break;
		}
	}
}

