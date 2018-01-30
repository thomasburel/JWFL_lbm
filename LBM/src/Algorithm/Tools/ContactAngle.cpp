/*
 * ============================================================================
 * ContactAngle.cpp
 *
 *  Created on: 10 Jan 2017
 *      Author: Thomas Burel
 *      Version     :
 *      Copyright   : Your copyright notice
 *      Description : LBM Code Developing by JWFL at University of Strathclyde
 * ============================================================================
 */

#include "ContactAngle.h"

ContactAngle::ContactAngle() {
	PtrNodeCa=0;
	PtrParamCa=0;
	PtrOppositeCa=0;
	Switch2D=0;Model=0;
	I_tmp=0;
	n1=0;n2=0;
	D1=0;D2=0;r=0;rMinus1=0;
	costeta=0;sinteta=0;
	teta=0;///< Contact angle
	D_tmp=0;epsilon=1.e-10;
	PtrOppositeCa=0;
	MapWallId=0;
	MapCornerConcaveId=0;
	MapCornerConvexId=0;
	InvSqrt2=1.0/std::sqrt(2.0);
	ExtrapolNormalInSolid=0;
	ImposeNormalOnWall=0;
}

ContactAngle::~ContactAngle() {
	// TODO Auto-generated destructor stub
}
void ContactAngle::AllocateTeta(NodeArrays2D *PtrNode,Parameters *PtrParam,double * &tetaIn){
	if(teta!=0)
		delete teta;

	switch(PtrParam->Get_ContactAngleType())
		{
		case ContactAngleEnum::NoTeta:
			teta=0;
			break;
		case ContactAngleEnum::UniformTeta:
			teta=new double[1];
			teta[0]=0;
			break;
		case ContactAngleEnum::UserTeta:
			teta=new double[PtrNode->CornerConcave.size()+PtrNode->CornerConvex.size()+PtrNode->NodeWall.size()+PtrNode->NodeSpecialWall.size()];
	/*		int nodeWallIdx=0;
			for (int j=0;j<PtrNode->CornerConcave.size();j++)
			{
				teta[nodeWallIdx]=tetaIn[PtrNode->NodeCorner[PtrNode->CornerConcave[j]].Get_index()];
				nodeWallIdx++;
			}
			for (int j=0;j<PtrNode->NodeWall.size();j++)
			{
				teta[nodeWallIdx]=tetaIn[PtrNode->NodeWall[j].Get_index()];
				nodeWallIdx++;
			}
			for (int j=0;j<PtrNode->CornerConvex.size();j++)
			{
				teta[nodeWallIdx]=tetaIn[PtrNode->NodeCorner[PtrNode->CornerConvex[j]].Get_index()];
				nodeWallIdx++;
			}
			for (int j=0;j<PtrNode->NodeSpecialWall.size();j++)
			{
				teta[nodeWallIdx]=tetaIn[PtrNode->NodeSpecialWall[j].Get_index()];
				nodeWallIdx++;
			}*/

			for(unsigned int i=0;i<PtrNode->CornerConcave.size()+PtrNode->CornerConvex.size()+PtrNode->NodeWall.size()+PtrNode->NodeSpecialWall.size();i++)
				teta[i]=0;

			break;
		default:
			std::cerr<<" Contact angle type not found."<<std::endl;
			break;
		}
	tetaIn=teta;
}
void ContactAngle::InitContactAngle(NodeArrays2D *PtrNode,Parameters *PtrParam,unsigned int *PtrOpposite, double * &tetaIn) {
	PtrNodeCa=PtrNode;
	PtrParamCa=PtrParam;
	PtrOppositeCa=new unsigned int[9];
	for(int i=0;i<9;i++)
		PtrOppositeCa[i]=PtrOpposite[i];
	epsilon=1.e-10;
	int nodeWallIdx=0;
	//teta=tetaIn;
//switch How to select the theoretical contact angle
	switch(PtrParamCa->Get_SwitchSelectTeta())
	{
		case ContactAngleEnum::Binary:
			Switch2D=&ContactAngle::BinarySwitchContactAngle2D;
			break;
		case ContactAngleEnum::Linear:
			Switch2D=&ContactAngle::LinearSwitchContactAngle2D;
			break;
		default:
			std::cerr<<"Method for switching of the contact angle depending of the side of the droplet is not found."<<std::endl;
	}
// Switch between models for contact angles
	if(PtrParam->Get_ContactAngleType()==ContactAngleEnum::NoTeta)
	{
		Model=&ContactAngle::ApplyNoTetaMethod;
		if(PtrParamCa->Get_ColourExtrapolNoramlInterface())
		{

			if(PtrParamCa->Get_Dimension()==SolverEnum::D2)
			{
				ExtrapolNormalInSolid=&ContactAngle::Extrapolate_NormalInSolid2D;
				if(PtrParamCa->Get_NormalExtrapolType()==ContactAngleEnum::TailorExtrapol)
					ExtraPolNormalCa.initExtrapolation(2,9,ModelEnum::TailorExtrapol);
				else
					ExtraPolNormalCa.initExtrapolation(2,9,ModelEnum::WeightDistanceExtrapol);
			}
			else
				std::cerr<<" 3D Extrapolation not yet implemented."<<std::endl;
		}
		else
			ExtrapolNormalInSolid=&ContactAngle::NoExtrapolate_NormalInSolid;
	}
	else
	{
		switch(PtrParamCa->Get_ContactAngleModel())
		{
		case ContactAngleEnum::Standard:
			Model=&ContactAngle::ApplyStandardMethod;
			if(PtrParamCa->Get_Dimension()==SolverEnum::D2)
			{
				if(PtrParam->Get_ContactAngleType()==ContactAngleEnum::UniformTeta)
					ImposeNormalOnWall=&ContactAngle::Impose_ContactAngleOnWall2DFixTeta;
				else
					ImposeNormalOnWall=&ContactAngle::Impose_ContactAngleOnWall2DNonCstTeta;
			}
			if(PtrParamCa->Get_ColourExtrapolNoramlInterface())
			{
				ExtrapolNormalInSolid=&ContactAngle::Extrapolate_NormalInSolid2D;
				if(PtrParamCa->Get_Dimension()==SolverEnum::D2)
				{
					if(PtrParamCa->Get_NormalExtrapolType()==ContactAngleEnum::TailorExtrapol)
						ExtraPolNormalCa.initExtrapolation(2,9,ModelEnum::TailorExtrapol);
					else
						ExtraPolNormalCa.initExtrapolation(2,9,ModelEnum::WeightDistanceExtrapol);
				}
				else
					std::cerr<<" 3D Extrapolation not yet implemented."<<std::endl;

			}
			else
				ExtrapolNormalInSolid=&ContactAngle::NoExtrapolate_NormalInSolid;
			break;
		case ContactAngleEnum::Interpol:
			Model=&ContactAngle::ApplyInterpolMethod;
			if(PtrParamCa->Get_Dimension()==SolverEnum::D2)
			{
				if(PtrParam->Get_ContactAngleType()==ContactAngleEnum::UniformTeta)
					ImposeNormalOnWall=&ContactAngle::Impose_ContactAngleInSolidAndInterpol2DFixTeta;
				else
					ImposeNormalOnWall=&ContactAngle::Impose_ContactAngleInSolidAndInterpol2DNonCstTeta;
			}
			else
				std::cerr<<" 3D Interpolation method is not yet implemented."<<std::endl;
			ExtrapolNormalInSolid=&ContactAngle::Extrapolate_NormalInSolid2D;
			if(PtrParamCa->Get_Dimension()==SolverEnum::D2)
			{
				if(PtrParamCa->Get_NormalExtrapolType()==ContactAngleEnum::TailorExtrapol)
					ExtraPolNormalCa.initExtrapolation(2,9,ModelEnum::TailorExtrapol);
				else
					ExtraPolNormalCa.initExtrapolation(2,9,ModelEnum::WeightDistanceExtrapol);
			}
			else
				std::cerr<<" 3D Extrapolation is not yet implemented."<<std::endl;
			PtrParamCa->Set_ColourExtrapolNoramlInterface(true);
			if(PtrParamCa->Get_Dimension()==SolverEnum::D2)
			{
				if(PtrParamCa->Get_NormalInterpolType()==ContactAngleEnum::LinearLeastSquareInterpol)
					InterPolNormalCa.initInterpolation(2,9,ModelEnum::LinearLeastSquareInterpol,PtrOppositeCa,PtrNodeCa,PtrParamCa);
				else
					InterPolNormalCa.initInterpolation(2,9,ModelEnum::LinearInterpol,PtrOppositeCa,PtrNodeCa,PtrParamCa);
			}
			else
				std::cerr<<" 3D Interpolation is not yet implemented."<<std::endl;
			break;
		default:
			std::cerr<<"Method for switching between models is not found."<<std::endl;
		}
		if(PtrParamCa->Get_Dimension()==SolverEnum::D2)
		{
			if(PtrParam->Get_ContactAngleType()==ContactAngleEnum::UniformTeta)
			{
				nodeWallIdx=0;
				n1=new double**[1];n2=new double**[1];
				n1[0]=new double*[9];n2[0]=new double*[9];
				for(int i=0;i<9;i++)
				{
					n1[0][i]=new double[2];n2[0][i]=new double[2];
				}
				teta[nodeWallIdx]=PtrParam->Get_ContactAngle();
				Set_TwoChoiceOfContactAngle2D(nodeWallIdx);
			}
			else
			{
				int nbwalls=PtrNodeCa->CornerConcave.size()+PtrNodeCa->NodeWall.size()+PtrNodeCa->CornerConvex.size()+PtrNodeCa->NodeSpecialWall.size();
				n1=new double**[nbwalls];n2=new double**[nbwalls];
				for(int j=0;j<nbwalls;j++)
				{
					n1[j]=new double*[9];n2[j]=new double*[9];
					for(int i=0;i<9;i++)
					{
						n1[j][i]=new double[2];n2[j][i]=new double[2];
					}
				}

				for(nodeWallIdx=0;nodeWallIdx<nbwalls;nodeWallIdx++)
				{
					Set_TwoChoiceOfContactAngle2D(nodeWallIdx);
				}
				/*for (int j=0;j<PtrNodeCa->CornerConcave.size();j++){

				}*/
			}
		}
		else
			std::cerr<<" 3D contact angle is not yet implemented."<<std::endl;
	}
}
void ContactAngle::ApplyContactAngle2D(double **&Normal){
	(this->*Model)(Normal);
}
void ContactAngle::ApplyNoTetaMethod(double **&Normal){
	(this->*ExtrapolNormalInSolid)(Normal);
}
void ContactAngle::ApplyStandardMethod(double **&Normal){
	(this->*ExtrapolNormalInSolid)(Normal);
	(this->*ImposeNormalOnWall)(Normal);
}
void ContactAngle::ApplyInterpolMethod(double **&Normal){
	(this->*ExtrapolNormalInSolid)(Normal);
	(this->*ImposeNormalOnWall)(Normal);
}
void ContactAngle::Set_TwoChoiceOfContactAngle2D(int &nodeWallIdx){
	costeta=std::cos(teta[nodeWallIdx]);
	sinteta=std::sin(teta[nodeWallIdx]);
	if (std::abs(costeta)<epsilon) {
		costeta=0.0;
		sinteta=copysign(1.0,sinteta);
	}
	if (std::abs(costeta)>1.0-epsilon) {
		costeta=copysign(1.0,costeta);
		sinteta=0.0;
	}

	n1[nodeWallIdx][0][0]=0.0;
	n2[nodeWallIdx][0][0]=n1[nodeWallIdx][0][0];
	n1[nodeWallIdx][0][1]=0.0;
	n2[nodeWallIdx][0][1]=n1[nodeWallIdx][0][1];
	n1[nodeWallIdx][1][0]=costeta;
	n2[nodeWallIdx][1][0]=n1[nodeWallIdx][1][0];
	n1[nodeWallIdx][1][1]=sinteta;
	n2[nodeWallIdx][1][1]=-n1[nodeWallIdx][1][1];
	n1[nodeWallIdx][2][0]=sinteta;
	n2[nodeWallIdx][2][0]=-n1[nodeWallIdx][2][0];
	n1[nodeWallIdx][2][1]=costeta;
	n2[nodeWallIdx][2][1]=n1[nodeWallIdx][2][1];
	n1[nodeWallIdx][3][0]=-costeta;
	n2[nodeWallIdx][3][0]=n1[nodeWallIdx][3][0];
	n1[nodeWallIdx][3][1]=sinteta;
	n2[nodeWallIdx][3][1]=-n1[nodeWallIdx][3][1];
	n1[nodeWallIdx][4][0]=sinteta;
	n2[nodeWallIdx][4][0]=-n1[nodeWallIdx][4][0];
	n1[nodeWallIdx][4][1]=-costeta;
	n2[nodeWallIdx][4][1]=n1[nodeWallIdx][4][1];

	n1[nodeWallIdx][5][0]=InvSqrt2*(costeta-sinteta);
	n1[nodeWallIdx][5][1]=InvSqrt2*(costeta+sinteta);
	n2[nodeWallIdx][5][0]=InvSqrt2*(costeta+sinteta);
	n2[nodeWallIdx][5][1]=InvSqrt2*(costeta-sinteta);

	n1[nodeWallIdx][6][0]=InvSqrt2*(-costeta+sinteta);
	n1[nodeWallIdx][6][1]=InvSqrt2*(costeta+sinteta);
	n2[nodeWallIdx][6][0]=InvSqrt2*(-costeta-sinteta);
	n2[nodeWallIdx][6][1]=InvSqrt2*(costeta-sinteta);

	n1[nodeWallIdx][7][0]=InvSqrt2*(-costeta-sinteta);
	n1[nodeWallIdx][7][1]=InvSqrt2*(-costeta+sinteta);
	n2[nodeWallIdx][7][0]=InvSqrt2*(-costeta+sinteta);
	n2[nodeWallIdx][7][1]=InvSqrt2*(-costeta-sinteta);

	n1[nodeWallIdx][8][0]=InvSqrt2*(costeta+sinteta);
	n1[nodeWallIdx][8][1]=InvSqrt2*(-costeta+sinteta);
	n2[nodeWallIdx][8][0]=InvSqrt2*(costeta-sinteta);
	n2[nodeWallIdx][8][1]=InvSqrt2*(-costeta-sinteta);
}

void ContactAngle::Extrapolate_NormalInSolid2D(double **&Normal){

	for (unsigned int j=0;j<PtrNodeCa->CornerConcave.size();j++)
	{
		ExtraPolNormalCa.ExtrapolationCornerConcaveToSolid(Normal[0],PtrNodeCa->NodeCorner[PtrNodeCa->CornerConcave[j]].Get_connect(),PtrNodeCa->NodeCorner[PtrNodeCa->CornerConcave[j]].Get_BcNormal());
		ExtraPolNormalCa.ExtrapolationCornerConcaveToSolid(Normal[1],PtrNodeCa->NodeCorner[PtrNodeCa->CornerConcave[j]].Get_connect(),PtrNodeCa->NodeCorner[PtrNodeCa->CornerConcave[j]].Get_BcNormal());
		//Normalise(Normal[0],Normal[1],PtrNodeCa->NodeCorner[PtrNodeCa->CornerConcave[j]].Get_connect()[PtrOppositeCa[PtrNodeCa->NodeCorner[PtrNodeCa->CornerConcave[j]].Get_BcNormal()]]);
		NormaliseConcaveCornerExtraPol(Normal[0],Normal[1],PtrNodeCa->NodeCorner[PtrNodeCa->CornerConcave[j]].Get_connect(),PtrNodeCa->NodeCorner[PtrNodeCa->CornerConcave[j]].Get_BcNormal());
	}
	for (unsigned int j=0;j<PtrNodeCa->NodeWall.size();j++)
	{
		ExtraPolNormalCa.ExtrapolationWallToSolid(Normal[0],PtrNodeCa->NodeWall[j].Get_connect(),PtrNodeCa->NodeWall[j].Get_BcNormal());
		ExtraPolNormalCa.ExtrapolationWallToSolid(Normal[1],PtrNodeCa->NodeWall[j].Get_connect(),PtrNodeCa->NodeWall[j].Get_BcNormal());
		Normalise(Normal[0],Normal[1],PtrNodeCa->NodeWall[j].Get_connect()[PtrOppositeCa[PtrNodeCa->NodeWall[j].Get_BcNormal()]]);
	}
	for (unsigned int j=0;j<PtrNodeCa->CornerConvex.size();j++)
	{
		ExtraPolNormalCa.ExtrapolationCornerConvexToSolid(Normal[0],PtrNodeCa->NodeCorner[PtrNodeCa->CornerConvex[j]].Get_connect(),PtrNodeCa->NodeCorner[PtrNodeCa->CornerConvex[j]].Get_BcNormal());
		ExtraPolNormalCa.ExtrapolationCornerConvexToSolid(Normal[1],PtrNodeCa->NodeCorner[PtrNodeCa->CornerConvex[j]].Get_connect(),PtrNodeCa->NodeCorner[PtrNodeCa->CornerConvex[j]].Get_BcNormal());
		Normalise(Normal[0],Normal[1],PtrNodeCa->NodeCorner[PtrNodeCa->CornerConvex[j]].Get_connect()[PtrOppositeCa[PtrNodeCa->NodeCorner[PtrNodeCa->CornerConvex[j]].Get_BcNormal()]]);
	}
	for (unsigned int j=0;j<PtrNodeCa->NodeSpecialWall.size();j++)
	{
		ExtraPolNormalCa.ExtrapolationCornerConcaveToSolid(Normal[0],PtrNodeCa->NodeSpecialWall[j].Get_connect(),PtrNodeCa->NodeSpecialWall[j].Get_BcNormal());
		ExtraPolNormalCa.ExtrapolationCornerConcaveToSolid(Normal[1],PtrNodeCa->NodeSpecialWall[j].Get_connect(),PtrNodeCa->NodeSpecialWall[j].Get_BcNormal());
		NormaliseConcaveCornerExtraPol(Normal[0],Normal[1],PtrNodeCa->NodeSpecialWall[j].Get_connect(),PtrNodeCa->NodeSpecialWall[j].Get_BcNormal());
	}
}
void ContactAngle::Impose_ContactAngleInSolidAndInterpol2DFixTeta(double **&Normal){
	int nodeWallIdx=0;
	//define the contact angle in the solid
	for (unsigned int j=0;j<PtrNodeCa->CornerConcave.size();j++)
	{
		ContactAngleConcaveCornerInSolid2D(nodeWallIdx,Normal,PtrNodeCa->NodeCorner[PtrNodeCa->CornerConcave[j]].Get_connect(),PtrNodeCa->NodeCorner[PtrNodeCa->CornerConcave[j]].Get_BcNormal());
	}
	for (unsigned int j=0;j<PtrNodeCa->NodeWall.size();j++)
	{
		ContactAngleWallInSolid2D(nodeWallIdx,Normal,PtrNodeCa->NodeWall[j].Get_connect(),PtrNodeCa->NodeWall[j].Get_BcNormal());
	}
	for (unsigned int j=0;j<PtrNodeCa->CornerConvex.size();j++)
	{
		ContactAngleConvexCornerInSolid2D(nodeWallIdx,Normal,PtrNodeCa->NodeCorner[PtrNodeCa->CornerConvex[j]].Get_connect(),PtrNodeCa->NodeCorner[PtrNodeCa->CornerConvex[j]].Get_BcNormal());
	}
	for (unsigned int j=0;j<PtrNodeCa->NodeSpecialWall.size();j++)
	{
		ContactAngleConcaveCornerInSolid2D(nodeWallIdx,Normal,PtrNodeCa->NodeSpecialWall[j].Get_connect(),PtrNodeCa->NodeSpecialWall[j].Get_BcNormal());
	}
	//interpolate the contact angle on the boundary from in the solid and in the interior fluid
	for (unsigned int j=0;j<PtrNodeCa->CornerConcave.size();j++)
	{
		InterPolNormalCa.InterpolationOnCornerConcave(Normal[0],Normal[1],PtrNodeCa->NodeCorner[PtrNodeCa->CornerConcave[j]].Get_connect(),PtrNodeCa->NodeCorner[PtrNodeCa->CornerConcave[j]].Get_BcNormal());
		Normalise(Normal[0],Normal[1],PtrNodeCa->NodeCorner[PtrNodeCa->CornerConcave[j]].Get_index());
	}
	for (unsigned int j=0;j<PtrNodeCa->NodeWall.size();j++)
	{
		InterPolNormalCa.InterpolationOnWall(Normal[0],Normal[1],PtrNodeCa->NodeWall[j].Get_connect(),PtrNodeCa->NodeWall[j].Get_BcNormal());
		Normalise(Normal[0],Normal[1],PtrNodeCa->NodeWall[j].Get_index());
	}
	for (unsigned int j=0;j<PtrNodeCa->CornerConvex.size();j++)
	{
		InterPolNormalCa.InterpolationOnCornerConvex(Normal[0],Normal[1],PtrNodeCa->NodeCorner[PtrNodeCa->CornerConvex[j]].Get_connect(),PtrNodeCa->NodeCorner[PtrNodeCa->CornerConvex[j]].Get_BcNormal());
		Normalise(Normal[0],Normal[1],PtrNodeCa->NodeCorner[PtrNodeCa->CornerConvex[j]].Get_index());
	}
	for (unsigned int j=0;j<PtrNodeCa->NodeSpecialWall.size();j++)
	{
	//	ExtraPolNormalCa.ExtrapolationOnCornerConcave(Normal[0],PtrNodeCa->NodeSpecialWall[j].Get_connect(),PtrNodeCa->NodeSpecialWall[j].Get_BcNormal());
	//	ExtraPolNormalCa.ExtrapolationOnCornerConcave(Normal[1],PtrNodeCa->NodeSpecialWall[j].Get_connect(),PtrNodeCa->NodeSpecialWall[j].Get_BcNormal());
		//InterPolNormalCa.InterpolationOnCornerConcave(Normal[0],Normal[1],PtrNodeCa->NodeSpecialWall[j].Get_connect(),PtrNodeCa->NodeSpecialWall[j].Get_BcNormal());
	//	Normalise(Normal[0],Normal[1],PtrNodeCa->NodeSpecialWall[j].Get_index());
/*		if(PtrNodeCa->TypeOfNode[PtrNodeCa->NodeSpecialWall[j].Get_connect()[1]]==Wall)
		{
			Normal[0][PtrNodeCa->NodeSpecialWall[j].Get_index()]=Normal[0][PtrNodeCa->NodeSpecialWall[j].Get_connect()[1]];
			Normal[1][PtrNodeCa->NodeSpecialWall[j].Get_index()]=Normal[1][PtrNodeCa->NodeSpecialWall[j].Get_connect()[1]];
		}
		else if(PtrNodeCa->TypeOfNode[PtrNodeCa->NodeSpecialWall[j].Get_connect()[2]]==Wall)
		{
			Normal[0][PtrNodeCa->NodeSpecialWall[j].Get_index()]=Normal[0][PtrNodeCa->NodeSpecialWall[j].Get_connect()[2]];
			Normal[1][PtrNodeCa->NodeSpecialWall[j].Get_index()]=Normal[1][PtrNodeCa->NodeSpecialWall[j].Get_connect()[2]];
		}
		else if(PtrNodeCa->TypeOfNode[PtrNodeCa->NodeSpecialWall[j].Get_connect()[3]]==Wall)
		{
			Normal[0][PtrNodeCa->NodeSpecialWall[j].Get_index()]=Normal[0][PtrNodeCa->NodeSpecialWall[j].Get_connect()[3]];
			Normal[1][PtrNodeCa->NodeSpecialWall[j].Get_index()]=Normal[1][PtrNodeCa->NodeSpecialWall[j].Get_connect()[3]];
		}
		else if(PtrNodeCa->TypeOfNode[PtrNodeCa->NodeSpecialWall[j].Get_connect()[4]]==Wall)
		{
			Normal[0][PtrNodeCa->NodeSpecialWall[j].Get_index()]=Normal[0][PtrNodeCa->NodeSpecialWall[j].Get_connect()[4]];
			Normal[1][PtrNodeCa->NodeSpecialWall[j].Get_index()]=Normal[1][PtrNodeCa->NodeSpecialWall[j].Get_connect()[4]];
		}
		*/
//		Normalise(Normal[0],Normal[1],PtrNodeCa->NodeSpecialWall[j].Get_index());
	}
}
void ContactAngle::Impose_ContactAngleInSolidAndInterpol2DNonCstTeta(double **&Normal){
	int nodeWallIdx=0;
	//define the contact angle in the solid
	for (unsigned int j=0;j<PtrNodeCa->CornerConcave.size();j++)
	{
		ContactAngleConcaveCornerInSolid2D(nodeWallIdx,Normal,PtrNodeCa->NodeCorner[PtrNodeCa->CornerConcave[j]].Get_connect(),PtrNodeCa->NodeCorner[PtrNodeCa->CornerConcave[j]].Get_BcNormal());
		nodeWallIdx++;
	}
	for (unsigned int j=0;j<PtrNodeCa->NodeWall.size();j++)
	{
		ContactAngleWallInSolid2D(nodeWallIdx,Normal,PtrNodeCa->NodeWall[j].Get_connect(),PtrNodeCa->NodeWall[j].Get_BcNormal());
		nodeWallIdx++;
	}
	for (unsigned int j=0;j<PtrNodeCa->CornerConvex.size();j++)
	{
		ContactAngleConvexCornerInSolid2D(nodeWallIdx,Normal,PtrNodeCa->NodeCorner[PtrNodeCa->CornerConvex[j]].Get_connect(),PtrNodeCa->NodeCorner[PtrNodeCa->CornerConvex[j]].Get_BcNormal());
		nodeWallIdx++;
	}
	for (unsigned int j=0;j<PtrNodeCa->NodeSpecialWall.size();j++)
	{
		ContactAngleConcaveCornerInSolid2D(nodeWallIdx,Normal,PtrNodeCa->NodeSpecialWall[j].Get_connect(),PtrNodeCa->NodeSpecialWall[j].Get_BcNormal());
		nodeWallIdx++;
	}
	//interpolate the contact angle on the boundary from in the solid and in the interior fluid
	for (unsigned int j=0;j<PtrNodeCa->CornerConcave.size();j++)
	{
		InterPolNormalCa.InterpolationOnCornerConcave(Normal[0],Normal[1],PtrNodeCa->NodeCorner[PtrNodeCa->CornerConcave[j]].Get_connect(),PtrNodeCa->NodeCorner[PtrNodeCa->CornerConcave[j]].Get_BcNormal());
		Normalise(Normal[0],Normal[1],PtrNodeCa->NodeCorner[PtrNodeCa->CornerConcave[j]].Get_index());
	}
	for (unsigned int j=0;j<PtrNodeCa->NodeWall.size();j++)
	{
		InterPolNormalCa.InterpolationOnWall(Normal[0],Normal[1],PtrNodeCa->NodeWall[j].Get_connect(),PtrNodeCa->NodeWall[j].Get_BcNormal());
		Normalise(Normal[0],Normal[1],PtrNodeCa->NodeWall[j].Get_index());
	}
	for (unsigned int j=0;j<PtrNodeCa->CornerConvex.size();j++)
	{
		InterPolNormalCa.InterpolationOnCornerConvex(Normal[0],Normal[1],PtrNodeCa->NodeCorner[PtrNodeCa->CornerConvex[j]].Get_connect(),PtrNodeCa->NodeCorner[PtrNodeCa->CornerConvex[j]].Get_BcNormal());
		Normalise(Normal[0],Normal[1],PtrNodeCa->NodeCorner[PtrNodeCa->CornerConvex[j]].Get_index());
	}
	for (unsigned int j=0;j<PtrNodeCa->NodeSpecialWall.size();j++)
	{
		InterPolNormalCa.InterpolationOnCornerConcave(Normal[0],Normal[1],PtrNodeCa->NodeSpecialWall[j].Get_connect(),PtrNodeCa->NodeSpecialWall[j].Get_BcNormal());
		Normalise(Normal[0],Normal[1],PtrNodeCa->NodeSpecialWall[j].Get_index());
	}
}
void ContactAngle::Select_ContactAngle2D(int &nodeWallIdx, int & Bcnormal,double & Nx, double & Ny){
	//Calcul the distance between the estimate normal vector and the two choices of impose normal vector
	D1=std::sqrt((Nx-n1[nodeWallIdx][Bcnormal][0])*(Nx-n1[nodeWallIdx][Bcnormal][0])+(Ny-n1[nodeWallIdx][Bcnormal][1])*(Ny-n1[nodeWallIdx][Bcnormal][1]));
	D2=std::sqrt((Nx-n2[nodeWallIdx][Bcnormal][0])*(Nx-n2[nodeWallIdx][Bcnormal][0])+(Ny-n2[nodeWallIdx][Bcnormal][1])*(Ny-n2[nodeWallIdx][Bcnormal][1]));
	//calculate the ratio
	if(D1+D2>0)
		r=D1/(D1+D2);
	else
		r=0.5;
	//Select binary switch or linear switch (quicker but less precise)
	(this->*Switch2D)(nodeWallIdx,r,Bcnormal,Nx,Ny);

}
void ContactAngle::LinearSwitchContactAngle2D(int &nodeWallIdx, double & rin,int & Bcnormal,double & Nx, double & Ny){
	rMinus1=1.0-rin;
	Nx=rMinus1*n1[nodeWallIdx][Bcnormal][0]+r*n2[nodeWallIdx][Bcnormal][0];
	Ny=rMinus1*n1[nodeWallIdx][Bcnormal][1]+r*n2[nodeWallIdx][Bcnormal][1];
}
void ContactAngle::BinarySwitchContactAngle2D(int &nodeWallIdx, double & rin,int & Bcnormal,double & Nx, double & Ny){
	if(rin<0.5)
	{
		Nx=n1[nodeWallIdx][Bcnormal][0];
		Ny=n1[nodeWallIdx][Bcnormal][1];
	}
	else
	{
		Nx=n2[nodeWallIdx][Bcnormal][0];
		Ny=n2[nodeWallIdx][Bcnormal][1];
	}
}
void ContactAngle::Impose_ContactAngleOnWall2DFixTeta(double **&Normal){
	int nodeWallIdx=0;
	//define the contact angle on walls
	for (unsigned int j=0;j<PtrNodeCa->CornerConcave.size();j++)
	{
		ContactAngleonWalls2D(nodeWallIdx,Normal,PtrNodeCa->NodeCorner[PtrNodeCa->CornerConcave[j]].Get_connect(),PtrNodeCa->NodeCorner[PtrNodeCa->CornerConcave[j]].Get_BcNormal());
	}
	for (unsigned int j=0;j<PtrNodeCa->NodeWall.size();j++)
	{
		ContactAngleonWalls2D(nodeWallIdx,Normal,PtrNodeCa->NodeWall[j].Get_connect(),PtrNodeCa->NodeWall[j].Get_BcNormal());
	}
	for (unsigned int j=0;j<PtrNodeCa->CornerConvex.size();j++)
	{
		ContactAngleonWalls2D(nodeWallIdx,Normal,PtrNodeCa->NodeCorner[PtrNodeCa->CornerConvex[j]].Get_connect(),PtrNodeCa->NodeCorner[PtrNodeCa->CornerConvex[j]].Get_BcNormal());
	}
	for (unsigned int j=0;j<PtrNodeCa->NodeSpecialWall.size();j++)
	{
		ContactAngleonWalls2D(nodeWallIdx,Normal,PtrNodeCa->NodeSpecialWall[j].Get_connect(),PtrNodeCa->NodeSpecialWall[j].Get_BcNormal());
	}
}
void ContactAngle::Impose_ContactAngleOnWall2DNonCstTeta(double **&Normal){
	int nodeWallIdx=0;
	//define the contact angle on walls
	for (unsigned int j=0;j<PtrNodeCa->CornerConcave.size();j++)
	{
		ContactAngleonWalls2D(nodeWallIdx,Normal,PtrNodeCa->NodeCorner[PtrNodeCa->CornerConcave[j]].Get_connect(),PtrNodeCa->NodeCorner[PtrNodeCa->CornerConcave[j]].Get_BcNormal());
		nodeWallIdx++;
	}
	for (unsigned int j=0;j<PtrNodeCa->NodeWall.size();j++)
	{
		ContactAngleonWalls2D(nodeWallIdx,Normal,PtrNodeCa->NodeWall[j].Get_connect(),PtrNodeCa->NodeWall[j].Get_BcNormal());
		nodeWallIdx++;
	}
	for (unsigned int j=0;j<PtrNodeCa->CornerConvex.size();j++)
	{
		ContactAngleonWalls2D(nodeWallIdx,Normal,PtrNodeCa->NodeCorner[PtrNodeCa->CornerConvex[j]].Get_connect(),PtrNodeCa->NodeCorner[PtrNodeCa->CornerConvex[j]].Get_BcNormal());
		nodeWallIdx++;
	}
	for (unsigned int j=0;j<PtrNodeCa->NodeSpecialWall.size();j++)
	{
		ContactAngleonWalls2D(nodeWallIdx,Normal,PtrNodeCa->NodeSpecialWall[j].Get_connect(),PtrNodeCa->NodeSpecialWall[j].Get_BcNormal());
		nodeWallIdx++;
	}
}
void ContactAngle::ContactAngleonWalls2D(int &nodeWallIdx, double **&Normal,int* connect, int & Bcnormal){
	//Select contact angle for the node in the opposite of the normal direction
	Select_ContactAngle2D(nodeWallIdx,Bcnormal,Normal[0][connect[0]],Normal[1][connect[0]]);
}
void ContactAngle::ContactAngleConcaveCornerInSolid2D(int &nodeWallIdx, double **&Normal,int* connect, int & Bcnormal){
	//Select contact angle for the node in the opposite of the normal direction
	Select_ContactAngle2D(nodeWallIdx,Bcnormal,Normal[0][connect[PtrOppositeCa[Bcnormal]]],Normal[1][connect[PtrOppositeCa[Bcnormal]]]);
	//Select contact angle for the other solid nodes
	switch(Bcnormal)
	{
	case 5:
		Select_ContactAngle2D(nodeWallIdx,IntRef(1),Normal[0][connect[3]],Normal[1][connect[3]]);
		Select_ContactAngle2D(nodeWallIdx,IntRef(2),Normal[0][connect[4]],Normal[1][connect[4]]);

		break;
	case 6:
		Select_ContactAngle2D(nodeWallIdx,IntRef(3),Normal[0][connect[1]],Normal[1][connect[1]]);
		Select_ContactAngle2D(nodeWallIdx,IntRef(2),Normal[0][connect[4]],Normal[1][connect[4]]);
		break;
	case 7:
		Select_ContactAngle2D(nodeWallIdx,IntRef(3),Normal[0][connect[1]],Normal[1][connect[1]]);
		Select_ContactAngle2D(nodeWallIdx,IntRef(4),Normal[0][connect[2]],Normal[1][connect[2]]);
		break;
	case 8:
		Select_ContactAngle2D(nodeWallIdx,IntRef(1),Normal[0][connect[3]],Normal[1][connect[3]]);
		Select_ContactAngle2D(nodeWallIdx,IntRef(4),Normal[0][connect[2]],Normal[1][connect[2]]);
		break;
	}
}
void ContactAngle::ContactAngleConvexCornerInSolid2D(int &nodeWallIdx, double **&Normal,int* connect, int & Bcnormal){
	Select_ContactAngle2D(nodeWallIdx,Bcnormal,Normal[0][connect[PtrOppositeCa[Bcnormal]]],Normal[1][connect[PtrOppositeCa[Bcnormal]]]);
}
void ContactAngle::ContactAngleWallInSolid2D(int &nodeWallIdx, double **&Normal,int* Connect, int & Bcnormal){
	Select_ContactAngle2D(nodeWallIdx,Bcnormal,Normal[0][Connect[PtrOppositeCa[Bcnormal]]],Normal[1][Connect[PtrOppositeCa[Bcnormal]]]);
}
