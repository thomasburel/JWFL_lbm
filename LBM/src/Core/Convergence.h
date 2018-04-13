/*
 * ============================================================================
 * Convergence.h
 *
 *  Created on: 30 Sep 2016
 *      Author: Thomas Burel
 *      Version     :
 *      Copyright   : Your copyright notice
 *      Description : LBM Code Developing by JWFL at University of Strathclyde
 * ============================================================================
 */

#ifndef SRC_CORE_CONVERGENCE_H_
#define SRC_CORE_CONVERGENCE_H_
#include "../Mesh/MultiBlock.h"
#include "Dictionary.h"
#include "../Mesh/SingleBlock/Patch/PatchBc.h"
#include "Viscosity.h"
#include "../Algorithm/Tools/Gradients.h"
#include "../User/UserConvergence.h"
enum TypeConverge{none,GlobalConvergence,FieldConvergence};
enum TypeConvergeScalar{SacalarGlobal, ScalarField};
enum TypeConvergeVector{VectorGlobal,VectorField};
class Convergence: public UserConvergence {
public:
	Convergence();
	virtual ~Convergence();
	void Set_Convergence();
	void Calcul_Error(int &Time);
	double Get_Error(){return Error;};

private:
	void Calcul_LpError();
	void CalCul_L1Error();

	void Calcul_Error_ScalarField();
	void Calcul_Error_VectorField();
	void NoCalcul_Error(){};
	double Calcul_Error_ScalarFieldInOneProc();
	void Set_ConvergencePatchBc();

	void Set_MarkFluidNodesError();
	bool IsBoundaryError(int idx);
	bool IsWrongSideDomainError(int idx);
	bool IsNormalLimitDomainError(int idx);
	bool IsInDomainError(int idx);
	bool IsInsideDomainError(int idx);
	bool IsGlobalCornerASpecialWallError(int idx);
	bool IsConvexCornerNormalOutsideError(int idx);


	void Calcul_PorousMediaConvergence(int &Time);
	double Calcul_ProductionRate(int &Time);
	double Calcul_Permeability_SinglePhase(int &Time);
	double Calcul_DarcyPermeability_SinglePhase(int &Time);
	double Calcul_Permeability_TwoPhases(int &Time);
	double Calcul_DarcyPermeability_TwoPhases(int &Time);
	void Calcul_localDeltaP();


	void Sum_ScalarNodeArray(std::vector<int>& NodeArray,std::vector<int>& NodeArraySpecialWall,std::vector<int>& NodeArrayGlobalCorner,double *&Var1,double &sum);
	void Sum_VectorNodeArray(std::vector<int>& NodeArray,std::vector<int>& NodeArraySpecialWall,std::vector<int>& NodeArrayGlobalCorner,double *&Var1,double *&Var2,double &sum1,double &sum2);
	void Sum_ScalarNodeArray(std::vector<int>& NodeArray,std::vector<int>& NodeArraySpecialWall,std::vector<int>& NodeArrayGlobalCorner,double *&Var1,double *&weight,bool phase1,double &sum);
	void Sum_VectorNodeArray(std::vector<int>& NodeArray,std::vector<int>& NodeArraySpecialWall,std::vector<int>& NodeArrayGlobalCorner,double *&Var1,double *&Var2,double *&weight,bool phase1,double &sum1,double &sum2);
	void Sum_WeightNodeArray(std::vector<int>& NodeArray,std::vector<int>& NodeArraySpecialWall,std::vector<int>& NodeArrayGlobalCorner,double *&weight,bool phase1,double &sum);
	void Sum_ViscosityNodeArray(std::vector<int>& NodeArray,std::vector<int>& NodeArraySpecialWall,std::vector<int>& NodeArrayGlobalCorner,double *&Rho,double &sum);
	void Sum_ViscosityNodeArray(std::vector<int>& NodeArray,std::vector<int>& NodeArraySpecialWall,std::vector<int>& NodeArrayGlobalCorner,double *&Rho,double *&weight,bool phase1,double &sum);


	double Porosity();
	void Set_MarkFluidNodes();
	//void SumFluidVolume(double & sumfluidvolume);
	void SumSolidVolume(double & sumsolidvolume);
	bool IsBoundary(int idx);
	bool IsWrongSideDomain(int idx);
	bool IsNormalLimitDomain(int idx);
	bool IsInDomain(int idx);
	bool IsInsideDomain(int idx);
	bool IsGlobalCornerASpecialWall(int idx);
	bool IsConvexCornerNormalOutside(int idx);

	//convert RhoN to alpha or 1-alpha for phase 2
	double Convert_RhoNToAlpha(double const &RhoN,bool const &phase1){if(phase1)return (RhoN+1.0)*0.5;else return (1.0-RhoN)*0.5;};
	double Convert_AlphaToRhoN(double const &alpha,bool const &phase1){if(phase1)return 2.0*alpha-1.0;else return 2.0*(1.0-alpha)-1.0;};
	double NoConvert_RhoNToAlpha(double const &alpha){return alpha;};
	double Get_Weigth(double const &RhoN,bool const &phase1){ return Convert_RhoNToAlpha(RhoN,phase1);};
	double Get_WeightInv(double const &Weight,bool const &phase1){ return Convert_AlphaToRhoN(Weight,phase1);};

	double RhoToP(double const Rho){return Rho/3.0;};
	double Permeability(double const &deltaP,double const &u,double const mu){return u*mu/deltaP;};

protected:
	MultiBlock *PtrMultiBlockConv;
	Dictionary *PtrDicConv;
	Parameters *PtrParmConv;
	PatchBc *PtrPatchBcConv;
	NodeArrays* PtrNodeArraysConv;
	Viscosity *PtrViscosityConv;

private:
	double Sum_Current;
	double Error;///< Save error and it used inside the sum if needed
	double Error_sum, Error_avg;
	double Error_tmp;
	std::vector<int> MarkFluidNodeError_V1,MarkFluidNodeError_V075,MarkFluidNodeError_V05,MarkFluidNodeError_V025;///< Mark fluid node for error over the domain

	double *RhoError,**UError;
	int nbnodesErrorLp;

	std::vector<int> InletPatchId,OutletPatchId;
	double *RhoNProductionRate,**UPerm,*RhoPerm,*RhoNPerm;
	std::vector<int> *NodeId,*NodeIdSpeWall,*NodeIdGloCorner;
	std::vector<int> MarkFluidNode_V1,MarkFluidNode_V075,MarkFluidNode_V05,MarkFluidNode_V025;///< Mark fluid node for sum over the porous domain

	double avg;
	double porosity,LuToPhy2;
	double *PressureConv;
	double **DeltaP;///< Local Pressure difference
	//double *UMag;///< Calculate Velocity magnitude
	Gradients DensityGradient;
	double fluidVolumesum;double solidVolumesum;
	double LengthMedium,SectionMedium;
	unsigned int minXMedia,maxXMedia,minYMedia,maxYMedia,minZMedia,maxZMedia;

	int NbNodes;

	double *Scalar_CurrentTime;///< Pointer on the scalar of the Current Time
	double **Vector_CurrentTime;///< Pointer on the scalar of the Current Time
	//Save previous time step
	double *Scalar_last;///< Save previous time step for scalars
	double **Vector_last;///< Save previous time step for vectors

	// Pointers on function
	///Simplify notation for pointer on a member function of D2Q9ColourFluid class for Colour Gradient methods
		typedef double(Convergence::*PtrFctPermeability)(int & time);
	//Define name for pointers on functions
		PtrFctPermeability PtrPermeabilityDarcy;///< Permeability Darcy
		PtrFctPermeability PtrPermeability;///< Permeability of the medium
};

#endif /* SRC_CORE_CONVERGENCE_H_ */
