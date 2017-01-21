/*
 * D2Q9TwoPhases.h
 *
 *  Created on: 9 Jun 2015
 *      Author: thomas
 */

#ifndef ALGORITHM_LOWORDER_D2Q9TwoPhases_H_
#define ALGORITHM_LOWORDER_D2Q9TwoPhases_H_
#include "D2Q9CommonVar.h"
#include "../../Core/Parameters.h"
#include "../../Core/Solver.h"
#include "../../Core/GlobalDef.h"
#include <iostream>
#include <cmath>
#include "Boundaries/D2Q9Bc.h"
#include "../Tools/ContactAngle.h"
class D2Q9TwoPhases: public SolverTwoPhasesLowOrder2D, public D2Q9Bc, protected D2Q9CommonVar{
public:
	D2Q9TwoPhases();
	virtual ~D2Q9TwoPhases();


protected:

	void InitD2Q9TwoPhases(MultiBlock* MultiBlock__,ParallelManager* parallel__,WriterManager* Writer__, Parameters* Parameters__,InitLBM& ini);
	void init(InitLBM& ini);
	void InitAllDomain(InitLBM& ini);
	void InitDomainBc(InitLBM& ini);
	void InitWall(InitLBM& ini);
	void InitInterior(InitLBM& ini);
	void StreamD2Q9();


protected:
	// Multiphase member functions
	double Convert_Alpha_To_Rho(double alpha);
	double Convert_Rho_To_Alpha(double Rho);

	//Streaming by type of node
	void SelectStream(int & nodenumber, unsigned int& direction);
	void TmptoDistri(unsigned int& direction, int& IdDistri);
	void InteriorStream(int & nodenumber, unsigned int& direction);
	void GhostStream(int & nodenumber, unsigned int& direction);
	void WallStream(int & nodenumber, unsigned int& direction);
	void CornerStream(int & nodenumber, unsigned int& direction);
	void PeriodicStream(int & nodenumber, unsigned int& direction);
	void VelocityStream(int & nodenumber, unsigned int& direction);
	void PressureStream(int & nodenumber, unsigned int& direction);
	void SolidStream(int & nodenumber, unsigned int& direction);

	//Set streaming for Boundaries conditions (Stream only what we are able to stream)
	void Set_BcType(int & nodenumber);
	void Set_GhostType(int & nodenumber);
	void Set_WallType(int & nodenumber);
	void Set_CornerType(int & nodenumber);
	void Set_VelocityType(int & nodenumber);
	void Set_PressureType(int & nodenumber);
	void Set_BcType();
	void Set_GhostType(NodeGhost2D& Node);
	void Set_WallType(NodeWall2D& Node);
	void Set_CornerType(NodeCorner2D& Node);
	void Set_VelocityType(NodeVelocity2D& Node);
	void Set_PressureType(NodePressure2D& Node);
	void Set_SymmetryType(NodeSymmetry2D& NodeIn);
	void Set_PeriodicType(NodePeriodic2D& NodeIn);
	//Set Orientation of Boundaries

	void StreamingOrientation(int & nodenumber, bool GhostStreaming[9]);
	void StreamingOrientation(NodeGhost2D& Node, bool GhostStreaming[9]);
	void StreamingOrientation(NodeWall2D& Node, bool WallStreaming[9]);
	void StreamingOrientation(NodeCorner2D& Node, bool CornerStreaming[9]);
	void StreamingOrientation(NodeVelocity2D& Node, bool VelocityStreaming[9]);
	void StreamingOrientation(NodePressure2D& Node, bool PressureStreaming[9]);
	void StreamingOrientation(NodeSymmetry2D& Node, bool SymmetryStreaming[9]);
	void StreamingOrientation(NodePeriodic2D& Node, bool PeriodicStreaming[9]);

	//Apply boundary conditions

	void ApplyHeZou_U(NodeVelocity2D& Node, int distID, double &U, double &V);
	void ApplyHeZou_U(NodeCorner2D& Node,int normal, int distID, double &U, double &V);//Use for Global Corners
	void ApplyHeZou_P(NodePressure2D& Node, int distID, double Mass, double &U, double &V);
	void ApplyHeZou_P(NodeCorner2D& Node,int normal, int distID, double Mass, double &U, double &V);//Use for Global Corners
	void ApplyBounceBackSymmetry(NodeWall2D& Node);//Use for special wall (Symmetry-Wall connection)
	void ApplyDiffuseWall(NodeWall2D& Node);
	void ApplyDiffuseWallSymmetry(NodeWall2D& Node);
	void ApplyDiffuseWall(NodeCorner2D& Node);
	void ApplySymmetryOnNode(NodeSymmetry2D& Node);
	double Cal_RhoCorner(NodeCorner2D& Node);

	//Communications
	void IniComVariables();
	void SyncFromGhost();
	void SyncToGhost();
	void GhostNodesSyncFromGhost();
	void GhostNodesSyncToGhost();
	void CornerNodesSyncFromGhost();
	void CornerNodesSyncToGhost();
	void SyncMacroVarToGhost();

	//Tools
inline	void Normalise(double* Var_x,double* Var_y, int nodenumber){
	// Normalise
	D_tmp=sqrt(Var_x[nodenumber]*Var_x[nodenumber]+Var_y[nodenumber]*Var_y[nodenumber]);
	if(D_tmp>0)
		{Var_x[nodenumber]/=D_tmp;Var_y[nodenumber]/=D_tmp;}
	else
		{Var_x[nodenumber]=0.0; Var_y[nodenumber]=0.0;}

}


protected:
//Multiphase variables
	double Rho1, Rho2;


//Common variables
	double* tmpDistribution;// variable to move the temporary distribution function to the distribution function without copy (memory moving)
	double* f_tmp;// variable to temporary copy the local distribution (for one node)

	unsigned int Opposite[9]; //opposite direction in the distribution function

	int** DiagConnect; //not use
	unsigned int tmpreturn; //return for connection between nodes
	int intTmpReturn;
	double doubleTmpReturn;

	int Nd_variables_sync;//number of variable has to be synchronise
	double ***buf_send, ***buf_recv; //buffers to send and receive
	int *size_buf; // size of buffers

	int Nd_MacroVariables_sync;//number of variable has to be synchronise
	std::vector<double*> SyncVar;
	double ***buf_MacroSend, ***buf_MacroRecv; //buffers to send and receive
	int *size_MacroBuf; // size of buffers



private:
	friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
       ar & boost::serialization::base_object<SolverTwoPhasesLowOrder2D>(*this);
  //     ar & boost::serialization::base_object<Boundaries>(*this);
       ar &  Nd_variables_sync;
    }

};

/*inline void D2Q9TwoPhases::CollideD2Q9() {
	for (int i=0;i<nbvelo;i++)
	{
		for (int j=0;j<nbnodes_real;j++)
			f[0]->f[i][j]=f[0]->f[i][j]-InvTau*(f[0]->f[i][j]-CollideLowOrder::EquiDistriFunct2D(Rho[j], U[0][j], U[1][j], &Ei[i][0], omega[i]));
	}
}*/
/// Calculate \f$\rho\f$ and \f$\vec{U}\f$ in the local domain
/*inline void D2Q9TwoPhases::UpdateMacroVariables(){
	for (int i=0;i<nbnodes_real;i++)
	{
		U[0][i]=0;
		U[1][i]=0;
		Rho[i]=0;
		for (int k=0; k<nbvelo;k++)
		{
			//Rho[i]=1;
			Rho[i]+=f[0]->f[k][i];
			for (int j=0;j<2;j++)
			{
				U[j][i]+=f[0]->f[k][i]*Ei[k][j];
			}
		}
		U[0][i]=U[0][i]/Rho[i];
		U[1][i]=U[1][i]/Rho[i];
	}
}*/
//int D2Q9TwoPhases::Get_BcNormal2(T & nodeIn);
//template <typename T>


#endif /* ALGORITHM_LOWORDER_D2Q9TwoPhases_H_ */
