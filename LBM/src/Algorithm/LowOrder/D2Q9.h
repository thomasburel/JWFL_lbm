/*
 * D2Q9.h
 *
 *  Created on: 9 Jun 2015
 *      Author: thomas
 */

#ifndef ALGORITHM_LOWORDER_D2Q9_H_
#define ALGORITHM_LOWORDER_D2Q9_H_
#include "../../Core/Parameters.h"
#include "../../Core/Solver.h"
#include "../../Core/GlobalDef.h"
#include "Boundaries.h"
#include <iostream>
#include <cmath>

class D2Q9: public SolverSinglePhaseLowOrder2D, public Boundaries {
public:
	D2Q9();
	D2Q9(MultiBlock* MultiBlock__,ParallelManager* parallel__,WriterManager* Writer__, Parameters* Parameters__,InitLBM& ini);
	virtual ~D2Q9();
	virtual void init(InitLBM& ini);
	virtual void run();

	void check_fi();

private:
	template <typename T>
	void Collide(T node);
	void CollideD2Q9();
	void CollideD2Q9_node();
	void StreamD2Q9();
	void UpdateMacroVariables();
	void UpdateMacroVariables_node();
	void MacroVariables(int& idx);
	bool checknode(int x, int y);
	void checkcomm();

private:
	//Streaming by type of node
	void SelectStream(int & nodenumber, unsigned int& direction);
	void TmptoDistri(unsigned int& direction);
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

	//Set Orientation of Boundaries

	void StreamingOrientation(int & nodenumber, bool GhostStreaming[9]);
	void StreamingOrientation(NodeGhost2D& Node, bool GhostStreaming[9]);
	void StreamingOrientation(NodeWall2D& Node, bool WallStreaming[9]);
	void StreamingOrientation(NodeCorner2D& Node, bool CornerStreaming[9]);
	void StreamingOrientation(NodeVelocity2D& Node, bool VelocityStreaming[9]);
	void StreamingOrientation(NodePressure2D& Node, bool PressureStreaming[9]);
	void StreamingOrientation(NodeSymmetry2D& Node, bool SymmetryStreaming[9]);

	//Apply boundary conditions
	void ApplyBc();
	void ApplyHeZou_U(int & nodenumber);
	void ApplyHeZou_P(int & nodenumber);
	void ApplyBounceBack(int & nodenumber);
	void ApplyCorner(int & nodenumber);
	void ApplyHeZou_U(NodeVelocity2D& Node);
	void ApplyHeZou_U(NodeCorner2D& Node,int normal);//Use for Global Corners
	void ApplyHeZou_P(NodePressure2D& Node);
	void ApplyHeZou_P(NodeCorner2D& Node,int normal);//Use for Global Corners
	void ApplyBounceBack(NodeWall2D& Node);
	void ApplyBounceBackSymmetry(NodeWall2D& Node);//Use for special wall
	void ApplyBounceBack(NodeCorner2D& Node);
	void ApplyDiffuseWall(NodeWall2D& Node);
	void ApplyDiffuseWallSymmetry(NodeWall2D& Node);
	void ApplyDiffuseWall(NodeCorner2D& Node);
	void ApplyCorner(NodeCorner2D& Node);
	void ApplyGlobalCorner(NodeCorner2D& Node);
	void ApplySymmetryOnNode(NodeSymmetry2D& Node);
	void ApplySymmetryPressureOnNode(NodeSymmetry2D& Node);
	unsigned int& Connect (int &NodeNumber,unsigned int& direction);
	//void Set_Connect(int &NodeNumber,unsigned int& direction);
	double Cal_RhoCorner(int &normalBc, int &nodenumber);
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

private:
	double* tmp;// variable to copy tmp to distribution function
	double Ei[9][2]; //Velocity in the distribution function
	double omega[9];//Weight in the distribution function
	unsigned int Opposite[9]; //opposite direction in the distribution function

	int** DiagConnect; //not use
	unsigned int tmpreturn; //return for connection between nodes
	int intTmpReturn;
	double doubleTmpReturn;

	int Nd_variables_sync;//number of variable has to be synchronise
	double ***buf_send, ***buf_recv; //buffers to send and receive
	int *size_buf; // size of buffers

	int Nd_MacroVariables_sync;//number of variable has to be synchronise
	double ***buf_MacroSend, ***buf_MacroRecv; //buffers to send and receive
	int *size_MacroBuf; // size of buffers

private:
	friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
       ar & boost::serialization::base_object<SolverSinglePhaseLowOrder2D>(*this);
       ar & boost::serialization::base_object<Boundaries>(*this);
       ar & Ei & omega & Opposite & Nd_variables_sync;
    }

};

inline void D2Q9::CollideD2Q9() {
	for (int i=0;i<nbvelo;i++)
	{
		for (int j=0;j<nbnodes_real;j++)
			f->f[i][j]=f->f[i][j]-InvTau*(f->f[i][j]-CollideLowOrder::EquiDistriFunct2D(Rho[j], U[0][j], U[1][j], &Ei[i][0], omega[i]));
	}
}
/// Calculate \f$\rho\f$ and \f$\vec{U}\f$ in the local domain
inline void D2Q9::UpdateMacroVariables(){
	for (int i=0;i<nbnodes_real;i++)
	{
		U[0][i]=0;
		U[1][i]=0;
		Rho[i]=0;
		for (int k=0; k<nbvelo;k++)
		{
			//Rho[i]=1;
			Rho[i]+=f->f[k][i];
			for (int j=0;j<2;j++)
			{
				U[j][i]+=f->f[k][i]*Ei[k][j];
			}
		}
		U[0][i]=U[0][i]/Rho[i];
		U[1][i]=U[1][i]/Rho[i];
	}
}
//int D2Q9::Get_BcNormal2(T & nodeIn);
//template <typename T>


#endif /* ALGORITHM_LOWORDER_D2Q9_H_ */
