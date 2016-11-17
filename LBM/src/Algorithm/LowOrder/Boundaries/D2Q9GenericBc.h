/*
 * ============================================================================
 * D2Q9GenericBc.h
 *
 *  Created on: 6 Sep 2016
 *      Author: Thomas Burel
 *      Version     :
 *      Copyright   : Your copyright notice
 *      Description : LBM Code Developing by JWFL at University of Strathclyde
 * ============================================================================
 */

#ifndef SRC_ALGORITHM_LOWORDER_BOUNDARIES_D2Q9GENERICBC_H_
#define SRC_ALGORITHM_LOWORDER_BOUNDARIES_D2Q9GENERICBC_H_
#include "D2Q9BoundariesList.h"
#include "../../../Core/Dictionary.h"

class D2Q9GenericBc {
public:
	D2Q9GenericBc();
	virtual ~D2Q9GenericBc();
	void SetBcObjects( Parameters *Param);
	///Apply Pressure as one 0-moment variable (single phase = 1 density)
	void ApplyPressure(int const &BcNormal,int const *Connect, double const Rho_def, DistriFunct * & f_in);
	///Apply Pressure as  several 0-moment variables
	void ApplyPressure(int const &BcNormal,int const *Connect, double const Rho_def, DistriFunct * & f_in,double * & Rho, double * &U, double * &V);
	///Apply Velocity as one 0-moment variable (single phase = 1 density)
	void ApplyVelocity(int const &BcNormal,int const *Connect, double const *UDef, DistriFunct * & f_in);
	///Apply Velocity as  several 0-moment variables
	void ApplyVelocity(int const &BcNormal,int const *Connect, double const *UDef, DistriFunct * & f_in,double * & Rho, double * &U, double * &V);
	///Apply Pressure/Velocity Special Wall as one 0-moment variable (single phase = 1 density)
	void ApplyPreVelSpecialWall(NodeWall2D& Node, DistriFunct * & f_in,double const & RhoDef,double const & UDef,double const & VDef,double * & Rho, double * &U, double * &V);
	///Apply Pressure/Velocity Special Wall as  several 0-moment variables
//	void ApplyPreVelSpecialWall(NodeWall2D& Node, DistriFunct * & f_in,double const & RhoDef,double const & UDef,double const & VDef,double * & Rho, double * &U, double * &V);
	///Apply Corner as one 0-moment variable (single phase = 1 density)
	void ApplyCorner(NodeCorner2D& Node, DistriFunct * & f_in,double const & RhoDef,double const & UDef,double const & VDef);
	///Apply Corner as  several 0-moment variables
	void ApplyCorner(NodeCorner2D& Node, DistriFunct * & f_in,double const & RhoDef,double const & UDef,double const & VDef,double * & Rho, double * &U, double * &V);
	///Apply CornerWall as one 0-moment variable (single phase = 1 density)
	void ApplyCornerWall(NodeCorner2D& Node, DistriFunct * & f_in);
	///Apply CornerWall as  several 0-moment variables
	void ApplyCornerWall(NodeCorner2D& Node, DistriFunct * & f_in,double * & Rho, double * &U, double * &V);
	///Apply Corner for SpecialWall as one 0-moment variable (single phase = 1 density)
	void ApplyCornerSpecialWall(NodeWall2D& Node, DistriFunct * & f_in,double * & Rho, double * &U, double * &V);
	///Apply Corner for SpecialWall as  several 0-moment variables
//	void ApplyCornerSpecialWall(NodeWall2D& Node, DistriFunct * & f_in,double * & Rho, double * &U, double * &V);
	///Apply Wall as one 0-moment variable (single phase = 1 density)
	void ApplyWall(int const &BcNormal,int const *Connect, DistriFunct * & f_in);
	///Apply Wall as  several 0-moment variables
	void ApplyWall(int const &BcNormal,int const *Connect, DistriFunct * & f_in,double * & Rho, double * &U, double * &V);
	///Apply SpecialWall as one 0-moment variable (single phase = 1 density)
//	void ApplySpecialWall(NodeWall2D& Node, DistriFunct * & f_in, std::map<int,NodeType> TypeOfNode_);
	///Apply SpecialWall as  several 0-moment variables
//	void ApplySpecialWall(NodeWall2D& Node, DistriFunct * & f_in, std::map<int,NodeType> TypeOfNode_,double * & Rho, double * &U, double * &V);
	///Apply Symmetry as one 0-moment variable (single phase = 1 density)
	void ApplySymmetry(int const &BcNormal,int const *Connect, double const Rho_def, double const *UDef, DistriFunct * & f_in);
	///Apply Symmetry as  several 0-moment variables
	void ApplySymmetry(int const &BcNormal,int const *Connect, double const Rho_def, double const *UDef, DistriFunct * & f_in,double * & Rho, double * &U, double * &V);
	///Apply Periodic as one 0-moment variable (single phase = 1 density)
	void ApplyPeriodic(int const &BcNormal,int const *Connect, double const Rho_def, double const *UDef, DistriFunct * & f_in);
	///Apply Periodic as  several 0-moment variables
	void ApplyPeriodic(int const &BcNormal,int const *Connect, double const Rho_def, double const *UDef, DistriFunct * & f_in,double * & Rho, double * &U, double * &V);

	const double& Get_Rho(int index){return PtrD2Q9BcRho[index];};
	const double& Get_U(int index){return PtrD2Q9BcU[index];};
	const double& Get_V(int index){return PtrD2Q9BcV[index];};
	double*& Get_Rho(){return PtrD2Q9BcRho;};
	double*& Get_U(){return PtrD2Q9BcU;};
	double*& Get_V(){return PtrD2Q9BcV;};

	inline void Get_CalculU(int const &BcNormal,int const *Connect, double const *UDef, double *U, double *V, double & Ureturn,double & Vreturn){D2Q9VelocityBc.Get_CalculU(BcNormal,Connect, UDef, U, V,Ureturn,Vreturn);};
	inline void Get_CalculRho(int const &BcNormal,int const *Connect, double const Rho_def, double *Rho, double & Rhoreturn){D2Q9PressureBc.Get_CalculRho(BcNormal,Connect, Rho_def, Rho, Rhoreturn);};

protected:
	double *PtrD2Q9BcRho, *PtrD2Q9BcU, *PtrD2Q9BcV;
	Dictionary *PtrD2Q9BcDic;
	//Boundary conditions objects
	D2Q9Pressure D2Q9PressureBc;
	D2Q9Velocity D2Q9VelocityBc;
	D2Q9Corner D2Q9CornerBc;
	D2Q9Wall D2Q9WallBc;
	D2Q9Symmetry D2Q9SymmetryBc;
	D2Q9Periodic D2Q9PeriodicBc;
};

#endif /* SRC_ALGORITHM_LOWORDER_BOUNDARIES_D2Q9GENERICBC_H_ */
