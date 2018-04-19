/*
 * ============================================================================
 * D2Q9SpecialWall.h
 *
 *  Created on: 2 Sep 2016
 *      Author: Thomas Burel
 *      Version     :
 *      Copyright   : Your copyright notice
 *      Description : LBM Code Developing by JWFL at University of Strathclyde
 * ============================================================================
 */

#ifndef SRC_ALGORITHM_LOWORDER_BOUNDARIES_D2Q9SPECIALWALL_H_
#define SRC_ALGORITHM_LOWORDER_BOUNDARIES_D2Q9SPECIALWALL_H_

#include "D2Q9BcVar.h"
#include "D2Q9GenericBc.h"


class D2Q9SpecialWall: public D2Q9BcVar {
public:
	D2Q9SpecialWall();
	virtual ~D2Q9SpecialWall();

	void Set_SpecialWall(Dictionary *PtrDic,NodeArrays2D* NodeArrays, Parameters *Param,D2Q9GenericBc* D2Q9GenericBc);
//	void ApplySpecialWall(NodeWall2D& Node, std::map<int,NodeType> TypeOfNode_, DistriFunct* f_in);
	//Specify in the solver: set values (Rho, U) and pointers on macroscopic variables
	void ApplySpecialWall(NodeWall2D& Node, double const Rho_def, double const UDef, double const VDef, std::map<int,NodeType> TypeOfNode_, DistriFunct* f_in,double * & Rho, double * &U, double * &V, unsigned int idxDistribution=0);
	void ApplyPeriodicWall(NodeWall2D& Node, double const Rho_def, double const UDef, double const VDef, std::map<int,NodeType> TypeOfNode_, DistriFunct* f_in,double * & Rho, double * &U, double * &V, unsigned int idxDistribution=0);
	void ApplySymmetryWall(NodeWall2D& Node, double const Rho_def, double const UDef, double const VDef, std::map<int,NodeType> TypeOfNode_, DistriFunct* f_in,double * & Rho, double * &U, double * &V, unsigned int idxDistribution=0);
	void ApplyPressureWall(NodeWall2D& Node, double const Rho_def, std::map<int,NodeType> TypeOfNode_, DistriFunct* f_in,double * & Rho, double * &U, double * &V, unsigned int idxDistribution=0);
	void ApplyVelocityWall(NodeWall2D& Node, double const UDef, double const VDef, std::map<int,NodeType> TypeOfNode_, DistriFunct* f_in,double * & Rho, double * &U, double * &V, unsigned int idxDistribution=0);

	void ApplyPeriodicWallPreStream(NodeWall2D& Node, double const Rho_def, double const UDef, double const VDef, std::map<int,NodeType> TypeOfNode_, DistriFunct* f_in,double * & Rho, double * &U, double * &V, unsigned int idxDistribution=0);
	void ApplySymmetryWallPreStream(NodeWall2D& Node, double const Rho_def, double const UDef, double const VDef, std::map<int,NodeType> TypeOfNode_, DistriFunct* f_in,double * & Rho, double * &U, double * &V, unsigned int idxDistribution=0);
	void ApplyPressureWallPreStream(NodeWall2D& Node, double const Rho_def, std::map<int,NodeType> TypeOfNode_, DistriFunct* f_in,double * & Rho, double * &U, double * &V, unsigned int idxDistribution=0);
	void ApplyVelocityWallPreStream(NodeWall2D& Node, double const UDef, double const VDef, std::map<int,NodeType> TypeOfNode_, DistriFunct* f_in,double * & Rho, double * &U, double * &V, unsigned int idxDistribution=0);
	
private:
	void FunctionSpecialWall(NodeWall2D& Node, double const Rho_def, double const *UDef, std::map<int,NodeType> &TypeOfNode_, DistriFunct* &f_in,double * & Rho, double * &U, double * &V, unsigned int idxDistribution=0);
	void FunctionPeriodicWall(NodeWall2D& Node, double const Rho_def, double const *UDef, std::map<int,NodeType> &TypeOfNode_, DistriFunct* &f_in,double * & Rho, double * &U, double * &V, unsigned int idxDistribution=0);
	void FunctionSymmetryWall(NodeWall2D& Node, double const Rho_def, double const *UDef, std::map<int,NodeType> &TypeOfNode_, DistriFunct* &f_in,double * & Rho, double * &U, double * &V, unsigned int idxDistribution=0);
	void FunctionPressureWall(NodeWall2D& Node, double const Rho_def, std::map<int,NodeType> &TypeOfNode_, DistriFunct* &f_in,double * & Rho, double * &U, double * &V, unsigned int idxDistribution=0);
	void FunctionVelocityWall(NodeWall2D& Node, double const *UDef, std::map<int,NodeType> &TypeOfNode_, DistriFunct* &f_in,double * & Rho, double * &U, double * &V, unsigned int idxDistribution=0);

	void FunctionSpecialWallPreStream(NodeWall2D& Node, double const Rho_def, double const *UDef, std::map<int,NodeType> &TypeOfNode_, DistriFunct* &f_in,double * & Rho, double * &U, double * &V, unsigned int idxDistribution=0);
	void FunctionPeriodicWallPreStream(NodeWall2D& Node, double const Rho_def, double const *UDef, std::map<int,NodeType> &TypeOfNode_, DistriFunct* &f_in,double * & Rho, double * &U, double * &V, unsigned int idxDistribution=0);
	void FunctionSymmetryWallPreStream(NodeWall2D& Node, double const Rho_def, double const *UDef, std::map<int,NodeType> &TypeOfNode_, DistriFunct* &f_in,double * & Rho, double * &U, double * &V, unsigned int idxDistribution=0);
	void FunctionPressureWallPreStream(NodeWall2D& Node, double const Rho_def, std::map<int,NodeType> &TypeOfNode_, DistriFunct* &f_in,double * & Rho, double * &U, double * &V, unsigned int idxDistribution=0);
	void FunctionVelocityWallPreStream(NodeWall2D& Node, double const *UDef, std::map<int,NodeType> &TypeOfNode_, DistriFunct* &f_in,double * & Rho, double * &U, double * &V, unsigned int idxDistribution=0);

	D2Q9GenericBc* BcMethods;
	Extrapolation Extrapol;
	double RhoDef_tmp,UDef_tmp,VDef_tmp;
	double U_tmp[2];
};

#endif /* SRC_ALGORITHM_LOWORDER_BOUNDARIES_D2Q9SPECIALWALL_H_ */
