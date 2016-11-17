/*
 * ============================================================================
 * D2Q9GenericBc.cpp
 *
 *  Created on: 6 Sep 2016
 *      Author: Thomas Burel
 *      Version     :
 *      Copyright   : Your copyright notice
 *      Description : LBM Code Developing by JWFL at University of Strathclyde
 * ============================================================================
 */

#include "D2Q9GenericBc.h"

D2Q9GenericBc::D2Q9GenericBc() {
	PtrD2Q9BcRho=0;PtrD2Q9BcU=0; PtrD2Q9BcV=0;
	PtrD2Q9BcDic=0;

}

D2Q9GenericBc::~D2Q9GenericBc() {
	// TODO Auto-generated destructor stub
}

void D2Q9GenericBc::SetBcObjects( Parameters *Param){
	D2Q9PressureBc.Set_PressureBcs(Param);
	D2Q9VelocityBc.SetVelocity(Param);
	D2Q9CornerBc.Set_Corner(Param);
	D2Q9WallBc.SetWall(Param);
	D2Q9SymmetryBc.Set_Symmetry(Param);
	D2Q9PeriodicBc.Set_Periodic(Param);
}
///Apply pressure boundary conditions. Can be used on pressure nodes or global corners)
void D2Q9GenericBc::ApplyPressure(int const &BcNormal,int const *Connect, double const Rho_def, DistriFunct * & f_in){
	D2Q9PressureBc.ApplyPressure(BcNormal,Connect, Rho_def,f_in,PtrD2Q9BcRho,PtrD2Q9BcU,PtrD2Q9BcV);
}
void D2Q9GenericBc::ApplyPressure(int const &BcNormal,int const *Connect, double const Rho_def, DistriFunct * & f_in,double * & Rho, double * &U, double * &V){
	D2Q9PressureBc.ApplyPressure(BcNormal,Connect, Rho_def,f_in,Rho,U,V);
}
//Apply velocity boundary conditions. Can be used on velocity nodes or global corners)
void D2Q9GenericBc::ApplyVelocity(int const &BcNormal,int const *Connect, double const *UDef, DistriFunct * & f_in){
	D2Q9VelocityBc.ApplyVelocity(BcNormal,Connect, UDef,f_in,PtrD2Q9BcRho,PtrD2Q9BcU,PtrD2Q9BcV);
}
void D2Q9GenericBc::ApplyVelocity(int const &BcNormal,int const *Connect, double const *UDef, DistriFunct * & f_in,double * & Rho, double * &U, double * &V){
	D2Q9VelocityBc.ApplyVelocity(BcNormal,Connect, UDef,f_in,Rho,U,V);
}
//Apply Symmetry boundary conditions. Can be used on Symmetry nodes or global corners)
void D2Q9GenericBc::ApplySymmetry(int const &BcNormal,int const *Connect, double const Rho_def, double const *UDef, DistriFunct * & f_in){
	D2Q9SymmetryBc.ApplySymmetry(BcNormal,Connect,Rho_def,UDef,f_in,PtrD2Q9BcRho,PtrD2Q9BcU,PtrD2Q9BcV);
}
void D2Q9GenericBc::ApplySymmetry(int const &BcNormal,int const *Connect, double const Rho_def, double const *UDef, DistriFunct * & f_in,double * & Rho, double * &U, double * &V){
	D2Q9SymmetryBc.ApplySymmetry(BcNormal,Connect,Rho_def,UDef,f_in,Rho,U,V);
}
void D2Q9GenericBc::ApplyPeriodic(int const &BcNormal,int const *Connect, double const Rho_def, double const *UDef, DistriFunct * & f_in){
	D2Q9PeriodicBc.ApplyPeriodic(BcNormal,Connect,Rho_def,UDef,f_in,PtrD2Q9BcRho,PtrD2Q9BcU,PtrD2Q9BcV);
}
void D2Q9GenericBc::ApplyPeriodic(int const &BcNormal,int const *Connect, double const Rho_def, double const *UDef, DistriFunct * & f_in,double * & Rho, double * &U, double * &V){
	D2Q9PeriodicBc.ApplyPeriodic(BcNormal,Connect,Rho_def,UDef,f_in,Rho,U,V);
}
void D2Q9GenericBc::ApplyWall(int const &BcNormal,int const *Connect, DistriFunct * & f_in){
	D2Q9WallBc.ApplyWall(BcNormal,Connect,f_in,PtrD2Q9BcRho,PtrD2Q9BcU,PtrD2Q9BcV);
}
void D2Q9GenericBc::ApplyWall(int const &BcNormal,int const *Connect, DistriFunct * & f_in,double * & Rho, double * &U, double * &V){
	D2Q9WallBc.ApplyWall(BcNormal,Connect,f_in,Rho,U,V);
}
/*void D2Q9GenericBc::ApplySpecialWall(NodeWall2D& Node, DistriFunct * & f_in, std::map<int,NodeType> TypeOfNode_){
	D2Q9WallBc.ApplySpecialWall(Node,f_in,TypeOfNode_,PtrD2Q9BcRho,PtrD2Q9BcU,PtrD2Q9BcV);
}
void D2Q9GenericBc::ApplySpecialWall(NodeWall2D& Node, DistriFunct * & f_in, std::map<int,NodeType> TypeOfNode_,double * & Rho, double * &U, double * &V){
	D2Q9WallBc.ApplySpecialWall(Node,f_in,TypeOfNode_,Rho,U,V);
}*/
void D2Q9GenericBc::ApplyCorner(NodeCorner2D& Node, DistriFunct * & f_in,double const & RhoDef,double const & UDef,double const & VDef){
	D2Q9CornerBc.ApplyCorner(Node,f_in,RhoDef,UDef,VDef,PtrD2Q9BcRho,PtrD2Q9BcU,PtrD2Q9BcV);
}
void D2Q9GenericBc::ApplyCorner(NodeCorner2D& Node, DistriFunct * & f_in,double const & RhoDef,double const & UDef,double const & VDef,double * & Rho, double * &U, double * &V){
	D2Q9CornerBc.ApplyCorner(Node,f_in,RhoDef,UDef,VDef,Rho,U,V);
}
void D2Q9GenericBc::ApplyCornerWall(NodeCorner2D& Node, DistriFunct * & f_in){
	D2Q9CornerBc.ApplyCornerWall(Node,f_in,PtrD2Q9BcRho,PtrD2Q9BcU,PtrD2Q9BcV);
}
void D2Q9GenericBc::ApplyCornerWall(NodeCorner2D& Node, DistriFunct * & f_in,double * & Rho, double * &U, double * &V){
	D2Q9CornerBc.ApplyCornerWall(Node,f_in,Rho,U,V);
}
void D2Q9GenericBc::ApplyPreVelSpecialWall(NodeWall2D& Node, DistriFunct * & f_in,double const & RhoDef,double const & UDef,double const & VDef,double * & Rho, double * &U, double * &V){
	D2Q9CornerBc.ApplyPreVelSpecialWall(Node,f_in,RhoDef,UDef,VDef);
}
void D2Q9GenericBc::ApplyCornerSpecialWall(NodeWall2D& Node, DistriFunct * & f_in,double * & Rho, double * &U, double * &V){
	D2Q9CornerBc.ApplyCornerSpecialWall(Node,f_in,Rho,U,V);
}
