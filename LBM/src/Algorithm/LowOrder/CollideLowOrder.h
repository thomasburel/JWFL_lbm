/*
 * Collide.h
 *
 *  Created on: 15 Jun 2015
 *      Author: Thomas Burel
 *
 *  This file is used to calculate the value of the distribution after collision for a "single" phase (one distribution function)
 *  The Local force is an additional force gives by the User and directly added to the distribution.
 *  The body force is a the body force in Navier-Stokes as in the Z. Guo
 *  Non Constant Tau is used for non Newtonian fluid or tow phases flow with different viscosity
 */

#ifndef ALGORITHM_LOWORDER_COLLIDELOWORDER_H_
#define ALGORITHM_LOWORDER_COLLIDELOWORDER_H_

#include "../../Core/GlobalDef.h"
#include "../../Core/Parameters.h"
#include "../../Core/Dictionary.h"
#include "../../Mesh/SingleBlock.h"
#include "../../User/UserForce.h"
#include "../../SpecialTreatment/PorousMedia.h"
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

enum CollideType{Std2D,Std2DLocal,Std2DBody,Std2DLocalBody,Std2DNonCstTau,Std2DNonCstTauLocal,Std2DNonCstTauBody,Std2DNonCstTauLocalBody};

class CollideLowOrder: public UserForce, public PorousMedia
 {
public:
	CollideLowOrder();
	virtual ~CollideLowOrder();

	void Collide_2D(int & i, double &fi,double &rho, double &u, double &v, double & Fx, double & Fy, double InvTau_tmp);
	void Collide_2D_V2(double * &fi,double &rho, double &u, double &v, double *&Fi,double & Fx, double & Fy, double InvTau_tmp, double KineViscosity);
	void Select_Collide_2D(CollideType Type,double Cs2,double referenceDensity=1,SolverEnum::fluidmodel model=SolverEnum::Compressible);
	void Select_Collide_2D_V2(Parameters const &Param,CollideType Type,double Cs2,double referenceDensity=1,SolverEnum::fluidmodel model=SolverEnum::Compressible);

	double CollideEquillibrium(double &rho_macro, double &u_macro, double &v_macro, double *u_i, double &omega);

	double CompressibleEquiDistriFunct2D(double &rho_macro, double &u_macro, double &v_macro, double *u_i, double &omega);
	double IncompressibleEquiDistriFunct2D(double &rho_macro, double &u_macro, double &v_macro, double *u_i, double &omega);
	void Collide_2D_SinglePhase(int & i, double &fi,double &rho, double &u, double &v, double & Fx, double & Fy, double & InvTau_tmp);
	void Collide_2D_SinglePhase_V2(int & i, double &fi,double &rho, double &u, double &v, double & Fi, double & Fx, double & Fy, double & InvTau_tmp);
	void Collide_2D_SinglePhase_With_LocalForce(int & i, double &fi, double &rho, double &u, double &v, double & Fx, double & Fy, double & InvTau_tmp);
	void Collide_2D_SinglePhase_With_LocalForce_V2(int & i, double &fi, double &rho, double &u, double &v, double & Fi, double & Fx, double & Fy, double & InvTau_tmp);
	void Collide_2D_SinglePhase_With_BodyForce(int & i, double &fi, double &rho, double &u, double &v, double & Fx, double & Fy, double & InvTau_tmp);
	void Collide_2D_SinglePhase_With_BodyForce_V2(int & i, double &fi, double &rho, double &u, double &v, double & Fi, double & Fx, double & Fy, double & InvTau_tmp);
	void Collide_2D_SinglePhase_With_LocalForceAndUserLocalForce(int & i, double &fi, double &rho, double &u, double &v, double & Fx, double & Fy, double & InvTau_tmp);
	void Collide_2D_SinglePhase_With_BodyForceAndUserLocalForce(int & i, double &fi, double &rho, double &u, double &v, double & Fx, double & Fy, double & InvTau_tmp);
	void Collide_2D_SinglePhase_With_LocalForceAndUserBodyForce(int & i, double &fi, double &rho, double &u, double &v, double & Fx, double & Fy, double & InvTau_tmp);
	void Collide_2D_SinglePhase_With_LocalForceAndBodyForce_V2(int & i, double &fi, double &rho, double &u, double &v, double & Fi, double & Fx, double & Fy, double & InvTau_tmp);
	void Collide_2D_SinglePhase_With_BodyForceAndUserBodyForce(int & i, double &fi, double &rho, double &u, double &v, double & Fx, double & Fy, double & InvTau_tmp);
	double Collide_2D_BodyForce(int & i, double &u, double &v, double Fx, double Fy);

	void Collide_2D_SinglePhase_Non_Constant_Tau(int & i, double &fi,double &rho, double &u, double &v, double & Fx, double & Fy, double & InvTau_tmp);
	void Collide_2D_SinglePhase_Non_Constant_Tau_V2(int & i, double &fi,double &rho, double &u, double &v, double & Fi, double & Fx, double & Fy, double & InvTau_tmp);
	void Collide_2D_SinglePhase_Non_Constant_Tau_With_LocalForce(int & i, double &fi, double &rho, double &u, double &v, double & Fx, double & Fy, double & InvTau_tmp);
	void Collide_2D_SinglePhase_Non_Constant_Tau_With_LocalForce_V2(int & i, double &fi, double &rho, double &u, double &v, double & Fi, double & Fx, double & Fy, double & InvTau_tmp);
	void Collide_2D_SinglePhase_Non_Constant_Tau_With_BodyForce(int & i, double &fi, double &rho, double &u, double &v, double & Fx, double & Fy, double & InvTau_tmp);
	void Collide_2D_SinglePhase_Non_Constant_Tau_With_BodyForce_V2(int & i, double &fi, double &rho, double &u, double &v, double & Fi, double & Fx, double & Fy, double & InvTau_tmp);
	void Collide_2D_SinglePhase_Non_Constant_Tau_With_LocalForceAndUserLocalForce(int & i, double &fi, double &rho, double &u, double &v, double & Fx, double & Fy, double & InvTau_tmp);
	void Collide_2D_SinglePhase_Non_Constant_Tau_With_BodyForceAndUserLocalForce(int & i, double &fi, double &rho, double &u, double &v, double & Fx, double & Fy, double & InvTau_tmp);
	void Collide_2D_SinglePhase_Non_Constant_Tau_With_LocalForceAndUserBodyForce(int & i, double &fi, double &rho, double &u, double &v, double & Fx, double & Fy, double & InvTau_tmp);
	void Collide_2D_SinglePhase_Non_Constant_Tau_With_LocalForceAndBodyForce_V2(int & i, double &fi, double &rho, double &u, double &v, double & Fi, double & Fx, double & Fy, double & InvTau_tmp);
	void Collide_2D_SinglePhase_Non_Constant_Tau_With_BodyForceAndUserBodyForce(int & i, double &fi, double &rho, double &u, double &v, double & Fx, double & Fy, double & InvTau_tmp);

	double Collide_2D_BodyForce_Non_Constant_Tau(int & i, double &u, double &v, double Fx, double Fy, double & InvTau_tmp);
	void NoHeleShawDrag(double const &u, double const &v,double const &mu,double &Fx,double &Fy,double const InterfaceFx=0,double const InterfaceFy=0){};
	double NoUserLocalForce(int const direction, double const Rho, double const U, double const V, double const W=0){return 0;};
	double NoUserBodyForce(int const direction_xyz, double const Rho, double const U, double const  V, double const W=0){return 0;};
protected:
	double InvTau;
	Block* PtrBlockCollide;
	DistriFunct* PtrFiCollide;
	double **EiCollide;
	double *omegaCollide;
	double InvCs2Collide,InvCs2_2Collide,InvCs4_2Collide;
	double RefDensity;
	double dot_E_U;
	int Nb_VelocityCollide;
//	double Ei[9][2];
//	double omega[9];
	double D_tmp;// Temporary double
	double* PtrD_tmp;// Temporary pointer for a double
	double DVec_2D_tmp[2];// Temporary vector 2D for a double

// Pointers on function
	//Simplify notation for pointer on member functions
	typedef void(CollideLowOrder::*Collide_2D_TypeDef)(int & i, double &fi,double &rho, double &u, double &v, double & Fx, double & Fy, double & InvTau_tmp);
	typedef void(CollideLowOrder::*Collide_2D_V2_TypeDef)(int & i, double &fi,double &rho, double &u, double &v, double & Fi,double & Fx, double & Fy, double & InvTau_tmp);
	typedef void(CollideLowOrder::*AddHeleShaw)(double const &u, double const &v,double const &mu,double &Fx,double &Fy,double const InterfaceFx,double const InterfaceFy);
	typedef double(CollideLowOrder::*AddUserLocalForce)(int const direction, double const Rho, double const U, double const V, double const W);
	typedef double(CollideLowOrder::*AddUserBodyForce)(int const direction_xyz, double const Rho, double const U, double const  V, double const W);

	Collide_2D_TypeDef PtrCollide_2D;
	Collide_2D_V2_TypeDef PtrCollide_2D_V2;
	AddHeleShaw PtrAddHeleShaw;
	AddUserLocalForce PtrAddUserLocalForce;
	AddUserBodyForce PtrAddUserBodyForce;
	typedef double(CollideLowOrder::*EquiDistri)(double &rho_macro, double &u_macro, double &v_macro, double *u_i, double &omega);
	EquiDistri PtrEquiDistri;

	Dictionary *PtrDicCollide;//For debugging
	double **ColSingle,**ColTwoPhase;
	double ColSingle_tmp,ColTwoPhase_tmp;
private:
	friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
       ar & InvTau;
    }
};
/**
 * Equilibrium Distribution Function for D2Q9 scheme with BKG approximation (Qian and al 1992)
\f[
  {f_i}^{eq}=\rho\omega_i
  [1+3(\vec{e_i}\cdot\vec{u})
  +\frac{9}{2}{(\vec{e_i}\cdot\vec{u})}^2
  -\frac{3}{2}{\vec{u}}^2]
\f]
*/
inline double CollideLowOrder::CompressibleEquiDistriFunct2D(double &rho_macro, double &u_macro, double &v_macro, double *u_i, double &omega){
	dot_E_U=(u_i[0]*u_macro+u_i[1]*v_macro);
	return rho_macro*omega*(1.0+InvCs2Collide*dot_E_U+InvCs4_2Collide*dot_E_U*dot_E_U-InvCs2_2Collide*(u_macro*u_macro+v_macro*v_macro));
}
inline double CollideLowOrder::IncompressibleEquiDistriFunct2D(double &rho_macro, double &u_macro, double &v_macro, double *u_i, double &omega){
	dot_E_U=(u_i[0]*u_macro+u_i[1]*v_macro);
	return omega*(rho_macro+RefDensity*(InvCs2Collide*dot_E_U+InvCs4_2Collide*dot_E_U*dot_E_U-InvCs2_2Collide*(u_macro*u_macro+v_macro*v_macro)));
}
inline double CollideLowOrder::CollideEquillibrium(double &rho_macro, double &u_macro, double &v_macro, double *u_i, double &omega){
	return (this->*PtrEquiDistri)(rho_macro, u_macro, v_macro,u_i, omega);
}
#endif /* ALGORITHM_LOWORDER_COLLIDELOWORDER_H_ */
