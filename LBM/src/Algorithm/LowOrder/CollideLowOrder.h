/*
 * Collide.h
 *
 *  Created on: 15 Jun 2015
 *      Author: thomas
 */

#ifndef ALGORITHM_LOWORDER_COLLIDELOWORDER_H_
#define ALGORITHM_LOWORDER_COLLIDELOWORDER_H_
#include "../../Core/GlobalDef.h"
#include "../../Mesh/SingleBlock.h"
#include "../../User/UserForce.h"
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
class CollideLowOrder: public UserForce {
public:
	CollideLowOrder();
	virtual ~CollideLowOrder();

	double EquiDistriFunct2D(double &rho_macro, double &u_macro, double &v_macro, double *u_i, double &omega);
	void Collide_SinglePhase(double &fi,double &rho, double &u, double &v, double *e_i, double &omega);

protected:
	double InvTau;
	Block* PtrBlockCollide;
	DistriFunct* PtrFiCollide;
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
inline double CollideLowOrder::EquiDistriFunct2D(double &rho_macro, double &u_macro, double &v_macro, double *u_i, double &omega){
	double dot_E_U=(u_i[0]*u_macro+u_i[1]*v_macro);
	return rho_macro*omega*(1+3*dot_E_U+4.5*dot_E_U*dot_E_U-1.5*(u_macro*u_macro+v_macro*v_macro));
}
class CollideD2Q9Colour: public CollideLowOrder {
public:
	CollideD2Q9Colour();
	virtual ~CollideD2Q9Colour();
	void Collide_ColorFluid(int & direction, double & fi,double &rho,double*  F,double & F_Norm, double & InvTau_, double &u, double &v);
	void Recoloring(double & f, double & fr, double & fb, double & Rho, double & Rho_r, double & Rho_b);
protected:
	double Ei[9][2];
	double omega[9];
	double Ak;

private:
	double TwoPhase_Collision_operator(int & i, double* F, double & F_Norm);
};
inline 	double CollideD2Q9Colour::TwoPhase_Collision_operator(int & i, double* F, double & F_Norm){
	 return Ak*0.5*F_Norm*(((F[0]* Ei[i][0]+F[1]* Ei[i][1])/F_Norm)*((F[0]* Ei[i][0]+F[1]* Ei[i][1])/F_Norm)-3/4);
}
inline void CollideD2Q9Colour::Recoloring(double & f, double & fr, double & fb, double & Rho, double & Rho_r, double & Rho_b)
{
		fr=Rho_r*f/Rho;
		fb=f-fr;

}
#endif /* ALGORITHM_LOWORDER_COLLIDELOWORDER_H_ */
