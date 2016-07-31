/*
 * Collide.cpp
 *
 *  Created on: 15 Jun 2015
 *      Author: thomas
 */

#include "CollideLowOrder.h"

CollideLowOrder::CollideLowOrder() {
	// TODO Auto-generated constructor stub
	PtrBlockCollide=0;
	PtrFiCollide=0;
	InvTau=1;
}

CollideLowOrder::~CollideLowOrder() {
	// TODO Auto-generated destructor stub
}


void CollideLowOrder::Collide_SinglePhase(double &fi,double &rho, double &u, double &v, double *e_i, double &omega){
	fi= fi-InvTau*(fi-EquiDistriFunct2D(rho, u, v,e_i, omega));
}

CollideD2Q9Colour::CollideD2Q9Colour(){
//	Ei=new double [9][2];
	Ak=0;
}
CollideD2Q9Colour::~CollideD2Q9Colour(){

}
void CollideD2Q9Colour::Collide_ColorFluid(int & direction, double & fi,double &rho,double*  F,double & F_Norm, double & InvTau_, double &u, double &v)
{
	InvTau=InvTau_;
//	double F_Norm=std::sqrt(F[0]*F[0]+F[1]+F[1]);
	//f(0) => Ei=0
//	f[0]= f[0]-InvTau*(f[0]-rho*omega*(1+-1.5*(u*u+v*v)));
//	f[0]+=-Ak*0.5*F_Norm*3/4;

		CollideLowOrder::Collide_SinglePhase(fi,rho, u, v, Ei[direction], omega[direction]);
		fi+=CollideD2Q9Colour::TwoPhase_Collision_operator(direction, F, F_Norm);

}
/*void CollideD2Q9Colour::Collide_ColorFluid(double &fr,double &fb,double &rho_r,double &rho_b, double &u, double &v, double *e_i, double &omega){
	fr= fr-InvTau*(fr-EquiDistriFunct2D(rho_r, u, v,e_i, omega));
	fb= fb-InvTau*(fb-EquiDistriFunct2D(rho_r, u, v,e_i, omega));
}*/
