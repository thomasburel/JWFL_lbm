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
	PtrD_tmp=0;
	PtrCollide_2D=0;
	D_tmp=0;
	InvTau=1;
	Select_Collide_2D(Std2D);
}

CollideLowOrder::~CollideLowOrder() {
	// TODO Auto-generated destructor stub
}

void CollideLowOrder::Collide_2D(int & i, double &fi,double &rho, double &u, double &v, double* F_tmp, double & InvTau_tmp){
	(this->*PtrCollide_2D)(i,fi,rho,u,v,F_tmp,InvTau_tmp);
}
void CollideLowOrder::Select_Collide_2D(CollideType type){
	switch (type)
	{
		case Std2D:
			PtrCollide_2D=&CollideLowOrder::Collide_2D_SinglePhase; break;
		case Std2DLocal:
			PtrCollide_2D=&CollideLowOrder::Collide_2D_SinglePhase_With_LocalForce; break;
		case Std2DBody:
			PtrCollide_2D=&CollideLowOrder::Collide_2D_SinglePhase_With_BodyForce; break;
		case Std2DNonCstTau:
			PtrCollide_2D=&CollideLowOrder::Collide_2D_SinglePhase_Non_Constant_Tau; break;
		case Std2DNonCstTauLocal:
			PtrCollide_2D=&CollideLowOrder::Collide_2D_SinglePhase_Non_Constant_Tau_With_LocalForce; break;
		case Std2DNonCstTauBody:
			PtrCollide_2D=&CollideLowOrder::Collide_2D_SinglePhase_Non_Constant_Tau_With_BodyForce; break;
		default:
			std::cerr<< "Collide type not found"<<std::endl;
	}
}

void CollideLowOrder::Collide_2D_SinglePhase(int & i, double &fi,double &rho, double &u, double &v, double* F_tmp, double & InvTau_tmp){
	fi= fi-InvTau*(fi-EquiDistriFunct2D(rho, u, v,Ei[i], omega[i]));
}
void CollideLowOrder::Collide_2D_SinglePhase_With_LocalForce(int & i, double &fi,double &rho, double &u, double &v, double* F_tmp, double & InvTau_tmp){
	fi= fi-InvTau*(fi-EquiDistriFunct2D(rho, u, v,Ei[i], omega[i]))+LocalForce(i, rho, u, v);
}
void CollideLowOrder::Collide_2D_SinglePhase_With_BodyForce(int & i, double &fi,double &rho, double &u, double &v, double* F_tmp, double & InvTau_tmp){
	F_tmp[0]=F_tmp[0]+BodyForce(0,rho,u,v);
	F_tmp[1]=F_tmp[1]+BodyForce(1,rho,u,v);
	fi= fi-InvTau*(fi-EquiDistriFunct2D(rho, u, v,Ei[i], omega[i]))+
			Collide_2D_BodyForce(i, u, v, F_tmp);
}
double CollideLowOrder::Collide_2D_BodyForce(int & i, double &u, double &v, double* F_tmp){
	return F_tmp[0]*(1-InvTau*0.5)*omega[i]*3*((Ei[i][0]-u)+3*Ei[i][0]*(Ei[i][0]*u+Ei[i][1]*u))+
	F_tmp[1]*(1-InvTau*0.5)*omega[i]*3*((Ei[i][1]-v)+3*Ei[i][1]*(Ei[i][0]*v+Ei[i][1]*v));

}

void CollideLowOrder::Collide_2D_SinglePhase_Non_Constant_Tau(int & i, double &fi,double &rho, double &u, double &v, double* F_tmp, double & InvTau_tmp){
	fi= fi-InvTau_tmp*(fi-EquiDistriFunct2D(rho, u, v,Ei[i], omega[i]));
}
void CollideLowOrder::Collide_2D_SinglePhase_Non_Constant_Tau_With_LocalForce(int & i, double &fi,double &rho, double &u, double &v, double* F_tmp, double & InvTau_tmp){
	fi= fi-InvTau_tmp*(fi-EquiDistriFunct2D(rho, u, v,Ei[i], omega[i]))+LocalForce(i, rho, u, v);;
}
void CollideLowOrder::Collide_2D_SinglePhase_Non_Constant_Tau_With_BodyForce(int & i, double &fi,double &rho, double &u, double &v, double* F_tmp, double & InvTau_tmp){
	F_tmp[0]=F_tmp[0]+BodyForce(0,rho,u,v);
	F_tmp[1]=F_tmp[1]+BodyForce(1,rho,u,v);
	fi= fi-InvTau_tmp*(fi-EquiDistriFunct2D(rho, u, v,Ei[i], omega[i]))+
			Collide_2D_BodyForce_Non_Constant_Tau(i, u, v, F_tmp, InvTau_tmp);
}
double CollideLowOrder::Collide_2D_BodyForce_Non_Constant_Tau(int & i, double &u, double &v, double* F_tmp, double & InvTau_tmp){
	return F_tmp[0]*(1-InvTau_tmp*0.5)*omega[i]*3*((Ei[i][0]-u)+3*Ei[i][0]*(Ei[i][0]*u+Ei[i][1]*u))+
	F_tmp[1]*(1-InvTau_tmp*0.5)*omega[i]*3*((Ei[i][1]-v)+3*Ei[i][1]*(Ei[i][0]*v+Ei[i][1]*v));

}
/*CollideD2Q9Colour::CollideD2Q9Colour(){
//	Ei=new double [9][2];
	Ak=0.01;
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

}*/
/*void CollideD2Q9Colour::Collide_ColorFluid(double &fr,double &fb,double &rho_r,double &rho_b, double &u, double &v, double *e_i, double &omega){
	fr= fr-InvTau*(fr-EquiDistriFunct2D(rho_r, u, v,e_i, omega));
	fb= fb-InvTau*(fb-EquiDistriFunct2D(rho_r, u, v,e_i, omega));
}*/
