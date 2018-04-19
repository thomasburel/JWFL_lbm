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
	PtrCorrectVelocity=0;
	D_tmp=0;
	InvTau=1;
	Select_Collide_2D(Std2D,1.0/3.0);
	InvCs2Collide=3.0;
	InvCs2_2Collide=3.0/2.0;

	InvCs4_2Collide=9.0/2.0;
	RefDensity=1.0;
	PtrEquiDistri=0;
}

CollideLowOrder::~CollideLowOrder() {
	// TODO Auto-generated destructor stub
}

void CollideLowOrder::Collide_2D(int & i, double &fi,double &rho, double &u, double &v, double & Fx, double & Fy, double InvTau_tmp){
	(this->*PtrCollide_2D)(i,fi,rho,u,v,Fx,Fy,InvTau_tmp);
}
void CollideLowOrder::Collide_2D_V2(double * &fi,double &rho, double &u, double &v, double *&Fi,double & Fx, double & Fy, double InvTau_tmp, double KineViscosity){
	//add Hele-Shaw Body Force
	(this->*PtrAddHeleShaw)(u, v,KineViscosity,Fx,Fy,0,0);//Need to update for adding the contribution of the interface forces
	//add user local force
	for(int i=0;i<Nb_VelocityCollide;i++)
	{
		Fi[i]+=(this->*PtrAddUserLocalForce)(i, rho, u, v,0);
	}
	//add User body force
	Fx+=(this->*PtrAddUserBodyForce)(0, rho, u, v,0);
	Fy+=(this->*PtrAddUserBodyForce)(1, rho, u, v,0);

	//add velocity from body force if needed

	(this->*PtrCorrectVelocity)(u,v,rho,Fx,Fy);
	for(int i=0;i<Nb_VelocityCollide;i++)
	{
		(this->*PtrCollide_2D_V2)(i,fi[i],rho,u,v,Fi[i],Fx,Fy,InvTau_tmp);
	}
}
///Collision and add source term on walls
void CollideLowOrder::Collide_2D_V2_Wall(double * &fi,double &rho, double &u, double &v, double *&Fi,double & Fx, double & Fy, double InvTau_tmp, double KineViscosity){
	//add Hele-Shaw Body Force
	(this->*PtrAddHeleShaw)(u, v,KineViscosity,Fx,Fy,0,0);//Need to update for adding the contribution of the interface forces
	//add user local force
	for(int i=0;i<Nb_VelocityCollide;i++)
	{
		Fi[i]+=(this->*PtrAddUserLocalForce)(i, rho, u, v,0);
	}
	//add User body force
	Fx+=(this->*PtrAddUserBodyForce)(0, rho, u, v,0);
	Fy+=(this->*PtrAddUserBodyForce)(1, rho, u, v,0);

	//add velocity from body force if needed

	(this->*PtrCorrectVelocity)(u,v,rho,Fx,Fy);
	for(int i=0;i<Nb_VelocityCollide;i++)
	{
		(this->*PtrCollide_2D_V2_Wall)(i,fi[i],rho,u,v,Fi[i],Fx,Fy,InvTau_tmp);
	}
}
void CollideLowOrder::Select_Collide_2D(CollideType type,double Cs2,double referenceDensity,SolverEnum::fluidmodel model){
	switch (type)
	{
		case Std2D:
			PtrCollide_2D=&CollideLowOrder::Collide_2D_SinglePhase;
			PtrCollide_2D_V2=&CollideLowOrder::Collide_2D_SinglePhase_V2;
			break;
		case Std2DLocal:
			PtrCollide_2D=&CollideLowOrder::Collide_2D_SinglePhase_With_LocalForce;
			PtrCollide_2D_V2=&CollideLowOrder::Collide_2D_SinglePhase_With_LocalForce_V2;
			break;
		case Std2DBody:
			PtrCollide_2D=&CollideLowOrder::Collide_2D_SinglePhase_With_BodyForce;
			PtrCollide_2D_V2=&CollideLowOrder::Collide_2D_SinglePhase_With_BodyForce_V2;
			break;
		case Std2DNonCstTau:
			PtrCollide_2D=&CollideLowOrder::Collide_2D_SinglePhase_Non_Constant_Tau;
			PtrCollide_2D_V2=&CollideLowOrder::Collide_2D_SinglePhase_Non_Constant_Tau_V2;
			break;
		case Std2DNonCstTauLocal:
			PtrCollide_2D=&CollideLowOrder::Collide_2D_SinglePhase_Non_Constant_Tau_With_LocalForce;
			PtrCollide_2D_V2=&CollideLowOrder::Collide_2D_SinglePhase_Non_Constant_Tau_With_LocalForce_V2;
			break;
		case Std2DNonCstTauBody:
			PtrCollide_2D=&CollideLowOrder::Collide_2D_SinglePhase_Non_Constant_Tau_With_BodyForce;
			PtrCollide_2D_V2=&CollideLowOrder::Collide_2D_SinglePhase_Non_Constant_Tau_With_BodyForce_V2;
			break;
		default:
			std::cerr<< "Collide type not found"<<std::endl;
	}
	InvCs2Collide=1.0/Cs2;
	InvCs2_2Collide=1.0/(2.0*Cs2);
	InvCs4_2Collide=1.0/(2.0*Cs2*Cs2);
	RefDensity=referenceDensity;
	if(model==SolverEnum::Incompressible)
		PtrEquiDistri=&CollideLowOrder::IncompressibleEquiDistriFunct2D;
	else
		PtrEquiDistri=&CollideLowOrder::CompressibleEquiDistriFunct2D;
}
void CollideLowOrder::Select_Collide_2D_V2(Parameters const &Param,CollideType type,double Cs2,double referenceDensity,SolverEnum::fluidmodel model){
//	std::cout<<"type: "<<type;

	switch (type)
	{
		case Std2D:
			PtrCollide_2D=&CollideLowOrder::Collide_2D_SinglePhase;
			if(Param.Get_UserForceType()==ModelEnum::LocalForce)
			{
				PtrCollide_2D_V2=&CollideLowOrder::Collide_2D_SinglePhase_With_LocalForce_V2;
				PtrAddUserLocalForce=&UserForce::LocalForce;
				PtrAddUserBodyForce=&CollideLowOrder::NoUserBodyForce;
				PtrCorrectVelocity=&CollideLowOrder::NoVelocityCorrection;
				if(Param.Get_CollisionOnWalls()) PtrCollide_2D_V2_Wall=PtrCollide_2D_V2;
				else PtrCollide_2D_V2_Wall=&CollideLowOrder::No_Collide_2D_SinglePhase_With_LocalForce_V2;
			}
			else if(Param.Get_UserForceType()==ModelEnum::BodyForce)
			{
				PtrCollide_2D_V2=&CollideLowOrder::Collide_2D_SinglePhase_With_BodyForce_V2;
				PtrAddUserLocalForce=&CollideLowOrder::NoUserLocalForce;
				PtrAddUserBodyForce=&UserForce::BodyForce;
				PtrCorrectVelocity=&CollideLowOrder::VelocityCorrection;
				if(Param.Get_CollisionOnWalls()) PtrCollide_2D_V2_Wall=PtrCollide_2D_V2;
				else PtrCollide_2D_V2_Wall=&CollideLowOrder::No_Collide_2D_SinglePhase_With_BodyForce_V2;
			}
			else
			{
				PtrCollide_2D_V2=&CollideLowOrder::Collide_2D_SinglePhase_V2;
				PtrAddUserLocalForce=&CollideLowOrder::NoUserLocalForce;
				PtrAddUserBodyForce=&CollideLowOrder::NoUserBodyForce;
				PtrCorrectVelocity=&CollideLowOrder::NoVelocityCorrection;
				if(Param.Get_CollisionOnWalls()) PtrCollide_2D_V2_Wall=PtrCollide_2D_V2;
				else PtrCollide_2D_V2_Wall=&CollideLowOrder::No_Collide_2D_SinglePhase_V2;
			}
			break;
		case Std2DLocal:
			PtrCollide_2D=&CollideLowOrder::Collide_2D_SinglePhase_With_LocalForce;
			if(Param.Get_UserForceType()==ModelEnum::LocalForce)
			{
				PtrCollide_2D_V2=&CollideLowOrder::Collide_2D_SinglePhase_With_LocalForce_V2;
				PtrAddUserLocalForce=&UserForce::LocalForce;
				PtrAddUserBodyForce=&CollideLowOrder::NoUserBodyForce;
				PtrCorrectVelocity=&CollideLowOrder::NoVelocityCorrection;
				if(Param.Get_CollisionOnWalls()) PtrCollide_2D_V2_Wall=PtrCollide_2D_V2;
				else PtrCollide_2D_V2_Wall=&CollideLowOrder::No_Collide_2D_SinglePhase_With_LocalForce_V2;
			}
			else if(Param.Get_UserForceType()==ModelEnum::BodyForce)
			{
				PtrCollide_2D_V2=&CollideLowOrder::Collide_2D_SinglePhase_With_LocalForceAndBodyForce_V2;
				PtrAddUserLocalForce=&CollideLowOrder::NoUserLocalForce;
				PtrAddUserBodyForce=&UserForce::BodyForce;
				PtrCorrectVelocity=&CollideLowOrder::VelocityCorrection;
				if(Param.Get_CollisionOnWalls()) PtrCollide_2D_V2_Wall=PtrCollide_2D_V2;
				else PtrCollide_2D_V2_Wall=&CollideLowOrder::No_Collide_2D_SinglePhase_With_LocalForceAndBodyForce_V2;
			}
			else
			{
				PtrCollide_2D_V2=&CollideLowOrder::Collide_2D_SinglePhase_With_LocalForce_V2;
				PtrAddUserLocalForce=&CollideLowOrder::NoUserLocalForce;
				PtrAddUserBodyForce=&CollideLowOrder::NoUserBodyForce;
				PtrCorrectVelocity=&CollideLowOrder::NoVelocityCorrection;
				if(Param.Get_CollisionOnWalls()) PtrCollide_2D_V2_Wall=PtrCollide_2D_V2;
				else PtrCollide_2D_V2_Wall=&CollideLowOrder::No_Collide_2D_SinglePhase_With_LocalForce_V2;
			}
			break;
		case Std2DBody:
			PtrCollide_2D=&CollideLowOrder::Collide_2D_SinglePhase_With_BodyForce;
			if(Param.Get_UserForceType()==ModelEnum::LocalForce)
			{
				PtrCollide_2D_V2=&CollideLowOrder::Collide_2D_SinglePhase_With_LocalForceAndBodyForce_V2;
				PtrAddUserLocalForce=&UserForce::LocalForce;
				PtrAddUserBodyForce=&CollideLowOrder::NoUserBodyForce;
				PtrCorrectVelocity=&CollideLowOrder::VelocityCorrection;
				if(Param.Get_CollisionOnWalls()) PtrCollide_2D_V2_Wall=PtrCollide_2D_V2;
				else PtrCollide_2D_V2_Wall=&CollideLowOrder::No_Collide_2D_SinglePhase_With_LocalForceAndBodyForce_V2;
			}
			else if(Param.Get_UserForceType()==ModelEnum::BodyForce)
			{
				PtrCollide_2D_V2=&CollideLowOrder::Collide_2D_SinglePhase_With_BodyForce_V2;
				PtrAddUserLocalForce=&CollideLowOrder::NoUserLocalForce;
				PtrAddUserBodyForce=&UserForce::BodyForce;
				PtrCorrectVelocity=&CollideLowOrder::VelocityCorrection;
				if(Param.Get_CollisionOnWalls()) PtrCollide_2D_V2_Wall=PtrCollide_2D_V2;
				else PtrCollide_2D_V2_Wall=&CollideLowOrder::No_Collide_2D_SinglePhase_With_BodyForce_V2;
			}
			else
			{
				PtrCollide_2D_V2=&CollideLowOrder::Collide_2D_SinglePhase_With_BodyForce_V2;
				PtrAddUserLocalForce=&CollideLowOrder::NoUserLocalForce;
				PtrAddUserBodyForce=&CollideLowOrder::NoUserBodyForce;
				PtrCorrectVelocity=&CollideLowOrder::VelocityCorrection;
				if(Param.Get_CollisionOnWalls()) PtrCollide_2D_V2_Wall=PtrCollide_2D_V2;
				else PtrCollide_2D_V2_Wall=&CollideLowOrder::No_Collide_2D_SinglePhase_With_BodyForce_V2;
			}
			break;
		case Std2DNonCstTau:
			PtrCollide_2D=&CollideLowOrder::Collide_2D_SinglePhase_Non_Constant_Tau;
			if(Param.Get_UserForceType()==ModelEnum::LocalForce)
			{
				PtrCollide_2D_V2=&CollideLowOrder::Collide_2D_SinglePhase_Non_Constant_Tau_With_LocalForce_V2;
				PtrAddUserLocalForce=&UserForce::LocalForce;
				PtrAddUserBodyForce=&CollideLowOrder::NoUserBodyForce;
				PtrCorrectVelocity=&CollideLowOrder::NoVelocityCorrection;
				if(Param.Get_CollisionOnWalls()) PtrCollide_2D_V2_Wall=PtrCollide_2D_V2;
				else PtrCollide_2D_V2_Wall=&CollideLowOrder::No_Collide_2D_SinglePhase_Non_Constant_Tau_With_LocalForce_V2;
			}
			else if(Param.Get_UserForceType()==ModelEnum::BodyForce)
			{
				PtrCollide_2D_V2=&CollideLowOrder::Collide_2D_SinglePhase_Non_Constant_Tau_With_BodyForce_V2;
				PtrAddUserLocalForce=&CollideLowOrder::NoUserLocalForce;
				PtrAddUserBodyForce=&UserForce::BodyForce;
				PtrCorrectVelocity=&CollideLowOrder::VelocityCorrection;
				if(Param.Get_CollisionOnWalls()) PtrCollide_2D_V2_Wall=PtrCollide_2D_V2;
				else PtrCollide_2D_V2_Wall=&CollideLowOrder::No_Collide_2D_SinglePhase_Non_Constant_Tau_With_BodyForce_V2;
			}
			else
			{
				PtrCollide_2D_V2=&CollideLowOrder::Collide_2D_SinglePhase_Non_Constant_Tau_V2;
				PtrAddUserLocalForce=&CollideLowOrder::NoUserLocalForce;
				PtrAddUserBodyForce=&CollideLowOrder::NoUserBodyForce;
				PtrCorrectVelocity=&CollideLowOrder::NoVelocityCorrection;
				if(Param.Get_CollisionOnWalls()) PtrCollide_2D_V2_Wall=PtrCollide_2D_V2;
				else PtrCollide_2D_V2_Wall=&CollideLowOrder::No_Collide_2D_SinglePhase_Non_Constant_Tau_V2;
			}
			break;
		case Std2DNonCstTauLocal:
			PtrCollide_2D=&CollideLowOrder::Collide_2D_SinglePhase_Non_Constant_Tau_With_LocalForce;
			if(Param.Get_UserForceType()==ModelEnum::LocalForce)
			{
				PtrCollide_2D_V2=&CollideLowOrder::Collide_2D_SinglePhase_Non_Constant_Tau_With_LocalForce_V2;
				PtrAddUserLocalForce=&UserForce::LocalForce;
				PtrAddUserBodyForce=&CollideLowOrder::NoUserBodyForce;
				PtrCorrectVelocity=&CollideLowOrder::NoVelocityCorrection;
				if(Param.Get_CollisionOnWalls()) PtrCollide_2D_V2_Wall=PtrCollide_2D_V2;
				else PtrCollide_2D_V2_Wall=&CollideLowOrder::No_Collide_2D_SinglePhase_Non_Constant_Tau_With_LocalForce_V2;
			}
			else if(Param.Get_UserForceType()==ModelEnum::BodyForce)
			{
				PtrCollide_2D_V2=&CollideLowOrder::Collide_2D_SinglePhase_Non_Constant_Tau_With_LocalForceAndBodyForce_V2;
				PtrAddUserLocalForce=&CollideLowOrder::NoUserLocalForce;
				PtrAddUserBodyForce=&UserForce::BodyForce;
				PtrCorrectVelocity=&CollideLowOrder::VelocityCorrection;
				if(Param.Get_CollisionOnWalls()) PtrCollide_2D_V2_Wall=PtrCollide_2D_V2;
				else PtrCollide_2D_V2_Wall=&CollideLowOrder::No_Collide_2D_SinglePhase_Non_Constant_Tau_With_LocalForceAndBodyForce_V2;
			}
			else
			{
				PtrCollide_2D_V2=&CollideLowOrder::Collide_2D_SinglePhase_Non_Constant_Tau_With_LocalForce_V2;
				PtrAddUserLocalForce=&CollideLowOrder::NoUserLocalForce;
				PtrAddUserBodyForce=&CollideLowOrder::NoUserBodyForce;
				PtrCorrectVelocity=&CollideLowOrder::NoVelocityCorrection;
				if(Param.Get_CollisionOnWalls()) PtrCollide_2D_V2_Wall=PtrCollide_2D_V2;
				else PtrCollide_2D_V2_Wall=&CollideLowOrder::No_Collide_2D_SinglePhase_Non_Constant_Tau_With_LocalForce_V2;
			}
			break;
		case Std2DNonCstTauBody:
			PtrCollide_2D=&CollideLowOrder::Collide_2D_SinglePhase_Non_Constant_Tau_With_BodyForce;
			if(Param.Get_UserForceType()==ModelEnum::LocalForce)
			{
				PtrCollide_2D_V2=&CollideLowOrder::Collide_2D_SinglePhase_Non_Constant_Tau_With_LocalForceAndBodyForce_V2;
				PtrAddUserLocalForce=&UserForce::LocalForce;
				PtrAddUserBodyForce=&CollideLowOrder::NoUserBodyForce;
				PtrCorrectVelocity=&CollideLowOrder::VelocityCorrection;
				if(Param.Get_CollisionOnWalls()) PtrCollide_2D_V2_Wall=PtrCollide_2D_V2;
				else PtrCollide_2D_V2_Wall=&CollideLowOrder::No_Collide_2D_SinglePhase_Non_Constant_Tau_With_LocalForceAndBodyForce_V2;
			}
			else if(Param.Get_UserForceType()==ModelEnum::BodyForce)
			{
				PtrCollide_2D_V2=&CollideLowOrder::Collide_2D_SinglePhase_Non_Constant_Tau_With_BodyForce_V2;
				PtrAddUserLocalForce=&CollideLowOrder::NoUserLocalForce;
				PtrAddUserBodyForce=&UserForce::BodyForce;
				PtrCorrectVelocity=&CollideLowOrder::VelocityCorrection;
				if(Param.Get_CollisionOnWalls()) PtrCollide_2D_V2_Wall=PtrCollide_2D_V2;
				else PtrCollide_2D_V2_Wall=&CollideLowOrder::No_Collide_2D_SinglePhase_Non_Constant_Tau_With_BodyForce_V2;
			}
			else
			{
				PtrCollide_2D_V2=&CollideLowOrder::Collide_2D_SinglePhase_Non_Constant_Tau_With_BodyForce_V2;
				PtrAddUserLocalForce=&CollideLowOrder::NoUserLocalForce;
				PtrAddUserBodyForce=&CollideLowOrder::NoUserBodyForce;
				PtrCorrectVelocity=&CollideLowOrder::VelocityCorrection;
				if(Param.Get_CollisionOnWalls()) PtrCollide_2D_V2_Wall=PtrCollide_2D_V2;
				else PtrCollide_2D_V2_Wall=&CollideLowOrder::No_Collide_2D_SinglePhase_Non_Constant_Tau_With_BodyForce_V2;
			}
			break;
		default:
			std::cerr<< "Collide type not found"<<std::endl;
	}

	if(Param.Get_Dimension()==SolverEnum::D2 && Param.Get_HeleShawBodyForce() != PorousMediaEnum::no)
	{
		depth=Param.Get_Depth();
		depth2=depth*depth;
		if(Param.Get_HeleShawBodyForce() == PorousMediaEnum::SinglePhase)
			PtrAddHeleShaw=&PorousMedia::AddHeleShawDragSinglePhase;
		else
			PtrAddHeleShaw=&PorousMedia::AddHeleShawDragTwoPhases;
	}
	else
		PtrAddHeleShaw=&CollideLowOrder::NoHeleShawDrag;

	InvCs2Collide=1.0/Cs2;
	InvCs2_2Collide=1.0/(2.0*Cs2);
	InvCs4_2Collide=1.0/(2.0*Cs2*Cs2);
	RefDensity=referenceDensity;
	if(model==SolverEnum::Incompressible)
		PtrEquiDistri=&CollideLowOrder::IncompressibleEquiDistriFunct2D;
	else
		PtrEquiDistri=&CollideLowOrder::CompressibleEquiDistriFunct2D;


}

///Standard Collision Single phase step
void CollideLowOrder::Collide_2D_SinglePhase(int & i, double &fi,double &rho, double &u, double &v, double & Fx, double & Fy, double & InvTau_tmp){
	fi= fi-InvTau*(fi-CollideEquillibrium(rho, u, v,EiCollide[i], omegaCollide[i]));
}
void CollideLowOrder::Collide_2D_SinglePhase_V2(int & i, double &fi,double &rho, double &u, double &v, double & Fi, double & Fx, double & Fy, double & InvTau_tmp){
	fi= fi-InvTau*(fi-CollideEquillibrium(rho, u, v,EiCollide[i], omegaCollide[i]));
}
///Collision Single phase step with a local (directional) force calculated from User force or model or both
void CollideLowOrder::Collide_2D_SinglePhase_With_LocalForce(int & i, double &fi,double &rho, double &u, double &v, double & Fx, double & Fy, double & InvTau_tmp){
	Fx+=LocalForce(i, rho, u, v);
	fi= fi-InvTau*(fi-CollideEquillibrium(rho, u, v,EiCollide[i], omegaCollide[i]))+Fx;
}
void CollideLowOrder::Collide_2D_SinglePhase_With_LocalForce_V2(int & i, double &fi,double &rho, double &u, double &v, double & Fi, double & Fx, double & Fy, double & InvTau_tmp){
	fi= fi-InvTau*(fi-CollideEquillibrium(rho, u, v,EiCollide[i], omegaCollide[i]))+Fi;
}
///Collision Single phase step with a body force calculated from User force or model or both
void CollideLowOrder::Collide_2D_SinglePhase_With_BodyForce(int & i, double &fi,double &rho, double &u, double &v, double & Fx, double & Fy, double & InvTau_tmp){

//	Fx=Fx+BodyUserForce_X(0,rho,u,v);
//	Fy=Fy+BodyForce(1,rho,u,v);
	fi=fi-InvTau*(fi-CollideEquillibrium(rho, u, v,EiCollide[i], omegaCollide[i]))+Collide_2D_BodyForce(i, u, v, Fx,Fy);
}
void CollideLowOrder::Collide_2D_SinglePhase_With_BodyForce_V2(int & i, double &fi,double &rho, double &u, double &v, double & Fi, double & Fx, double & Fy, double & InvTau_tmp){
	fi=fi-InvTau*(fi-CollideEquillibrium(rho, u, v,EiCollide[i], omegaCollide[i]))+Collide_2D_BodyForce(i, u, v, Fx,Fy);
}
///Collision Single phase step with a local (directional) force calculated from User force or model or both
void CollideLowOrder::Collide_2D_SinglePhase_With_LocalForceAndUserLocalForce(int & i, double &fi,double &rho, double &u, double &v, double & Fx, double & Fy, double & InvTau_tmp){
	Fx+=LocalForce(i, rho, u, v);
	fi= fi-InvTau*(fi-CollideEquillibrium(rho, u, v,EiCollide[i], omegaCollide[i]))+Fx;
}
///Collision Single phase step with a body force calculated from User force or model or both
void CollideLowOrder::Collide_2D_SinglePhase_With_BodyForceAndUserLocalForce(int & i, double &fi,double &rho, double &u, double &v, double & Fx, double & Fy, double & InvTau_tmp){

	fi=fi-InvTau*(fi-CollideEquillibrium(rho, u, v,EiCollide[i], omegaCollide[i]))+Collide_2D_BodyForce(i, u, v, Fx,Fy);
}
///Collision Single phase step with a local (directional) force calculated from User force or model or both
void CollideLowOrder::Collide_2D_SinglePhase_With_LocalForceAndUserBodyForce(int & i, double &fi,double &rho, double &u, double &v, double & Fx, double & Fy, double & InvTau_tmp){
	Fx+=LocalForce(i, rho, u, v);
	fi= fi-InvTau*(fi-CollideEquillibrium(rho, u, v,EiCollide[i], omegaCollide[i]))+Fx;
}
void CollideLowOrder::Collide_2D_SinglePhase_With_LocalForceAndBodyForce_V2(int & i, double &fi,double &rho, double &u, double &v, double & Fi, double & Fx, double & Fy, double & InvTau_tmp){

	fi= fi-InvTau*(fi-CollideEquillibrium(rho, u, v,EiCollide[i], omegaCollide[i]))+Fi+Collide_2D_BodyForce(i, u, v, Fx,Fy);
}
void CollideLowOrder::No_Collide_2D_SinglePhase_With_LocalForce_V2(int & i, double &fi, double &rho, double &u, double &v, double & Fi, double & Fx, double & Fy, double & InvTau_tmp){
	fi= fi+Fi;
}
void CollideLowOrder::No_Collide_2D_SinglePhase_With_BodyForce_V2(int & i, double &fi, double &rho, double &u, double &v, double & Fi, double & Fx, double & Fy, double & InvTau_tmp){
	fi= fi+Collide_2D_BodyForce(i, u, v, Fx,Fy);
}
void CollideLowOrder::No_Collide_2D_SinglePhase_With_LocalForceAndBodyForce_V2(int & i, double &fi, double &rho, double &u, double &v, double & Fi, double & Fx, double & Fy, double & InvTau_tmp){
	fi= fi+Fi+Collide_2D_BodyForce(i, u, v, Fx,Fy);
}
///Collision Single phase step with a body force calculated from User force or model or both
void CollideLowOrder::Collide_2D_SinglePhase_With_BodyForceAndUserBodyForce(int & i, double &fi,double &rho, double &u, double &v, double & Fx, double & Fy, double & InvTau_tmp){

//	Fx=Fx+BodyUserForce_X(0,rho,u,v);
//	Fy=Fy+BodyForce(1,rho,u,v);
	fi=fi-InvTau*(fi-CollideEquillibrium(rho, u, v,EiCollide[i], omegaCollide[i]))+Collide_2D_BodyForce(i, u, v, Fx,Fy);
}

///Body force function
double CollideLowOrder::Collide_2D_BodyForce(int & i, double &u, double &v, double Fx, double Fy){
	return
			Fx*(1.0-InvTau*0.5)*omegaCollide[i]*InvCs2Collide*((EiCollide[i][0]-u)+InvCs2Collide*EiCollide[i][0]*(EiCollide[i][0]*u+EiCollide[i][1]*v))+
			Fy*(1.0-InvTau*0.5)*omegaCollide[i]*InvCs2Collide*((EiCollide[i][1]-v)+InvCs2Collide*EiCollide[i][1]*(EiCollide[i][0]*u+EiCollide[i][1]*v));

}
///Standard Collision Single phase step with a local Tau (relaxation time which is related to viscosity and/or Knudsen)
void CollideLowOrder::Collide_2D_SinglePhase_Non_Constant_Tau(int & i, double &fi,double &rho, double &u, double &v, double & Fx, double & Fy, double & InvTau_tmp){
	fi= fi-InvTau_tmp*(fi-CollideEquillibrium(rho, u, v,EiCollide[i], omegaCollide[i]));
}
void CollideLowOrder::Collide_2D_SinglePhase_Non_Constant_Tau_V2(int & i, double &fi,double &rho, double &u, double &v, double & Fi, double & Fx, double & Fy, double & InvTau_tmp){
	fi= fi-InvTau_tmp*(fi-CollideEquillibrium(rho, u, v,EiCollide[i], omegaCollide[i]));
}
///Collision Single phase step with a local (directional) force calculated from User force or model or both and with a local Tau (relaxation time which is related to viscosity and/or Knudsen)
void CollideLowOrder::Collide_2D_SinglePhase_Non_Constant_Tau_With_LocalForce(int & i, double &fi,double &rho, double &u, double &v, double & Fx, double & Fy, double & InvTau_tmp){
	fi= fi-InvTau_tmp*(fi-CollideEquillibrium(rho, u, v,EiCollide[i], omegaCollide[i]))+Fx;
}
void CollideLowOrder::Collide_2D_SinglePhase_Non_Constant_Tau_With_LocalForce_V2(int & i, double &fi,double &rho, double &u, double &v, double & Fi, double & Fx, double & Fy, double & InvTau_tmp){
	fi= fi-InvTau_tmp*(fi-CollideEquillibrium(rho, u, v,EiCollide[i], omegaCollide[i]))+Fi;
}
///Collision Single phase step with a body force calculated from User force or model or both and with a local Tau (relaxation time which is related to viscosity and/or Knudsen)
void CollideLowOrder::Collide_2D_SinglePhase_Non_Constant_Tau_With_BodyForce(int & i, double &fi,double &rho, double &u, double &v, double & Fx, double & Fy, double & InvTau_tmp){
	fi=fi-InvTau_tmp*(fi-CollideEquillibrium(rho, u, v,EiCollide[i], omegaCollide[i]))+Collide_2D_BodyForce_Non_Constant_Tau(i, u, v, Fx,Fy, InvTau_tmp);
}
void CollideLowOrder::Collide_2D_SinglePhase_Non_Constant_Tau_With_BodyForce_V2(int & i, double &fi,double &rho, double &u, double &v, double & Fi, double & Fx, double & Fy, double & InvTau_tmp){
	fi=fi-InvTau_tmp*(fi-CollideEquillibrium(rho, u, v,EiCollide[i], omegaCollide[i]))+Collide_2D_BodyForce_Non_Constant_Tau(i, u, v, Fx,Fy, InvTau_tmp);
}
///Collision Single phase step with a local (directional) force calculated from User force or model or both and with a local Tau (relaxation time which is related to viscosity and/or Knudsen)
void CollideLowOrder::Collide_2D_SinglePhase_Non_Constant_Tau_With_LocalForceAndUserLocalForce(int & i, double &fi,double &rho, double &u, double &v, double & Fx, double & Fy, double & InvTau_tmp){
	fi= fi-InvTau_tmp*(fi-CollideEquillibrium(rho, u, v,EiCollide[i], omegaCollide[i]))+LocalForce(i, rho, u, v)+Fx;
}
///Collision Single phase step with a body force calculated from User force or model or both and with a local Tau (relaxation time which is related to viscosity and/or Knudsen)
void CollideLowOrder::Collide_2D_SinglePhase_Non_Constant_Tau_With_BodyForceAndUserLocalForce(int & i, double &fi,double &rho, double &u, double &v, double & Fx, double & Fy, double & InvTau_tmp){
	fi=fi-InvTau_tmp*(fi-CollideEquillibrium(rho, u, v,EiCollide[i], omegaCollide[i]))+Collide_2D_BodyForce_Non_Constant_Tau(i, u, v, Fx,Fy, InvTau_tmp)+LocalForce(i, rho, u, v);
}
///Collision Single phase step with a local (directional) force calculated from User force or model or both and with a local Tau (relaxation time which is related to viscosity and/or Knudsen)
void CollideLowOrder::Collide_2D_SinglePhase_Non_Constant_Tau_With_LocalForceAndUserBodyForce(int & i, double &fi,double &rho, double &u, double &v, double & Fx, double & Fy, double & InvTau_tmp){
	DVec_2D_tmp[0]=BodyForce(0,rho,u,v);
	DVec_2D_tmp[1]=BodyForce(1,rho,u,v);
	fi= fi-InvTau_tmp*(fi-CollideEquillibrium(rho, u, v,EiCollide[i], omegaCollide[i]))+LocalForce(i, rho, u, v)+Collide_2D_BodyForce_Non_Constant_Tau(i, u, v, DVec_2D_tmp[0],DVec_2D_tmp[1], InvTau_tmp);
}
void CollideLowOrder::Collide_2D_SinglePhase_Non_Constant_Tau_With_LocalForceAndBodyForce_V2(int & i, double &fi,double &rho, double &u, double &v, double & Fi, double & Fx, double & Fy, double & InvTau_tmp){
	fi= fi-InvTau_tmp*(fi-CollideEquillibrium(rho, u, v,EiCollide[i], omegaCollide[i]))+Fi+Collide_2D_BodyForce_Non_Constant_Tau(i, u, v, Fx,Fy, InvTau_tmp);
}
///Collision Single phase step with a body force calculated from User force or model or both and with a local Tau (relaxation time which is related to viscosity and/or Knudsen)
void CollideLowOrder::Collide_2D_SinglePhase_Non_Constant_Tau_With_BodyForceAndUserBodyForce(int & i, double &fi,double &rho, double &u, double &v, double & Fx, double & Fy, double & InvTau_tmp){
	fi=fi-InvTau_tmp*(fi-CollideEquillibrium(rho, u, v,EiCollide[i], omegaCollide[i]))+Collide_2D_BodyForce_Non_Constant_Tau(i, u, v, Fx,Fy, InvTau_tmp);
}
void CollideLowOrder::No_Collide_2D_SinglePhase_Non_Constant_Tau_With_LocalForce_V2(int & i, double &fi, double &rho, double &u, double &v, double & Fi, double & Fx, double & Fy, double & InvTau_tmp){
	fi=fi+Fi;
}
void CollideLowOrder::No_Collide_2D_SinglePhase_Non_Constant_Tau_With_BodyForce_V2(int & i, double &fi, double &rho, double &u, double &v, double & Fi, double & Fx, double & Fy, double & InvTau_tmp){
	fi=fi+Collide_2D_BodyForce_Non_Constant_Tau(i, u, v, Fx,Fy, InvTau_tmp);
}
void CollideLowOrder::No_Collide_2D_SinglePhase_Non_Constant_Tau_With_LocalForceAndBodyForce_V2(int & i, double &fi, double &rho, double &u, double &v, double & Fi, double & Fx, double & Fy, double & InvTau_tmp){
	fi=fi+Fi+Collide_2D_BodyForce_Non_Constant_Tau(i, u, v, Fx,Fy, InvTau_tmp);
}
///Body force function with a local Tau (relaxation time which is related to viscosity and/or Knudsen)
double CollideLowOrder::Collide_2D_BodyForce_Non_Constant_Tau(int & i, double &u, double &v, double Fx, double Fy, double & InvTau_tmp){
	return
			Fx*(1.0-InvTau_tmp*0.5)*omegaCollide[i]*InvCs2Collide*((EiCollide[i][0]-u)+InvCs2Collide*EiCollide[i][0]*(EiCollide[i][0]*u+EiCollide[i][1]*v))+
			Fy*(1.0-InvTau_tmp*0.5)*omegaCollide[i]*InvCs2Collide*((EiCollide[i][1]-v)+InvCs2Collide*EiCollide[i][1]*(EiCollide[i][0]*u+EiCollide[i][1]*v));

}
/*CollideD2Q9Colour::CollideD2Q9Colour(){
//	EiCollide=new double [9][2];
	Ak=0.01;
}
CollideD2Q9Colour::~CollideD2Q9Colour(){

}
void CollideD2Q9Colour::Collide_ColorFluid(int & direction, double & fi,double &rho,double*  F,double & F_Norm, double & InvTau_, double &u, double &v)
{
	InvTau=InvTau_;
//	double F_Norm=std::sqrt(F[0]*F[0]+F[1]+F[1]);
	//f(0) => EiCollide=0
//	f[0]= f[0]-InvTau*(f[0]-rho*omegaCollide*(1+-1.5*(u*u+v*v)));
//	f[0]+=-Ak*0.5*F_Norm*3/4;

		CollideLowOrder::Collide_SinglePhase(fi,rho, u, v, EiCollide[direction], omegaCollide[direction]);
		fi+=CollideD2Q9Colour::TwoPhase_Collision_operator(direction, F, F_Norm);

}*/
/*void CollideD2Q9Colour::Collide_ColorFluid(double &fr,double &fb,double &rho_r,double &rho_b, double &u, double &v, double *e_i, double &omegaCollide){
	fr= fr-InvTau*(fr-EquiDistriFunct2D(rho_r, u, v,e_i, omegaCollide));
	fb= fb-InvTau*(fb-EquiDistriFunct2D(rho_r, u, v,e_i, omegaCollide));
}*/
