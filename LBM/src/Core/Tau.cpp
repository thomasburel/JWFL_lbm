/*
 * Tau.cpp
 *
 *  Created on: 15 Aug 2016
 *      Author: Thomas Burel
 */

#include "Tau.h"

Tau::Tau() {
	// TODO Auto-generated constructor stub
	tau=0;
	RhoTau=1;RhoNTau=0;
	TTau=20;
}
void Tau::IniTau(Parameters *Param) {
	// TODO Auto-generated constructor stub
	if(tau!=0)
		delete tau;
	tau=new BGKCollision();
	tau->Init_Tau(Param->Get_Tau(),Param->Get_cs2(),Param->Get_deltaT());
	tau->Init_Viscosity(Param);
	/*
	if(Param->Get_Model()==SolverEnum::SinglePhase ||Param->Get_ViscosityType()==ConstViscosity)
	{
		tau=new ConstTau();
		tau->Init_Tau(Param->Get_Tau(),Param->Get_cs2(),Param->Get_deltaT());
	}
	else if(Param->Get_ViscosityType()==LinearViscosity)
	{
		tau=new LinearTau();
		tau->Init_Tau(Param->Get_Tau(),Param->Get_cs2(),Param->Get_deltaT());
		tau->Set_Tau(Param->Get_Tau_1(),Param->Get_Tau_2());
	}
	else
	{
		tau=new HarmonicTau();
		tau->Init_Tau(Param->Get_Tau(),Param->Get_cs2(),Param->Get_deltaT());
		tau->Set_mu(Param->Convert_NutoMu(Param->Convert_TauToNu(Param->Get_Tau_1()),Param->Get_Rho_1()),Param->Convert_NutoMu(Param->Convert_TauToNu(Param->Get_Tau_2()),Param->Get_Rho_2()));
	}*/
}
Tau::~Tau() {
	delete tau;
}
double& Tau::Get_InvTau(double const Rho, double const RhoN){
return tau->Get_InvTau(Rho,RhoN);
}
/*double*& Tau::Get_InvTau(){
	return tau->Get_InvTau(RhoTau,RhoNTau);
}*/
TauDEF::TauDEF() {
	InvTau=1;
	cs2=1;
	deltaT=1;
//	visco=0;
}

TauDEF::~TauDEF() {
}

void TauDEF::Init_Viscosity(Parameters *Param) {

	visco.Set_Viscosity(Param);
}
BGKCollision::BGKCollision(){

}
BGKCollision::~BGKCollision(){

}
double& BGKCollision::Get_InvTau(double const &Rho, double const &RhoN){
	InvTau=1.0/Convert_NuToTau(visco.Get_Nu(Rho,RhoN));
	return InvTau;
}
/*
ConstTau::ConstTau() {
	InvTau=1;
}

ConstTau::~ConstTau() {
}
HarmonicTau::HarmonicTau(){
	mu2_1=0.01;
	mu2_2=mu2_1;
}

HarmonicTau::~HarmonicTau(){

}
double& HarmonicTau::Get_InvTau(double const &Rho, double const &RhoN){
	InvTau=1.0/CalculTau(Rho,RhoN);
	return InvTau;}

double HarmonicTau::CalculTau(double const &Rho, double const &RhoN){
	return Convert_NuToTau(Convert_MutoNu(CalculHarmonicVisco(RhoN),Rho));
}
double HarmonicTau::CalculHarmonicVisco(double const &RhoN){
	return 1.0/((1.0+RhoN)/mu2_1+(1.0-RhoN)/mu2_2);
}

LinearTau::LinearTau(){
	Tau_1=1;
	Tau_2=Tau_1;
}

LinearTau::~LinearTau(){

}
double& LinearTau::Get_InvTau(double const &Rho, double const &RhoN){
	InvTau=1.0/CalculTau(RhoN);
	return InvTau;}

double LinearTau::CalculTau(double const &RhoN){
	return 0.5*((1.0+RhoN)*Tau_1+(1.0-RhoN)*Tau_2);
}
*/
