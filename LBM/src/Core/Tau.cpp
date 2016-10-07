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
}
void Tau::IniTau(Parameters *Param) {
	// TODO Auto-generated constructor stub
	if(Param->Get_Model()==SinglePhase ||Param->Get_ViscosityType()==ConstViscosity)
	{
		tau=new ConstTau();
		tau->Init_Tau(Param->Get_Tau(),Param->Get_cs2(),Param->Get_deltaT());
	}
	else if(Param->Get_ViscosityType()==LinearViscosity)
	{
		tau=new LinearVisco();
		tau->Init_Tau(Param->Get_Tau(),Param->Get_cs2(),Param->Get_deltaT());
		tau->Set_Tau(Param->Get_Tau_1(),Param->Get_Tau_2());
	}
	else
	{
		tau=new HarmonicVisco();
		tau->Init_Tau(Param->Get_Tau(),Param->Get_cs2(),Param->Get_deltaT());
		tau->Set_mu(Param->Convert_NutoMu(Param->Convert_TauToNu(Param->Get_Tau_1()),Param->Get_Rho_1()),Param->Convert_NutoMu(Param->Convert_TauToNu(Param->Get_Tau_2()),Param->Get_Rho_2()));
	}
}
Tau::~Tau() {
	delete tau;
}
double& Tau::Get_InvTau(double const Rho, double const RhoN){
return tau->Get_InvTau(Rho,RhoN);
}
TauDEF::TauDEF() {
	InvTau=1;
	cs2=1;
	deltaT=1;
}

TauDEF::~TauDEF() {
}

ConstTau::ConstTau() {
	InvTau=1;
}

ConstTau::~ConstTau() {
}
HarmonicVisco::HarmonicVisco(){
	mu2_1=0.01;
	mu2_2=mu2_1;
}

HarmonicVisco::~HarmonicVisco(){

}
double& HarmonicVisco::Get_InvTau(double const &Rho, double const &RhoN){
	InvTau=1.0/CalculTau(Rho,RhoN);
	return InvTau;}

double HarmonicVisco::CalculTau(double const &Rho, double const &RhoN){
	return Convert_NuToTau(Convert_MutoNu(CalculHarmonicVisco(RhoN),Rho));
}
double HarmonicVisco::CalculHarmonicVisco(double const &RhoN){
	return 1.0/((1.0+RhoN)/mu2_1+(1.0-RhoN)/mu2_2);
}

LinearVisco::LinearVisco(){
	Tau_1=1;
	Tau_2=Tau_1;
}

LinearVisco::~LinearVisco(){

}
double& LinearVisco::Get_InvTau(double const &Rho, double const &RhoN){
	InvTau=1.0/CalculTau(RhoN);
	return InvTau;}

double LinearVisco::CalculTau(double const &RhoN){
	return 0.5*((1.0+RhoN)*Tau_1+(1.0-RhoN)*Tau_2);
}

