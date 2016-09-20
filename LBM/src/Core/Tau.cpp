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
Tau::Tau(Parameters *Param) {
	// TODO Auto-generated constructor stub
	if(Param->Get_Model()==SinglePhase)
	{
		tau=new ConstTau();
		tau->Init_Tau(Param->Get_Tau(),Param->Get_cs2(),Param->Get_deltaT());
	}
	else
	{
		tau=new HarmonicViscosity();
		tau->Init_Tau(Param->Get_Tau(),Param->Get_cs2(),Param->Get_deltaT());
		tau->Set_mu(Param->Convert_NutoMu(Param->Convert_TauToNu(Param->Get_Tau_1()),Param->Get_Rho_1()),Param->Convert_NutoMu(Param->Convert_TauToNu(Param->Get_Tau_2()),Param->Get_Rho_2()));
	}
}
Tau::~Tau() {
	delete tau;
}
double Tau::Get_InvTau(double const Rho, double const RhoN){
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
HarmonicViscosity::HarmonicViscosity(){
	mu1=0.01;
	mu2=mu1;
}

HarmonicViscosity::~HarmonicViscosity(){

}
double HarmonicViscosity::Get_InvTau(double const &Rho, double const &RhoN){
	return 1/CalculTau(Rho,RhoN);}

double HarmonicViscosity::CalculTau(double const &Rho, double const &RhoN){
	Convert_NuToTau(Convert_MutoNu(CalculHarmonicViscosity(RhoN),Rho));
	return Convert_NuToTau(Convert_MutoNu(CalculHarmonicViscosity(RhoN),Rho));
}
double HarmonicViscosity::CalculHarmonicViscosity(double const &RhoN){
	return 1/((1+RhoN)/mu1+(1-RhoN)/mu2);
}
