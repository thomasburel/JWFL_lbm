/*
 * ============================================================================
 * Viscosity.cpp
 *
 *  Created on: 5 Apr 2017
 *      Author: Thomas Burel
 *      Version     :
 *      Copyright   : Your copyright notice
 *      Description : LBM Code Developing by JWFL at University of Strathclyde
 * ============================================================================
 */

#include "Viscosity.h"

Viscosity::Viscosity() {
	// TODO Auto-generated constructor stub
	visco=0;
}

Viscosity::~Viscosity() {
	// TODO Auto-generated destructor stub
}
void Viscosity::Set_Viscosity(Parameters *Param){
	if(visco!=0)
		delete visco;
	if(Param->Get_Model()==SolverEnum::SinglePhase ||Param->Get_ViscosityType()==ConstViscosity)
	{
		visco=new ConstVisco();
		visco->Set_mu(Param->Convert_NutoMu(Param->Convert_TauToNu(Param->Get_Tau_1()),Param->Get_Rho_1()));
	}
	else if(Param->Get_ViscosityType()==LinearViscosity)
	{
		visco=new LinearVisco();
		visco->Set_mu(Param->Convert_NutoMu(Param->Convert_TauToNu(Param->Get_Tau_1()),Param->Get_Rho_1()),Param->Convert_NutoMu(Param->Convert_TauToNu(Param->Get_Tau_2()),Param->Get_Rho_2()));
	}
	else
	{
		visco=new HarmonicVisco();
		visco->Set_mu(Param->Convert_NutoMu(Param->Convert_TauToNu(Param->Get_Tau_1()),Param->Get_Rho_1()),Param->Convert_NutoMu(Param->Convert_TauToNu(Param->Get_Tau_2()),Param->Get_Rho_2()));
	}
}
double Viscosity::Get_Nu(double const Rho, double const RhoN){
	return visco->Get_Nu(Rho,RhoN);
}
double Viscosity::Get_Mu(double const Rho, double const RhoN){
	return visco->Get_Mu(Rho,RhoN);
}
double Viscosity::Get_Mu_1(){
	return visco->Get_Mu_1();
}
double Viscosity::Get_Mu_2(){
	return visco->Get_Mu_2();
}
ViscosityDEF::ViscosityDEF() {

}

ViscosityDEF::~ViscosityDEF() {
}

ConstVisco::ConstVisco() {
	mu_1=0.03;
}

ConstVisco::~ConstVisco() {
}
HarmonicVisco::HarmonicVisco(){
	mu2_1=0.03;
	mu2_2=mu2_1;
}

HarmonicVisco::~HarmonicVisco(){

}

double HarmonicVisco::CalculHarmonicVisco(double const &RhoN){
	return 1.0/((1.0+RhoN)/mu2_1+(1.0-RhoN)/mu2_2);
}

LinearVisco::LinearVisco(){
	mu_1=0.03;
	mu_2=mu_1;
}

LinearVisco::~LinearVisco(){

}
double LinearVisco::CalculLinearVisco(double const &RhoN){
	return 0.5*((1.0+RhoN)*mu_1+(1.0-RhoN)*mu_2);
}

