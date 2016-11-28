/*
 * Tau.h
 *
 *  Created on: 15 Aug 2016
 *      Author: Thomas Burel
 */

#ifndef SRC_ALGORITHM_LOWORDER_TAU_H_
#define SRC_ALGORITHM_LOWORDER_TAU_H_

#include "Parameters.h"

class TauDEF {
public:
	TauDEF();
	virtual ~TauDEF();
	//For constant viscosity -> Update the local variable in the collision class
	double& Get_InvTau(){return InvTau;}
	virtual double& Get_InvTau(double const &Rho, double const &RhoN)=0;
	void Init_Tau(double Tau_,double cs2_,double deltaT_){InvTau=1.0/Tau_;cs2=cs2_;deltaT=deltaT_;};
	virtual void Set_mu(double mu1_,double mu2_=0)=0;
	virtual void Set_Tau(double Tau1_,double Tau2_=1)=0;
inline	double Convert_TauToNu(double Tau_){return (Tau_-0.5)*cs2*deltaT;};
inline	double Convert_MutoNu(double mu_,double rho){return mu_/rho;};
inline	double Convert_NutoMu(double nu_,double rho){return nu_*rho;};
inline	double Convert_NuToTau(double nu_){return nu_/(cs2*deltaT)+0.5;};
protected:
	double InvTau;
	double cs2,deltaT;
};

class ConstTau: public TauDEF {
public:
	ConstTau();
	virtual ~ConstTau();
	virtual double& Get_InvTau(double const &Rho, double const &RhoN){return InvTau;};
	virtual void Set_mu(double mu1_,double mu2_=0){};
	virtual void Set_Tau(double Tau1_,double Tau2_=1){InvTau=1.0/Tau1_;};
};

class HarmonicVisco: public TauDEF {
public:
	HarmonicVisco();
	virtual ~HarmonicVisco();
	virtual double& Get_InvTau(double const &Rho, double const &RhoN);
	virtual void Set_mu(double mu1_,double mu2_=0){mu2_1=2.0*mu1_;mu2_2=2.0*mu2_;};
	virtual void Set_Tau(double Tau1_,double Tau2_=1){};
private:
	double CalculTau(double const &Rho, double const &RhoN);
	double CalculHarmonicVisco(double const &RhoN);
	double mu2_1,mu2_2;//two times kinematic viscosity

};

class LinearVisco: public TauDEF {
public:
	LinearVisco();
	virtual ~LinearVisco();
	virtual double& Get_InvTau(double const &Rho, double const &RhoN);
	virtual void Set_mu(double mu1_,double mu2_=0){};
	virtual void Set_Tau(double Tau1_,double Tau2_=1){InvTau=1.0/Tau1_;Tau_1=Tau1_;Tau_2=Tau2_;};
private:
	double CalculTau(double const &RhoN);
	double Tau_1,Tau_2;//two times kinematic viscosity

};


class Tau {
public:
	Tau();
	void IniTau(Parameters *Param);
	virtual ~Tau();
	double& Get_InvTau(double const Rho=1, double const RhoN=0);
private:
	TauDEF* tau;
};




#endif /* SRC_ALGORITHM_LOWORDER_TAU_H_ */
