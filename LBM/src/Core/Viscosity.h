/*
 * ============================================================================
 * Viscosity.h
 *
 *  Created on: 5 Apr 2017
 *      Author: Thomas Burel
 *      Version     :
 *      Copyright   : Your copyright notice
 *      Description : LBM Code Developing by JWFL at University of Strathclyde
 * ============================================================================
 */

#ifndef SRC_CORE_VISCOSITY_H_
#define SRC_CORE_VISCOSITY_H_
#include "Parameters.h"

class ViscosityDEF {
public:
	ViscosityDEF();
	virtual ~ViscosityDEF();
	//For constant viscosity -> Update the local variable in the collision class
	virtual void Set_mu(double mu1_,double mu2_=0)=0;
	virtual double Get_Nu(double const &Rho, double const &RhoN)=0;
	virtual double Get_Mu(double const &Rho, double const &RhoN)=0;
inline	double Convert_MuToNu(double mu_,double rho){return mu_/rho;};
inline	double Convert_NuToMu(double nu_,double rho){return nu_*rho;};
	virtual double Get_Mu_1()=0;
	virtual double Get_Mu_2()=0;
};

class ConstVisco: public ViscosityDEF {
public:
	ConstVisco();
	virtual ~ConstVisco();
	virtual void Set_mu(double mu1_,double mu2_=0){mu_1=mu1_;};
	virtual double Get_Nu(double const &Rho, double const &RhoN){return Convert_MuToNu(mu_1,Rho);};
	virtual double Get_Mu(double const &Rho, double const &RhoN){return mu_1;};
	virtual double Get_Mu_1(){return mu_1;};
	virtual double Get_Mu_2(){return 0;};
private:
	double mu_1;//kinematic viscosity
};

class HarmonicVisco: public ViscosityDEF {
public:
	HarmonicVisco();
	virtual ~HarmonicVisco();
	virtual void Set_mu(double mu1_,double mu2_=0){mu2_1=2.0*mu1_;mu2_2=2.0*mu2_;};
	virtual double Get_Nu(double const &Rho, double const &RhoN){return Convert_MuToNu(CalculHarmonicVisco(RhoN),Rho);};
	virtual double Get_Mu(double const &Rho, double const &RhoN){return CalculHarmonicVisco(RhoN);};
	virtual double Get_Mu_1(){return mu2_1*0.5;};
	virtual double Get_Mu_2(){return mu2_2*0.5;};
private:
	double CalculHarmonicVisco(double const &RhoN);
	double mu2_1,mu2_2;//two times kinematic viscosity

};

class LinearVisco: public ViscosityDEF {
public:
	LinearVisco();
	virtual ~LinearVisco();
	virtual void Set_mu(double mu1_,double mu2_=0){mu_1=mu1_;mu_2=mu2_;};
	virtual double Get_Nu(double const &Rho, double const &RhoN){return Convert_MuToNu(CalculLinearVisco(RhoN),Rho);};
	virtual double Get_Mu(double const &Rho, double const &RhoN){return CalculLinearVisco(RhoN);};
	virtual double Get_Mu_1(){return mu_1;};
	virtual double Get_Mu_2(){return mu_2;};
private:
	double CalculLinearVisco(double const &RhoN);
	double mu_1,mu_2;//kinematic viscosities


};

class Viscosity {
public:
	Viscosity();
	virtual ~Viscosity();

	void Set_Viscosity(Parameters *Param);
	double Get_Nu(double const Rho=1, double const RhoN=0);
	double Get_Mu(double const Rho=1, double const RhoN=0);
	double Get_Mu_1();
	double Get_Mu_2();
private:
	ViscosityDEF *visco;
};

#endif /* SRC_CORE_VISCOSITY_H_ */
