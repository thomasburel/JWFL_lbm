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
inline	double Convert_MuToNu(double mu_,double rho){return mu_/rho;};
inline	double Convert_NuToMu(double nu_,double rho){return nu_*rho;};

};

class ConstVisco: public ViscosityDEF {
public:
	ConstVisco();
	virtual ~ConstVisco();
	virtual void Set_mu(double mu1_,double mu2_=0){mu_1=mu1_;};
	virtual double Get_Nu(double const &Rho, double const &RhoN){return Convert_MuToNu(mu_1,Rho);};
	double mu_1;//kinematic viscosity
};

class HarmonicVisco: public ViscosityDEF {
public:
	HarmonicVisco();
	virtual ~HarmonicVisco();
	virtual void Set_mu(double mu1_,double mu2_=0){mu2_1=2.0*mu1_;mu2_2=2.0*mu2_;};
	virtual double Get_Nu(double const &Rho, double const &RhoN){return Convert_MuToNu(CalculHarmonicVisco(RhoN),Rho);};
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
private:
	double CalculLinearVisco(double const &RhoN);
	double mu_1,mu_2;//kinematic viscosities


};

class Viscosity {
public:
	Viscosity();
	virtual ~Viscosity();

	void Set_Viscosity(Parameters *Param);
	double Get_Nu(double const &Rho, double const &RhoN);

private:
	ViscosityDEF *visco;
};

#endif /* SRC_CORE_VISCOSITY_H_ */
