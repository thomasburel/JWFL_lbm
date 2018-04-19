/*
 * UserParameters.h
 *
 *  Created on: 10 Jun 2015
 *      Author: thomas
 */

#ifndef USER_USERPARAMETERS_H_
#define USER_USERPARAMETERS_H_
#include <boost/archive/xml_oarchive.hpp>
#include <boost/serialization/nvp.hpp>
class UserParameters {
	private:
	friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_NVP(Umax);
        ar & BOOST_SERIALIZATION_NVP(H);
        ar & BOOST_SERIALIZATION_NVP(L);
        ar & BOOST_SERIALIZATION_NVP(Pmax);
        ar & BOOST_SERIALIZATION_NVP(Pmin);
    }
public:
	UserParameters();
	virtual ~UserParameters();
	double Get_UserH()const {return H;};
	double Get_UserL()const {return L;};
	void Get_UserParameters(double & Umax_,double & H_,double & L_,double & Pmax_,double & Pmin_);
	void Set_UserParameters(double  Umax_,double  H_,double  L_,double  Pmax_,double  Pmin_);
	void Set_UserDroplets(double angle, double sigma, double diameter, double Re, double Ca);
	void Get_UserDroplets(double & angle, double & sigma, double & diameter, double & Re, double & Ca);
	void Set_ReUser(double Re){ReUser=Re;};
	double Get_ReUser(){return ReUser;};
	void Set_HchannelUser(double H_channel){Hchannel=H_channel;};
	double Get_HchannelUser(){return Hchannel;};
	void Set_ViscoUser(double viscosity){visco=viscosity;};
	double Get_ViscoUser(){return visco;};
private:
	double Umax,H,L,Pmax,Pmin;
	double ReUser,CaUser,DiameterUser,sigmaUser,ContactAngleUser;
	double Hchannel;
	double visco;

};

#endif /* USER_USERPARAMETERS_H_ */
