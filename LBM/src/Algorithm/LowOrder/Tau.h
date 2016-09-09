/*
 * Tau.h
 *
 *  Created on: 15 Aug 2016
 *      Author: Thomas Burel
 */

#ifndef SRC_ALGORITHM_LOWORDER_TAU_H_
#define SRC_ALGORITHM_LOWORDER_TAU_H_

class TauDEF {
public:
	TauDEF();
	TauDEF(double InvTau_);
	virtual ~TauDEF();
	virtual double Get_Tau()=0;

protected:
	double InvTau;
};

class ConstTau: public TauDEF {
public:
	ConstTau();
	ConstTau(double InvTau_);
	virtual ~ConstTau();
	virtual double Get_Tau(){return InvTau;}
};

class Tau {
public:
	Tau();
	Tau(double InvTau_);
	virtual ~Tau();
	double Get_Tau();

private:
	TauDEF* tau;
};




#endif /* SRC_ALGORITHM_LOWORDER_TAU_H_ */
