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

Tau::Tau(double InvTau_) {
	// TODO Auto-generated constructor stub
	tau=new ConstTau(InvTau_);
}

Tau::~Tau() {
	delete tau;
}

TauDEF::TauDEF() {
	InvTau=1;
}
TauDEF::TauDEF(double InvTau_) {
	InvTau=InvTau_;
}
TauDEF::~TauDEF() {
}

ConstTau::ConstTau() {
	InvTau=1;
}

ConstTau::ConstTau(double InvTau_) {
	InvTau=InvTau_;
}
ConstTau::~ConstTau() {
}
