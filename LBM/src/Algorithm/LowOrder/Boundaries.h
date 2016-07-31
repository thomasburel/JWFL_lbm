/*
 * Boundaries.h
 *
 *  Created on: 17 Jul 2015
 *      Author: thomas
 */

#ifndef ALGORITHM_LOWORDER_BOUNDARIES_H_
#define ALGORITHM_LOWORDER_BOUNDARIES_H_
//#include <boost/math/special_functions/sign.hpp>
#include <cmath>        // std::abs
#include <iostream>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
/*
#include <fstream>
#include "mpi.h"*/
class Boundaries {
public:
	Boundaries();
	virtual ~Boundaries();
	void BC_HeZou_U(int & NormalBc, double* fi, double U, double V);
	void BC_HeZou_P(int & NormalBc, double* fi, double Rho, double & U, double & V);
	void BC_corner(int & NormalCorner, double* fi, double Rho, double U, double V);
	void BC_corner_no_vel_concave(int & NormalCorner, double* fi, double Rho);
	void BC_corner_no_vel_convex(int & NormalCorner, double* fi);


private:
	void FUNC_HeZou_U (double & a,double & b,double & c,double & d,double & e,double & f,double & g,double & h,double & i,double & U,double & V);
	void FUNC_HeZou_P (double & a,double & b,double & c,double & d,double & e,double & f,double & g,double & h,double & i,double & V,double & Rho);
	void FUNC_corner (double & a,double & b,double & c,double & d,double & e,double & f,double & g,double & h,double & i,double & U,double & V,double & Rho);
	void FUNC_corner_no_vel_concave (double & a,double & b,double & c,double & d,double & e,double & f,double & g,double & h,double & i,double & Rho);
	void FUNC_corner_no_vel_convex (double & a,double & b,double & c,double & d,double & e,double & f);

protected:
	double rhodiff;
	double SumWeightS,SumWeightE,SumWeightN,SumWeightW;
	double SumWeightConcaveSE,SumWeightConcaveNE,SumWeightConcaveNW,SumWeightConcaveSW;
	double SumWeightConvexSE,SumWeightConvexNE,SumWeightConvexNW,SumWeightConvexSW;

private:
		double ru1,ru2;
		double q23,q16,q13;

		/*int rank;
		std::ofstream myFlux;*/

private:
	friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
       ar & ru1 & ru2;
       ar & q23 & q16 & q13;
    }
};

#endif /* ALGORITHM_LOWORDER_BOUNDARIES_H_ */
