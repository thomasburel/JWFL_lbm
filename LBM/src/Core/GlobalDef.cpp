/*
 * GlobalDef.cpp
 *
 *  Created on: 12 Jun 2015
 *      Author: thomas
 */

#include "GlobalDef.h"

DistriFunct::DistriFunct()
:f(0),NbNodes(1),NbVelocities(1)
{
	// TODO Auto-generated constructor stub
	f=new double* [NbVelocities];
	for (int i=0;i<NbVelocities;i++)
	{
		f[i]=new double [NbNodes];
		for (int j=0;j<NbNodes;j++)
			f[i][j]=0;
	}
}

DistriFunct::DistriFunct(int NbNodes_, int NbVelocities_)
:f(0),NbNodes(NbNodes_),NbVelocities(NbVelocities_)
{
	f=new double* [NbVelocities];
	for (int i=0;i<NbVelocities;i++)
	{
		f[i]=new double [NbNodes];
		for (int j=0;j<NbNodes;j++)
			f[i][j]=0;
	}
}

DistriFunct::~DistriFunct() {
	for (int i=0;i<NbVelocities;i++)
	{
			delete [] f[i];
	}
	delete [] f;
}

