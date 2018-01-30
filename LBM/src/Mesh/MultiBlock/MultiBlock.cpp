/*
 * MultiBlock.cpp
 *
 *  Created on: 5 May 2015
 *      Author: thomas
 */

#include "MultiBlock.h"

MultiBlock::MultiBlock() :
PtrParameters(0),parallel(0)
{
	ndims=2;
	// TODO Auto-generated constructor stub

}

MultiBlock::~MultiBlock() {
	// TODO Auto-generated destructor stub
}
bool MultiBlock::Check_CleaningMesh(int nbTotalSolidRemoved, int nbTotalSolidadded){
	int maxnbTotalSolidRemoved=0,maxnbTotalSolidadded=0;
	MPI_Allreduce(&nbTotalSolidRemoved,&maxnbTotalSolidRemoved,1, MPI_INT , MPI_SUM ,parallel->getGlobalCommunicator());
	MPI_Allreduce(&nbTotalSolidadded,&maxnbTotalSolidadded,1, MPI_INT , MPI_SUM ,parallel->getGlobalCommunicator());
//	std::cout<<"processor: "<<parallel->getRank()<<" sum of solid removed: "<<maxnbTotalSolidRemoved<<" sum of solid added: "<<maxnbTotalSolidadded<<std::endl;
	if(maxnbTotalSolidRemoved>0||maxnbTotalSolidadded>0)
		return true;
	else
		return false;

}
