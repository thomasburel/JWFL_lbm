/*
 * TecplotIO.h
 *
 *  Created on: 26 Jul 2016
 *      Author: thomas
 */

#ifndef IO_TECPLOTIO_H_
#define IO_TECPLOTIO_H_

#include "WriterManager.h"
//#include "TECIO.h"
#include "mpi.h"
#include <cstdio>
#include <cstring>
#include "../Core/Parameters.h" //type dimension
#include <sstream>      // std::stringstream
class TecplotIO: public WriterManager {
public:
	TecplotIO();
	virtual ~TecplotIO();
	virtual void Write_Output(int& Nstep);
	virtual void Read_data(double * &d_, std::string variablename, std::string filename);
	virtual void Write_breakpoint(Parameters &Parameters);
	virtual void Set_solution(double **d_, std::string *str, int nbvar);
	virtual void Set_breakpoint(double **b_, std::string *str, int nbvar);
};

#endif /* IO_TECPLOTIO_H_ */
