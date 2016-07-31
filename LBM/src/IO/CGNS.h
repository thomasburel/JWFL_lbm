/*
 * CGNS.h
 *
 *  Created on: 8 May 2015
 *      Author: thomas
 */

#ifndef IO_CGNS_H_
#define IO_CGNS_H_
#include "WriterManager.h"
#include "pcgnslib.h"
#include "mpi.h"
#include <cstdio>
#include <cstring>
#include "../Core/Parameters.h" //type dimension
#include <sstream>      // std::stringstream


/*
 *
 */
class CGNS : public WriterManager{
public:
	CGNS();
	CGNS(dimension dimension_,std::string outputfilename_,int tot_nodes_, int tot_elems_, int start_nodes_, int end_nodes_, int start_elems_, int end_elems_, const double *x_, const double *y_, const double *z_,int *e_);
	virtual ~CGNS();
	virtual CGNS& Get_class(){return *this;};
	virtual void Write_Output(int& Nstep);
	virtual void Write_breakpoint(Parameters &Parameters);
	virtual void Set_solution(double **d_, std::string *str, int nbvar);
	virtual void Set_breakpoint(double **b_, std::string *str, int nbvar);
	void Write_nodeproperties(Parameters &Parameters);///To keep the connectivity, the type of node with their specific property...
	void Write_ParametersProperties(Parameters &Parameters);
private:
    int F, B, Z, E, S, Fs, A, Cx, Cy, Cz;
    const double* x, *y, *z;
    double  **d,**b;
    cgsize_t sizes[3], *e, start_nodes, end_nodes, start_elems, end_elems, ncells;


private:
	friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
    	ar & boost::serialization::base_object<WriterManager>(*this);
       ar & F & B & Z & E & S & Fs & Cx & Cy & Cz ;
   //    ar & x & y & z;

      ar & sizes & start_nodes & end_nodes & start_elems & end_elems & ncells;
/*       ar & outfile;
       ar & NbVariableOutput & NbVariableBreakpoint;
       ar & VariableOutput & VariableBreakpoint;
       ar & Dimension;
       ar & outputfilename & Breakpointfilename;*/
    }

};

#endif /* IO_CGNS_H_ */
