/*
 * WriterManager.h
 *
 *  Created on: 10 May 2015
 *      Author: thomas
 */

#ifndef IO_WRITERMANAGER_H_
#define IO_WRITERMANAGER_H_
/*
 *
 */
//#include <string>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/list.hpp>
#include <boost/serialization/assume_abstract.hpp>
#include "../Core/Parameters.h"
class WriterManager {
public:
	WriterManager();
	std::string& Get_classtype() {return classtype;};
	void UpdateFileNames(std::string NewFileName);
//	virtual CGNS& Get_class(){return *this;}=0;
	virtual ~WriterManager();
	virtual void Write_Output(int& Nstep)=0;
	virtual void Write_breakpoint(Parameters &Parameters)=0;
	virtual void Set_solution(double **d_, std::string *str, int nbvar)=0;
	virtual void Set_breakpoint(double **d_, std::string *str, int nbvar)=0;

protected:
	int tot_nnodes, tot_nelems, nnodes, nelems; //total in the domain and per process
    char *outfile; //temporary variable to write cgns => no archive
    int NbVariableOutput,NbVariableBreakpoint;
    std::string *VariableOutput,*VariableBreakpoint;
	bool Dimension;
    std::string outputfilename,Breakpointfilename;
    std::string classtype;

private:
	friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
      ar & tot_nnodes & tot_nelems & nnodes & nelems ;
      ar & NbVariableOutput & NbVariableBreakpoint;
      for (int i=0;i<NbVariableOutput;i++)
    	  ar & VariableOutput[i];
      for (int i=0;i<NbVariableBreakpoint;i++)
    	  ar & VariableBreakpoint[i];
      ar & Dimension;
      ar & outputfilename & Breakpointfilename;
      ar & classtype;
    }
};

#endif /* IO_WRITERMANAGER_H_ */
