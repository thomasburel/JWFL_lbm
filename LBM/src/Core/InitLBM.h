/*
 * InitLBM.h
 *
 *  Created on: 5 May 2015
 *      Author: thomas
 */

#ifndef CORE_INITLBM_H_
#define CORE_INITLBM_H_
#include "../Parallelism/MpiManager.h"
#include "mpi.h"
#include "Parameters.h"
#include "../User/UserInit.h"
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/list.hpp>
#include <boost/serialization/assume_abstract.hpp>
/*
 *
 */

class InitLBM:
		public MpiManager,
		public UserInit
		{
public:
	InitLBM();
	virtual ~InitLBM();
	void Set_Parameters(Parameters *Parameters_);
	void IniMPI(ParallelManager* parallel_,int *argc, char ***argv, bool verbous_=false);
	void IniDomain(int rank,Node2D & Node,int elem, int nodenumber, double* pos,double& Rho, double* U,double alpha=1);


private:
	Parameters *PtrParameters;

private:
	friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
       ar & boost::serialization::base_object<UserInit>(*this) & PtrParameters;
    }
};


namespace serialize{
namespace text{
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>


/*void Save_IniLBM(InitLBM & object){
	//std::string filename(boost::archive::tmpdir());
	std::string filename("Init_save.txt");
	//filename = "/Init_save.txt";
	std::cout<< "the object InitLBM is saving in the text file: "<<filename<<std::endl;
	// make an archive
	std::ofstream ofs(filename.c_str());
	//assert(ofs.good());
	boost::archive::text_oarchive oa(ofs);
	oa << object;
}*/
/*void Save_IniLBM(InitLBM & object,std::string &filename){
	//std::string filename(boost::archive::tmpdir());
	std::cout<< "the object InitLBM is saving in the text file: "<<filename<<std::endl;
	// make an archive
	std::ofstream ofs(filename.c_str());
	//assert(ofs.good());
	boost::archive::text_oarchive oa(ofs);
	oa << object;
}*/
}
}

#endif /* CORE_INITLBM_H_ */
