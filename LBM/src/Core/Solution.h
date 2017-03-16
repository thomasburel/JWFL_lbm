/*
 * Solution.h
 *
 *  Created on: 30 Apr 2015
 *      Author: thomas
 */

#ifndef CORE_SOLUTION_H_
#define CORE_SOLUTION_H_
//#include "../Mesh/SingleBlock/Node2D.h"
#include "Dictionary.h"
#include "../Mesh/MultiBlock.h"
#include "../Parallelism/ParallelManager.h"
#include "../IO/CGNS.h"
#include <iostream>

// include this header to serialize vectors
#include <boost/serialization/vector.hpp>
// include this header to serialize string
#include <boost/serialization/string.hpp>
/*
 *
 */

class Solution{
public:
	Solution();
	virtual ~Solution();
	virtual void Set_Solution(Parameters *Param)=0;
	virtual void Set_output()=0;
	virtual void Set_breakpoint()=0;

protected:
//Pointer to general objects
	MultiBlock* MultiBlock_;
	ParallelManager* parallel;
	WriterManager* Writer;

//Solution variables
	Dictionary* Dic;
	double **U, *Rho;
	int nbnodes_real, nbnodes_total; //total include ghost nodes, real without ghost nodes
};
class Solution2D:public Solution {
public:
	Solution2D();
	virtual ~Solution2D();
	virtual void Set_Solution(Parameters *Param);

protected:
	void Set_output();
	void Set_breakpoint();
	//Read a variable from file
	void Read_Variable(std::string variablename, std::string filename);


protected:
//Mesh variables
	//std::vector<Node2D*> *Node;
	NodeArrays2D* NodeArrays;
	std::vector<int> IdBoundaries;
	PatchBc* PatchBc;



private:
	friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
    	ar & nbnodes_real & nbnodes_total & IdBoundaries;
    	for (int i=0;i<nbnodes_total;i++)
    		ar & Rho[i];
    	for (int j=0;j<nbnodes_total;j++)
    	{
    		ar & U[0][j];
    		ar & U[1][j];
    	}
    }
};

#endif /* CORE_SOLUTION_H_ */
