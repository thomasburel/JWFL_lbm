/*
 * GlobalDef.h
 *
 *  Created on: 20 Apr 2015
 *      Author: thomas
 */

#ifndef CORE_GLOBALDEF_H_
#define CORE_GLOBALDEF_H_
#include <boost/serialization/split_member.hpp>
#include <boost/serialization/split_free.hpp>
///Class to manage distribution function to destroyed arrays correctly after used
class DistriFunct {
private:
	friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
    	boost::serialization::split_member(ar, *this, version);
    }
    template<class Archive>
    void save(Archive & ar, const unsigned int version) const
    {
        ar  & NbNodes;
        ar  & NbVelocities;
    	for (int i=0;i<NbVelocities;i++)
    	{
    		for (int j=0;j<NbNodes;j++)
    			ar  & f[i][j];
    	}
    }
    template<class Archive>
    void load(Archive & ar, const unsigned int version)

    {
        ar  & NbNodes ;
        ar  & NbVelocities;
        f=new double* [NbVelocities];
    	for (int i=0;i<NbVelocities;i++)
    	{
    		f[i]=new double [NbNodes];
    		for (int j=0;j<NbNodes;j++)
    			ar  & f[i][j];
    	}
    }
public:
	DistriFunct();
	DistriFunct(int NbNodes_, int NbVelocities_);
	virtual ~DistriFunct();

public: //variables in public to easy access
	double **f;

private:
	int NbNodes, NbVelocities;

};


  //  BOOST_SERIALIZATION_SPLIT_MEMBER(DistriFunct)


#endif /* CORE_GLOBALDEF_H_ */
