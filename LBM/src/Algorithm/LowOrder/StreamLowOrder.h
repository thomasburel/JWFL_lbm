/*
 * Stream.h
 *
 *  Created on: 15 Jun 2015
 *      Author: thomas
 */

#ifndef ALGORITHM_STREAM_H_
#define ALGORITHM_STREAM_H_
#include "../../Core/GlobalDef.h"
#include "../../Mesh/SingleBlock.h"
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
class StreamLowOrder {
public:
	StreamLowOrder();
	virtual ~StreamLowOrder();
	void Set_StreamLowOrder();

protected:
	Block* PtrBlockStream;
	DistriFunct* PtrFiStream;
private:
	friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive &ar, const unsigned int version)
    {

    }
};

#endif /* ALGORITHM_STREAM_H_ */
