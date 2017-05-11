/*
 * PorousMedia.h
 *
 *  Created on: 17 Oct 2015
 *      Author: thomas
 */

#ifndef SPECIALTREATMENT_POROUSMEDIA_H_
#define SPECIALTREATMENT_POROUSMEDIA_H_
#include <iostream>
#include <fstream>
#include <ios>
#include <string>
// include this header to serialize vectors
#include <boost/serialization/vector.hpp>

class PorousMedia {
public:
	PorousMedia();
	virtual ~PorousMedia();
	void AddHeleShawDragSinglePhase(double const  &u,double const  &v,double const  &mu,double &Fx,double &Fy,double const InterfaceFx=0,double const InterfaceFy=0);
	void AddHeleShawDragTwoPhases(double const  &u,double const  &v,double const  &mu,double &Fx,double &Fy,double const InterfaceFx=0,double const InterfaceFy=0);
protected:
	void ReadBinaryImage(std::string filename,int nx, int ny, int nz);
	void AdaptSolid();

private:
	void DomainInteriorTreatment2D();
	void DomainBoundaryTreatment2D();
	void DomainCornerTreatment2D();
	void CorrectPoreConnections2D(int & i, int & j,bool & allneighbours);
	void RemoveWrongCorners2D(int i, int j,bool allneighbours);
	bool CheckWrongCorners2D(int i, int j,bool & allneighbours, int & diagonal);
protected:
	bool ***image;
//	std::vector<std::vector<std::vector<bool>>> imagevect;
	int image_nx,image_ny,image_nz;
	int xstartmedia,ystartmedia,zstartmedia;
	double depth,depth2;
};

#endif /* SPECIALTREATMENT_POROUSMEDIA_H_ */
