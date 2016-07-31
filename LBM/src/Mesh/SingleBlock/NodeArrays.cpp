/*
 * NodeArray.cpp
 *
 *  Created on: 9 Dec 2015
 *      Author: thomas
 */

#include "NodeArrays.h"

NodeArrays::NodeArrays() {
	// TODO Auto-generated constructor stub

}

NodeArrays::~NodeArrays() {
	// TODO Auto-generated destructor stub
}
NodeArrays2D::NodeArrays2D() {
	// TODO Auto-generated constructor stub

}

NodeArrays2D::~NodeArrays2D() {
	NodeInterior.clear();
	NodeSolid.clear();
	NodeGhost.clear();
	NodeCorner.clear();
	NodeGlobalCorner.clear();
	NodeWall.clear();
	NodeSpecialWall.clear();
	NodePeriodic.clear();
	NodeVelocity.clear();
	NodePressure.clear();
	NodeSymmetry.clear();
}
