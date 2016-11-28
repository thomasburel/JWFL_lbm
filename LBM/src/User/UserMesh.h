/*
 * UserMesh.h
 *
 *  Created on: 17 Oct 2015
 *      Author: thomas
 */

#ifndef USER_USERMESH_H_
#define USER_USERMESH_H_
#include "../SpecialTreatment/PorousMedia.h"
#include "../Mesh/SingleBlock/Node2D.h"
#include "../Core/Parameters.h"
#include <cmath>
class UserMesh :public PorousMedia {
public:
	UserMesh();
	virtual ~UserMesh();

	void SetUserMeshVariables();
	void ChangeNode(Node2D &Node, bool &solid );
	void SetSymmetryType(SymmetryType &Type, double x, double y);
	//void ChangeNode(Node3D &Node, bool solid );
protected:
double a,b,epsilon;
double R,H,L,teta,pi;
double beta;
Parameters *PtrParametersUserMesh;
};

#endif /* USER_USERMESH_H_ */
