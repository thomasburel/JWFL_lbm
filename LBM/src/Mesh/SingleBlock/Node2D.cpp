/*
 * Node2D.cpp
 *
 *  Created on: 21 Apr 2015
 *      Author: thomas
 */

#include "Node2D.h"
#include <boost/serialization/export.hpp>
//BOOST_CLASS_EXPORT_IMPLEMENT(Node2D);
BOOST_CLASS_EXPORT_IMPLEMENT(NodeInterior2D);
BOOST_CLASS_EXPORT_IMPLEMENT(NodeSolid2D);
BOOST_CLASS_EXPORT_IMPLEMENT(NodeGhost2D);
BOOST_CLASS_EXPORT_IMPLEMENT(NodeCorner2D);
BOOST_CLASS_EXPORT_IMPLEMENT(NodeWall2D);
BOOST_CLASS_EXPORT_IMPLEMENT(NodePeriodic2D);
BOOST_CLASS_EXPORT_IMPLEMENT(NodeVelocity2D);
BOOST_CLASS_EXPORT_IMPLEMENT(NodePressure2D);
//**************General Methods for Nodes in 2D***************
///Constructor for General Nodes in 2D
Node2D::Node2D(signed short int x_,signed short int y_): x(x_), y(y_),NodeType_(Ghost),Connect_N(0),Connect_S(0),Connect_W(0),Connect_E(0),NbVelocity(0),Connect(0)
{index=0;}
Node2D::Node2D(Node2D const& other){
	x=other.x;
	y=other.y;
	NodeType_=other.NodeType_;
	index=other.index;
	Connect_N=other.Connect_N;Connect_S=other.Connect_S;Connect_W=other.Connect_W;Connect_E=other.Connect_E;
	NbVelocity=other.NbVelocity;
	if(other.Connect!=0)
	{
		Connect=new int[NbVelocity];
		for(int i=0;i<NbVelocity;i++)
			Connect[i]=other.Connect[i];
	}
	else
		Connect=0;
}
///Destructor for General Nodes in 2D
Node2D::~Node2D()
{
	if(Connect!=0)
	delete [] Connect;
}

double Node2D::get_x() const
{
	return x;
}

double Node2D::get_y() const
{
	return y;
}
void Node2D::set_x(signed short int x_)
{
	x=x_;
}

void Node2D::set_y(signed short int y_)
{
	y=y_;
}
NodeType& Node2D::get_NodeType()
{
	return NodeType_;
}
/*void Node2D::Set_SolidType(Node2D* Node)
{
	unsigned int x_tmp=x;
	unsigned int y_tmp=y;
	unsigned int Connect_N_tmp=Connect_N;
	unsigned int Connect_S_tmp=Connect_S;
	unsigned int Connect_W_tmp=Connect_W;
	unsigned int Connect_E_tmp=Connect_E;
	unsigned int NbVelocity_tmp=NbVelocity;
	delete Node;
	Node=new NodeSolid2D(x,y);
	Connect_N=Connect_N_tmp;
	Connect_S=Connect_S_tmp;
	Connect_W=Connect_W_tmp;
	Connect_E=Connect_E_tmp;
	NodeType_=Solid;
	NbVelocity=NbVelocity_tmp;

}
void Node2D::Set_WallType(Node2D* Node)
{
	unsigned int x_tmp=x;
	unsigned int y_tmp=y;
	unsigned int Connect_N_tmp=Connect_N;
	unsigned int Connect_S_tmp=Connect_S;
	unsigned int Connect_W_tmp=Connect_W;
	unsigned int Connect_E_tmp=Connect_E;
	unsigned int NbVelocity_tmp=NbVelocity;
	delete Node;
	Node=new NodeWall2D(x,y);
	Connect_N=Connect_N_tmp;
	Connect_S=Connect_S_tmp;
	Connect_W=Connect_W_tmp;
	Connect_E=Connect_E_tmp;
	NodeType_=Wall;
	NbVelocity=NbVelocity_tmp;

}*/
/*void Node2D::Set_CornerType(CornerType CornerType_)
{
	switch(CornerType_)
	{
	case ConcaveCornerNoVel:
		NodeType_==ConvexCorner;
		break;
	case ConvexCornerNoVel:
		NodeType_=ConvexCorner;
		break;
	case ConcaveCornerWithVel:
		NodeType_=ConcaveCorner;
		break;
	case ConvexCornerWithVel:
		NodeType_=ConcaveCorner;
		break;
	}
}*/
void Node2D::Set_CornerType(NodeType CornerType_)
{
	if((CornerType_==ConvexCorner)||(CornerType_==ConcaveCorner))
	{
		NodeType_=CornerType_;
	}
	else
		std::cerr<<" Error corner type"<<std::endl;
}
void Node2D::Set_Connect(unsigned int Connect_N_,unsigned int Connect_S_,unsigned int Connect_W_,unsigned int Connect_E_){
	Connect_N=Connect_N_-1;
	Connect_S=Connect_S_-1;
	Connect_W=Connect_W_-1;
	Connect_E=Connect_E_-1;
}
void Node2D::Remove_Connect(int direction){
	switch (direction)
	{
	case 0:
		Connect_S=index;
		break;
	case 1:
		Connect_E=index;
		break;
	case 2:
		Connect_N=index;
		break;
	case 3:
		Connect_W=index;
		break;
	default:
		std::cerr<<"Wrong Connection ask. Default Connection in South node "<<std::endl;
		Connect_S=index;
		break;
	}
}
void Node2D::Set_Connect(int* Connect_,int NbVelocity)
{
	if(Connect!=0)
		delete []Connect;
	Connect=new int[NbVelocity];
	for(int i=0; i<NbVelocity;i++)
		Connect[i]=Connect_[i];
}

unsigned int&  Node2D::Get_connect(unsigned int direction){

	switch (direction)
	{
	case 0:
		return Connect_S;
		break;
	case 1:
		return Connect_E;
		break;
	case 2:
		return Connect_N;
		break;
	case 3:
		return Connect_W;
		break;
	default:
		std::cerr<<"Wrong Connection ask. Default Connection in South node "<<std::endl;
		return Connect_S;
		break;
	}
}
void Node2D::Set_Index(int idx){index=idx;};
int& Node2D::Get_index(){return index;};
//**************Interior Node Methods***************
///Constructor for Interior Nodes
NodeInterior2D::NodeInterior2D(signed short int x_,signed short int y_)
{
	x=x_;
	y=y_;
	NodeType_=Interior;
}
///Destructor for Interior Nodes
NodeInterior2D::~NodeInterior2D()
{}

bool* NodeInterior2D::stream(){return 0;} //Don't need
double* NodeInterior2D::Get_UDef(){return 0;}//Don't need
double NodeInterior2D::Get_RhoDef(){return 0;}//Don't need

void NodeInterior2D::Set_stream(bool* streaming,int NbVelocity){}//Don't need
void NodeInterior2D::Set_UDef(double UDef, double VDef){}//Don't need
void NodeInterior2D::Set_RhoDef(double RhoDef){}//Don't need

//**************Solid Node Methods***************
///Constructor for Solid Nodes
NodeSolid2D::NodeSolid2D(signed short int x_,signed short int y_)
{
	x=x_;
	y=y_;
	NodeType_=Solid;
	firstlayer=false;
}
NodeSolid2D::NodeSolid2D(NodeSolid2D const& other): Node2D( other ){
	firstlayer=other.firstlayer;
}
///Destructor for Solid Nodes
NodeSolid2D::~NodeSolid2D()
{}
bool* NodeSolid2D::stream(){return 0;}//Don't need
double* NodeSolid2D::Get_UDef(){return 0;}//Don't need
double NodeSolid2D::Get_RhoDef(){return 0;}//Don't need

void NodeSolid2D::Set_stream(bool* streaming,int NbVelocity){}//Don't need
void NodeSolid2D::Set_UDef(double UDef, double VDef){}//Don't need
void NodeSolid2D::Set_RhoDef(double RhoDef){}//Don't need

//**************Ghost Node Methods***************
///Constructor for Ghost Nodes
NodeGhost2D::NodeGhost2D(signed short int x_,signed short int y_, unsigned int Connect2Rank_,unsigned int Connect2Cell_, unsigned int Connect2Node_,NodeType GhostType)
{
	x=x_;
	y=y_;
	NodeType_=GhostType;
	Connect2Rank=Connect2Rank_;
	Connect2Cell=Connect2Cell_;
	Connect2Node=Connect2Node_;
	GhostStream=0;
	NodeType_=Ghost;
	GhostType_=InteriorGhostType;
}
///Destructor for GhostNodes
NodeGhost2D::~NodeGhost2D()
{delete [] GhostStream;}
/// Set the rank of the processor which is connected to this ghost node
void NodeGhost2D::SetConnect(int Connect2Rank_,int Connect2Cell_, int Connect2Node_)
{
	Connect2Rank=Connect2Rank_;
	Connect2Cell=Connect2Cell_;
	Connect2Node=Connect2Node_;
}
/// Get the rank of the processor which is connected to this ghost node
int* NodeGhost2D::GetConnect() {
	int* Connect= new int[3];
	Connect[0]=Connect2Rank;
	Connect[1]=Connect2Cell;
	Connect[2]=Connect2Node;
	return Connect;
}
bool* NodeGhost2D::stream(){return GhostStream;}
double* NodeGhost2D::Get_UDef(){return 0;}//Don't need
double NodeGhost2D::Get_RhoDef(){return 0;}//Don't need

void NodeGhost2D::Set_stream(bool* streaming,int NbVelocity_){
	NbVelocity=NbVelocity_;
	GhostStream=new bool [NbVelocity];
	for (int i=0;i<NbVelocity;i++)
		GhostStream[i]=streaming[i];
}
void NodeGhost2D::Set_UDef(double UDef, double VDef){}//Don't need
void NodeGhost2D::Set_RhoDef(double RhoDef){}//Don't need

//**************Velocity Node Methods***************
///Constructor for Velocity Nodes
NodeVelocity2D::NodeVelocity2D(signed short int x_,signed short int y_,double* UDef_){
	x=x_;
	y=y_;
	NodeType_=Velocity;
	VelocityStream=0;
	BcNormal=0;
	if(UDef_!=0)
	{
		UDef[0]=UDef_[0];
		UDef[1]=UDef_[1];
	}
}
///Destructor for Velocity Nodes
NodeVelocity2D::~NodeVelocity2D() {
	delete [] VelocityStream;
}

/// Set the Velocity profile at the node
void NodeVelocity2D::SetU(double Ux, double Uy) {
	UDef[0]=Ux;
	UDef[1]=Uy;
}
/// Get the Velocity profile at the node
double* NodeVelocity2D::GetU() {
	return UDef;
}


bool* NodeVelocity2D::stream(){return VelocityStream;}
double* NodeVelocity2D::Get_UDef(){return UDef;}
double NodeVelocity2D::Get_RhoDef(){return 0;}//Don't need

void NodeVelocity2D::Set_stream(bool* streaming,int NbVelocity_){
	NbVelocity=NbVelocity_;
	VelocityStream=new bool [NbVelocity];
	for (int i=0;i<NbVelocity;i++)
		VelocityStream[i]=streaming[i];
}
void NodeVelocity2D::Set_UDef(double U_Def, double V_Def){UDef[0]=U_Def;UDef[1]=V_Def;}
void NodeVelocity2D::Set_RhoDef(double RhoDef){}

//**************Pressure Node Methods***************
///Constructor for Pressure Nodes
NodePressure2D::NodePressure2D(signed short int x_,signed short int y_,double PDef_) {
	x=x_;
	y=y_;
	NodeType_=Pressure;
	RhoDef=PDef_;
	PressureStream=0;
	BcNormal=0;
}
///Destructor for Pressure Nodes
NodePressure2D::~NodePressure2D() {
	delete [] PressureStream;
}
/// Set the Pressure at the node
void NodePressure2D::SetP(double PDef_) {
	RhoDef=PDef_;

}
/// Get the Pressure at the node
double NodePressure2D::GetP() {
	return RhoDef;
}
/// Get the Density at the node
double NodePressure2D::GetRho() {
	return RhoDef;
}


bool* NodePressure2D::stream(){return PressureStream;}//Don't need
double* NodePressure2D::Get_UDef(){return 0;}
double NodePressure2D::Get_RhoDef(){return RhoDef;}

void NodePressure2D::Set_stream(bool* streaming,int NbVelocity_){
	NbVelocity=NbVelocity_;
	PressureStream=new bool [NbVelocity];
	for (int i=0;i<NbVelocity;i++)
		PressureStream[i]=streaming[i];
}
void NodePressure2D::Set_UDef(double UDef, double VDef){}
void NodePressure2D::Set_RhoDef(double RhoDef_){RhoDef=RhoDef_;}


//**************Periodic Node Methods***************
///Constructor for Periodic Nodes
NodePeriodic2D::NodePeriodic2D(signed short int x_,signed short int y_) {
	x=x_;
	y=y_;
	NodeType_=Periodic;
	RhoDef=1;
	UDef[0]=0;
	UDef[1]=0;
	PeriodicStream=0;
	BcNormal=0;
}
///Destructor for Periodic Nodes
NodePeriodic2D::~NodePeriodic2D() {

}


bool* NodePeriodic2D::stream(){return PeriodicStream;}
double* NodePeriodic2D::Get_UDef(){return UDef;}
double NodePeriodic2D::Get_RhoDef(){return RhoDef;}

void NodePeriodic2D::Set_stream(bool* streaming,int NbVelocity_){
	NbVelocity=NbVelocity_;
	PeriodicStream=new bool [NbVelocity];
	for (int i=0;i<NbVelocity;i++)
		PeriodicStream[i]=streaming[i];
}
void NodePeriodic2D::Set_UDef(double U_Def, double V_Def){UDef[0]=U_Def;UDef[1]=V_Def;}
void NodePeriodic2D::Set_RhoDef(double RhoDef_){RhoDef=RhoDef_;}


//**************Corner Node Methods***************
///Constructor for Corner Nodes
NodeCorner2D::NodeCorner2D(signed short int x_,signed short int y_,double PDef_,double* UDef_) {
	x=x_;
	y=y_;
	NodeType_=Corner;
	RhoDef=PDef_;
	CornerStream=0;
	BcNormal=0;
	CornerType_=Concave;
	if(UDef_!=0)
	{
		UDef[0]=UDef_[0];
		UDef[1]=UDef_[1];
	}

}
///Destructor for Corner Nodes
NodeCorner2D::~NodeCorner2D() {
	delete [] CornerStream;
}
/// Set the Pressure at the corner
void NodeCorner2D::SetP(double PDef_) {
	RhoDef=PDef_;

}
/// Get the Pressure at the corner
double NodeCorner2D::GetP() {
	return RhoDef;
}
/// Set the Density at the corner
void NodeCorner2D::SetRho(double RhoDef_){
	RhoDef=RhoDef_;
}
/// Get the Density at the corner
double NodeCorner2D::GetRho() {
	return RhoDef;
}
/// Set the Velocity profile at the corner
void NodeCorner2D::SetU(double Ux, double Uy) {
	UDef[0]=Ux;
	UDef[1]=Uy;
}
/// Get the Velocity profile at the corner
double* NodeCorner2D::GetU() {
	return UDef;
}



bool* NodeCorner2D::stream(){return CornerStream;}//Don't need
double* NodeCorner2D::Get_UDef(){return UDef;}
double NodeCorner2D::Get_RhoDef(){return RhoDef;}

void NodeCorner2D::Set_stream(bool* streaming,int NbVelocity_){
	NbVelocity=NbVelocity_;
	CornerStream=new bool [NbVelocity];
	for (int i=0;i<NbVelocity;i++)
		CornerStream[i]=streaming[i];
}
void NodeCorner2D::Set_UDef(double U_Def, double V_Def){UDef[0]=U_Def;UDef[1]=V_Def;}
void NodeCorner2D::Set_RhoDef(double Rho_Def){RhoDef=Rho_Def;}


//**************General Wall Node Methods no distinction between kind of wall nodes as bounce back, half way bounce back***************
///Constructor for General Wall Nodes
NodeWall2D::NodeWall2D(signed short int x_,signed short int y_) {
	x=x_;
	y=y_;
	NodeType_=Wall;
	WallStream=0;
	BcNormal=0;
	WallSpecialType_=Standard;
	alphaDef=0;
	RhoDef=1;
	UDef[0]=0;UDef[1]=0;
}
///Destructor for Wall Nodes
NodeWall2D::~NodeWall2D() {
	delete [] WallStream;
}


bool* NodeWall2D::stream(){return WallStream;}
double* NodeWall2D::Get_UDef(){return UDef;}
double NodeWall2D::Get_RhoDef(){return RhoDef;}

void NodeWall2D::Set_stream(bool* streaming,int NbVelocity_){
	NbVelocity=NbVelocity_;
	WallStream=new bool [NbVelocity];
	for (int i=0;i<NbVelocity;i++)
		WallStream[i]=streaming[i];
}
void NodeWall2D::Set_UDef(double UDef_, double VDef){UDef[0]=UDef_;UDef[1]=VDef;}
void NodeWall2D::Set_RhoDef(double RhoDef_){RhoDef=RhoDef_;}

//**************General Symmetry Node Methods no distinction between kind of Symmetry nodes as on node or half way ***************
///Constructor for General Symmetry Nodes
NodeSymmetry2D::NodeSymmetry2D(signed short int x_,signed short int y_) {
	x=x_;
	y=y_;
	NodeType_=Symmetry;
	SymmetryStream=0;
	BcNormal=0;
	SymmetryType_=OnNode;
	UDef[0]=0;UDef[1]=0;
	RhoDef=1;
}
///Destructor for Symmetry Nodes
NodeSymmetry2D::~NodeSymmetry2D() {
	delete [] SymmetryStream;
}


bool* NodeSymmetry2D::stream(){return SymmetryStream;}
double* NodeSymmetry2D::Get_UDef(){return UDef;}
double NodeSymmetry2D::Get_RhoDef(){return RhoDef;}

void NodeSymmetry2D::Set_stream(bool* streaming,int NbVelocity_){
	NbVelocity=NbVelocity_;
	SymmetryStream=new bool [NbVelocity];
	for (int i=0;i<NbVelocity;i++)
		SymmetryStream[i]=streaming[i];
}
void NodeSymmetry2D::Set_UDef(double U_Def, double V_Def){UDef[0]=U_Def;UDef[1]=V_Def;}
void NodeSymmetry2D::Set_RhoDef(double Rho_Def){RhoDef=Rho_Def;}
