/*
 * Node2D.h
 *
 *  Created on: 21 Apr 2015
 *      Author: thomas
 */

#ifndef MESH_SINGLEBLOCK_NODE2D_H_
#define MESH_SINGLEBLOCK_NODE2D_H_
#include <iostream>
// include this header to serialize vectors
#include <boost/serialization/vector.hpp>
#include <boost/serialization/export.hpp>
#include <boost/serialization/assume_abstract.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <boost/serialization/utility.hpp>
#include <boost/serialization/list.hpp>

#include <boost/serialization/string.hpp>
#include <boost/serialization/nvp.hpp>
#include <boost/serialization/version.hpp>
#include <boost/archive/tmpdir.hpp>
/*
 *
 */

enum NodeType {Interior, Solid, Ghost, Corner, Wall, Periodic, Velocity, Symmetry, Pressure,ConcaveCorner,ConvexCorner, SolidGhost,GlobalCorner, SpecialWall};
enum CornerType {Concave,Convex};
enum GhostType{InteriorGhostType,SolidGhostType};
enum SymmetryType{OnNode,HalfWay,SymmetryPressureOnNode};
enum WallSpecialType{Standard,SymmetryWall, PressureWall, VelocityWall};

// Class Abstract
class Node2D {
public:
	Node2D(Node2D const& other);
	Node2D(signed short int x_=0, signed short int y_=0);
	virtual ~Node2D();
	double get_x() const;
	double get_y() const;
	void set_x(signed short int x_=0) ;
	void set_y(signed short int y_=0) ;
	NodeType& get_NodeType();
	void Set_NodeType(NodeType NodeType){NodeType_=NodeType;};
	//void Set_CornerType(CornerType CornerType_);
	void Set_CornerType(NodeType CornerType_);
	//void Set_WallType(Node2D* Node);
	//void Set_CornerType(Node2D* Node);
	void Set_Connect(unsigned int Connect_N_,unsigned int Connect_S_,unsigned int Connect_W_,unsigned int Connect_E_);
	void Remove_Connect(int direction);
	void Set_Connect(int* Connect_=0,int NbVelocity=0);
	unsigned int& Get_connect(unsigned int direction);
	int* Get_connect(){return Connect;};
	unsigned int& Get_NbVelocity(){return NbVelocity;};
	void Set_NbVelocity(unsigned int NbVelocitytmp){NbVelocity=NbVelocitytmp;};

	virtual bool* stream()=0;
	virtual double* Get_UDef()=0;
	virtual double Get_RhoDef()=0;

	virtual void Set_stream(bool* streaming=0,int NbVelocity=0)=0;
	virtual void Set_UDef(double UDef, double VDef)=0;
	virtual void Set_RhoDef(double RhoDef)=0;

	void Set_Index(int idx);
	int& Get_index();

protected:
	signed short int x,y;
	int index;
	unsigned int Connect_N,Connect_S,Connect_W,Connect_E;
	NodeType NodeType_;
	unsigned int NbVelocity;
	int* Connect;
private:
	friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
    	ar & x & y & Connect_N & Connect_S & Connect_W & Connect_E;
    	ar & NodeType_;
    	ar & NbVelocity;
    }

};
//BOOST_SERIALIZATION_ASSUME_ABSTRACT(Node2D);
class NodeInterior2D : public Node2D {
public:
	NodeInterior2D(signed short int x_=0,signed short int y_=0);
//	virtual void Set_Connect(unsigned int Connect_N_,unsigned int Connect_S_,unsigned int Connect_W_,unsigned int Connect_E_);
	virtual ~NodeInterior2D();


	virtual bool* stream();
	virtual double* Get_UDef();
	virtual double Get_RhoDef();

	virtual void Set_stream(bool* streaming=0,int NbVelocity=0);
	virtual void Set_UDef(double UDef, double VDef);
	virtual void Set_RhoDef(double RhoDef);
private:
	friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
    	ar & boost::serialization::base_object<Node2D>(*this);
    }
};

class NodeSolid2D : public Node2D {
public:
	NodeSolid2D(NodeSolid2D const& other);
	NodeSolid2D(signed short int x_=0, signed short int y_=0);
//	virtual void Set_Connect(unsigned int Connect_N_,unsigned int Connect_S_,unsigned int Connect_W_,unsigned int Connect_E_);
	virtual ~NodeSolid2D();

	virtual bool* stream();
	virtual double* Get_UDef();
	virtual double Get_RhoDef();

	virtual void Set_stream(bool* streaming=0,int NbVelocity=0);
	virtual void Set_UDef(double UDef, double VDef);
	virtual void Set_RhoDef(double RhoDef);

	void Set_FirstLayer(bool IsFirstLayer=false){firstlayer=IsFirstLayer;};
	bool IsFirstLayer(){return firstlayer;}
	bool operator()(NodeSolid2D& obj){return obj.Get_index() == index;}
	bool operator()(int& idxtest){return idxtest == index;}
private:
	bool firstlayer;

private:
	friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
    	ar & boost::serialization::base_object<Node2D>(*this);
    }
};

class NodeGhost2D : public Node2D {
public:
	NodeGhost2D(signed short int x_=0, signed short int y_=0, unsigned int Connect2Rank_=0,unsigned int Connect2Cell_=0,unsigned int Connect2Node_=0,NodeType GhostType=Ghost);
//	virtual void Set_Connect(unsigned int Connect_N_,unsigned int Connect_S_,unsigned int Connect_W_,unsigned int Connect_E_);
	virtual ~NodeGhost2D();

	void SetConnect(int Connect2Rank_=0,int Connect2Cell_=0,int Connect2Node_=0);
	int* GetConnect();

	virtual bool* stream();
	virtual double* Get_UDef();
	virtual double Get_RhoDef();

	virtual void Set_stream(bool* streaming=0,int NbVelocity=0);
	virtual void Set_UDef(double UDef, double VDef);
	virtual void Set_RhoDef(double RhoDef);

	void Set_GhostType(GhostType GhostType__){GhostType_=GhostType__;};
	GhostType& Get_GhostType(){return GhostType_;};

private:
	int Connect2Rank;
	int Connect2Cell;
	int Connect2Node;
	bool* GhostStream;
	GhostType GhostType_;
private:
	friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
    	ar & boost::serialization::base_object<Node2D>(*this);
    	ar & Connect2Rank & Connect2Cell & Connect2Node;
    	for (int i=0;i<NbVelocity;i++)
    		ar & GhostStream[i];
    }
};

class NodeVelocity2D : public Node2D {
public:
	NodeVelocity2D(signed short int x_=0, signed short int y_=0,double* UDef_=0 );
//	virtual void Set_Connect(unsigned int Connect_N_,unsigned int Connect_S_,unsigned int Connect_W_,unsigned int Connect_E_);
	void SetU(double Ux, double Uy);
	double* GetU();
	virtual ~NodeVelocity2D();


	virtual bool* stream();
	virtual double* Get_UDef();
	virtual double Get_RhoDef();

	virtual void Set_stream(bool* streaming=0,int NbVelocity=0);
	virtual void Set_UDef(double UDef, double VDef);
	virtual void Set_RhoDef(double RhoDef);

	virtual void Set_AlphaDef(double Alpha){alphaDef=Alpha;};
	virtual double Get_AlphaDef(){return alphaDef;};

	void Set_BcNormal(int BcNormal_){BcNormal=BcNormal_;};
	int& Get_BcNormal(){return BcNormal;};

private:
	double alphaDef;
	double UDef[2];
	bool* VelocityStream;
	int BcNormal;
private:
	friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
    	ar & boost::serialization::base_object<Node2D>(*this);
    	ar & UDef;
    	for (int i=0;i<NbVelocity;i++)
    		ar & VelocityStream[i];
    }
};

class NodePressure2D : public Node2D {
public:
	NodePressure2D(signed short int x_=0, signed short int y_=0,double PDef_=1.0);
//	virtual void Set_Connect(unsigned int Connect_N_,unsigned int Connect_S_,unsigned int Connect_W_,unsigned int Connect_E_);
	virtual ~NodePressure2D();

	void SetP(double PressureDefine_);
	double GetP();
	void SetRho(double RhoDef_);
	double GetRho();

	virtual bool* stream();
	virtual double* Get_UDef();
	virtual double Get_RhoDef();

	virtual void Set_stream(bool* streaming=0,int NbVelocity=0);
	virtual void Set_UDef(double UDef, double VDef);
	virtual void Set_RhoDef(double RhoDef);
	virtual void Set_AlphaDef(double Alpha){alphaDef=Alpha;};
	virtual double Get_AlphaDef(){return alphaDef;};

	void Set_BcNormal(int BcNormal_){BcNormal=BcNormal_;};
	int& Get_BcNormal(){return BcNormal;};
private:
	double RhoDef, alphaDef;//alpha is for two-phases flows
	bool* PressureStream;
	int BcNormal;
private:
	friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
    	ar & boost::serialization::base_object<Node2D>(*this);
    	ar & RhoDef;
    	for (int i=0;i<NbVelocity;i++)
    		ar & PressureStream[i];
    }
};

class NodePeriodic2D : public Node2D {
public:
	NodePeriodic2D(signed short int x_=0, signed short int y_=0);
//	virtual void Set_Connect(unsigned int Connect_N_,unsigned int Connect_S_,unsigned int Connect_W_,unsigned int Connect_E_);
	virtual ~NodePeriodic2D();
	virtual bool* stream();
	virtual double* Get_UDef();
	virtual double Get_RhoDef();

	virtual void Set_stream(bool* streaming=0,int NbVelocity=0);
	virtual void Set_UDef(double UDef, double VDef);
	virtual void Set_RhoDef(double RhoDef);

	void Set_BcNormal(int BcNormal_){BcNormal=BcNormal_;};
	int& Get_BcNormal(){return BcNormal;};
private:
	double RhoDef;
	double UDef[2];
	bool* PeriodicStream;
	int BcNormal;
private:
	friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
    	ar & boost::serialization::base_object<Node2D>(*this);
    }

};

class NodeCorner2D : public Node2D {
public:
	NodeCorner2D(signed short int x_=0, signed short int y_=0,double PDef_=1.0,double* UDef_=0);
//	virtual void Set_Connect(unsigned int Connect_N_,unsigned int Connect_S_,unsigned int Connect_W_,unsigned int Connect_E_);
	virtual ~NodeCorner2D();
	void SetP(double PDef_);
	double GetP();
	void SetRho(double RhoDef_);
	double GetRho();
	void SetU(double Ux, double Uy);
	double* GetU();

	virtual bool* stream();
	virtual double* Get_UDef();
	virtual double Get_RhoDef();

	virtual void Set_stream(bool* streaming=0,int NbVelocity=0);
	virtual void Set_UDef(double UDef, double VDef);
	virtual void Set_RhoDef(double RhoDef);
	virtual void Set_AlphaDef(double Alpha){alphaDef=Alpha;};
	virtual double Get_AlphaDef(){return alphaDef;};

	void Set_BcNormal(int BcNormal_){BcNormal=BcNormal_;};
	int& Get_BcNormal(){return BcNormal;};

	void Set_CornerType(CornerType CornerType__){CornerType_=CornerType__;};
	CornerType& Get_CornerType(){return CornerType_;};

	//Save distribution
	void Ini_SaveData(unsigned int NbSaveData);
	double Get_SaveData(int idx){return SaveData[idx];};
	void Set_SaveData(int idx,double value){SaveData[idx]=value;}
private:
	double RhoDef,alphaDef; //alpha is for two-phases flows
	double UDef[2];
	bool* CornerStream;
	int BcNormal;
	CornerType CornerType_;
	std::vector<double> SaveData;
private:
	friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
    	ar & boost::serialization::base_object<Node2D>(*this);
    	ar & RhoDef;
    	ar & UDef;
    	for (int i=0;i<NbVelocity;i++)
    		ar & CornerStream[i];
    }
};

class NodeWall2D : public Node2D {
public:
	NodeWall2D(signed short int x_=0, signed short int y_=0);
//	virtual void Set_Connect(unsigned int Connect_N_,unsigned int Connect_S_,unsigned int Connect_W_,unsigned int Connect_E_);
	virtual ~NodeWall2D();

	virtual bool* stream();
	virtual double* Get_UDef();
	virtual double Get_RhoDef();

	virtual void Set_stream(bool* streaming=0,int NbVelocity=0);
	virtual void Set_UDef(double UDef, double VDef);
	virtual void Set_RhoDef(double RhoDef);

	void Set_BcNormal(int BcNormal_){BcNormal=BcNormal_;};
	int& Get_BcNormal(){return BcNormal;};

	void Set_WallType(WallSpecialType WallSpecialType__){WallSpecialType_=WallSpecialType__;};
	WallSpecialType& Get_WallType(){return WallSpecialType_;};

	//for special walls
	virtual void Set_AlphaDef(double Alpha){alphaDef=Alpha;};
	virtual double Get_AlphaDef(){return alphaDef;};
	//Save distribution
	void Ini_SaveData(unsigned int NbSaveData);
	double Get_SaveData(int idx){return SaveData[idx];};
	void Set_SaveData(int idx,double value){SaveData[idx]=value;}
private:
	double RhoDef, alphaDef;
	double UDef[2];
	bool* WallStream;
	int BcNormal;
	WallSpecialType WallSpecialType_;
	std::vector<double> SaveData;
private:
	friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
    	ar & boost::serialization::base_object<Node2D>(*this);
    	for (int i=0;i<NbVelocity;i++)
    		ar & WallStream[i];
    }

};
class NodeSymmetry2D : public Node2D {
public:
	NodeSymmetry2D(signed short int x_=0, signed short int y_=0);
//	virtual void Set_Connect(unsigned int Connect_N_,unsigned int Connect_S_,unsigned int Connect_W_,unsigned int Connect_E_);
	virtual ~NodeSymmetry2D();

	virtual bool* stream();
	virtual double* Get_UDef();
	virtual double Get_RhoDef();

	virtual void Set_stream(bool* streaming=0,int NbVelocity=0);
	virtual void Set_UDef(double UDef, double VDef);
	virtual void Set_RhoDef(double RhoDef);

	void Set_BcNormal(int BcNormal_){BcNormal=BcNormal_;};
	int& Get_BcNormal(){return BcNormal;};

	void Set_SymmetryType(SymmetryType SymmetryType__){SymmetryType_=SymmetryType__;};
	SymmetryType& Get_SymmetryType(){return SymmetryType_;};

private:
	double RhoDef;
	double UDef[2];
	bool* SymmetryStream;
	int BcNormal;
	SymmetryType SymmetryType_;

private:
	friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
    	ar & boost::serialization::base_object<Node2D>(*this);
    	for (int i=0;i<NbVelocity;i++)
    		ar & SymmetryStream[i];
    }

};
BOOST_CLASS_EXPORT_KEY(Node2D);
BOOST_CLASS_EXPORT_KEY(NodeInterior2D);
BOOST_CLASS_EXPORT_KEY(NodeSolid2D);
BOOST_CLASS_EXPORT_KEY(NodeGhost2D);
BOOST_CLASS_EXPORT_KEY(NodeCorner2D);
BOOST_CLASS_EXPORT_KEY(NodeWall2D);
BOOST_CLASS_EXPORT_KEY(NodePeriodic2D);
BOOST_CLASS_EXPORT_KEY(NodeVelocity2D);
BOOST_CLASS_EXPORT_KEY(NodePressure2D);
BOOST_CLASS_EXPORT_KEY(NodeSymmetry2D);
#endif /* MESH_SINGLEBLOCK_NODE2D_H_ */
