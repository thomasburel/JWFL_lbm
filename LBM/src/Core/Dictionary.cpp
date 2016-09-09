/*
 * Dictionary.cpp
 *
 *  Created on: 24 Aug 2016
 *      Author: Thomas Burel
 */

#include "Dictionary.h"

Dictionary::Dictionary() {
	dimension=2;
	numberNodes=0;
	PtrExportVar=0;
	PtrStringExportVar=0;
	PtrExportVarBreakpoint=0;
	PtrStringExportVarBreakpoint=0;
}
Dictionary::Dictionary(int dim, int nbNodes) {
	dimension=dim;
	numberNodes=nbNodes;
	PtrExportVar=0;
	PtrStringExportVar=0;
	PtrExportVarBreakpoint=0;
	PtrStringExportVarBreakpoint=0;
}
Dictionary::~Dictionary() {
	for(int i=0;i<Var.size();i++)
		delete [] Var[i];
	StringExportVar.clear();
	ExportVar.clear();
	StringExportVarBreakpoint.clear();
	ExportVarBreakpoint.clear();
	SyncVar.clear();
	Var.clear();
	if(PtrStringExportVar!=0)
		delete PtrStringExportVar;
	if(PtrExportVar!=0)
			delete PtrExportVar;
	if(PtrStringExportVarBreakpoint!=0)
		delete PtrStringExportVarBreakpoint;
	if(PtrExportVarBreakpoint!=0)
			delete PtrExportVarBreakpoint;
}

void Dictionary::AddDistributionBreakpoint(std::string VarName,double* & Var_x)
{
	ExportVarBreakpoint.push_back(Var_x);
	StringExportVarBreakpoint.push_back(VarName);
}
void Dictionary::AddSync(std::string VarName,double* & Var_x){
	SyncVar.push_back(Var_x);
	StringSyncVar.push_back(VarName);
}
void Dictionary::AddVar(dataType type,std::string VarName,bool BoolExportVar, bool BoolExportBreakpoint, bool BoolSynchronise, double* & Var_x){

	if(type==Scalar)
	{


		Var_x=new double [numberNodes];
		mapVar[VarName]=Var.size();
		Var.push_back(Var_x);


		if(BoolExportVar)
		{
			ExportVar.push_back(Var_x);
			StringExportVar.push_back(VarName);
		}

		if(BoolExportBreakpoint)
		{
			ExportVarBreakpoint.push_back(Var_x);
			StringExportVarBreakpoint.push_back(VarName);
		}

		if(BoolSynchronise)
		{
			SyncVar.push_back(Var_x);
			StringSyncVar.push_back(VarName);
		}
	}
	else
	{
		std::cerr<<"Variable named: "<<VarName<<" is not a scalar. Adding in the dictionary failed."<<std::endl;
	}

}


void Dictionary::AddVar(dataType type,std::string VarName,bool BoolExportVar, bool BoolExportBreakpoint, bool BoolSynchronise, double* & Var_x, double* & Var_y){
	if(dimension==2 &&type==Vector)
	{
		Var_x=new double [numberNodes];
		mapVar[VarName+"X"]=Var.size();
		Var.push_back(Var_x);
		Var_y=new double [numberNodes];
		mapVar[VarName+"Y"]=Var.size();
		Var.push_back(Var_y);

		if(BoolExportVar)
		{
			ExportVar.push_back(Var_x);
			ExportVar.push_back(Var_y);

			StringExportVar.push_back(VarName+"X");
			StringExportVar.push_back(VarName+"Y");
		}

		if(BoolExportBreakpoint)
		{
			ExportVarBreakpoint.push_back(Var_x);
			ExportVarBreakpoint.push_back(Var_y);

			StringExportVarBreakpoint.push_back(VarName+"X");
			StringExportVarBreakpoint.push_back(VarName+"Y");
		}

		if(BoolSynchronise)
		{
			SyncVar.push_back(Var_x);
			SyncVar.push_back(Var_y);

			StringSyncVar.push_back(VarName+"X");
			StringSyncVar.push_back(VarName+"Y");
		}
	}
	else
	{
		std::cerr<<"Variable named: "<<VarName<<" is not a vector in 2D. Adding in the dictionary failed."<<std::endl;
	}
}


void Dictionary::AddVar(dataType type,std::string VarName,bool BoolExportVar, bool BoolExportBreakpoint, bool BoolSynchronise, double* & Var_x, double* & Var_y, double* & Var_z){
	if(dimension==3 && type==Vector)
	{
		Var_x=new double [numberNodes];
		mapVar[VarName+"X"]=Var.size();
		Var.push_back(Var_x);
		Var_y=new double [numberNodes];
		mapVar[VarName+"Y"]=Var.size();
		Var.push_back(Var_y);
		Var_z=new double [numberNodes];
		mapVar[VarName+"Z"]=Var.size();
		Var.push_back(Var_z);

		if(BoolExportVar)
		{
			ExportVar.push_back(Var_x);
			ExportVar.push_back(Var_y);
			ExportVar.push_back(Var_z);

			StringExportVar.push_back(VarName+"X");
			StringExportVar.push_back(VarName+"Y");
			StringExportVar.push_back(VarName+"Z");
		}

		if(BoolExportBreakpoint)
		{
			ExportVarBreakpoint.push_back(Var_x);
			ExportVarBreakpoint.push_back(Var_y);
			ExportVarBreakpoint.push_back(Var_z);

			StringExportVarBreakpoint.push_back(VarName+"X");
			StringExportVarBreakpoint.push_back(VarName+"Y");
			StringExportVarBreakpoint.push_back(VarName+"Z");
		}

		if(BoolSynchronise)
		{
			SyncVar.push_back(Var_x);
			SyncVar.push_back(Var_y);
			SyncVar.push_back(Var_z);

			StringSyncVar.push_back(VarName+"X");
			StringSyncVar.push_back(VarName+"Y");
			StringSyncVar.push_back(VarName+"Z");
		}
	}
	else
	{
		std::cerr<<"Variable named: "<<VarName<<" is not a vector in 3D. Adding in the dictionary failed."<<std::endl;
	}
}
double** Dictionary::Get_PtrExportVar()
{
	if(PtrExportVar!=0)
		delete []PtrExportVar;
	PtrExportVar=new double*[Get_NbExportVar()];
	for(int i=0;i<Get_NbExportVar();i++)
		PtrExportVar[i]=ExportVar[i];
	return PtrExportVar;
}
std::string* Dictionary::Get_PtrExportVarName()
{
	if(PtrStringExportVar!=0)
		delete [] PtrStringExportVar;
	PtrStringExportVar=new std::string[Get_NbExportVar()];
	for(int i=0;i<Get_NbExportVar();i++)
		PtrStringExportVar[i]=StringExportVar[i];
	return PtrStringExportVar;
}
double** Dictionary::Get_PtrExportVarBreakpoint()
{
	if(PtrExportVarBreakpoint!=0)
		delete PtrExportVarBreakpoint;
	PtrExportVarBreakpoint=new double*[Get_NbExportVarBreakpoint()];
	for(int i=0;i<Get_NbExportVarBreakpoint();i++)
		PtrExportVarBreakpoint[i]=ExportVarBreakpoint[i];
	return PtrExportVarBreakpoint;
}
std::string* Dictionary::Get_PtrExportVarBreakpointName()
{
	if(PtrStringExportVarBreakpoint!=0)
		delete [] PtrStringExportVarBreakpoint;
	PtrStringExportVarBreakpoint=new std::string[Get_NbExportVarBreakpoint()];
	for(int i=0;i<Get_NbExportVarBreakpoint();i++)
		PtrStringExportVarBreakpoint[i]=StringExportVarBreakpoint[i];
	return PtrStringExportVarBreakpoint;
}
int Dictionary::Get_Id_Var(std::string VarName){
	 std::map<std::string,int>::iterator it;
	 it=mapVar.find(VarName);
	 if(it==mapVar.end())
	 {
		 std::cerr<<" Variable name ( "<< VarName<<" )in the dictionary is not found"<<std::endl;
		 std::cerr<<"Variable available are: " ;
		 for (it = mapVar.begin() ; it != mapVar.end() ; ++it)
			 std::cerr<<it->first<<" ; ";
		 std::cerr<<std::endl;
	 }
	return mapVar[VarName];
}
