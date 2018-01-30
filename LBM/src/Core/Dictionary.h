/*
 * Dictionary.h
 *
 *  Created on: 24 Aug 2016
 *      Author: Thomas Burel
 */

#ifndef SRC_CORE_DICTIONARY_H_
#define SRC_CORE_DICTIONARY_H_

#include <iostream>
#include <vector>
#include <map>

enum dataType{Scalar,Vector};
class Dictionary {
public:
	Dictionary();
	Dictionary(int dim, int nbNodes);
	virtual ~Dictionary();

	void AddDistributionBreakpoint(std::string VarName,double* & Var_x);
	void UpdatedDistributionBreakpoint(int idx,double* & Var_x);
	void AddSync(std::string VarName,double* & Var_x);

	void AddVar(dataType type,std::string VarName,bool ExportVar, bool ExportBreakpoint, bool Synchronise, double* & Var_x);
	void AddVar(dataType type,std::string VarName,bool ExportVar, bool ExportBreakpoint, bool Synchronise, double* & Var_x, double* & Var_y);
	void AddVar(dataType type,std::string VarName,bool ExportVar, bool ExportBreakpoint, bool Synchronise, double* & Var_x, double* & Var_y, double* & Var_z);

	std::vector<std::string> Get_ExportVarName(){return StringExportVar;};
	std::string* Get_PtrExportVarName();//{return PtrStringExportVar;};
	std::vector<double*> Get_ExportVar(){return ExportVar;};
	double** Get_PtrExportVar();//{return PtrExportVar;};
	int Get_NbExportVar(){return ExportVar.size();};

	std::string *Get_PtrExportVarBreakpointName();
	double** Get_PtrExportVarBreakpoint();
	int Get_NbExportVarBreakpoint();

	std::vector<double*> Get_SyncVar(){return SyncVar;};
	std::vector<std::string> Get_SyncVarName(){return StringSyncVar;};
	int Get_NbSyncVar(){return SyncVar.size();};

	/// Get the ID of the variable from its name. It is a slow access.
	int Get_Id_Var(std::string VarName, bool &Var_found);
	/// Get the pointer of the variable from its Id. It is a quick access.
	double* &Get_PtrVar(int Id_Var){return Var[Id_Var];};
	/// Get the pointer of the variable from its name. It is a slow access.
	void Get_PtrVar(std::string VarName, double* & Var_out, bool &Var_found){Var_out=Var[Get_Id_Var(VarName,Var_found)];};

	/// Get number of nodes per variable
	int Get_NbNodes(){return numberNodes;};


private:
	int dimension;
	int numberNodes;//with ghost

//Save Variables (save the full array). Distribution functions are saved outside.
	std::vector<double*> Var;
	std::map<std::string,int> mapVar;

//Save export variables (array of pointers)
	std::vector<std::string> StringExportVar;
	std::string* PtrStringExportVar;
	std::vector<double*> ExportVar;
	double** PtrExportVar;


//Save breakpoint variables (array of pointers)
	std::vector<std::string> StringExportVarBreakpoint;
	std::vector<std::string> StringExportVarBreakpointTmp;
	std::vector<std::string> StringExportVarBreakpointDistriTmp;
	std::string* PtrStringExportVarBreakpoint;
	std::vector<double*> ExportVarBreakpoint;
	std::vector<double*> ExportVarBreakpointTmp;
	std::vector<double*> ExportVarBreakpointDistriTmp;
	double** PtrExportVarBreakpoint;


//Save which variables need to be synchronised (array of pointers)
	std::vector<double*> SyncVar;
	std::vector<std::string> StringSyncVar;

};

#endif /* SRC_CORE_DICTIONARY_H_ */
