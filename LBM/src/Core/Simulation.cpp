/*
 * Simulation.cpp
 *
 *  Created on: 17 Apr 2015
 *      Author: thomas
 */

#include "Simulation.h"

Simulation::Simulation()
:time(0),Solver_(0),MultiBlock_(0),parallel(0),Writer(0),PtrParameters(0),WriterCGNS(0),WriterTecplot(0),MultiBlock2D_(0),SolverD2Q9(0),SolverD2Q9TwoPhases(0)
{}
void Simulation::InitSimu(Parameters &Parameters_, bool create_mesh)
{
	PtrParameters=&Parameters_;
	PtrParameters->CheckParameters();
	parallel=new MpiManager;
	ini.Set_Parameters(PtrParameters);
	ini.IniMPI(parallel,PtrParameters->Get_Argc(), PtrParameters->Get_Argv(), PtrParameters->Get_Verbous());
	time=get_time();

	if(create_mesh)
	{
		Create_Mesh();
	}
	else
	{
		Import_Mesh(PtrParameters->Get_MeshFile());
	}
	MultiBlock_->GeneratePatchBc();
	MultiBlock_->Modify_Block();
	MultiBlock_->ConvertToPhysicalUnit();
	MultiBlock_->reorganizeNodeByType();
	SelectSolver();
	//SolverD2Q9=new D2Q9(MultiBlock_,parrallel,Writer,PtrParameters,InitSimu());
	//Solver_= SolverD2Q9;
	//ini.IniDomain(Parameters_.Get_Nx(),Parameters_.Get_Ny(),Parameters_.Get_Nz());
	if (parallel->isMainProcessor())
	cout<< "Time to initialise the simulation:    "<<get_time()-time << endl;
	time=get_time();
}


Simulation::~Simulation() {
/*	Simulation::FinalizeSimu();
	delete MultiBlock_;
	//delete Writer;
	delete SolverD2Q9;
	delete SolverD2Q9TwoPhases;
	delete Solver_;*/
}
/*
void Simulation::Test_ReadData(double * &d_, std::string variablename, std::string filename){
	Writer->Read_data(d_,variablename,filename);
}*/
double Simulation::get_time() {
		return parallel->getTime();
}
InitLBM& Simulation::InitSimu(){
	return ini;
}
void Simulation::Create_Mesh(){ //(dimension dim, int Nx, int Ny, int Nz) {
	if(PtrParameters->Get_Dimension()==SolverEnum::D2)
	{

		//MultiBlock_=new MultiBlock2D(parrallel,PtrParameters->Get_Nx(),PtrParameters->Get_Ny());
		MultiBlock2D_=new MultiBlock2D(parallel,PtrParameters);
		MultiBlock_=MultiBlock2D_;
		MultiBlock_->Partitioning();

		int tot_nodes_=(PtrParameters->Get_Nx()+1)*(PtrParameters->Get_Ny()+1);//MultiBlock_->Get_End_Nodes()-MultiBlock_->Get_Start_Nodes()+1;//(Nx+1)*(Ny+1);
		int tot_elems_=PtrParameters->Get_Nx()*PtrParameters->Get_Ny();//MultiBlock_->Get_End_Elems()-MultiBlock_->Get_Start_Elems()+1;//Nx*Ny;

		SelectOutputFormat();

		WriterCGNS=new CGNS(PtrParameters->Get_Dimension(),PtrParameters->Get_OutputFileName(),tot_nodes_,tot_elems_,MultiBlock_->Get_Start_Nodes(),MultiBlock_->Get_End_Nodes(),MultiBlock_->Get_Start_Elems(),MultiBlock_->Get_End_Elems(),MultiBlock_->Get_X0(),MultiBlock_->Get_Y0(),0,MultiBlock_->Get_Elems0());
		Writer=WriterCGNS;
		int it=0;
		if(PtrParameters->Get_Verbous())
		{
			Writer->Write_Output(it);
			if(parallel->isMainProcessor())
				std::cout<<"End of Writing of the mesh"<<std::endl;
		}
	}
	else
	{
//		MultiBlock_=new MultiBlock3D;
	}
	Simulation::barrier();
}

void Simulation::Import_Mesh(string MeshFile_) {

}
void Simulation::FinalizeSimu() {
	if (parallel->isMainProcessor())
		cout<<"Deleting Solvers"<<endl;
	Solver_->~Solver();
	if (parallel->isMainProcessor())
		cout<<"Deleting Writers"<<endl;
	Writer->~WriterManager();
	if (parallel->isMainProcessor())
		cout<<"Deleting Multiblocks"<<endl;
	MultiBlock_->~MultiBlock();
	if (parallel->isMainProcessor())
		cout<<"Deleting Parallel manager"<<endl;
	parallel->~ParallelManager();

}
void Simulation::barrier() {
	parallel->barrier();
}
void Simulation::RunSimu()
{
	Solver_->Set_UserConvergence(PtrParameters);
	Solver_->SetUserForce(PtrParameters);
	Solver_->run();
	if (parallel->isMainProcessor())
	cout<< "Time to run the simulation:    "<<get_time()-time << endl;
	time=get_time();
	Save_Simulation();
	cout<< " End Simulation " << endl;
}
void Simulation::RunSimu(Parameters &UpdatedParam)
{
	PtrParameters=&UpdatedParam;
	Solver_->Set_UserConvergence(PtrParameters);
	Solver_->SetUserForce(PtrParameters);
	Writer->UpdateFileNames(PtrParameters->Get_OutputFileName());
	Solver_->run(PtrParameters);
	if (parallel->isMainProcessor())
	cout<< "Time to run the simulation:    "<<get_time()-time << endl;
	time=get_time();

}
void Simulation::UpdateAllDomainFromFile(Parameters &UpdatedParam){
	ini.Set_Parameters(&UpdatedParam);
	Solver_->UpdateAllDomainFromFile(&UpdatedParam,ini);
}
void Simulation::UpdateAllDomain(Parameters &UpdatedParam){
	ini.Set_Parameters(&UpdatedParam);
	Solver_->UpdateAllDomain(&UpdatedParam,ini);
}
void Simulation::UpdateDomainBc(Parameters &UpdatedParam){
	ini.Set_Parameters(&UpdatedParam);
	Solver_->UpdateDomainBc(&UpdatedParam,ini);
}
void Simulation::UpdateWall(Parameters &UpdatedParam){
	ini.Set_Parameters(&UpdatedParam);
	Solver_->UpdateWall(&UpdatedParam,ini);
}
void Simulation::UpdateInterior(Parameters &UpdatedParam){
	ini.Set_Parameters(&UpdatedParam);
	Solver_->UpdateInterior(&UpdatedParam,ini);
}
void Simulation::SelectSolver()
{
	if(PtrParameters->Get_Dimension()==SolverEnum::D2)
	{
		switch(PtrParameters->Get_Model())
		{
		case SolverEnum::SinglePhase:
			SolverD2Q9=new D2Q9(MultiBlock_,parallel,Writer,PtrParameters,InitSimu());
			Solver_= SolverD2Q9;
			break;
		case SolverEnum::ColourFluid:
			//SolverD2Q9=new D2Q9(MultiBlock_,parrallel,Writer,PtrParameters,InitSimu());
			SolverD2Q9TwoPhases=new D2Q9ColourFluid(MultiBlock_,parallel,Writer,PtrParameters,InitSimu());
			Solver_= SolverD2Q9TwoPhases;
			break;
		default:
			std::cerr<<" Model not found"<<std::endl;

			break;
		}
	}
	else
		std::cerr<<" Dimension not found"<<std::endl;
}
void Simulation::SelectOutputFormat()
{
	switch(PtrParameters->Get_OutputFormat())
	{
	case CGNSFormat:
		WriterCGNS=new CGNS();//new CGNS(PtrParameters->Get_Dimension(),PtrParameters->Get_OutputFileName(),tot_nodes_,tot_elems_,MultiBlock_->Get_Start_Nodes(),MultiBlock_->Get_End_Nodes(),MultiBlock_->Get_Start_Elems(),MultiBlock_->Get_End_Elems(),MultiBlock_->Get_X0(),MultiBlock_->Get_Y0(),0,MultiBlock_->Get_Elems0());
		Writer=WriterCGNS;
		break;
	case TecplotFormat:
		WriterTecplot=new TecplotIO();
		Writer=WriterTecplot;
		break;
	default:
		std::cerr<<" Output Format not found"<<std::endl;

		break;
	}
}
void Simulation::Save_Simulation(){
	if (parallel->isMainProcessor())
	{
	std::string filename("Save_Simulation.txt");
	std::ofstream ofs(filename.c_str()); assert(ofs.good());
	boost::archive::text_oarchive oa(ofs);
	std::string filenameParam("Save_Parameters.xml");
	std::ofstream ofsParam(filenameParam.c_str()); assert(ofsParam.good());
	boost::archive::xml_oarchive oaParam(ofsParam);
	Save_Parameters(*PtrParameters,filenameParam);
	Save_InitLBM(ini,filename);
	Save_Writer(filename,ios::app);
	Save_MultiBlock(filename,ios::app);
	Save_Solver(filename,ios::app);
	}
}

void Simulation::Save_Parameters(Parameters & object){

	std::string filename("Parameters_save.xml");
	std::cout<< "Parameters are saving in the xml file: "<<filename<<std::endl;
	// make an archive
	std::ofstream ofs(filename.c_str());
	assert(ofs.good());
	boost::archive::xml_oarchive oa(ofs);
	oa << BOOST_SERIALIZATION_NVP(object);
}
void Simulation::Save_Parameters(Parameters & object,std::string filename){
	std::cout<< "Parameters are saving in the xml file: "<<filename<<std::endl;
	// make an archive
	std::ofstream ofs(filename.c_str());
	assert(ofs.good());
	boost::archive::xml_oarchive oa(ofs);
	oa << BOOST_SERIALIZATION_NVP(object);
}

void Simulation::Load_Parameters(Parameters & object){
	std::string filename("Parameters_save.xml");
	std::cout<< "Parameters are loading from the xml file: "<<filename<<std::endl;
	// open an archive
	std::ifstream ifs(filename.c_str());
	assert(ifs.good());
	boost::archive::xml_iarchive ia(ifs);
	ia >> BOOST_SERIALIZATION_NVP(object);
	object.Set_VariablesOutput(object.Get_output_density(),object.Get_output_velocity(),object.Get_output_pressure());
}





void Simulation::Save_InitLBM(InitLBM & object){
	std::string filename("Init_save.txt");
	std::cout<< "the object InitLBM is saving in the text file: "<<filename<<std::endl;
	// make an archive
	std::ofstream ofs(filename.c_str());
	boost::archive::text_oarchive oa(ofs);
	oa << object;
}
void Simulation::Save_InitLBM(InitLBM & object,std::string &filename){
	std::cout<< "the object InitLBM is saving in the text file: "<<filename<<std::endl;
	// make an archive
	std::ofstream ofs(filename.c_str(),ios::trunc);
	boost::archive::text_oarchive oa(ofs);
	oa << object;
}

void Simulation::Save_Writer(std::string &filename,ios_base::openmode mode){
	std::cout<< "the object Writer is saving in the text file: "<<filename<<std::endl;
	// make an archive
	std::ofstream ofs(filename.c_str(),mode);
	boost::archive::text_oarchive oa(ofs);
	if (Writer->Get_classtype()=="CGNS")
		oa << *WriterCGNS;
	else
		oa << *Writer;
}
void Simulation::Save_MultiBlock(std::string &filename,ios_base::openmode mode){
	std::cout<< "the object MultiBlock is saving in the text file: "<<filename<<std::endl;
	// make an archive
	std::ofstream ofs(filename.c_str(),mode);
	boost::archive::text_oarchive oa(ofs);
	if (PtrParameters->Get_Dimension()==SolverEnum::D2)
		oa << *MultiBlock2D_;
	else
		oa << *MultiBlock_;
}
void Simulation::Save_Solver(std::string &filename,ios_base::openmode mode){
	std::cout<< "the object Solver is saving in the text file: "<<filename<<std::endl;
	// make an archive
	std::ofstream ofs(filename.c_str(),mode);
	boost::archive::text_oarchive oa(ofs);
	if (PtrParameters->Get_Dimension()==SolverEnum::D2){
		if(PtrParameters->Get_Scheme()==SolverEnum::Q9){
			if(PtrParameters->Get_Model()==SolverEnum::SinglePhase)
				oa << *SolverD2Q9;
			else
				oa << *SolverD2Q9TwoPhases;
		}
	}
	else
		oa << *Solver_;
}
