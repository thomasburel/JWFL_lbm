/*
 * CGNS.cpp
 *
 *  Created on: 8 May 2015
 *      Author: thomas
 */

#include "CGNS.h"
#include <fstream>
#include <iostream>

using namespace std;
CGNS::CGNS() :
F(0),B(0),Z(0),E(0),S(0), A(0), Cx(0), Cy(0), Cz(0), Fs(0),
x(0),y(0),z(0),d(0),b(0),e(0),
start_nodes(1),end_nodes(0),start_elems(1),end_elems(0),
ncells(0)
{
	outfile=0;
    /* total number of nodes and hex elements */
    tot_nnodes = (50+1)*(50+1);
    tot_nelems = 50*50;
    NbVariableOutput=0;
    VariableOutput=0;
    VariableOutput=0;
    Dimension=true; // True=2D ; False 3D
    outputfilename="LBM_Output";
}
CGNS::CGNS(SolverEnum::dimension dimension_,std::string outputfilename_, int tot_nodes_, int tot_elems_, int start_nodes_, int end_nodes_, int start_elems_, int end_elems_, const double *x_, const double *y_, const double *z_, int *e_) :
		F(0),B(0),Z(0),E(0),S(0), A(0), Cx(0), Cy(0), Cz(0), Fs(0),x(x_),y(y_),z(z_),d(0),b(0), 	start_nodes(start_nodes_),end_nodes(end_nodes_),start_elems(start_elems_),end_elems(end_elems_)
{
	e=e_;
	outfile=0;
    /* total number of nodes and hex elements */
    tot_nnodes = tot_nodes_;
    tot_nelems = tot_elems_;
    ncells = tot_nelems;
    NbVariableOutput=0;
    VariableOutput=0;
    if(dimension_==SolverEnum::D2)
    	Dimension=true; // True=2D ; False 3D
    else
    	Dimension=false;
    outputfilename=outputfilename_;
	outputfilename+="_";
	Breakpointfilename=outputfilename;
	Breakpointfilename+="Breakpoint.cgns";
	classtype="CGNS";
}
CGNS::~CGNS() {
	//for (int i=0;i<NbVariableOutput;i++)
		delete [] d;
}
void CGNS::Set_solution(double **d_, std::string *str, int nbvar){
	NbVariableOutput=nbvar;
	d=new double* [NbVariableOutput];
	Fs=new int [NbVariableOutput];
	VariableOutput=new std::string [NbVariableOutput];
	for (int i=0;i<NbVariableOutput;i++)
	{
		VariableOutput[i]=str[i];
		d[i]=d_[i];
	}

}
#pragma optimize("", off)
void CGNS::Read_data(double * &d_, std::string variablename, std::string filename){

	int Ftmp=0; int Btmp=0; int Ztmp=0;int Stmp=0; int Fstmp=0;
	char basenametmp[50];
	int cell_dimtmp;
	int phys_dimtmp;
	cgsize_t sizetmp;
	char zonenametmp[50];
	char solnametmp[50];
	char fieldname[50];
	GridLocation_t locationtmp;
	DataType_t datatypetmp;

	char *ary = new char[filename.size()+1];
	std::strcpy ( ary, filename.c_str() );
	char *ary2 = new char[variablename.size()+1];
	std::strcpy ( ary2, variablename.c_str() );

	//ary=filename.c_str();
	int nbases=0;
	int nzones=0;
	int nsols=0;
	int nfields=0;

	if(d_!=0)
	//	d_=new double[end_nodes-start_nodes+1];
	{
//open the cgns file in reading mode
	if (cgp_open(ary, CG_MODE_READ, &Ftmp))
		cgp_error_exit();

//Find the base called "Base"
	if (cg_nbases(Ftmp, &nbases))
		cgp_error_exit();
	if(nbases>0)
	{
		for(int i=1;i<nbases+1;i++)
		{
			if (cg_base_read(Ftmp, i, &basenametmp[0], &cell_dimtmp,&phys_dimtmp))
				cgp_error_exit();
			if(strcmp (basenametmp,"Base")==0)
				Btmp=i;
		}
		if(Btmp>0)
		{
			//std::cout<<"Base is: "<<Btmp<<std::endl;
	//Get number of zones
			if (cg_nzones(Ftmp, Btmp, &nzones))
				cgp_error_exit();
			if(nzones>0)
			{
		//Find the zone called "Zone"
				for(int i=1;i<nzones+1;i++)
				{
					if (cg_zone_read(Ftmp, Btmp, i, &zonenametmp[0], &sizetmp))
						cgp_error_exit();
					if(strcmp (zonenametmp,"Zone")==0)
						Ztmp=i;
				}
				if(Ztmp>0)
				{
					//std::cout<<"Zone is: "<<Ztmp<<std::endl;
		//Get number of solutions
					if(cg_nsols(Ftmp, Btmp, Ztmp, &nsols))
						cgp_error_exit();
					if(nsols>0)
					{
						//ary=variablename.c_str();
			//Find the variable in the solution list
						for(int i=1;i<nsols+1;i++)
						{
							if (cg_sol_info(Ftmp, Btmp, Ztmp, i, &solnametmp[0], &locationtmp))
								cg_error_exit();
							if(strcmp ("Solution",solnametmp)==0)
								Stmp=i;
						}
			//Read the data
						if(Stmp>0)
						{
							//std::cout<<"Solution is: "<<Stmp<<std::endl;
							cg_nfields(Ftmp, Btmp, Ztmp, Stmp, &nfields);
							if(nfields>0)
							{
								for(int i=1;i<nfields+1;i++)
								{
									if (cg_field_info(Ftmp, Btmp, Ztmp, Stmp, i, &datatypetmp, &fieldname[0]))
										cg_error_exit();
									if(strcmp (ary2,fieldname)==0)
										Fstmp=i;
								}
								//std::cout<<"Field is: "<<Fstmp<<std::endl;
								if(Fstmp>0)
								{
								//	std::cout<<"Start node: "<<start_nodes<<" End node: "<<end_nodes<<std::endl;
									if (cgp_field_read_data(Ftmp, Btmp, Ztmp, Stmp, Fstmp, &start_nodes, &end_nodes, d_))
										cgp_error_exit();
								}
							}

						}
					}
				}
				else
					std::cerr<<"No Zone found in the file"<<std::endl;
			}
		}
		else
			std::cerr<<"Base is unfound in the file"<<std::endl;
	}
	else
		std::cerr<<"No Base found in the file"<<std::endl;

    /* close the file and terminate MPI */
    cgp_close(Ftmp);
	}
	else
		std::cerr<<"Pointer not initialised before passing in the function Read_data"<<std::endl;

	delete ary;delete ary2;
}
#pragma optimize("", on)
void CGNS::Write_Output(int& Noutput) {
   /*  Tmp data   */

	//end=1;

    /* Create the File name   */
	std::string str =outputfilename;
	std::stringstream convert; // stringstream used for the conversion
	convert << Noutput;//add the value of Number to the characters in the stream
	str+= convert.str();//set Result to the content of the stream
	str+=".cgns";
	char *ary = new char[str.size()+1];
	std::strcpy ( ary, str.c_str() );
	outfile=ary;



    /* open the file and create base and zone */
    sizes[0] = tot_nnodes; //Number of Vertex
    sizes[1] = tot_nelems; //Number of Cell
    sizes[2] = 0; //Number of Bound Vertex

    /* the default here is to use MPI_COMM_WORLD,
       but this allows assigning of another communicator
    cgp_mpi_comm(MPI_COMM_WORLD); */

    if (cgp_open(outfile, CG_MODE_WRITE, &F)|| cg_base_write(F, "Base", 2, 2, &B))
        cgp_error_exit();
    if (cg_zone_write(F, B, "Zone", sizes, Unstructured, &Z))
    	cgp_error_exit();
    /* print info */

//        printf("writing %d coordinates and %d hex elements to %s\n start is %d and end is %d\n",
//            tot_nnodes, tot_nelems, outfile,start_nodes,end_nodes);
//        std::cout<<x<<std::endl;

        /* create data nodes for coordinates */
        if (cgp_coord_write(F, B, Z, RealDouble, "CoordinateX", &Cx) ||
            cgp_coord_write(F, B, Z, RealDouble, "CoordinateY", &Cy)/* ||
            cgp_coord_write(F, B, Z, RealSingle, "CoordinateZ", &Cz)*/)
            cgp_error_exit();

        /* write the coordinate data in parallel to the queue */
          if (cgp_queue_set(1) ||
            cgp_coord_write_data(F, B, Z, Cx, &start_nodes, &end_nodes, x) ||
            cgp_coord_write_data(F, B, Z, Cy, &start_nodes, &end_nodes, y) /*||
            cgp_coord_write_data(F, B, Z, Cz, &start, &end, z)*/)
            cgp_error_exit();

        /* write out the queued coordinate data */
       if (cgp_queue_flush()) cgp_error_exit();
        cgp_queue_set(0);

        /* create data node for elements */
        if (cgp_section_write(F, B, Z, "Quad", QUAD_4, 1, tot_nelems, 0, &E)) // Hexa "Hex", HEXA_8
            cgp_error_exit();

        /* write the element connectivity in parallel */
      if (cgp_elements_write_data(F, B, Z, E, start_elems, end_elems, e))
           cgp_error_exit();


      /* create a vertex solution */
      if (cg_sol_write(F, B, Z, "Solution", Vertex, &S))
          cgp_error_exit();


      /* write the solution field data in parallel */
      if(d!=0)
      {
    	  for (int i=0;i<NbVariableOutput;i++)
    	  {
			  char buffer[50];
			  std::strcpy(buffer,VariableOutput[i].c_str());
			  if (cgp_field_write(F, B, Z, S, RealDouble, buffer , &Fs[i]))
				  cgp_error_exit();
			  if (cgp_field_write_data(F, B, Z, S, Fs[i], &start_nodes, &end_nodes, d[i]))
				  cgp_error_exit();
    	  }
      }
        /* close the file and terminate MPI */
        cgp_close(F);
        outfile=0;
        delete ary;
}
void CGNS::Write_breakpoint(Parameters &Parameters){
    /* Create the File name   */
	char *ary = new char[Breakpointfilename.size()+1];
	std::strcpy ( ary, Breakpointfilename.c_str() );
	outfile=ary;
	//delete ary;



    /* open the file and create base and zone */
    sizes[0] = tot_nnodes; //Number of Vertex
    sizes[1] = tot_nelems; //Number of Cell
    sizes[2] = 0; //Number of Bound Vertex

    /* the default here is to use MPI_COMM_WORLD,
       but this allows assigning of another communicator
    cgp_mpi_comm(MPI_COMM_WORLD); */

    if (cgp_open(outfile, CG_MODE_WRITE, &F))
        cgp_error_exit();
    if(cg_base_write(F, "Base", 2, 2, &B))
    	cgp_error_exit();
    CGNS::Write_ParametersProperties(Parameters);
    if (cg_zone_write(F, B, "Zone", sizes, Unstructured, &Z))
    	cgp_error_exit();
    /* print info */

//        printf("writing %d coordinates and %d hex elements to %s\n start is %d and end is %d\n",
//            tot_nnodes, tot_nelems, outfile,start_nodes,end_nodes);
//        std::cout<<x<<std::endl;

        /* create data nodes for coordinates */
        if (cgp_coord_write(F, B, Z, RealDouble, "CoordinateX", &Cx) ||
            cgp_coord_write(F, B, Z, RealDouble, "CoordinateY", &Cy)/* ||
            cgp_coord_write(F, B, Z, RealSingle, "CoordinateZ", &Cz)*/)
            cgp_error_exit();

        /* write the coordinate data in parallel to the queue */
          if (cgp_queue_set(1) ||
            cgp_coord_write_data(F, B, Z, Cx, &start_nodes, &end_nodes, x) ||
            cgp_coord_write_data(F, B, Z, Cy, &start_nodes, &end_nodes, y) /*||
            cgp_coord_write_data(F, B, Z, Cz, &start, &end, z)*/)
            cgp_error_exit();

        /* write out the queued coordinate data */
       if (cgp_queue_flush()) cgp_error_exit();
        cgp_queue_set(0);

        /* create data node for elements */
        if (cgp_section_write(F, B, Z, "Quad", QUAD_4, 1, tot_nelems, 0, &E)) // Hexa "Hex", HEXA_8
            cgp_error_exit();

        /* write the element connectivity in parallel */
      if (cgp_elements_write_data(F, B, Z, E, start_elems, end_elems, e))
           cgp_error_exit();


      /* create a centered solution */
      if (cg_sol_write(F, B, Z, "Solution", Vertex, &S))
          cgp_error_exit();


      /* write the solution field data in parallel */
      if(b!=0)
      {
    	  for (int i=0;i<NbVariableBreakpoint;i++)
    	  {
			  char buffer[50];
			  std::strcpy(buffer,VariableBreakpoint[i].c_str());
			  if (cgp_field_write(F, B, Z, S, RealDouble, buffer , &Fs[i]))
				  cgp_error_exit();
			  if (cgp_field_write_data(F, B, Z, S, Fs[i], &start_nodes, &end_nodes, b[i]))
				  cgp_error_exit();
    	  }
      }
        /* close the file and terminate MPI */
        cgp_close(F);

        outfile=0;
        delete ary;
}
void CGNS::Set_breakpoint(double **b_, std::string *str, int nbvar){
	NbVariableBreakpoint=nbvar;
	b=new double* [NbVariableBreakpoint];
	VariableBreakpoint=new std::string [NbVariableBreakpoint];
	for (int i=0;i<NbVariableBreakpoint;i++)
	{
		VariableBreakpoint[i]=str[i];
		b[i]=b_[i];
	}
}
void CGNS::Write_ParametersProperties(Parameters &Param)
{
int rank,size;
MPI_Comm_size(MPI_COMM_WORLD,&size);
if (MPI_Comm_rank(MPI_COMM_WORLD,&rank)==0)
{
	//cg_open(outfile, CG_MODE_MODIFY, &F);
	cg_goto(F, B, "end");
	cg_descriptor_write("About", "This is a breakpoint created by the in-house LBM code written at JWFL and by Thomas BUREL");
	cg_descriptor_write("Generating by", "LBM JWFL");
	std::string str;
	std::stringstream convert; // stringstream used for the conversion
	convert << size;//add the value of Number to the characters in the stream
	str+= convert.str();//set Result to the content of the stream
	char *ary = new char[str.size()+1];
	std::strcpy ( ary, str.c_str() );
	cg_descriptor_write("NbPartition", ary);
}
}
