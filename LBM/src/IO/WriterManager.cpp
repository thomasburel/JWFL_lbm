/*
 * WriterManager.cpp
 *
 *  Created on: 10 May 2015
 *      Author: thomas
 */

#include "WriterManager.h"

WriterManager::WriterManager() :
tot_nnodes(0), tot_nelems(0), nnodes(0), nelems(0),outfile(0),VariableOutput(0),NbVariableOutput(0),VariableBreakpoint(0),NbVariableBreakpoint(0)
{
    Dimension=true; // True=2D ; False 3D
    outputfilename="LBM_Output";
}

WriterManager::~WriterManager() {
	delete outfile,VariableOutput;
}

