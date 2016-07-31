/*
 * ParrallelManager.h
 *
 *  Created on: 15 Apr 2015
 *      Author: thomas
 */

#ifndef PARRALLELISM_PARRALLELMANAGER_H_
#define PARRALLELISM_PARRALLELMANAGER_H_
#include "mpi.h"
// Abstract class
class ParallelManager {
public:
	ParallelManager();
	virtual ~ParallelManager();
    virtual void init(int *argc, char ***argv, bool verbous=false)=0;
    /// Initializes the MPI manager, but assumes that the MPI
    ///   machine is handled by another instance.
    virtual void init(MPI_Comm globalCommunicator_)=0;
    /// Initializes the MPI manager, but assumes that the MPI
    ///   machine is handled by another instance.
    virtual void init()=0;
    /// Returns the number of processes (pure virtual method)
    virtual int getSize() const = 0;
    /// Returns the process ID (pure virtual method)
    virtual int getRank() const = 0;

    virtual int bossId() const=0;
    /// Tells whether current processor is main processor
    virtual bool isMainProcessor() const=0;
    /// Returns universal MPI-time in seconds
    virtual double getTime() const=0;
    /// Returns the global communicator for this program or library instance.
    virtual MPI_Comm getGlobalCommunicator() const=0;

    /// Synchronizes the processes
    virtual void barrier()=0;

    virtual void sendRecv(double *sendBuf, double *recvBuf, int count_dest, int dest, int count_source, int source, int tag=0)=0;

protected:
    int numTasks, taskId;
};

#endif /* PARRALLELISM_PARRALLELMANAGER_H_ */
