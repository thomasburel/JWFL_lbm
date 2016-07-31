/*
 * MPIManager.h
 *
 *  Created on: 15 Apr 2015
 *      Author: thomas
 *      Modified the Palabos file to remove complex number
 */
/* This file is part of the Palabos library.
 *
 * Copyright (C) 2011-2015 FlowKit Sarl
 * Route d'Oron 2
 * 1010 Lausanne, Switzerland
 * E-mail contact: contact@flowkit.com
 *
 * The most recent release of Palabos can be downloaded at
 * <http://www.palabos.org/>
 *
 * The library Palabos is free software: you can redistribute it and/or
 * modify it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * The library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

/** \file
 * Wrapper functions that simplify the use of MPI
 */

#ifndef PARRALLELISM_LBMMPI_H_
#define PARRALLELISM_LBMMPI_H_

#include <vector>
#include <string>
#include <iostream>

#include "ParallelManager.h"
class MpiManager: public ParallelManager {
public:
    /// Initializes the MPI manager and the MPI machine.
    virtual void init(int *argc, char ***argv, bool verbous=false);
    /// Initializes the MPI manager, but assumes that the MPI
    ///   machine is handled by another instance.
    virtual void init(MPI_Comm globalCommunicator_);
    /// Initializes the MPI manager, but assumes that the MPI
    ///   machine is handled by another instance.
    virtual void init();
    /// Returns the number of processes
    virtual int getSize() const;
    /// Returns the process ID
    virtual int getRank() const;
    /// Returns process ID of main processor
    virtual int bossId() const;
    /// Tells whether current processor is main processor
    virtual bool isMainProcessor() const;
    /// Returns universal MPI-time in seconds
    virtual double getTime() const;
    /// Returns the global communicator for this program or library instance.
    virtual MPI_Comm getGlobalCommunicator() const;

    /// Synchronizes the processes
    virtual void barrier();

    /// Sends data at *buf, blocking
    template <typename T>
    void send( T *buf, int count, int dest, int tag = 0 );

    /// Sends data at *buf, non blocking
    template <typename T>
    void iSend( T *buf, int count, int dest, MPI_Request* request, int tag = 0 );

    /// Sends data at *buf, assuming that receiver is ready.
    template <typename T>
    void rSend( T *buf, int count, int dest, int tag = 0 );

    /// Sends data at *buf, non blocking and request free
    template <typename T>
    void iSendRequestFree( T *buf, int count, int dest, int tag = 0 );

    /// Receives data at *buf, blocking
    template <typename T>
    void receive( T *buf, int count, int source, int tag = 0 );

    /// Receives data at *buf, non blocking
    template <typename T>
    void iRecv( T *buf, int count, int source, MPI_Request* request, int tag = 0 );

    /// Send and receive data between two partners
    void sendRecv( char *sendBuf, char *recvBuf, int count, int dest,
                   int source, int tag = 0 );
    void sendRecv( int *sendBuf, int *recvBuf, int count, int dest,
                      int source, int tag = 0 );
    void sendRecv( bool *sendBuf, bool *recvBuf, int count, int dest,
                      int source, int tag = 0 );
    void sendRecv( long *sendBuf, long *recvBuf, int count, int dest,
                      int source, int tag = 0 );
    void sendRecv( float *sendBuf, float *recvBuf, int count, int dest,
                      int source, int tag = 0 );
    void sendRecv( double *sendBuf, double *recvBuf, int count_dest, int dest,
    		int count_source,int source, int tag = 0 );
    void sendRecv( long double *sendBuf, long double *recvBuf, int count, int dest,
                      int source, int tag = 0 );


    /// Send and receive vector between two partners
    template <typename T>
    void sendRecvVect(std::vector<T>& sendBuf, std::vector<T>& recvBuf, int count, int dest,
                   int source, int tag = 0 );

    /// Sends data to master processor
    template <typename T>
    void sendToMaster( T* sendBuf, int sendCount, bool iAmRoot );

    /// Scatter data from one processor over multiple processors
    template <typename T>
    void scatterV( T *sendBuf, T *recvBuf, int* sendCounts, int root = 0 );

    /// Gather data from multiple processors to one processor
    template <typename T>
    void gatherV( T* sendBuf, T* recvBuf, int *recvCounts, int root = 0 );


    /// Broadcast data from one processor to multiple processors
    template <typename T>
    void bCast( T* sendBuf, int sendCount, int root = 0 );

    /// Special case for broadcasting strings. Memory handling is automatic.
    void bCast( std::string& message, int root = 0 );

    /// Broadcast data when root is unknown to other processors
    template <typename T>
    void bCastThroughMaster( T* sendBuf, int sendCount, bool iAmRoot );

    /// Reduction operation toward one processor
    template <typename T>
    void reduce( T sendVal, T& recvVal, MPI_Op op, int root = 0 );

    /// Element-per-element reduction of a vector of data
    template <typename T>
    void reduceVect( std::vector<T>& sendVal, std::vector<T>& recvVal,
                     MPI_Op op, int root = 0 );

    /// Inplace element-per-element reduction of a vector of data; result
    ///   available on all MPI threads.
    template <typename T>
    void allReduceVect( std::vector<T>& sendRecvVal, MPI_Op op );

    /// Reduction operation, followed by a broadcast
    template <typename T>
    void reduceAndBcast(T& reductVal, MPI_Op op, int root = 0 );

    /// Complete a non-blocking MPI operation
    void wait(MPI_Request* request, MPI_Status* status);

private:
    /// Implementation code for Scatter
    template <typename T>
    void scatterv_impl(T *sendBuf, int* sendCounts, int* displs,
                       T* recvBuf, int recvCount, int root);

    /// Implementation code for Gather
    template <typename T>
    void gatherv_impl(T* sendBuf, int sendCount, T* recvBuf, int* recvCounts,
                      int* displs, int root);
//private:
public:
    MpiManager();
    ~MpiManager();
private:
//    int numTasks, taskId;
    bool ok;
    bool responsibleForMpiMachine;
    MPI_Comm globalCommunicator;

friend MpiManager& mpi();
};

inline MpiManager& mpi() {
    static MpiManager instance;
    return instance;
}

#endif /* PARRALLELISM_MPIMANAGER_H_ */
