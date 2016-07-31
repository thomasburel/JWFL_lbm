/*
 * MpiManager.cpp
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
 * Wrapper functions that simplify the use of MPI, template instatiations
 */

#include "../Parallelism/MpiManager.h"


MpiManager::MpiManager()
    : ok(false),
      responsibleForMpiMachine(false),
	  globalCommunicator(0)
{ }

MpiManager::~MpiManager() {
    if (responsibleForMpiMachine) {
        MPI_Finalize();
        ok = false;
        responsibleForMpiMachine = false;
    }
}

void MpiManager::init(int *argc, char ***argv, bool verbous) {
    if (verbous) {
        std::cerr << "Constructing an MPI thread" << std::endl;
    }
    int ok1 = MPI_Init(argc, argv);
    // If I'm the one who calls MPI_Init, then I need to be
    // the one who calls MPI_Finalize.
    responsibleForMpiMachine = true;
    globalCommunicator = MPI_COMM_WORLD;
    int ok2 = MPI_Comm_rank(getGlobalCommunicator(),&taskId);
    int ok3 = MPI_Comm_size(getGlobalCommunicator(),&numTasks);
    ok = (ok1==0 && ok2==0 && ok3==0);
}

void MpiManager::init(MPI_Comm globalCommunicator_) {
    globalCommunicator = globalCommunicator_;
    int ok1 = MPI_Comm_rank(getGlobalCommunicator(),&taskId);
    int ok2 = MPI_Comm_size(getGlobalCommunicator(),&numTasks);
    ok = (ok1==0 && ok2==0);
}

void MpiManager::init() {
    init(MPI_COMM_WORLD);
}


int MpiManager::getSize() const {
    return numTasks;
}

int MpiManager::getRank() const {
    return taskId;
}

int MpiManager::bossId() const {
    return 0;
}

bool MpiManager::isMainProcessor() const {
    return bossId() == getRank();
}

double MpiManager::getTime() const {
    if (!ok) return 0.;
    return MPI_Wtime();
}

MPI_Comm MpiManager::getGlobalCommunicator() const {
    return globalCommunicator;
}

void MpiManager::barrier() {
    if (!ok) return;
    MPI_Barrier(getGlobalCommunicator());
}

template <>
void MpiManager::send<char>(char *buf, int count, int dest, int tag) {
    if (!ok) return;
    MPI_Send(static_cast<void*>(buf), count, MPI_CHAR, dest, tag, getGlobalCommunicator());
}

template <>
void MpiManager::send<int>(int *buf, int count, int dest, int tag) {
    if (!ok) return;
    MPI_Send(static_cast<void*>(buf), count, MPI_INT, dest, tag, getGlobalCommunicator());
}

template <>
void MpiManager::send<bool>(bool *buf, int count, int dest, int tag) {
    if (!ok) return;
    MPI_Send(static_cast<void*>(buf), count*sizeof(bool), MPI_CHAR, dest, tag, getGlobalCommunicator());
}


template <>
void MpiManager::send<long>(long *buf, int count, int dest, int tag) {
    if (!ok) return;
    MPI_Send(static_cast<void*>(buf), count, MPI_LONG, dest, tag, getGlobalCommunicator());
}

template <>
void MpiManager::send<float>(float *buf, int count, int dest, int tag) {
    if (!ok) return;
    MPI_Send(static_cast<void*>(buf), count, MPI_FLOAT, dest, tag, getGlobalCommunicator());
}

template <>
void MpiManager::send<double>(double *buf, int count, int dest, int tag) {
    if (!ok) return;
    MPI_Send(static_cast<void*>(buf), count, MPI_DOUBLE, dest, tag, getGlobalCommunicator());
}

template <>
void MpiManager::send<long double>(long double *buf, int count, int dest, int tag) {
    if (!ok) return;
    MPI_Send(static_cast<void*>(buf), count, MPI_LONG_DOUBLE, dest, tag, getGlobalCommunicator());
}

template <>
void MpiManager::iSend<char>
    (char *buf, int count, int dest, MPI_Request* request, int tag)
{
    if (ok) {
        MPI_Isend(static_cast<void*>(buf), count, MPI_CHAR, dest, tag, getGlobalCommunicator(), request);
    }
}

template <>
void MpiManager::iSend<int>
    (int *buf, int count, int dest, MPI_Request* request, int tag)
{
    if (ok) {
        MPI_Isend(static_cast<void*>(buf), count, MPI_INT, dest, tag, getGlobalCommunicator(), request);
    }
}

template <>
void MpiManager::iSend<bool>
    (bool *buf, int count, int dest, MPI_Request* request, int tag)
{
    if (ok) {
        MPI_Isend(static_cast<void*>(buf), count*sizeof(bool), MPI_CHAR, dest, tag, getGlobalCommunicator(), request);
    }
}


template <>
void MpiManager::iSend<long>
    (long *buf, int count, int dest, MPI_Request* request, int tag)
{
    if (ok) {
        MPI_Isend(static_cast<void*>(buf), count, MPI_LONG, dest, tag, getGlobalCommunicator(), request);
    }
}

template <>
void MpiManager::iSend<float>
    (float *buf, int count, int dest, MPI_Request* request, int tag)
{
    if (ok) {
        MPI_Isend(static_cast<void*>(buf), count, MPI_FLOAT, dest, tag, getGlobalCommunicator(), request);
    }
}

template <>
void MpiManager::iSend<double>
    (double *buf, int count, int dest, MPI_Request* request, int tag)
{
    if (ok) {
        MPI_Isend(static_cast<void*>(buf), count, MPI_DOUBLE, dest, tag, getGlobalCommunicator(), request);
    }
}

template <>
void MpiManager::iSend<long double>
    (long double *buf, int count, int dest, MPI_Request* request, int tag)
{
    if (ok) {
        MPI_Isend(static_cast<void*>(buf), count, MPI_LONG_DOUBLE, dest, tag, getGlobalCommunicator(), request);
    }
}

template <>
void MpiManager::rSend<char>(char *buf, int count, int dest, int tag) {
    if (!ok) return;
    MPI_Rsend(static_cast<void*>(buf), count, MPI_CHAR, dest, tag, getGlobalCommunicator());
}

template <>
void MpiManager::rSend<int>(int *buf, int count, int dest, int tag) {
    if (!ok) return;
    MPI_Rsend(static_cast<void*>(buf), count, MPI_INT, dest, tag, getGlobalCommunicator());
}

template <>
void MpiManager::rSend<bool>(bool *buf, int count, int dest, int tag) {
    if (!ok) return;
    MPI_Rsend(static_cast<void*>(buf), count*sizeof(bool), MPI_CHAR, dest, tag, getGlobalCommunicator());
}


template <>
void MpiManager::rSend<long>(long *buf, int count, int dest, int tag) {
    if (!ok) return;
    MPI_Rsend(static_cast<void*>(buf), count, MPI_LONG, dest, tag, getGlobalCommunicator());
}

template <>
void MpiManager::rSend<float>(float *buf, int count, int dest, int tag) {
    if (!ok) return;
    MPI_Rsend(static_cast<void*>(buf), count, MPI_FLOAT, dest, tag, getGlobalCommunicator());
}

template <>
void MpiManager::rSend<double>(double *buf, int count, int dest, int tag) {
    if (!ok) return;
    MPI_Rsend(static_cast<void*>(buf), count, MPI_DOUBLE, dest, tag, getGlobalCommunicator());
}

template <>
void MpiManager::rSend<long double>(long double *buf, int count, int dest, int tag) {
    if (!ok) return;
    MPI_Rsend(static_cast<void*>(buf), count, MPI_LONG_DOUBLE, dest, tag, getGlobalCommunicator());
}


template <>
void MpiManager::iSendRequestFree<char>
    (char *buf, int count, int dest, int tag)
{
    MPI_Request request;
    if (ok) {
        MPI_Isend(static_cast<void*>(buf), count, MPI_CHAR, dest, tag, getGlobalCommunicator(), &request);
    }
    MPI_Request_free(&request);
}

template <>
void MpiManager::iSendRequestFree<int>
    (int *buf, int count, int dest, int tag)
{
    MPI_Request request;
    if (ok) {
        MPI_Isend(static_cast<void*>(buf), count, MPI_INT, dest, tag, getGlobalCommunicator(), &request);
    }
    MPI_Request_free(&request);
}

template <>
void MpiManager::iSendRequestFree<bool>
    (bool *buf, int count, int dest, int tag)
{
    MPI_Request request;
    if (ok) {
        MPI_Isend(static_cast<void*>(buf), count*sizeof(bool), MPI_CHAR, dest, tag, getGlobalCommunicator(), &request);
    }
    MPI_Request_free(&request);
}


template <>
void MpiManager::iSendRequestFree<long>
    (long *buf, int count, int dest, int tag)
{
    MPI_Request request;
    if (ok) {
        MPI_Isend(static_cast<void*>(buf), count, MPI_LONG, dest, tag, getGlobalCommunicator(), &request);
    }
    MPI_Request_free(&request);
}

template <>
void MpiManager::iSendRequestFree<float>
    (float *buf, int count, int dest, int tag)
{
    MPI_Request request;
    if (ok) {
        MPI_Isend(static_cast<void*>(buf), count, MPI_FLOAT, dest, tag, getGlobalCommunicator(), &request);
    }
    MPI_Request_free(&request);
}

template <>
void MpiManager::iSendRequestFree<double>
    (double *buf, int count, int dest, int tag)
{
    MPI_Request request;
    if (ok) {
        MPI_Isend(static_cast<void*>(buf), count, MPI_DOUBLE, dest, tag, getGlobalCommunicator(), &request);
    }
    MPI_Request_free(&request);
}

template <>
void MpiManager::iSendRequestFree<long double>
    (long double *buf, int count, int dest, int tag)
{
    MPI_Request request;
    if (ok) {
        MPI_Isend(static_cast<void*>(buf), count, MPI_LONG_DOUBLE, dest, tag, getGlobalCommunicator(), &request);
    }
    MPI_Request_free(&request);
}


template <>
void MpiManager::receive<char>(char *buf, int count, int source, int tag)
{
    if (!ok) return;
    MPI_Status status;
    MPI_Recv(static_cast<void*>(buf), count, MPI_CHAR, source, tag, getGlobalCommunicator(), &status);
}

template <>
void MpiManager::receive<int>(int *buf, int count, int source, int tag)
{
    if (!ok) return;
    MPI_Status status;
    MPI_Recv(static_cast<void*>(buf), count, MPI_INT, source, tag, getGlobalCommunicator(), &status);
}

template <>
void MpiManager::receive<bool>(bool *buf, int count, int source, int tag)
{
    if (!ok) return;
    MPI_Status status;
    MPI_Recv(static_cast<void*>(buf), count*sizeof(bool), MPI_CHAR, source, tag, getGlobalCommunicator(), &status);
}


template <>
void MpiManager::receive<long>(long *buf, int count, int source, int tag)
{
    if (!ok) return;
    MPI_Status status;
    MPI_Recv(static_cast<void*>(buf), count, MPI_LONG, source, tag, getGlobalCommunicator(), &status);
}

template <>
void MpiManager::receive<float>(float *buf, int count, int source, int tag)
{
    if (!ok) return;
    MPI_Status status;
    MPI_Recv(static_cast<void*>(buf), count, MPI_FLOAT, source, tag, getGlobalCommunicator(), &status);
}

template <>
void MpiManager::receive<double>(double *buf, int count, int source, int tag)
{
    if (!ok) return;
    MPI_Status status;
    MPI_Recv(static_cast<void*>(buf), count, MPI_DOUBLE, source, tag, getGlobalCommunicator(), &status);
}

template <>
void MpiManager::receive<long double>(long double *buf, int count, int source, int tag)
{
    if (!ok) return;
    MPI_Status status;
    MPI_Recv(static_cast<void*>(buf), count, MPI_LONG_DOUBLE, source, tag, getGlobalCommunicator(), &status);
}


template <>
void MpiManager::sendToMaster<char>(char* sendBuf, int sendCount, bool iAmRoot)
{
    if (!ok) return;
    if (iAmRoot && !isMainProcessor()) {
        send(sendBuf, sendCount, 0);
    }
    if (isMainProcessor() && !iAmRoot) {
        receive(sendBuf, sendCount, MPI_ANY_SOURCE);
    }
}

template <>
void MpiManager::sendToMaster<int>(int* sendBuf, int sendCount, bool iAmRoot)
{
    if (!ok) return;
    if (iAmRoot && !isMainProcessor()) {
        send(sendBuf, sendCount, 0);
    }
    if (isMainProcessor() && !iAmRoot) {
        receive(sendBuf, sendCount, MPI_ANY_SOURCE);
    }
}

template <>
void MpiManager::sendToMaster<bool>(bool* sendBuf, int sendCount, bool iAmRoot)
{
    if (!ok) return;
    if (iAmRoot && !isMainProcessor()) {
        send(sendBuf, sendCount, 0);
    }
    if (isMainProcessor() && !iAmRoot) {
        receive(sendBuf, sendCount, MPI_ANY_SOURCE);
    }
}


template <>
void MpiManager::sendToMaster<long>(long* sendBuf, int sendCount, bool iAmRoot)
{
    if (!ok) return;
    if (iAmRoot && !isMainProcessor()) {
        send(sendBuf, sendCount, 0);
    }
    if (isMainProcessor() && !iAmRoot) {
        receive(sendBuf, sendCount, MPI_ANY_SOURCE);
    }
}

template <>
void MpiManager::sendToMaster<float>(float* sendBuf, int sendCount, bool iAmRoot)
{
    if (!ok) return;
    if (iAmRoot && !isMainProcessor()) {
        send(sendBuf, sendCount, 0);
    }
    if (isMainProcessor() && !iAmRoot) {
        receive(sendBuf, sendCount, MPI_ANY_SOURCE);
    }
}

template <>
void MpiManager::sendToMaster<double>(double* sendBuf, int sendCount, bool iAmRoot)
{
    if (!ok) return;
    if (iAmRoot && !isMainProcessor()) {
        send(sendBuf, sendCount, 0);
    }
    if (isMainProcessor() && !iAmRoot) {
        receive(sendBuf, sendCount, MPI_ANY_SOURCE);
    }
}

template <>
void MpiManager::sendToMaster<long double>(long double* sendBuf, int sendCount, bool iAmRoot)
{
    if (!ok) return;
    if (iAmRoot && !isMainProcessor()) {
        send(sendBuf, sendCount, 0);
    }
    if (isMainProcessor() && !iAmRoot) {
        receive(sendBuf, sendCount, MPI_ANY_SOURCE);
    }
}



template <>
void MpiManager::iRecv<char>(char *buf, int count, int source, MPI_Request* request, int tag)
{
    if (ok) {
      MPI_Irecv(static_cast<void*>(buf), count, MPI_CHAR, source, tag, getGlobalCommunicator(), request);
    }
}

template <>
void MpiManager::iRecv<int>(int *buf, int count, int source, MPI_Request* request, int tag)
{
    if (ok) {
      MPI_Irecv(static_cast<void*>(buf), count, MPI_INT, source, tag, getGlobalCommunicator(), request);
    }
}

template <>
void MpiManager::iRecv<bool>(bool *buf, int count, int source, MPI_Request* request, int tag)
{
    if (ok) {
      MPI_Irecv(static_cast<void*>(buf), count*sizeof(bool), MPI_CHAR, source, tag, getGlobalCommunicator(), request);
    }
}


template <>
void MpiManager::iRecv<long>(long *buf, int count, int source, MPI_Request* request, int tag)
{
    if (ok) {
      MPI_Irecv(static_cast<void*>(buf), count, MPI_LONG, source, tag, getGlobalCommunicator(), request);
    }
}

template <>
void MpiManager::iRecv<float>(float *buf, int count, int source, MPI_Request* request, int tag)
{
    if (ok) {
      MPI_Irecv(static_cast<void*>(buf), count, MPI_FLOAT, source, tag, getGlobalCommunicator(), request);
    }
}

template <>
void MpiManager::iRecv<double>(double *buf, int count, int source, MPI_Request* request, int tag)
{
    if (ok) {
      MPI_Irecv(static_cast<void*>(buf), count, MPI_DOUBLE, source, tag, getGlobalCommunicator(), request);
    }
}

template <>
void MpiManager::iRecv<long double>(long double *buf, int count, int source, MPI_Request* request, int tag)
{
    if (ok) {
      MPI_Irecv(static_cast<void*>(buf), count, MPI_DOUBLE, source, tag, getGlobalCommunicator(), request);
    }
}



void MpiManager::sendRecv(char *sendBuf, char *recvBuf, int count, int dest, int source, int tag)
{
    if (!ok) return;
    MPI_Status status;
    MPI_Sendrecv(static_cast<void*>(sendBuf),
                 count,
                 MPI_CHAR, dest, tag,
                 static_cast<void*>(recvBuf),
                 count,
                 MPI_CHAR, source, tag, getGlobalCommunicator(), &status);
}


void MpiManager::sendRecv(int *sendBuf, int *recvBuf, int count, int dest, int source, int tag)
{
    if (!ok) return;
    MPI_Status status;
    MPI_Sendrecv(static_cast<void*>(sendBuf),
                 count,
                 MPI_INT, dest, tag,
                 static_cast<void*>(recvBuf),
                 count,
                 MPI_INT, source, tag, getGlobalCommunicator(), &status);
}


void MpiManager::sendRecv(bool *sendBuf, bool *recvBuf, int count, int dest, int source, int tag)
{
    if (!ok) return;
    MPI_Status status;
    MPI_Sendrecv(static_cast<void*>(sendBuf),
                 count*sizeof(bool),
                 MPI_CHAR, dest, tag,
                 static_cast<void*>(recvBuf),
                 count,
                 MPI_LONG, source, tag, getGlobalCommunicator(), &status);
}


void MpiManager::sendRecv(long *sendBuf, long *recvBuf, int count, int dest, int source, int tag)
{
    if (!ok) return;
    MPI_Status status;
    MPI_Sendrecv(static_cast<void*>(sendBuf),
                 count,
                 MPI_LONG, dest, tag,
                 static_cast<void*>(recvBuf),
                 count,
                 MPI_LONG, source, tag, getGlobalCommunicator(), &status);
}


void MpiManager::sendRecv(float *sendBuf, float *recvBuf, int count, int dest, int source, int tag)
{
    if (!ok) return;
    MPI_Status status;
    MPI_Sendrecv(static_cast<void*>(sendBuf),
                 count,
                 MPI_FLOAT, dest, tag,
                 static_cast<void*>(recvBuf),
                 count,
                 MPI_FLOAT, source, tag, getGlobalCommunicator(), &status);
}


void MpiManager::sendRecv(double *sendBuf, double *recvBuf, int count_send, int dest, int count_recv, int source, int tag)
{
    if (!ok) return;
    MPI_Status status;
    MPI_Sendrecv(static_cast<void*>(sendBuf),
                 count_send,
                 MPI_DOUBLE, source, tag,
                 static_cast<void*>(recvBuf),
                 count_recv,
                 MPI_DOUBLE, dest, tag, getGlobalCommunicator(), &status);

}


void MpiManager::sendRecv(long double *sendBuf, long double *recvBuf, int count, int dest, int source, int tag)
{
    if (!ok) return;
    MPI_Status status;
    MPI_Sendrecv(static_cast<void*>(sendBuf),
                 count,
                 MPI_LONG_DOUBLE, dest, tag,
                 static_cast<void*>(recvBuf),
                 count,
                 MPI_LONG_DOUBLE, source, tag, getGlobalCommunicator(), &status);

}

void MpiManager::wait(MPI_Request* request, MPI_Status* status)
{
    if (!ok) return;
    MPI_Wait(request, status);
}


