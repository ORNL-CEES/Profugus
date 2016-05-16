//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Utils/comm/MPI.t.hh
 * \author Thomas M. Evans
 * \date   Wed Jan  2 13:42:57 2008
 * \brief  MPI template function declarations (defined in Functions.hh).
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef Utils_comm_MPI_t_hh
#define Utils_comm_MPI_t_hh

#include <Utils/config.h>
#ifdef COMM_MPI

#include "MPI.hh"

#include <vector>
#include <mpi.h>

#include "harness/DBC.hh"
#include "MPI_Traits.hh"
#include "Request.hh"
#include "Definitions.hh"

namespace profugus
{
//---------------------------------------------------------------------------//
// BLOCKING SEND/RECEIVE OPERATIONS
//---------------------------------------------------------------------------//

template<class T>
int send(const T *buffer,
         int      size,
         int      destination,
         int      tag)
{
    MPI_Send(const_cast<T *>(buffer), size, MPI_Traits<T>::element_type(),
             destination, tag, communicator);
    return COMM_SUCCESS;
}

//---------------------------------------------------------------------------//

template<class T>
int send_comm(
        const T*              buffer,
        const Communicator_t& comm,
        int                   size,
        int                   destination,
        int                   tag)
{
    MPI_Send(const_cast<T *>(buffer), size, MPI_Traits<T>::element_type(),
             destination, tag, comm);
    return COMM_SUCCESS;
}

//---------------------------------------------------------------------------//

template<class T>
int receive(
        T*  buffer,
        int size,
        int source,
        int tag)
{
    int count = 0;

    // get a handle to the MPI_Status
    MPI_Status status;

    // do the blocking receive
    MPI_Recv(buffer, size, MPI_Traits<T>::element_type(), source, tag,
             communicator, &status);

    // get the count of received data
    MPI_Get_count(&status, MPI_Traits<T>::element_type(), &count);
    return count;
}

//---------------------------------------------------------------------------//

template<class T>
int receive_comm(
        T*                    buffer,
        const Communicator_t& comm,
        int                   size,
        int                   source,
        int                   tag)
{
    int count = 0;

    // get a handle to the MPI_Status
    MPI_Status status;

    // do the blocking receive
    MPI_Recv(buffer, size, MPI_Traits<T>::element_type(), source, tag,
             comm, &status);

    // get the count of received data
    MPI_Get_count(&status, MPI_Traits<T>::element_type(), &count);
    return count;
}

//---------------------------------------------------------------------------//
// NON-BLOCKING SEND/RECEIVE OPERATIONS
//---------------------------------------------------------------------------//

template<class T>
Request send_async(
        const T * buffer,
        int       size,
        int       destination,
        int       tag)
{
    // make a comm request handle
    Request request;

    // do an MPI_Isend (non-blocking send)
    MPI_Isend(const_cast<T *>(buffer), size, MPI_Traits<T>::element_type(),
              destination, tag, communicator, &request.r());

    // set the request to active
    request.set();
    return request;
}

//---------------------------------------------------------------------------//

template<class T>
void send_async(Request &request,
                const T *buffer,
                int      size,
                int      destination,
                int      tag)
{
    REQUIRE(!request.inuse());

    // set the request
    request.set();

    // post an MPI_Isend
    MPI_Isend(const_cast<T *>(buffer), size, MPI_Traits<T>::element_type(),
              destination, tag, communicator, &request.r());
}

//---------------------------------------------------------------------------//

template<class T>
void send_async_comm(
        Request &             request,
        const T *             buffer,
        const Communicator_t& comm,
        int                   size,
        int                   destination,
        int                   tag)
{
    REQUIRE(!request.inuse());

    // set the request
    request.set();

    // post an MPI_Isend
    MPI_Isend(const_cast<T *>(buffer), size, MPI_Traits<T>::element_type(),
              destination, tag, comm, &request.r());
}

//---------------------------------------------------------------------------//

template<class T>
Request receive_async(T* buffer,
                      int  size,
                      int  source,
                      int  tag)
{
    // make a comm request handle
    Request request;

    // post an MPI_Irecv (non-blocking receive)
    MPI_Irecv(buffer, size, MPI_Traits<T>::element_type(), source, tag,
              communicator, &request.r());

    // set the request to active
    request.set();
    return request;
}

//---------------------------------------------------------------------------//

template<class T>
void receive_async(Request &request,
                   T       *buffer,
                   int      size,
                   int      source,
                   int      tag)
{
    REQUIRE(!request.inuse());

    // set the request
    request.set();

    // post an MPI_Irecv
    MPI_Irecv(buffer, size, MPI_Traits<T>::element_type(), source, tag,
              communicator, &request.r());
}

//---------------------------------------------------------------------------//

template<class T>
void receive_async_comm(
        Request &             request,
        T       *             buffer,
        const Communicator_t& comm,
        int                   size,
        int                   source,
        int                   tag)
{
    REQUIRE(!request.inuse());

    // set the request
    request.set();

    // post an MPI_Irecv
    MPI_Irecv(buffer, size, MPI_Traits<T>::element_type(), source, tag,
              comm, &request.r());
}

//---------------------------------------------------------------------------//
// BROADCAST
//---------------------------------------------------------------------------//

template<class T>
int broadcast(T  *buffer,
              int size,
              int root)
{
    int r = MPI_Bcast(buffer, size, MPI_Traits<T>::element_type(),
                      root, communicator);
    return r;
}

//---------------------------------------------------------------------------//

template<class T>
int broadcast(
        T*                    buffer,
        const Communicator_t& comm,
        int                   size,
        int                   root)
{
    int r = MPI_Bcast(buffer, size, MPI_Traits<T>::element_type(),
                      root, comm);
    return r;
}

//---------------------------------------------------------------------------//
// GLOBAL REDUCTIONS
//---------------------------------------------------------------------------//

template<class T>
void global_sum(T &x)
{
    // copy data into send buffer
    T y = x;

    // do global MPI reduction (result is on all processors) into x
    MPI_Allreduce(&y, &x, 1, MPI_Traits<T>::element_type(), MPI_SUM,
                  communicator);
}

//---------------------------------------------------------------------------//

template<class T>
void global_prod(T &x)
{
     // copy data into send buffer
    T y = x;

    // do global MPI reduction (result is on all processors) into x
    MPI_Allreduce(&y, &x, 1, MPI_Traits<T>::element_type(), MPI_PROD,
                  communicator);
}

//---------------------------------------------------------------------------//

template<class T>
void global_min(T &x)
{
     // copy data into send buffer
    T y = x;

    // do global MPI reduction (result is on all processors) into x
    MPI_Allreduce(&y, &x, 1, MPI_Traits<T>::element_type(), MPI_MIN,
                  communicator);
}

//---------------------------------------------------------------------------//

template<class T>
void global_min(T &x, const Communicator_t& comm)
{
     // copy data into send buffer
    T y = x;

    // do global MPI reduction (result is on all processors) into x
    MPI_Allreduce(&y, &x, 1, MPI_Traits<T>::element_type(), MPI_MIN,
                  comm);
}

//---------------------------------------------------------------------------//

template<class T>
void global_max(T &x)
{
     // copy data into send buffer
    T y = x;

    // do global MPI reduction (result is on all processors) into x
    MPI_Allreduce(&y, &x, 1, MPI_Traits<T>::element_type(), MPI_MAX,
                  communicator);
}

//---------------------------------------------------------------------------//

template<class T>
void global_max(T &x, const Communicator_t& comm)
{
     // copy data into send buffer
    T y = x;

    // do global MPI reduction (result is on all processors) into x
    MPI_Allreduce(&y, &x, 1, MPI_Traits<T>::element_type(), MPI_MAX,
                  comm);
}

//---------------------------------------------------------------------------//

template<class T>
void global_sum(T  *x,
                int n)
{
    REQUIRE(x);

    // copy data into a send buffer
    std::vector<T> send_buffer(x, x + n);

    // do a element-wise global reduction (result is on all processors) into
    // x
    MPI_Allreduce(&send_buffer[0], x, n, MPI_Traits<T>::element_type(),
                  MPI_SUM, communicator);
}

//---------------------------------------------------------------------------//

template<class T>
void global_prod(T  *x,
                 int n)
{
    REQUIRE(x);

    // copy data into a send buffer
    std::vector<T> send_buffer(x, x + n);

    // do a element-wise global reduction (result is on all processors) into
    // x
    MPI_Allreduce(&send_buffer[0], x, n, MPI_Traits<T>::element_type(),
                  MPI_PROD, communicator);
}

//---------------------------------------------------------------------------//

template<class T>
void global_min(T  *x,
                int n)
{
    REQUIRE(x);

    // copy data into a send buffer
    std::vector<T> send_buffer(x, x + n);

    // do a element-wise global reduction (result is on all processors) into
    // x
    MPI_Allreduce(&send_buffer[0], x, n, MPI_Traits<T>::element_type(),
                  MPI_MIN, communicator);
}

//---------------------------------------------------------------------------//

template<class T>
void global_max(T  *x,
                int n)
{
    REQUIRE(x);

    // copy data into a send buffer
    std::vector<T> send_buffer(x, x + n);

    // do a element-wise global reduction (result is on all processors) into
    // x
    MPI_Allreduce(&send_buffer[0], x, n, MPI_Traits<T>::element_type(),
                  MPI_MAX, communicator);
}

//---------------------------------------------------------------------------//
// REDUCTIONS
//---------------------------------------------------------------------------//

template<class T>
void sum(T  &x,
         int to_node)
{
    // copy data into send buffer
    T y = x;

    // do MPI reduction (result is on to_node) into x
    MPI_Reduce(&y, &x, 1, MPI_Traits<T>::element_type(), MPI_SUM, to_node,
               communicator);
}

//---------------------------------------------------------------------------//

template<class T>
void prod(T  &x,
          int to_node)
{
    // copy data into send buffer
    T y = x;

    // do MPI reduction (result is on to_node) into x
    MPI_Reduce(&y, &x, 1, MPI_Traits<T>::element_type(), MPI_PROD, to_node,
               communicator);
}

//---------------------------------------------------------------------------//

template<class T>
void min(T  &x,
         int to_node)
{
    // copy data into send buffer
    T y = x;

    // do MPI reduction (result is on to_node) into x
    MPI_Reduce(&y, &x, 1, MPI_Traits<T>::element_type(), MPI_MIN, to_node,
               communicator);
}

//---------------------------------------------------------------------------//

template<class T>
void max(T  &x,
         int to_node)
{
    // copy data into send buffer
    T y = x;

    // do MPI reduction (result is on to_node) into x
    MPI_Reduce(&y, &x, 1, MPI_Traits<T>::element_type(), MPI_MAX, to_node,
               communicator);
}

//---------------------------------------------------------------------------//

template<class T>
void sum(T  *x,
         int n,
         int to_node)
{
    REQUIRE(x);

    // copy data into send buffer
    std::vector<T> send_buffer(x, x + n);

    // do MPI reduction (result is on to_node) into x
    MPI_Reduce(&send_buffer[0], x, n, MPI_Traits<T>::element_type(), MPI_SUM,
               to_node, communicator);
}

//---------------------------------------------------------------------------//

template<class T>
void prod(T  *x,
          int n,
          int to_node)
{
    REQUIRE(x);

    // copy data into send buffer
    std::vector<T> send_buffer(x, x + n);

    // do MPI reduction (result is on to_node) into x
    MPI_Reduce(&send_buffer[0], x, n, MPI_Traits<T>::element_type(), MPI_PROD,
               to_node, communicator);
}

//---------------------------------------------------------------------------//

template<class T>
void min(T  *x,
         int n,
         int to_node)
{
    REQUIRE(x);

    // copy data into send buffer
    std::vector<T> send_buffer(x, x + n);

    // do MPI reduction (result is on to_node) into x
    MPI_Reduce(&send_buffer[0], x, n, MPI_Traits<T>::element_type(), MPI_MIN,
               to_node, communicator);
}

//---------------------------------------------------------------------------//

template<class T>
void max(T  *x,
         int n,
         int to_node)
{
    REQUIRE(x);

    // copy data into send buffer
    std::vector<T> send_buffer(x, x + n);

    // do MPI reduction (result is on to_node) into x
    MPI_Reduce(&send_buffer[0], x, n, MPI_Traits<T>::element_type(), MPI_MAX,
               to_node, communicator);
}

//---------------------------------------------------------------------------//
// REDUCE-SCATTERS
//---------------------------------------------------------------------------//

template<class T>
void sum_scatter(const T   *x,
                 T         *y,
                 const int *steer)
{
    REQUIRE(x);
    REQUIRE(y);
    REQUIRE(steer);

    // do MPI Reduce-Scatter into y defined by steer
    MPI_Reduce_scatter(const_cast<T *>(x), y, const_cast<int *>(steer),
                       MPI_Traits<T>::element_type(), MPI_SUM, communicator);
}

//---------------------------------------------------------------------------//

template<class T>
void prod_scatter(const T   *x,
                  T         *y,
                  const int *steer)
{
    REQUIRE(x);
    REQUIRE(y);
    REQUIRE(steer);

    // do MPI Reduce-Scatter into y defined by steer
    MPI_Reduce_scatter(const_cast<T *>(x), y, const_cast<int *>(steer),
                       MPI_Traits<T>::element_type(), MPI_PROD, communicator);

}

//---------------------------------------------------------------------------//
// GATHERS
//---------------------------------------------------------------------------//

template<class T>
void all_gather(const T *sendbuf,
                T       *recvbuf,
                int      num_els)
{
    REQUIRE(sendbuf);
    REQUIRE(recvbuf);
    REQUIRE(num_els > 0);

    MPI_Allgather(
            const_cast<T*>(sendbuf), num_els, MPI_Traits<T>::element_type(),
            recvbuf, num_els, MPI_Traits<T>::element_type(),
            communicator);
}

//---------------------------------------------------------------------------//

template<class T>
void gather(const T *sendbuf,
            T       *recvbuf,
            int      num_els,
            int      root)
{
    REQUIRE(sendbuf);
    REQUIRE(root ? recvbuf != 0 : true);
    REQUIRE(num_els > 0);

#ifdef REQUIRE_ON
    int nodes = 0;
    MPI_Comm_size(communicator, &nodes);
    VALIDATE(root < nodes,
              "root node " << root << " > number of requested nodes");
#endif

    MPI_Gather(
        const_cast<T*>(sendbuf), num_els, MPI_Traits<T>::element_type(),
        recvbuf, num_els, MPI_Traits<T>::element_type(), root,
        communicator);
}

//---------------------------------------------------------------------------//
// ALL-TO-ALLS
//---------------------------------------------------------------------------//

template<class T>
void all_to_all(const T *sendbuf,
                T       *recvbuf,
                int      n)
{
    REQUIRE(sendbuf);
    REQUIRE(recvbuf);

    // Do all-to-all communication
    MPI_Alltoall(const_cast<T *>(sendbuf), n, MPI_Traits<T>::element_type(),
                 recvbuf, n, MPI_Traits<T>::element_type(), communicator);
}

//---------------------------------------------------------------------------//

template<class T>
void all_to_all(const T   *sendbuf,     const int *sendcounts,
                const int *sendoffsets, T         *recvbuf,
                const int *recvcounts,  const int *recvoffsets)
{
    REQUIRE(sendbuf);
    REQUIRE(sendcounts);
    REQUIRE(sendoffsets);
    REQUIRE(recvbuf);
    REQUIRE(recvcounts);
    REQUIRE(recvoffsets);

    // Do all-to-all communication
    MPI_Alltoallv(const_cast<T *>(sendbuf), const_cast<int *>(sendcounts),
                  const_cast<int *>(sendoffsets),
                  MPI_Traits<T>::element_type(),
                  recvbuf, const_cast<int *>(recvcounts),
                  const_cast<int *>(recvoffsets),
                  MPI_Traits<T>::element_type(), communicator);
}

} // end namespace profugus

#endif // COMM_MPI

#endif // Utils_comm_MPI_t_hh

//---------------------------------------------------------------------------//
//              end of comm/MPI.t.hh
//---------------------------------------------------------------------------//
