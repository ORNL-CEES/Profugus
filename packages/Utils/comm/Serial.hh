//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   comm/Serial.hh
 * \author Thomas M. Evans
 * \date   Wed Jan  2 14:32:30 2008
 * \brief
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef comm_Serial_hh
#define comm_Serial_hh

#include <map>
#include <cstring>
#include <algorithm>

#include <comm/config.h>
#ifdef COMM_SCALAR

#include "harness/DBC.hh"

#include "Functions.hh"
#include "Request.hh"
#include "Definitions.hh"

namespace profugus
{

//---------------------------------------------------------------------------//
// NULL COMMUNICATOR TYPE
//---------------------------------------------------------------------------//

extern Communicator_t default_communicator;
extern Communicator_t communicator;

//---------------------------------------------------------------------------//
// BLOCKING SEND/RECEIVE OPERATIONS
//---------------------------------------------------------------------------//

template<class T>
int send(
        const T * buffer,
        int       size,
        int       destination,
        int       tag)
{
    // blocking sends should never be used for send-to-self; althought the
    // spec doesn't specifically disallow it but says it is unsafe and
    // non-portable, so we disallow it here:
    throw profugus::assertion("Blocking self-communication not allowed.");
}

//---------------------------------------------------------------------------//

template<class T>
int send_comm(
        const T* const        buffer,
        const Communicator_t& comm,
        int                   size,
        int                   destination,
        int                   tag)
{
    // blocking sends should never be used for send-to-self; althought the
    // spec doesn't specifically disallow it but says it is unsafe and
    // non-portable, so we disallow it here:
    throw profugus::assertion("Blocking self-communication not allowed.");
}

//---------------------------------------------------------------------------//

template<class T>
int receive(
        T*  buffer,
        int size,
        int source,
        int tag)
{
    // blocking receives should never be used for receive-from-self; althought
    // the spec doesn't specifically disallow it but says it is unsafe and
    // non-portable, so we disallow it here:
    throw profugus::assertion("Blocking self-communication not allowed.");
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
    // blocking receives should never be used for receive-from-self; althought
    // the spec doesn't specifically disallow it but says it is unsafe and
    // non-portable, so we disallow it here:
    throw profugus::assertion("Blocking self-communication not allowed.");
}

//---------------------------------------------------------------------------//
// NON-BLOCKING SEND/RECEIVE OPERATIONS
//---------------------------------------------------------------------------//
// Receive buffers.
namespace internals
{
extern std::map<int, void*> buffers;
}

template<class T>
Request send_async(
        const T * buffer,
        int       size,
        int       destination,
        int       tag)
{
    REQUIRE(internals::buffers.count(tag));

    // make a comm request handle (no need to set it because we will do an
    // immediate copy into the receive buffer)
    Request request;

    // set the request
    request.set();

    // receive buffer
    T *recv_buffer = reinterpret_cast<T*>(internals::buffers[tag]);

    // memcopy into the tagged buffer
    std::memcpy(recv_buffer, buffer, sizeof(T) * size);

    return request;
}

//---------------------------------------------------------------------------//

template<class T>
void send_async(
        Request  & request,
        const T  * buffer,
        int        size,
        int        destination,
        int        tag)
{
    REQUIRE(internals::buffers.count(tag));

    // set the request
    request.set();

    // receive buffer
    T *recv_buffer = reinterpret_cast<T*>(internals::buffers[tag]);

    // memcopy into the tagged buffer
    std::memcpy(recv_buffer, buffer, sizeof(T) * size);
}

//---------------------------------------------------------------------------//

template<class T>
void send_async_comm(
        Request  &            request,
        const T  *            buffer,
        const Communicator_t& comm,
        int                   size,
        int                   destination,
        int                   tag)
{
    REQUIRE(internals::buffers.count(tag));

    // set the request
    request.set();

    // receive buffer
    T *recv_buffer = reinterpret_cast<T*>(internals::buffers[tag]);

    // memcopy into the tagged buffer
    std::memcpy(recv_buffer, buffer, sizeof(T) * size);
}

//---------------------------------------------------------------------------//

template<class T>
Request receive_async(
        T*  buffer,
        int size,
        int source,
        int tag)
{
    REQUIRE(source == 0);

    // NON-BLOCKING SEND-TO-SELF ALLOWED IN SPEC

    // make a comm request handle
    Request request;

    // set it
    request.set();

    // add buffer to internal storage
    internals::buffers[tag] = reinterpret_cast<void *>(buffer);

    // return the request
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
    REQUIRE(source == 0);
    REQUIRE(!request.inuse());

    // set it
    request.set();

    // add buffer to internal storage
    internals::buffers[tag] = reinterpret_cast<void *>(buffer);
}

//---------------------------------------------------------------------------//

template<class T>
void receive_async_comm(Request &request,
                        T       *buffer,
                        const Communicator_t& comm,
                        int      size,
                        int      source,
                        int      tag)
{
    REQUIRE(source == 0);
    REQUIRE(!request.inuse());

    // set it
    request.set();

    // add buffer to internal storage
    internals::buffers[tag] = reinterpret_cast<void *>(buffer);
}

//---------------------------------------------------------------------------//
// BROADCAST
//---------------------------------------------------------------------------//

template<class T>
int broadcast(
        T  * buffer,
        int  size,
        int  root)
{
    return COMM_SUCCESS;
}

//---------------------------------------------------------------------------//

template<class T>
int broadcast(
        T*                    buffer,
        const Communicator_t& comm,
        int                   size,
        int                   root)
{
    return COMM_SUCCESS;
}

//---------------------------------------------------------------------------//

template<class ForwardIterator, class OutputIterator>
void broadcast(
        ForwardIterator first,
        ForwardIterator last,
        OutputIterator  result)
{
    // No communication needed for Serial use.
    return;
}

//---------------------------------------------------------------------------//
// GLOBAL REDUCTIONS
//---------------------------------------------------------------------------//

template<class T>
void global_sum(T &x)
{
}

//---------------------------------------------------------------------------//

template<class T>
void global_prod(T &x)
{
}

//---------------------------------------------------------------------------//

template<class T>
void global_min(T &x)
{
}

//---------------------------------------------------------------------------//

template<class T>
void global_min(T &x, const Communicator_t& comm)
{
}

//---------------------------------------------------------------------------//

template<class T>
void global_max(T &x)
{
}

//---------------------------------------------------------------------------//

template<class T>
void global_max(T &x, const Communicator_t& comm)
{
}

//---------------------------------------------------------------------------//

template<class T>
void global_sum(T *x, int n)
{
}

//---------------------------------------------------------------------------//

template<class T>
void global_prod(T *x, int n)
{
}

//---------------------------------------------------------------------------//

template<class T>
void global_min(T *x, int n)
{
}

//---------------------------------------------------------------------------//

template<class T>
void global_max(T *x, int n)
{
}

//---------------------------------------------------------------------------//
// REDUCTIONS
//---------------------------------------------------------------------------//

template<class T>
void sum(T &x, int to_node)
{
}

//---------------------------------------------------------------------------//

template<class T>
void prod(T &x, int to_node)
{
}

//---------------------------------------------------------------------------//

template<class T>
void min(T &x, int to_node)
{
}

//---------------------------------------------------------------------------//

template<class T>
void max(T &x, int to_node)
{
}

//---------------------------------------------------------------------------//

template<class T>
void sum(T *x, int n, int to_node)
{
}

//---------------------------------------------------------------------------//

template<class T>
void prod(T *x, int n, int to_node)
{
}

//---------------------------------------------------------------------------//

template<class T>
void min(T *x, int n, int to_node)
{
}

//---------------------------------------------------------------------------//

template<class T>
void max(T *x, int n, int to_node)
{
}

//---------------------------------------------------------------------------//

template<class T>
void sum_scatter(const T   *x,
                 T         *y,
                 const int *steer)
{
    // Reduction is a no-op.  Scatter is simply a copy.
    std::copy(x, x+steer[0], y);
}

//---------------------------------------------------------------------------//

template<class T>
void prod_scatter(const T   *x,
                  T         *y,
                  const int *steer)
{
    // Reduction is a no-op.  Scatter is simply a copy.
    std::copy(x, x+steer[0], y);
}

//---------------------------------------------------------------------------//

template<class T>
void all_to_all(const T *sendbuf,
                T       *recvbuf,
                int      n)
{
    REQUIRE(sendbuf);
    REQUIRE(recvbuf);

    // all-to-all is copy
    std::copy(sendbuf, sendbuf+n, recvbuf);
}

//---------------------------------------------------------------------------//

template<class T>
void all_to_all(const T   *sendbuf,     const int *sendcounts,
                const int *sendoffsets, T         *recvbuf,
                const int *recvcounts,  const int *recvoffsets)
{
    REQUIRE(sendbuf);
    REQUIRE(recvbuf);
    REQUIRE(sendcounts[0]==recvcounts[0]);

    // all-to-all is copy
    std::copy(sendbuf, sendbuf+sendcounts[0], recvbuf);
}

//---------------------------------------------------------------------------//

template<class T>
void all_gather(const T* sendbuf, T* recvbuf, int num_els)
{
    REQUIRE(sendbuf);
    REQUIRE(recvbuf);
    REQUIRE(num_els > 0);

    // all_gather is copy
    std::copy(sendbuf, sendbuf + num_els, recvbuf);
}

//---------------------------------------------------------------------------//

template<class T>
void gather(const T* sendbuf, T* recvbuf, int num_els, int root)
{
    REQUIRE(sendbuf);
    REQUIRE(recvbuf);
    REQUIRE(num_els > 0);
    REQUIRE(root == 0);

    // gather is copy
    std::copy(sendbuf, sendbuf + num_els, recvbuf);
}

//---------------------------------------------------------------------------//
} // end namespace profugus

#endif // COMM_SCALAR

#endif // comm_Serial_hh

//---------------------------------------------------------------------------//
//              end of comm/Serial.hh
//---------------------------------------------------------------------------//
