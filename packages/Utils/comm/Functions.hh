//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   comm/Functions.hh
 * \author Thomas M. Evans
 * \date   Wed Jan  2 13:25:15 2008
 * \brief  Communication functions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 *
 * This file contains the declarations for communication functions provided
 * by comm.
 */
//---------------------------------------------------------------------------//

#ifndef comm_Functions_hh
#define comm_Functions_hh

#include <string>

#include <Utils/config.h>

#ifdef HAVE_SYS_TIMES_H
#include <sys/times.h> // defines the struct tms.
#endif

#include "Comm_Traits.hh"
#include "Definitions.hh"

namespace profugus
{

//---------------------------------------------------------------------------//
/*!
 * Comm unit tests.
 */
/*! \example comm/test/tstAbort.cc
 * Example of MPI abort functions.
 */
/*! \example comm/test/tstBroadcast.cc
 * Example of MPI broadcast-like functions
 */
/*! \example comm/test/tstComm_Dup.cc
 * Example of communicator duplication.
 */
/*! \example comm/test/tstPingPong.cc
 * Example of point-to-point communications
 */
/*! \example comm/test/tstReduction.cc
 * Example of many-to-one communications
 */
/*! \example comm/test/tstSplit.cc
 * Example of communicator splitting.
 */
/*! \example comm/test/tstSerial.cc
 * Example of serial build.
 */
//---------------------------------------------------------------------------//

// Forward declarations.
class Request;

//---------------------------------------------------------------------------//
// SETUP FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Initialize a parallel job.
 */
void initialize(int &argc, char **&argv);

//---------------------------------------------------------------------------//
/*!
 * \brief Finish a parallel job.
 */
void finalize();

//---------------------------------------------------------------------------//
/*!
 * \brief Changes the default internal communicator, and switches to the
 *        new default.
 */
void set_default(const Communicator_t&);

//---------------------------------------------------------------------------//
/*!
 * \brief Inherit a communicator from another application.
 */
void inherit(const Communicator_t&);

//---------------------------------------------------------------------------//
/*!
 * \brief Free an inherited communicator from another application.
 */
void free_inherited_comm();

//---------------------------------------------------------------------------//
/*!
 * \brief Split a communicator into groups specified by color and key.
 *
 * The new communicator will contain all processes with the same color.  Key
 * can be used to determine the rank for the color group.  Giving the same
 * rank to multiple processes in a color will preserve the same relative
 * rank-ordering of the processes in the parent-communicator.
 *
 * The communicator should be destroyed free_comm().
 */
void split(int color, int key, Communicator_t& new_comm);

//---------------------------------------------------------------------------//
/*!
 * \brief Split a communicator into groups specified by color and key
 *
 */
void split(
        int                   color,
        int                   key,
        Communicator_t&       new_comm,
        const Communicator_t& comm);

//---------------------------------------------------------------------------//
/*!
 * \brief Duplicate arbitrary communicator.
 */
void duplicate_comm(Communicator_t& new_comm, const Communicator_t& comm);

//---------------------------------------------------------------------------//
/*!
 * \brief Free a communicator.
 */
void free_comm(Communicator_t&);

//---------------------------------------------------------------------------//
/*!
 * \brief Set internal communicator.
 *
 * Allows client to change the underlying communicator used by the comm
 * package.
 */
void set_internal_comm(const Communicator_t& comm);

//---------------------------------------------------------------------------//
/*!
 * \brief Reset the internal communicator to the default communicator.
 *
 * For MPI the default is MPI_COMM_WORLD
 */
void reset_internal_comm();

//---------------------------------------------------------------------------//
// QUERY FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Get the node (rank) of the current processor.
 *
 * The rank is determined by the current communicator.
 */
int node();

//---------------------------------------------------------------------------//
/*!
 * \brief Get the node (rank) of the current processor.
 *
 * The rank is determined by the given communicator.
 */
int node(const Communicator_t& comm);

//---------------------------------------------------------------------------//
/*!
 * \brief Get the number of processors used for this job.
 *
 * The number of nodes is determined by the current communicator.
 */
int nodes();

//---------------------------------------------------------------------------//
/*!
 * \brief Get the number of processors used for this job.
 *
 * The number of nodes is determined by the given communicator.
 */
int nodes(const Communicator_t& comm);

//---------------------------------------------------------------------------//
/*!
 * \brief Get the processor name, which may contain the hardware node number.
 */
std::string proc_name();

//---------------------------------------------------------------------------//
// BARRIER FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Set a global barrier for the communicator.
 */
void global_barrier();

//---------------------------------------------------------------------------//
/*!
 * \brief Set a global barrier for arbitrary communicator.
 */
void barrier(const Communicator_t& comm);

//---------------------------------------------------------------------------//
// BLOCKING SEND/RECEIVE OPERATIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Do a point-to-point, blocking send.
 */
template<class T>
int send(const T *buffer, int size, int destination,
         int tag = Comm_Traits<T*>::tag);

//---------------------------------------------------------------------------//
/*!
 * \brief Do a point-to-point, blocking send for arbitrary communicator.
 */
template<class T>
int send_comm(const T* const buffer, const Communicator_t& comm, int size,
        int destination, int tag = Comm_Traits<T*>::tag);

//---------------------------------------------------------------------------//
/*!
 * \brief Do a point-to-point, blocking receive.
 */
template<class T>
int receive(T *buffer, int size, int source, int tag = Comm_Traits<T*>::tag);

//---------------------------------------------------------------------------//
/*!
 * \brief Do a point-to-point, blocking receive for arbitrary communicator.
 */
template<class T>
int receive_comm(T *buffer, const Communicator_t& comm, int size, int source,
            int tag = Comm_Traits<T*>::tag);

//---------------------------------------------------------------------------//
// NON-BLOCKING SEND/RECEIVE OPERATIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Do a point-to-point, non-blocking send.
 *
 * \return Request object to handle communciation requests
 */
template<class T>
Request send_async(const T *buffer, int size, int destination,
                  int tag = Comm_Traits<T*>::tag);

//---------------------------------------------------------------------------//
/*!
 * \brief Do a point-to-point, non-blocking send.
 */
template<class T>
void send_async(Request &request, const T *buffer, int size, int destination,
                int tag = Comm_Traits<T*>::tag);

//---------------------------------------------------------------------------//
/*!
 * \brief Do a point-to-point, non-blocking send for arbitrary communicator.
 */
template<class T>
void send_async_comm(Request &request, const T *buffer,
        const Communicator_t& comm, int size, int destination,
        int tag = Comm_Traits<T*>::tag);

//---------------------------------------------------------------------------//
/*!
 * \brief Do a point-to-point, non-blocking receive.
 *
 * \return Request object to handle communciation requests
 */
template<class T>
Request receive_async(T *buffer, int size, int source,
                     int tag = Comm_Traits<T*>::tag);

//---------------------------------------------------------------------------//
/*!
 * \brief Do a point-to-point, non-blocking receive.
 */
template<class T>
void receive_async(Request& request, T *buffer, int size, int source,
                   int tag = Comm_Traits<T*>::tag);

//---------------------------------------------------------------------------//
/*!
 * \brief Do a point-to-point, non-blocking receive for arbitrary communicator.
 */
template<class T>
void receive_async_comm(Request& request, T *buffer,
        const Communicator_t& comm, int size,
        int source, int tag = Comm_Traits<T*>::tag);

//---------------------------------------------------------------------------//
// BROADCAST
//---------------------------------------------------------------------------//

template<class T>
int broadcast(T *buffer, int size, int root);

/*---------------------------------------------------------------------------*/

template<class T>
int broadcast(T *buffer, const Communicator_t& comm, int size, int root);

/*---------------------------------------------------------------------------*/
/*!
 * \brief Send data from processor 0 to all other processors.
 */
template<class ForwardIterator, class OutputIterator>
void broadcast(ForwardIterator first, ForwardIterator last,
               OutputIterator result);

//---------------------------------------------------------------------------//
// GLOBAL REDUCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Do a global sum of a scalar variable.
 */
template<class T>
void global_sum(T &x);

//---------------------------------------------------------------------------//
/*!
 * \brief Do a global product of a scalar variable.
 */
template<class T>
void global_prod(T &x);

//---------------------------------------------------------------------------//
/*!
 * \brief Do a global minimum of a scalar variable.
 */
template<class T>
void global_min(T &x);

//---------------------------------------------------------------------------//
/*!
 * \brief Do a global minimum of a scalar variable for arbitrary communicator.
 */
template<class T>
void global_min(T &x, const Communicator_t& comm);

//---------------------------------------------------------------------------//
/*!
 * \brief Do a global maximum of a scalar variable.
 */
template<class T>
void global_max(T &x);

//---------------------------------------------------------------------------//
/*!
 * \brief Do a global maximum of a scalar variable for arbitrary communicator.
 */
template<class T>
void global_max(T &x, const Communicator_t& comm);

//---------------------------------------------------------------------------//
/*!
 * \brief Do an element-wise, global sum of an array.
 */
template<class T>
void global_sum(T *x, int n);

//---------------------------------------------------------------------------//
/*!
 * \brief Do an element-wise, global product of an array.
 */
template<class T>
void global_prod(T *x, int n);

//---------------------------------------------------------------------------//
/*!
 * \brief Do an element-wise, global minimum of an array.
 */
template<class T>
void global_min(T *x, int n);

//---------------------------------------------------------------------------//
/*!
 * \brief Do an element-wise, global maximum of an array.
 */
template<class T>
void global_max(T *x, int n);

//---------------------------------------------------------------------------//
// REDUCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Do a sum of a scalar variable.
 *
 * \param x scalar variable
 * \param to_node node that result is stored in
 */
template<class T>
void sum(T &x, int to_node);

//---------------------------------------------------------------------------//
/*!
 * \brief Do a product of a scalar variable.
 *
 * \param x scalar variable
 * \param to_node node that result is stored in
 */
template<class T>
void prod(T &x, int to_node);

//---------------------------------------------------------------------------//
/*!
 * \brief Find the minimum of a scalar variable.
 *
 * \param x scalar variable
 * \param to_node node that result is stored in
 */
template<class T>
void min(T &x, int to_node);

//---------------------------------------------------------------------------//
/*!
 * \brief Find the maximum of a scalar variable.
 *
 * \param x scalar variable
 * \param to_node node that result is stored in
 */
template<class T>
void max(T &x, int to_node);

//---------------------------------------------------------------------------//
/*!
 * \brief Do element-wise sum of an array.
 *
 * \param x array of length \c n
 * \param to_node node that result is stored in
 */
template<class T>
void sum(T *x, int n, int to_node);

//---------------------------------------------------------------------------//
/*!
 * \brief Do an element-wise  product of an array.
 *
 * \param x array of length \ n
 * \param to_node node that result is stored in
 */
template<class T>
void prod(T *x, int n, int to_node);

//---------------------------------------------------------------------------//
/*!
 * \brief Find the element-wise minimum of an array.
 *
 * \param x array of length \c n
 * \param to_node node that result is stored in
 */
template<class T>
void min(T *x, int n, int to_node);

//---------------------------------------------------------------------------//
/*!
 * \brief Find the element-wise maximum of an array.
 *
 * \param x array of length \c n
 * \param to_node node that result is stored in
 */
template<class T>
void max(T *x, int n, int to_node);

//---------------------------------------------------------------------------//
// REDUCE-SCATTERS
//---------------------------------------------------------------------------//
/*!
 * \brief Sum an array across all processors and scatter result.
 *
 * \param x array of length \f$ n = \sum_{i}^{N_p} recv_counts[i]\f$
 *
 * \param y array of length \c recv_counts[p] where \c p is the current
 *          processor id; result on \c p is written into this array
 *
 * \param steer steering array of length \c Np that gives the number of
 *              elements of \c x to write into \c y on a particular processor
 */
template<class T>
void sum_scatter(const T *x, T *y, const int *steer);

//---------------------------------------------------------------------------//
/*!
 * \brief Multiply an array across all processors and scatter result.
 *
 * \param x array of length \f$ n = \sum_{i}^{N_p} recv_counts[i]\f$
 *
 * \param y array of length \c recv_counts[p] where \c p is the current
 *          processor id; result on \c p is written into this array
 *
 * \param steer steering array of length \c Np that gives the number of
 *              elements of \c x to write into \c y on a particular processor
 */
template<class T>
void prod_scatter(const T *x, T *y, const int *steer);

//---------------------------------------------------------------------------//
// GATHERS
//---------------------------------------------------------------------------//
/*!
 * \brief Gather multiple elements to all processors
 *
 * The data in \c sendbuf[0:num_els] from node \c n will be placed in
 * \c recvbuf[num_els*n:num_els*(n+1)] on all processors.
 *
 * \param sendbuf Array of length num_els to send to other processors
 * \param recvbuf Array of length num_els * nodes to receive with
 * \param num_els Number of elements in sendbuf
 */
template<class T>
void all_gather(const T* sendbuf, T* recvbuf, int num_els);

//---------------------------------------------------------------------------//
/*!
 * \brief Gather from multiple elements to a root process.
 *
 * The data in \c sendbuf[0:num_els] from node \c n will be placed in
 * \c recvbuf[num_els*n:num_els*(n+1)] on the root process.
 *
 * \param sendbuf Array of length num_els to send to other processors
 * \param recvbuf Array of length num_els * nodes to receive with
 * \param num_els Number of elements in sendbuf
 * \param root    process that data is gathered on
 */
template<class T>
void gather(const T* sendbuf, T* recvbuf, int num_els, int root);

//---------------------------------------------------------------------------//
// ALL-TO-ALLS
//---------------------------------------------------------------------------//
/*!
 * \brief All-to-all communication with constant message size.
 *
 * \param sendbuf array of length \f$ n \times N_p \f$ where elements
                  \f$ n \times p \f$ will be send to processor \c p
 *
 * \param recvbuf array of length \f$ n \times N_p \f$ where elements
                  \f$ n \times p \f$ will be received from processor \c p
 *
 * \param n number of elements to be sent/received from each processor
 */
template<class T>
void all_to_all(const T *sendbuf, T *recvbuf, int n);

//---------------------------------------------------------------------------//
/*!
 * \brief All-to-all communication with variable message size.
 *
 * \param sendbuf array of length \f$ n = \sum_{i}^{N_p} sendcounts[i]\f$
 *                containing outgoing messages
 *
 * \param sendcounts array of length \f$ N_p \f$ containing number of
 *                   elements to be sent to each processor
 *
 * \param sendoffsets array of length \f$ N_p \f$ containing offset into
 *                    sendbuf array for each processor
 *
 * \param recvbuf array of length \f$ n = \sum_{i}^{N_p} recvcounts[i]\f$
 *                containing incoming messages
 *
 * \param recvcounts array of length \f$ N_p \f$ containing number of
 *                   elements to be received from each processor
 *
 * \param recvoffsets array of length \f$ N_p \f$ containing offset into
 *                    recvbuf array for each processor
 *
 */
template<class T>
void all_to_all(const T   *sendbuf,     const int *sendcounts,
                const int *sendoffsets, T         *recvbuf,
                const int *recvcounts,  const int *recvoffsets);

//---------------------------------------------------------------------------//
// TIMING FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Return the wall-clock time in seconds.
 */
double wall_clock_time();

//---------------------------------------------------------------------------//
#ifdef HAVE_SYS_TIMES_H
/*!
 * \brief Return the wall-clock time in seconds and initialize a POSIX-timer.
 */
double wall_clock_time( tms & now );
#endif

//---------------------------------------------------------------------------//
/*!
 * \brief Return the resolution of wall_clock_time.
 */
double wall_clock_resolution();

//---------------------------------------------------------------------------//
// PROBE/WAIT FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief See if a message is pending.
 *
 * \param source
 * Processor from which a message may be pending.
 * \param tag
 * Tag for pending message.
 * \param message_size
 * On return, size of the pending message in bytes.
 * \return \c true if a message from the specified processor with the
 * specified tag is pending; \c false otherwise.
 */
bool probe(int source, int tag, int &message_size);

//---------------------------------------------------------------------------//
/*!
 * \brief Wait until a message (of unknown size) is pending.
 *
 * \param source
 * Processor from which a message of unknown size is expected.
 * \param tag
 * Tag for pending message.
 * \param message_size
 * On return, size of the pending message in bytes.
 */
void blocking_probe(int source, int tag, int &message_size);

//---------------------------------------------------------------------------//
// ABORT
//---------------------------------------------------------------------------//
/*!
 * \brief Abort across all processors.
 *
 * \param error suggested return error, defaults to 1
 */
int abort(int error = 1);

//---------------------------------------------------------------------------//
// isScalar
//---------------------------------------------------------------------------//
/*!
 * \brief Is COMM executing in scalar-only mode?
 */
bool isScalar();

} // end namespace profugus

#endif // comm_Functions_hh

//---------------------------------------------------------------------------//
//              end of comm/Functions.hh
//---------------------------------------------------------------------------//
