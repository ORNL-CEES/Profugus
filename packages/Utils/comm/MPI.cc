//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Utils/comm/MPI.cc
 * \author Thomas M. Evans
 * \date   Wed Jan  2 13:33:33 2008
 * \brief  Comm MPI Implementation.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <Utils/config.h>
#ifdef COMM_MPI

#include <string>

#ifdef HAVE_SYS_TIMES_H
#include <sys/times.h>
#endif

#include "Functions.hh"
#include "MPI.hh"
#include "P_Stream.hh"
#include "OMP.hh"

namespace profugus
{

//---------------------------------------------------------------------------//
// MPI COMMUNICATOR
//---------------------------------------------------------------------------//

//! Default communicator to use for inheriting
Communicator_t default_communicator = MPI_COMM_WORLD;

//! Current communicator
Communicator_t communicator = default_communicator;

//! Whether a successful MPI_Init call has been performed
bool initialized = false;

//---------------------------------------------------------------------------//
// SETUP FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Changes the default internal communicator, and switches to the
 *        new default.
 */
void set_default(const Communicator_t& comm)
{
    int result = MPI_Comm_dup(comm, &default_communicator);
    ENSURE(result == MPI_SUCCESS);

    reset_internal_comm();
}

//---------------------------------------------------------------------------//

void inherit(const Communicator_t& comm)
{
    int result = MPI_Comm_dup(comm, &communicator);
    ENSURE(result == MPI_SUCCESS);
}

//---------------------------------------------------------------------------//

void split(
        int             color,
        int             key,
        Communicator_t& new_comm)
{
    int result = MPI_Comm_split(communicator, color, key, &new_comm);
    ENSURE(result == MPI_SUCCESS);
}

//---------------------------------------------------------------------------//

void split(
        int                   color,
        int                   key,
        Communicator_t&       new_comm,
        const Communicator_t& comm)
{
    int result = MPI_Comm_split(comm, color, key, &new_comm);
    ENSURE(result == MPI_SUCCESS);
}

//---------------------------------------------------------------------------//

void duplicate_comm(Communicator_t& new_comm, const Communicator_t& comm)
{
    int result = MPI_Comm_dup(comm, &new_comm);
    ENSURE(result == MPI_SUCCESS);
}

//---------------------------------------------------------------------------//

void free_comm(Communicator_t& comm)
{
    int result = MPI_Comm_free(&comm);
    ENSURE(result == MPI_SUCCESS);
}

//---------------------------------------------------------------------------//

void set_internal_comm(const Communicator_t& comm)
{
    communicator = comm;
}

//---------------------------------------------------------------------------//

void initialize(int &argc, char **&argv)
{
    int result = MPI_Init(&argc, &argv);
    initialized = (result == MPI_SUCCESS);
    CHECK( initialized );

    // Resync clocks for Darwin mpich
    MPI_Wtick();

    // Initialize globally defined parallel streams to the local node

    // conditionally set pnout
#ifdef UTILS_POUT
    pnout.set_master(node());
#endif

    // always set pncout
    pncout.set_master(node());

    // Set dynamic threading off by default
    turn_off_dynamic_threading();

    // Default to 1 thread
    set_num_threads(1);
}

//---------------------------------------------------------------------------//

void finalize()
{
    MPI_Finalize();
}

//---------------------------------------------------------------------------//

void free_inherited_comm()
{
    if (communicator != default_communicator)
    {
        MPI_Comm_free(&communicator);
        communicator = default_communicator;
    }
}

//---------------------------------------------------------------------------//

void reset_internal_comm()
{
    if (communicator != default_communicator)
    {
        communicator = default_communicator;
    }
}

//---------------------------------------------------------------------------//
// QUERY FUNCTIONS
//---------------------------------------------------------------------------//

int node()
{
    int node = 0;
    MPI_Comm_rank(communicator, &node);
    CHECK(node >= 0);
    return node;
}

//---------------------------------------------------------------------------//

int nodes()
{
    int nodes = 0;
    MPI_Comm_size(communicator, &nodes);
    CHECK(nodes > 0);
    return nodes;
}

//---------------------------------------------------------------------------//

int node(const Communicator_t& comm)
{
    int node = 0;
    MPI_Comm_rank(comm, &node);
    CHECK(node >= 0);
    return node;
}

//---------------------------------------------------------------------------//

int nodes(const Communicator_t& comm)
{
    int nodes = 0;
    MPI_Comm_size(comm, &nodes);
    CHECK(nodes > 0);
    return nodes;
}

//---------------------------------------------------------------------------//
// SETUP FUNCTIONS
//---------------------------------------------------------------------------//

std::string proc_name()
{
    std::string result = "";

    char name[MPI_MAX_PROCESSOR_NAME];
    int  len;

    MPI_Get_processor_name(name, &len);
    result.assign( name, len );

    return result;
}

//---------------------------------------------------------------------------//
// BARRIER FUNCTIONS
//---------------------------------------------------------------------------//

void global_barrier()
{
    MPI_Barrier(communicator);
}

//---------------------------------------------------------------------------//

void barrier(const Communicator_t& comm)
{
    MPI_Barrier(comm);
}

//---------------------------------------------------------------------------//
// TIMING FUNCTIONS
//---------------------------------------------------------------------------//
// overloaded function (no args)

double wall_clock_time()
{
    return MPI_Wtime();
}

//---------------------------------------------------------------------------//
// overloaded function (provide POSIX timer information).

#ifdef HAVE_SYS_TIMES_H
double wall_clock_time( tms & now )
{
    // obtain posix timer information and return it to the user via the
    // reference value argument "now".
    times( &now );
    // This funtion will return the MPI wall-clock time.
    return MPI_Wtime();
}
#endif

//---------------------------------------------------------------------------//

double wall_clock_resolution()
{
    return MPI_Wtick();
}

//---------------------------------------------------------------------------//
// PROBE/WAIT FUNCTIONS
//---------------------------------------------------------------------------//

bool probe(int  source,
           int  tag,
           int &message_size)
{
    REQUIRE(source>=0 && source<nodes());

    int flag;
    MPI_Status status;

    // post an MPI_Irecv (non-blocking receive)
    MPI_Iprobe(source, tag, communicator, &flag, &status);

    if (!flag) return false;

    MPI_Get_count(&status, MPI_CHAR, &message_size);

    return true;
}

//---------------------------------------------------------------------------//

void blocking_probe(int  source,
                    int  tag,
                    int &message_size)
{
    REQUIRE(source>=0 && source<nodes());

    MPI_Status status;
    MPI_Probe(source, tag, communicator, &status);
    MPI_Get_count(&status, MPI_CHAR, &message_size);
}

//---------------------------------------------------------------------------//
// ABORT
//---------------------------------------------------------------------------//

int abort(int error)
{
    // This test is not recorded as tested by BullseyeCoverage because abort
    // terminates the execution and BullseyeCoverage only reports coverage for
    // function that return control to main().

    int rerror = MPI_Abort(communicator, error);
    return rerror;
}

//---------------------------------------------------------------------------//

bool isScalar()
{
    return ! initialized;
}

} // end namespace profugus

#endif // COMM_MPI

//---------------------------------------------------------------------------//
//                 end of MPI.cc
//---------------------------------------------------------------------------//
