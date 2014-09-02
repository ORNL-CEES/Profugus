//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   comm/Serial.cc
 * \author Thomas M. Evans
 * \date   Wed Jan  2 14:32:30 2008
 * \brief  Serial communication function definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <Utils/config.h>
#ifdef COMM_SCALAR

#ifdef HAVE_SYS_TIMES_H
#include <unistd.h>    // defines sysconf
#include <sys/times.h> // defines tms
#endif

#include <cstdlib>
#include <string>

#include "Serial.hh"

namespace profugus
{

Communicator_t default_communicator = 0;
Communicator_t communicator = 0;

//---------------------------------------------------------------------------//
// CUSTOM COMMUNICATOR FUNCTIONS
//---------------------------------------------------------------------------//
void set_default(const Communicator_t&)
{
}

//---------------------------------------------------------------------------//

void inherit(const Communicator_t&)
{
}

//---------------------------------------------------------------------------//

void split(int, int, Communicator_t&)
{
}

//---------------------------------------------------------------------------//

void split( int, int, Communicator_t&, const Communicator_t&)
{
}

//---------------------------------------------------------------------------//

void duplicate_comm(Communicator_t&, const Communicator_t&)
{
}

//---------------------------------------------------------------------------//

void free_comm(Communicator_t&)
{
}

//---------------------------------------------------------------------------//

void set_internal_comm(const Communicator_t&)
{
}

//---------------------------------------------------------------------------//

int node(const Communicator_t&)
{
    return 0;
}

//---------------------------------------------------------------------------//

int nodes(const Communicator_t&)
{
    return 1;
}

//---------------------------------------------------------------------------//

void barrier(const Communicator_t&)
{
}

//---------------------------------------------------------------------------//
// SETUP FUNCTIONS
//---------------------------------------------------------------------------//

void initialize(int &argc, char **&argv)
{
}

//---------------------------------------------------------------------------//

void finalize()
{
}

//---------------------------------------------------------------------------//

void free_inherited_comm()
{
}

//---------------------------------------------------------------------------//

void reset_internal_comm()
{
}

//---------------------------------------------------------------------------//
// QUERY FUNCTIONS
//---------------------------------------------------------------------------//

int node()
{
    return 0;
}

//---------------------------------------------------------------------------//

int nodes()
{
    return 1;
}

//---------------------------------------------------------------------------//

std::string proc_name()
{
    std::string result;

    char name[] = "";
    int  len = 0;
    result.assign( name, len );

    return result;
}

//---------------------------------------------------------------------------//
// BARRIER FUNCTIONS
//---------------------------------------------------------------------------//

void global_barrier()
{
}

//---------------------------------------------------------------------------//
// TIMING FUNCTIONS
//---------------------------------------------------------------------------//

double wall_clock_time()
{
#ifdef HAVE_SYS_TIMES_H
    tms now;
    return times( &now ) / wall_clock_resolution();
#else
    return 0.0;
#endif
}

//---------------------------------------------------------------------------//

#ifdef HAVE_SYS_TIMES_H
double wall_clock_time( tms & now )
{
    return times( &now ) / wall_clock_resolution();
}
#endif

//---------------------------------------------------------------------------//

double wall_clock_resolution()
{
#ifdef HAVE_SYS_TIMES_H
    return static_cast<double>(sysconf(_SC_CLK_TCK));
#else
    return 1.0;
#endif
}

//---------------------------------------------------------------------------//
// PROBE/WAIT FUNCTIONS
//---------------------------------------------------------------------------//

bool probe(int  /* source */,
           int  /* tag */,
           int &/* message_size */)
{
    return false;
}

void blocking_probe(int  /* source */,
                    int  /* tag */,
                    int &/* message_size */)
{
    INSIST(false, "no messages expected in serial programs!");
}

//---------------------------------------------------------------------------//
// ABORT
//---------------------------------------------------------------------------//

int abort(int error)
{
    // This test is not recorded as tested by BullseyeCoverage because abort
    // terminates the execution and BullseyeCoverage only reports coverage for
    // function that return control to main().

    // call system exit
    std::abort();
    return error;
}

//---------------------------------------------------------------------------//

bool isScalar()
{
    return true;
}

//---------------------------------------------------------------------------//
// NON-BLOCKING BUFFER DATA
//---------------------------------------------------------------------------//

namespace internals
{
std::map<int, void*> buffers;
}

} // end namespace profugus

#endif // COMM_SCALAR

//---------------------------------------------------------------------------//
//                 end of Serial.cc
//---------------------------------------------------------------------------//
