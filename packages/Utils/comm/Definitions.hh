//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Utils/comm/Definitions.hh
 * \author Seth R Johnson
 * \date   Mon Oct 14 11:55:57 2013
 * \brief  Definitions used for MPI
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef Utils_comm_Definitions_hh
#define Utils_comm_Definitions_hh

#include <Utils/config.h>

#ifdef COMM_MPI
#include <mpi.h>
#endif

namespace profugus
{
//---------------------------------------------------------------------------//

#ifdef COMM_MPI
//! Use an MPI communicator
typedef MPI_Comm Communicator_t;
#else
//! Use a pretend communicator
typedef int Communicator_t;
#endif


//---------------------------------------------------------------------------//
/*!
 * \brief COMM Success tag.
 */
const int COMM_SUCCESS = 0;

//---------------------------------------------------------------------------//

} // end namespace profugus

#endif // Utils_comm_Definitions_hh

//---------------------------------------------------------------------------//
//                 end of Definitions.hh
//---------------------------------------------------------------------------//
