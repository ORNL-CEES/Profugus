//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   comm/Definitions.hh
 * \author Seth R Johnson
 * \date   Mon Oct 14 11:55:57 2013
 * \brief  Definitions used for MPI
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef comm_Definitions_hh
#define comm_Definitions_hh

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

#endif // comm_Definitions_hh

//---------------------------------------------------------------------------//
//                 end of Definitions.hh
//---------------------------------------------------------------------------//
