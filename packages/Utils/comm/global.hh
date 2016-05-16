//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Utils/comm/global.hh
 * \author Thomas M. Evans
 * \date   Wed Jan  2 15:39:55 2008
 * \brief  Include file for comm package.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 *
 * This file allows the client to include the message passing services
 * provided by comm.  The function declarations and class definitions are
 * contained in the profugus namespace.
 */
//---------------------------------------------------------------------------//

#ifndef Utils_comm_global_hh
#define Utils_comm_global_hh

// comm package configure
#include <Utils/config.h>

// comm Message Passing Functions
#include "Functions.hh"

// comm Request handler
#include "Request.hh"

//---------------------------------------------------------------------------//
// Load OpenMP API

#include "OMP.hh"

//---------------------------------------------------------------------------//
// Include the appropriate header for an underlying message passing
// implementation.  This allows the definition of inline functions declared
// in Functions.hh.

#ifdef COMM_SCALAR
#include "Serial.hh"
#endif

#ifdef COMM_MPI
#include "MPI.hh"
#endif

#endif // Utils_comm_global_hh

//---------------------------------------------------------------------------//
//              end of comm/global.hh
//---------------------------------------------------------------------------//
