//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   comm/global.hh
 * \author Thomas M. Evans
 * \date   Wed Jan  2 15:39:55 2008
 * \brief  Include file for comm package.
 * \note   Copyright (C) 2008 Oak Ridge National Laboratory, UT-Battelle, LLC.
 *
 * This file allows the client to include the message passing services
 * provided by comm.  The function declarations and class definitions are
 * contained in the nemesis namespace.
 */
//---------------------------------------------------------------------------//

#ifndef comm_global_hh
#define comm_global_hh

// comm package configure
#include <comm/config.h>

// comm Message Passing Functions
#include "Functions.hh"

// comm Request handler
#include "Request.hh"

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

#endif // comm_global_hh

//---------------------------------------------------------------------------//
//              end of comm/global.hh
//---------------------------------------------------------------------------//
