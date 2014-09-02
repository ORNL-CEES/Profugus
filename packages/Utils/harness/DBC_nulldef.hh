//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   harness/DBC_nulldef.hh
 * \author Seth R Johnson
 * \date   Thu Oct  3 10:51:42 2013
 * \brief  Null-op DBC Macro definitions
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 *
 * Turn all DBC macros into null-ops. This is primarily meant for use in CUDA
 * device code.
 *
 * \note There should be no include guards in this file.
 */
//---------------------------------------------------------------------------//

#define REQUIRE(c)
#define CHECK(c)
#define ENSURE(c)
#define Remember(c)
#define Assert(c)
#define Insist(c, m)
#define NOT_IMPLEMENTED(m)
#define VALIDATE(c, m)

//---------------------------------------------------------------------------//
//              end of harness/DBC_nulldef.hh
//---------------------------------------------------------------------------//
