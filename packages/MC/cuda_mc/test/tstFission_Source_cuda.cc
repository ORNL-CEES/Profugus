//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/test/tstFission_Source_cuda.cc
 * \author Stuart Slattery
 * \date   Tue May 06 11:54:26 2014
 * \brief  Fission_Source unit test.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <memory>

#include "Teuchos_Array.hpp"
#include "Teuchos_RCP.hpp"

#include "comm/P_Stream.hh"
#include "geometry/Cartesian_Mesh.hh"
#include "../Fission_Source.hh"

#include "gtest/utils_gtest.hh"
#include "SourceTestBase.hh"

//---------------------------------------------------------------------------//
//                 end of tstFission_Source.cc
//---------------------------------------------------------------------------//
