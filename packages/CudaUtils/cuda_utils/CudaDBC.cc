//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   CudaUtils/cuda_utils/CudaDBC.cc
 * \author Seth R Johnson
 * \date   Thu Jun 27 15:21:33 2013
 * \brief  DBC member definitions.
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "CudaDBC.hh"

#include "harness/DBC.hh"

namespace cuda
{

//---------------------------------------------------------------------------//
/*!
 * \brief Throw a profugus::assertion for CudaCall/CudaInsist macros.
 *
 * This function is not compiled if the CUDA runtime is not installed.
 */
void toss_cuda_cookies(
        const char* message,
        const char* detail,
        const char* file,
        int         line)
{
    std::ostringstream out;
    out << "CUDA runtime failure: " << message
        << " {" << detail << "}";
#if PROFUGUS_DBC > 0
    out << "\n ^^^ at " << file << ":" << line << "\n";
#else
    // Debug mode is off, so don't trouble the user with the particulars
    (void)sizeof(file);
    (void)sizeof(line);
#endif

    throw ::profugus::assertion(out.str());
}
//---------------------------------------------------------------------------//

} // end namespace cuda

//---------------------------------------------------------------------------//
//                 end of CudaDBC.cc
//---------------------------------------------------------------------------//
