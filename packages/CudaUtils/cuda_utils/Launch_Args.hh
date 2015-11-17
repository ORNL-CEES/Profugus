//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_utils/Launch_Args.hh
 * \author Seth R Johnson
 * \date   Wed Oct 02 13:16:37 2013
 * \brief  Launch class definition.
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_utils_Launch_Args_hh
#define cuda_utils_Launch_Args_hh

#include "harness/DBC.hh"
#include "Stream.hh"

namespace cuda
{
//===========================================================================//
#ifdef __NVCC__
/*! \brief For CUDA, use the kernel launch syntax
 *
 * This should be called like: \code
    CUDA_LAUNCH(kernel<Arch_t>, kd.launch_args)(kernel_args);
 * \endcode
 *
 * but note that if your kernel has multiple template parameters, you'll have
 * to launch with parentheses around the first macro argument for the sake of
 * the C preprocessor.
 */
#define CUDA_LAUNCH(GLOBAL_FUNCTION, LAUNCH_ARGS) \
    Require(LAUNCH_ARGS.is_valid()); \
    GLOBAL_FUNCTION<<< \
        LAUNCH_ARGS.grid_size, \
        LAUNCH_ARGS.block_size, \
        LAUNCH_ARGS.shared_mem, \
        LAUNCH_ARGS.stream_handle()>>>
#else
//! If not using CUDA, just call the function
#define CUDA_LAUNCH(GLOBAL_FUNCTION, LAUNCH_ARGS) \
    GLOBAL_FUNCTION
#endif

//===========================================================================//
/*!
 * \struct Launch_Args
 * \brief  Arguments used to launch a 1-D CUDA kernel.
 */
//===========================================================================//

template<class Arch_T>
struct Launch_Args
{
    // >>> TYPEDEFS
    typedef Stream<Arch_T> Stream_t;

    // >>> CONSTRUCTION

    //! Zero out launch arguments on construction
    Launch_Args()
    {
        grid_size  = 0;
        block_size = 0;
        shared_mem = 0;
    }

    // >>> ACCESSORS

    //! Stream handle for CUDA (defaults to zero)
    typename Stream_t::stream_t stream_handle() const
    {
        return stream.handle();
    }

    //! Check whether launch arguments are valid
    bool is_valid() const
    {
        return (grid_size > 0)
            && (block_size > 0);
        // could also add checks for being less than hardware capacities
    }

    // >>> PUBLIC DATA

    unsigned int grid_size;
    unsigned int block_size;
    size_t       shared_mem;
    Stream_t     stream;
};

} // end namespace cuda

#endif // cuda_utils_Launch_Args_hh

//---------------------------------------------------------------------------//
//                 end of Launch_Args.hh
//---------------------------------------------------------------------------//
