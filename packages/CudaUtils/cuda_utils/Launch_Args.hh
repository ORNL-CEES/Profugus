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
#include "Definitions.hh"
#include "Hardware.hh"

#include <type_traits>

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
    REQUIRE(LAUNCH_ARGS.is_valid()); \
    GLOBAL_FUNCTION<<< \
        LAUNCH_ARGS.grid_size(), \
        LAUNCH_ARGS.block_size(), \
        LAUNCH_ARGS.shared_mem(), \
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
class Launch_Args
{
  public:
    // >>> TYPEDEFS
    typedef Stream<Arch_T> Stream_t;

    // >>> CONSTRUCTION

    //! Zero out launch arguments on construction
    Launch_Args()
    {
        d_grid_size    = 0;
        d_block_size   = Hardware<Arch_T>::default_block_size();
        d_shared_mem   = 0;
        d_num_elements = 0;
    }

    //! Set block size
    void set_block_size(unsigned int block_size)
    {
        // Block size can only be modified before grid_size/num_elements
        REQUIRE( d_grid_size    == 0 );
        REQUIRE( d_num_elements == 0 );
        REQUIRE( block_size > 0 );
        d_block_size = block_size;
    }

    //! Set grid size
    void set_grid_size(unsigned int grid_size)
    {
        REQUIRE( grid_size > 0 );
        REQUIRE( d_block_size > 0 );
        REQUIRE( d_num_elements == 0 );
        d_grid_size = grid_size;
        d_num_elements = d_grid_size * d_block_size;
        ENSURE( d_num_elements > 0 );
    }

    //! Set number of elements for kernel
    void set_num_elements(std::size_t num_elements)
    {
        REQUIRE( num_elements > 0 );
        REQUIRE( d_block_size > 0 );
        REQUIRE( d_grid_size == 0 );
        d_num_elements = num_elements;
        d_grid_size = d_num_elements / d_block_size;
        if( d_num_elements % d_block_size > 0 )
            d_grid_size++;
    }

    // >>> ACCESSORS

    //! Get block size
    std::size_t block_size() const
    {
        REQUIRE( d_block_size > 0 );
        return d_block_size;
    }

    //! Get grid size
    std::size_t grid_size() const
    {
        REQUIRE( d_grid_size > 0 );
        return d_grid_size;
    }

    //! Get number of elements
    std::size_t num_elements() const
    {
        REQUIRE( d_num_elements > 0 );
        return d_num_elements;
    }

    //! Get shared memory size
    std::size_t shared_mem() const
    {
        return d_shared_mem;
    }

    //! Set stream
    void set_stream( Stream_t stream )
    {
        d_stream = stream;
    }

    //! Stream handle for CUDA (defaults to zero)
    typename Stream_t::stream_t stream_handle() const
    {
        return d_stream.handle();
    }

    //! Check whether launch arguments are valid
    bool is_valid() const
    {
        return (d_grid_size > 0)
            && (d_block_size > 0)
            && (d_num_elements <= d_grid_size*d_block_size);
        // could also add checks for being less than hardware capacities
    }

  private:

    // >>> PRIVATE DATA
    unsigned int d_grid_size;
    unsigned int d_block_size;
    std::size_t  d_shared_mem;
    Stream_t     d_stream;
    std::size_t  d_num_elements;
};

//---------------------------------------------------------------------------//
// PARALLEL LAUNCH
//---------------------------------------------------------------------------//
// Launch a function on the host or the device.

template <class Kernel>
void parallel_launch( Kernel& kernel,
                      const Launch_Args<cuda::arch::Host>& launch_args );

template <class Kernel>
void parallel_launch( Kernel& kernel,
                      const Launch_Args<cuda::arch::Device>& launch_args );
//---------------------------------------------------------------------------//

} // end namespace cuda

#endif // cuda_utils_Launch_Args_hh

//---------------------------------------------------------------------------//
//                 end of Launch_Args.hh
//---------------------------------------------------------------------------//
