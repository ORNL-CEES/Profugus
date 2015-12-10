//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_utils/Cuda_RNG.hh
 * \author Stuart Slattery
 * \date   Thu Dec 10 09:42:28 2015
 * \brief  Cuda RNG definition.
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_utils_RNG_hh
#define cuda_utils_RNG_hh

#include <memory>

#include "harness/DBC.hh"

#include "rng/RNG_Control.hh"

#include "Definitions.hh"
#include "Device_Vector.hh"

#include <cuda_runtime.h>
#include <curand_kernel.h>

namespace cuda
{
//---------------------------------------------------------------------------//
/*
 * \class Global CUDA RNG Controller for host-side implementations.
 */
class Cuda_Global_RNG_Control
{
  public:
    static profugus::RNG_Control d_rng_control;
    static void initialize( const int global_seed );
};

//---------------------------------------------------------------------------//
/*!
 * \class Cuda_RNG_State_Traits
 */
template<typename Arch_T>
class Cuda_RNG_State_Traits;

//---------------------------------------------------------------------------//
// Host specialization.
template<>
class Cuda_RNG_State_Traits<cuda::arch::Host>
{
  public:

    // Arch type.
    using Arch_t = cuda::arch::Host;

    // RNG state type.
    using RNG_State = profugus::RNG;

    // RNG state creator.
    static RNG_State create()
    { return profugus::RNG(); }

    // RNG state destructor.
    static void destroy( RNG_State& rng_state )
    { /* ... */ }

    // Initialize the state with a seed.
    static inline void initialize( RNG_State& rng_state, const int seed )
    {
	rng_state = Cuda_Global_RNG_Control::d_rng_control.rng( seed );
    }

    // Get a random number uniformly between 0 and 1 from the state.
    static inline double ran( RNG_State& rng_state )
    {
	return rng_state.ran();
    }
};

//---------------------------------------------------------------------------//
// Device specialization.
template<>
class Cuda_RNG_State_Traits<cuda::arch::Device>
{
  public:

    // Arch type.
    using Arch_t = cuda::arch::Device;

    // RNG state type.
    using RNG_State = curandState*;

    // RNG state creator.
    static RNG_State create()
    { 
	curandState* state;
	cudaMalloc( (void**) &state, sizeof(curandState) );
	return state;
    }

    // RNG state destructor.
    static void destroy( RNG_State& rng_state )
    { cudaFree( rng_state ); }

    // Initialize the state with a seed.
    __device__ static inline void 
    initialize( RNG_State& rng_state, const int seed )
    { curand_init( seed, 0, 0, rng_state ); }

    // Get a random number uniformly between 0 and 1 from the state.
    __device__ static inline double ran( RNG_State& rng_state )
    { return curand( rng_state ); }
};

//---------------------------------------------------------------------------//
/*!
 * \class Cuda_RNG
 */
template<class Arch_T>
class Cuda_RNG
{
  public:

    //! Type aliases.
    using Arch_t = Arch_T;
    using rng_state_traits = Cuda_RNG_State_Traits<Arch_t>; 
    using RNG_State = typename rng_state_traits::RNG_State;
    
  public:
    
    //! Constructor.
    Cuda_RNG() : d_rng_state( rng_state_traits::create() )
    { /* ... */ }

    //! Destructor.
    ~Cuda_RNG() { rng_state_traits::destroy(d_rng_state); }

    // Use default copy and assignment.
    Cuda_RNG( const Cuda_RNG& ) = default;
    Cuda_RNG& operator=( const Cuda_RNG& ) = default;

    // Initialize a random number generator with a seed.
    __device__ inline void initialize( const int seed )
    { rng_state_traits::initialize(d_rng_state,seed); }

    // Get a random number uniformly between 0 and 1.
    __device__ inline double ran()
    { return rng_state_traits::ran(d_rng_state); }

  private:

    // RNG state.
    RNG_State d_rng_state;
};

//---------------------------------------------------------------------------//
/*!
 * \class Cuda_RNG_Vector_Init_Kernel
 *
 * \brief Kernel for initializing a device vector of rngs.
 */
template<class Arch_T>
class Cuda_RNG_Vector_Init_Kernel
{
  public:
    
    // Type aliases.
    using Arch_t = Arch_T;
    using RNG = Cuda_RNG<Arch_t>;
    template<class T> using Device_Vector = cuda::Device_Vector<Arch_t,T>;
    template<class T> using SP_Device_Vector = std::shared_ptr<Device_Vector<T> >;

    // Constructor.
    Cuda_RNG_Vector_Init_Kernel( const SP_Device_Vector<RNG>& rng,
				 const SP_Device_Vector<int>& seeds )
	: d_rng( rng )
	, d_device_rng( rng.data() )
	, d_seeds( seeds )
	, d_device_seeds( seeds.data() )
    {
	REQUIRE( d_rng.size() == d_seeds.size() );
    }

    // Fill operator.
    __host__ __device__ void operator()( const std::size_t i )
    {
        d_device_rng[i].initialize( d_device_seeds[i] );
    }
    
  private:

    // RNG vector to fill the vector with.
    SP_Device_Vector<RNG> d_rng;

    // Raw pointer to rngs.
    RNG* d_device_rng;

    // Device vector to fill.
    SP_Device_Vector<double> d_seeds;

    // Raw pointer to vector seeds.
    double* d_device_seeds;
};

//---------------------------------------------------------------------------//
/*!
 * \class RNG_Vector_Fill_Kernel
 *
 * \brief Kernel for filling a device vector with RNG data.
 */
template<class Arch_T>
class Cuda_RNG_Vector_Fill_Kernel
{
  public:
    
    // Type aliases.
    using Arch_t = Arch_T;
    using RNG = Cuda_RNG<Arch_t>;
    template<class T> using Device_Vector = cuda::Device_Vector<Arch_t,T>;
    template<class T> using SP_Device_Vector = std::shared_ptr<Device_Vector<T> >;

    // Constructor.
    Cuda_RNG_Vector_Fill_Kernel( const SP_Device_Vector<RNG>& rng,
				 const SP_Device_Vector<double>& data )
	: d_rng( rng )
	, d_device_rng( rng->data() )
	, d_data( data )
	, d_device_data( data->data() )
    {
	REQUIRE( d_rng->size() == d_data->size() );
    }

    // Fill operator.
    __host__ __device__ void operator()( const std::size_t i )
    {
        d_device_data[i] = d_device_rng[i].ran();
    }
    
  private:

    // RNG vector to fill the vector with.
    SP_Device_Vector<RNG> d_rng;

    // Raw pointer to rngs.
    RNG* d_device_rng;

    // Device vector to fill.
    SP_Device_Vector<double> d_data;

    // Raw pointer to vector data.
    double* d_device_data;
};

//---------------------------------------------------------------------------//

} // end namespace cuda

#endif // cuda_utils_RNG_hh

//---------------------------------------------------------------------------//
//              end of cuda_utils/Cuda RNG.hh
//---------------------------------------------------------------------------//
