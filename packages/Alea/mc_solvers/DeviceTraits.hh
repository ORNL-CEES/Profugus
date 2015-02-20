//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   DeviceTraits.hh
 * \author Steven Hamilton
 * \brief  Templated interface for Kokkos devices.
 */
//---------------------------------------------------------------------------//

#ifndef mc_solvers_DeviceTraits_hh
#define mc_solvers_DeviceTraits_hh

#include "KokkosCore_config.h"

#ifdef KOKKOS_HAVE_OPENMP
#include <omp.h>
#endif

#include "Kokkos_hwloc.hpp"

// "New" Kokkos nodes
#include "Kokkos_Serial.hpp"
#include "Kokkos_OpenMP.hpp"
#include "Kokkos_Threads.hpp"

namespace alea
{

//---------------------------------------------------------------------------//
/*!
 * \class DeviceTraits
 * \brief Templated interface for initializing/finalizing Kokkos devices.
 *
 * This is a generic implementation for a single thread with trivial
 * initialization/finalization which should be suitable for Kokkos::Serial.
 */
//---------------------------------------------------------------------------//
template <class DeviceType>
class DeviceTraits
{
  public:

    //! \brief Initialize device on specified number of threads.
    static inline void initializeDevice(int threads)
    {
        TEUCHOS_ASSERT(threads==1);
    }

    //! \brief Finalize device.
    static inline void finalizeDevice(){}
};

//---------------------------------------------------------------------------//
/*!
 * \class DeviceTraits<Kokkos::OpenMP>
 * \brief Specialization of DeviceTraits for Kokkos::OpenMP.
 */
//---------------------------------------------------------------------------//
#ifdef KOKKOS_HAVE_OPENMP
template <>
class DeviceTraits<Kokkos::OpenMP>
{
  public:

    //! \brief Initialize device on specified number of threads.
    static inline void initializeDevice(int threads)
    {
        if( !Kokkos::OpenMP::is_initialized() )
            Kokkos::OpenMP::initialize( threads );
    }

    //! \brief Finalize device.
    static inline void finalizeDevice()
    {
        Kokkos::OpenMP::finalize();
    }
};
#endif

//---------------------------------------------------------------------------//
/*!
 * \class DeviceTraits<Kokkos::Threads>
 * \brief Specialization of DeviceTraits for Kokkos::Threads.
 */
//---------------------------------------------------------------------------//
template <>
class DeviceTraits<Kokkos::Threads>
{
  public:

    //! \brief Initialize device on specified number of threads.
    static inline void initializeDevice(int threads)
    {
        if( !Kokkos::Threads::is_initialized() )
            Kokkos::Threads::initialize( threads );
        //Kokkos::Threads::print_configuration(std::cout,true);
    }

    //! \brief Finalize device.
    static inline void finalizeDevice()
    {
        Kokkos::Threads::finalize();
    }
};

} // namespace alea

#endif // mc_solvers_DeviceTraits_hh

