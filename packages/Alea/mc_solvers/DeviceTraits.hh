//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Alea/mc_solvers/DeviceTraits.hh
 * \author Steven Hamilton
 * \brief  Templated interface for Kokkos devices.
 */
//---------------------------------------------------------------------------//

#ifndef Alea_mc_solvers_DeviceTraits_hh
#define Alea_mc_solvers_DeviceTraits_hh

#include "KokkosCore_config.h"

#ifdef KOKKOS_HAVE_OPENMP
#include <omp.h>
#endif

#include "Kokkos_hwloc.hpp"

// "New" Kokkos nodes
#include "Kokkos_Serial.hpp"
#include "Kokkos_OpenMP.hpp"
#include "Kokkos_Threads.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

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
    static inline void initialize(
        Teuchos::RCP<Teuchos::ParameterList> pl)
    {
    }

    //! \brief Finalize device.
    static inline void finalize(){}
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
    static inline void initialize(
        Teuchos::RCP<Teuchos::ParameterList> pl)
    {
        int threads = pl->get("num_threads",1);
        bool print_config = pl->get("print_config",false);
        if( !Kokkos::OpenMP::is_initialized() )
            Kokkos::OpenMP::initialize( threads );
        if( print_config )
            Kokkos::OpenMP::print_configuration(std::cout,true);
    }

    //! \brief Finalize device.
    static inline void finalize()
    {
        if( Kokkos::OpenMP::is_initialized() )
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
    static inline void initialize(
        Teuchos::RCP<Teuchos::ParameterList> pl)
    {
        int threads = pl->get("num_threads",1);
        bool print_config = pl->get("print_config",false);
        if( !Kokkos::Threads::is_initialized() )
            Kokkos::Threads::initialize( threads );
        if( print_config )
            Kokkos::Threads::print_configuration(std::cout,true);
    }

    //! \brief Finalize device.
    static inline void finalize()
    {
        if( Kokkos::Threads::is_initialized() )
            Kokkos::Threads::finalize();
    }
};

//---------------------------------------------------------------------------//
/*!
 * \class DeviceTraits<Kokkos::Cuda>
 * \brief Specialization of DeviceTraits for Kokkos::Cuda.
 */
//---------------------------------------------------------------------------//
#ifdef KOKKOS_HAVE_CUDA
template <>
class DeviceTraits<Kokkos::Cuda>
{
  public:

    //! \brief Initialize device
    static inline void initialize(
        Teuchos::RCP<Teuchos::ParameterList> pl)
    {
        int threads = pl->get("num_threads",1);
        bool print_config = pl->get("print_config",false);
        if( !HOST::is_initialized() )
        {
            HOST::initialize( threads );
        }

        if( !Kokkos::Cuda::is_initialized() )
        {
            Kokkos::Cuda::initialize();
        }

        if( print_config )
        {
            std::cout << "Host configuration:" << std::endl;
            HOST::print_configuration(std::cout,true);
            std::cout << "Device configuration:" << std::endl;
            Kokkos::Cuda::print_configuration(std::cout,true);
        }
    }

    //! \brief Finalize device.
    static inline void finalize()
    {
        if( HOST::is_initialized() )
        {
            HOST::finalize();
        }
        if( Kokkos::Cuda::is_initialized() )
        {
            Kokkos::Cuda::finalize();
        }
    }
};
#endif // KOKKOS_HAVE_CUDA

} // namespace alea

#endif // Alea_mc_solvers_DeviceTraits_hh

