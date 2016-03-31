//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   MC/cuda_mc/Manager_Builder_Cuda.hh
 * \author Steven Hamilton
 * \date   Wed Mar 30 15:02:55 2016
 * \brief  Manager_Builder_Cuda class declaration.
 * \note   Copyright (c) 2016 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef MC_cuda_mc_Manager_Builder_Cuda_hh
#define MC_cuda_mc_Manager_Builder_Cuda_hh

#include <memory>

#include "mc_driver/Manager_Base.hh"

namespace cuda_mc
{

//===========================================================================//
/*!
 * \class Manager_Builder_Cuda
 * \brief Build Cuda Manager.
 *
 * This class primarily exists to allow CPU-only driver-level classes to
 * construct the manager.  The actual Cuda manager has Cuda includes in the
 * header and this class provides a CPU-accessible interface.
 */
//===========================================================================//

class Manager_Builder_Cuda
{
  public:
    typedef std::shared_ptr<mc::Manager_Base>       SP_Manager_Base;
    typedef Teuchos::RCP<Teuchos::ParameterList>    RCP_ParameterList;

    // Build manager from PL
    static SP_Manager_Base build( const RCP_ParameterList &pl );
};

//---------------------------------------------------------------------------//
} // end namespace cuda_mc

//---------------------------------------------------------------------------//
#endif // MC_cuda_mc_Manager_Builder_Cuda_hh

//---------------------------------------------------------------------------//
// end of MC/cuda_mc/Manager_Builder_Cuda.hh
//---------------------------------------------------------------------------//
