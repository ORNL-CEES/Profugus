//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   MC/cuda_mc/test/Test_XS.hh
 * \author Steven Hamilton
 * \date   Wed May 04 16:31:30 2016
 * \brief  Test_XS class declaration.
 * \note   Copyright (c) 2016 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef MC_cuda_mc_test_Test_XS_hh
#define MC_cuda_mc_test_Test_XS_hh

#include "xs/XS.hh"
#include "Teuchos_RCP.hpp"

//===========================================================================//
/*!
 * \class Test_XS
 * \brief Helper class for building XS objects
 */
//===========================================================================//

class Test_XS
{

  public:

      static Teuchos::RCP<profugus::XS> build_xs(int num_groups);

  private:

      static Teuchos::RCP<profugus::XS> build_1grp_xs();
      static Teuchos::RCP<profugus::XS> build_2grp_xs();
      static Teuchos::RCP<profugus::XS> build_3grp_xs();
      static Teuchos::RCP<profugus::XS> build_5grp_xs();
};

#endif // MC_cuda_mc_test_Test_XS_hh

//---------------------------------------------------------------------------//
// end of MC/cuda_mc/test/Test_XS.hh
//---------------------------------------------------------------------------//
