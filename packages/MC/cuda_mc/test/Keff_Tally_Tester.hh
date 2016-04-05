//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/Keff_Tally_Tester.hh
 * \author Steven Hamilton
 * \date   Wed Jan 20 16:13:24 2016
 * \brief  Keff_Tally_Tester class definition.
 * \note   Copyright (C) 2016 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_mc_Keff_Tally_Tester_hh
#define cuda_mc_Keff_Tally_Tester_hh

#include <vector>
#include <memory>
#include "Matprop/xs/XS.hh"

namespace cuda_mc
{

//===========================================================================//
/*!
 * \class Keff_Tally_Tester
 * \brief Wrapper to kernel launches to test Keff_Tally class.
 */
//===========================================================================//

class Keff_Tally_Tester
{
  public:

      typedef std::vector<int>              Vec_Int;
      typedef std::vector<unsigned int>     Vec_UInt;
      typedef std::vector<double>           Vec_Dbl;
      typedef Teuchos::RCP<profugus::XS>    RCP_XS;

      static void test_tally( const Vec_Dbl  &x_edges,
                              const Vec_Dbl  &y_edges,
                              const Vec_Dbl  &z_edges,
                                    RCP_XS    xs,
                                    double   &keff,
                                    int       num_particles );
};

} // end namespace cuda_mc

#endif // cuda_mc_Keff_Tally_Tester_hh

//---------------------------------------------------------------------------//
//                 end of Keff_Tally_Tester.hh
//---------------------------------------------------------------------------//
