//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/Cell_Tally_Tester.hh
 * \author Steven Hamilton
 * \date   Wed Jan 20 16:13:24 2016
 * \brief  Cell_Tally_Tester class definition.
 * \note   Copyright (C) 2016 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_mc_Cell_Tally_Tester_hh
#define cuda_mc_Cell_Tally_Tester_hh

#include <vector>
#include <memory>
#include "Matprop/xs/XS.hh"

namespace cuda_mc
{

//===========================================================================//
/*!
 * \class Cell_Tally_Tester
 * \brief Wrapper to kernel launches to test Cell_Tally class.
 */
//===========================================================================//

class Cell_Tally_Tester
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
                              const Vec_Int  &cells,
                                    Vec_Dbl  &tally_result,
                                    int       num_particles );
};

} // end namespace cuda_mc

#endif // cuda_mc_Cell_Tally_Tester_hh

//---------------------------------------------------------------------------//
//                 end of Cell_Tally_Tester.hh
//---------------------------------------------------------------------------//
