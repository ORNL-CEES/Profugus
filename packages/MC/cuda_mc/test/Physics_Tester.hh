//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/Physics_Tester.hh
 * \author Steven Hamilton
 * \date   Wed Jan 20 16:13:24 2016
 * \brief  Physics_Tester class definition.
 * \note   Copyright (C) 2016 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_mc_Physics_Tester_hh
#define cuda_mc_Physics_Tester_hh

#include <vector>
#include <memory>
#include "Matprop/xs/XS.hh"

namespace cuda_mc
{

//===========================================================================//
/*!
 * \class Physics_Tester
 * \brief Wrapper to kernel launches to test Physics class.
 */
//===========================================================================//

class Physics_Tester
{
  public:

      typedef std::vector<int>              Vec_Int;
      typedef std::vector<double>           Vec_Dbl;
      typedef Teuchos::RCP<profugus::XS>    RCP_XS;

      static void test_total( const Vec_Dbl  &x_edges,
                              const Vec_Dbl  &y_edges,
                              const Vec_Dbl  &z_edges,
                              const Vec_Int  &matids,
                                    RCP_XS    xs,
                                    Vec_Dbl  &totals );

      static void test_collide( const Vec_Dbl  &x_edges,
                                const Vec_Dbl  &y_edges,
                                const Vec_Dbl  &z_edges,
                                const Vec_Int  &matids,
                                      RCP_XS    xs,
                                      int       num_particles );
};

} // end namespace cuda_mc

#endif // cuda_mc_Physics_Tester_hh

//---------------------------------------------------------------------------//
//                 end of Physics_Tester.hh
//---------------------------------------------------------------------------//
