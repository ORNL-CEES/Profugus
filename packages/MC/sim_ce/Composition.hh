//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   MC/sim_ce/Composition.hh
 * \author Thomas M Evans
 * \date   Fri Jan 08 15:23:48 2016
 * \brief  Composition class declaration.
 * \note   Copyright (c) 2016 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef MC_sim_ce_Composition_hh
#define MC_sim_ce_Composition_hh

#include <cstdint>
#include <vector>

#include "harness/DBC.hh"

namespace profugus
{

//===========================================================================//
/*!
 * \class Composition
 * \brief Simple composition definition for simulated CE data.
 */
/*!
 * \example sim_ce/test/tstComposition.cc
 *
 * Test of Composition.
 */
//===========================================================================//

class Composition
{
  public:
    //@{
    //! Typedefs.
    using size_t    = std::uint32_t;
    using Vec_Zaids = std::vector<size_t>;
    using Vec_Dbl   = std::vector<double>;
    //@}

  private:
    // ZAID.
    Vec_Zaids d_zaids;

    // Number density (/barn-cm).
    Vec_Dbl d_N;

  public:
    Composition() = default;

    //! Add zaids and number densities.
    void add(Vec_Zaids zaids, Vec_Dbl num_den)
    {
        d_zaids = std::move(zaids);
        d_N     = std::move(num_den);
    }

    //! Number of nuclides.
    size_t num_nuclides() const { return d_N.size(); }

    //@{
    //! ZAID access.
    const Vec_Zaids& zaids() const { return d_zaids; }
    size_t zaid(size_t n) const
    {
        REQUIRE(n < d_zaids.size());
        return d_zaids[n];
    }
    //@}

    //@{
    //! Number density access.
    const Vec_Dbl& number_dens() const { return d_N; }
    double number_den(size_t n) const
    {
        REQUIRE(n < d_N.size());
        return d_N[n];
    }
    //@}
};

} // end namespace profugus

#endif // MC_sim_ce_Composition_hh

//---------------------------------------------------------------------------//
// end of MC/sim_ce/Composition.hh
//---------------------------------------------------------------------------//
