//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   MC/mc/Multi_Pin_Subchannel.hh
 * \author Steven Hamilton
 * \date   Thu Aug 09 10:48:27 2018
 * \brief  Multi_Pin_Subchannel class declaration.
 * \note   Copyright (c) 2018 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef MC_mc_Multi_Pin_Subchannel_hh
#define MC_mc_Multi_Pin_Subchannel_hh

#include <memory>
#include <vector>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Assembly_Model.hh"
#include "Single_Pin_Subchannel.hh"

namespace mc
{

//===========================================================================//
/*!
 * \class Multi_Pin_Subchannel
 * \brief Solver channel-centered subchannel equations over multiple pins.
 */
/*!
 * \example mc/test/tstMulti_Pin_Subchannel.cc
 *
 * Test of Multi_Pin_Subchannel.
 */
//===========================================================================//

class Multi_Pin_Subchannel
{
  public:
    //@{
    //! Typedefs
    using SP_Assembly = std::shared_ptr<Assembly_Model>;
    using RCP_PL = Teuchos::RCP<Teuchos::ParameterList>;
    using SP_Pin_Subchannel = std::shared_ptr<Single_Pin_Subchannel>;
    //@}

  private:
    // >>> DATA
    SP_Assembly d_assembly;

    // Note that for channel-centered equations, Nx and Ny are one greater
    //  than number of pins in each direction.
    int d_Nx, d_Ny, d_Nz;

    // Cross sectional area of each channel
    std::vector<double> d_areas;

    // Mass flow rate (kg/s) in each channel
    std::vector<double> d_mdots;

    SP_Pin_Subchannel d_pin_subchannel;

  public:

    // Constructor
    Multi_Pin_Subchannel(SP_Assembly                assembly,
                         RCP_PL                     params,
                         const std::vector<double>& dz);

    // Solve
    void solve(const std::vector<double>& power,
                     std::vector<double>& channel_temp,
                     std::vector<double>& channel_density);

  private:

    int channel_index(int ix, int iy) const
    {
        REQUIRE(ix < d_Nx);
        REQUIRE(iy < d_Ny);
        return ix + d_Nx * iy;
    }
};

//---------------------------------------------------------------------------//
} // end namespace mc

//---------------------------------------------------------------------------//
#endif // MC_mc_Multi_Pin_Subchannel_hh

//---------------------------------------------------------------------------//
// end of MC/mc/Multi_Pin_Subchannel.hh
//---------------------------------------------------------------------------//
