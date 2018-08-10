//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   MC/mc_driver/Single_Pin_Subchannel.hh
 * \author Steven Hamilton
 * \date   Wed Jul 25 14:00:39 2018
 * \brief  Single_Pin_Subchannel class declaration.
 * \note   Copyright (c) 2018 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef MC_mc_driver_Single_Pin_Subchannel_hh
#define MC_mc_driver_Single_Pin_Subchannel_hh

#include <vector>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Utils/harness/DBC.hh"

namespace mc
{

//===========================================================================//
/*!
 * \class Single_Pin_Subchannel
 * \brief Solve subchannel flow in single channel.
 *
 * This class implements a two-equation subchannel model involving
 * conservation of energy and axial momentum.  No lateral flow is accounted
 * for.  The model excludes friction momentum losses and any spacer grid
 * effects.  It is useful for giving qualitatively-correct behavior but
 * should not be used for actual analysis.
 */
/*!
 * \example mc_driver/test/tstSingle_Pin_Subchannel.cc
 *
 * Test of Single_Pin_Subchannel.
 */
//===========================================================================//

class Single_Pin_Subchannel
{
  public:
    //@{
    //! Typedefs
    using RCP_PL = Teuchos::RCP<Teuchos::ParameterList>;
    //@}

    enum Verbosity {NONE, LOW, HIGH};

  private:
    // >>> DATA

    // Solve parameters
    double d_tol;
    int d_max_iters;
    Verbosity d_verbosity;

    // Axial grid
    std::vector<double> d_delta_z;

    // Channel geometry
    double d_area;

    // Flow conditions
    double d_mdot;
    double d_T_inlet;
    double d_p_exit;

  public:

    // Constructor
    Single_Pin_Subchannel(RCP_PL&                    parameters,
                          const std::vector<double>& delta_z);

    // Set channel area (cm^2)
    void set_channel_area(double area)
    {
        // Channel area internally is needed in m^2
        REQUIRE(area > 0.0);
        d_area = 1e-4 * area;
    }

    // Set mass flow rate (kg/s)
    void set_mass_flow_rate(double mdot)
    {
        REQUIRE(mdot > 0.0);
        REQUIRE(mdot < 2.0);
        d_mdot = mdot;
    }

    // Set inlet temperature (K)
    void set_inlet_temperature(double T)
    {
        REQUIRE(T > 0);
        REQUIRE(T < 1000);
        d_T_inlet = T;
    }

    // Set exit pressure (Pa)
    void set_exit_pressure(double p)
    {
        REQUIRE(p > 0);
        REQUIRE(p < 2.2e7);
        d_p_exit = p;
    }

    // Solve for single subchannel given power distribution
    void solve(const std::vector<double>& power,
               std::vector<double>&       temperature,
               std::vector<double>&       density);
};

//---------------------------------------------------------------------------//
} // end namespace mc

//---------------------------------------------------------------------------//
#endif // MC_mc_driver_Single_Pin_Subchannel_hh

//---------------------------------------------------------------------------//
// end of MC/mc_driver/Single_Pin_Subchannel.hh
//---------------------------------------------------------------------------//
