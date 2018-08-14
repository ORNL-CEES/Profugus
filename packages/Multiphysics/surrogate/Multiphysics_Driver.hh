//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   MC/mc/Multiphysics_Driver.hh
 * \author Steven Hamilton
 * \date   Fri Aug 10 09:00:11 2018
 * \brief  Multiphysics_Driver class declaration.
 * \note   Copyright (c) 2018 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef MC_mc_Multiphysics_Driver_hh
#define MC_mc_Multiphysics_Driver_hh

#include <memory>
#include <vector>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Assembly_Model.hh"
#include "Multi_Pin_Conduction.hh"
#include "Multi_Pin_Subchannel.hh"
#include "Two_Group_Diffusion.hh"

namespace mc
{

//===========================================================================//
/*!
 * \class Multiphysics_Driver
 * \brief Simple driver for multiphysics coupling.
 */
/*!
 * \example mc/test/tstMultiphysics_Driver.cc
 *
 * Test of Multiphysics_Driver.
 */
//===========================================================================//

class Multiphysics_Driver
{
  public:
    //@{
    //! Typedefs
    using SP_Assembly   = std::shared_ptr<Assembly_Model>;
    using SP_Subchannel = std::shared_ptr<Multi_Pin_Subchannel>;
    using SP_Conduction = std::shared_ptr<Multi_Pin_Conduction>;
    using SP_Diffusion  = std::shared_ptr<Two_Group_Diffusion>;
    using RCP_PL        = Teuchos::RCP<Teuchos::ParameterList>;
    using Vec_Dbl       = std::vector<double>;
    using Vec_BC        = std::vector<Two_Group_Diffusion::BC_TYPE>;
    //@}

  private:
    // >>> DATA
    SP_Assembly   d_assembly;
    SP_Subchannel d_subchannel;
    SP_Conduction d_conduction;
    SP_Diffusion  d_diffusion;

    int d_Nx, d_Ny, d_Nz, d_N;
    Vec_Dbl d_dz;

    double d_tol;
    int d_max_iters;
    double d_damping;

    double d_power_magnitude;
    Vec_Dbl d_power;
    Vec_Dbl d_fuel_temperature;
    Vec_Dbl d_coolant_temperature;
    Vec_Dbl d_coolant_density;

  public:

    // Constructor
    Multiphysics_Driver(SP_Assembly    assembly,
                        RCP_PL         params,
                        const Vec_Dbl& dz,
                        const Vec_BC&  bcs);

    // Solve problem
    void solve();

    // Accessors
    const Vec_Dbl& power()               const {return d_power;}
    const Vec_Dbl& fuel_temperature()    const {return d_fuel_temperature;}
    const Vec_Dbl& coolant_temperature() const {return d_coolant_temperature;}
    const Vec_Dbl& coolant_density()     const {return d_coolant_density;}
};

//---------------------------------------------------------------------------//
} // end namespace mc

//---------------------------------------------------------------------------//
#endif // MC_mc_Multiphysics_Driver_hh

//---------------------------------------------------------------------------//
// end of MC/mc/Multiphysics_Driver.hh
//---------------------------------------------------------------------------//
