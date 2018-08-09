//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   MC/mc/Multi_Pin_Conduction.hh
 * \author Steven Hamilton
 * \date   Thu Aug 09 09:12:06 2018
 * \brief  Multi_Pin_Conduction class declaration.
 * \note   Copyright (c) 2018 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef MC_mc_Multi_Pin_Conduction_hh
#define MC_mc_Multi_Pin_Conduction_hh

#include <memory>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Assembly_Model.hh"
#include "Single_Pin_Conduction.hh"

namespace mc
{

//===========================================================================//
/*!
 * \class Multi_Pin_Conduction
 * \brief Solve heat equation over multiple fuel pins
 */
/*!
 * \example mc/test/tstMulti_Pin_Conduction.cc
 *
 * Test of Multi_Pin_Conduction.
 */
//===========================================================================//

class Multi_Pin_Conduction
{
  public:
    //@{
    //! Typedefs
    using SP_Assembly = std::shared_ptr<Assembly_Model>;
    using SP_Pin_Conduction = std::shared_ptr<Single_Pin_Conduction>;
    using RCP_PL = Teuchos::RCP<Teuchos::ParameterList>;
    //@}

  private:
    // >>> DATA
    SP_Assembly d_assembly;
    int d_Nz;

    // Single pin solver
    SP_Pin_Conduction d_pin_conduction;

  public:

    // Constructor
    Multi_Pin_Conduction(SP_Assembly                assembly,
                         RCP_PL                     parameters,
                         const std::vector<double>& dz);

    // Solve for all pins
    void solve(const std::vector<double>& power,
               const std::vector<double>& channel_temp,
                     std::vector<double>& fuel_temp);
};

//---------------------------------------------------------------------------//
} // end namespace mc

//---------------------------------------------------------------------------//
#endif // MC_mc_Multi_Pin_Conduction_hh

//---------------------------------------------------------------------------//
// end of MC/mc/Multi_Pin_Conduction.hh
//---------------------------------------------------------------------------//
