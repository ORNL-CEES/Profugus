//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   MC/mc/Single_Pin_Conduction.hh
 * \author Steven Hamilton
 * \date   Thu Jul 26 09:06:36 2018
 * \brief  Single_Pin_Conduction class declaration.
 * \note   Copyright (c) 2018 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef MC_mc_Single_Pin_Conduction_hh
#define MC_mc_Single_Pin_Conduction_hh

#include <vector>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseSolver.hpp"
#include "Teuchos_SerialDenseVector.hpp"

#include "Utils/harness/DBC.hh"

namespace mc
{

//===========================================================================//
/*!
 * \class Single_Pin_Conduction
 * \brief Solve thermal conduction in single fuel pin.
 *
 * This class implements a simple 1-D radial thermal conduction equation
 * for each axial level of a fuel pin.  Axial conduction is ignored and there
 * is currently no gap conductance model.  This model should only be used for
 * very approximate qualitative behavior.
 */
/*!
 * \example mc/test/tstSingle_Pin_Conduction.cc
 *
 * Test of Single_Pin_Conduction.
 */
//===========================================================================//

class Single_Pin_Conduction
{
  public:
    //@{
    //! Typedefs
    using RCP_PL = Teuchos::RCP<Teuchos::ParameterList>;
    using Matrix = Teuchos::SerialDenseMatrix<int,double>;
    using Vector = Teuchos::SerialDenseVector<int,double>;
    using Solver = Teuchos::SerialDenseSolver<int,double>;
    //@}

  private:

    // >>> DATA

    // Solve parameters
    double d_tol;
    int d_max_iters;

    // Geometry
    double d_fuel_radius;
    double d_clad_radius;
    double d_delta_r_fuel;
    double d_delta_r_clad;
    std::vector<double> d_delta_z;

    // Thermal conductivities
    double d_k_fuel;
    double d_k_clad;

  public:

    // Constructor
    Single_Pin_Conduction(RCP_PL&                    parameters,
                          const std::vector<double>& delta_z);

    // Set fuel radius (cm)
    void set_fuel_radius(double r)
    {
        REQUIRE(r > 0.0);
        REQUIRE(r < 2.0);
        d_fuel_radius = r;
    }

    // Set clad radius (cm)
    void set_clad_radius(double r)
    {
        REQUIRE(r > 0.0);
        REQUIRE(r < 2.0);
        d_clad_radius = r;
    }

    // Solve conduction equation at each axial level
    void solve(const std::vector<double>& power,
               const std::vector<double>& channel_temperature,
               std::vector<double>&       fuel_temperature);
};

//---------------------------------------------------------------------------//
} // end namespace mc

//---------------------------------------------------------------------------//
#endif // MC_mc_Single_Pin_Conduction_hh

//---------------------------------------------------------------------------//
// end of MC/mc/Single_Pin_Conduction.hh
//---------------------------------------------------------------------------//
