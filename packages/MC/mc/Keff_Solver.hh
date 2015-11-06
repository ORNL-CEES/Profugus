//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   MC/mc/Keff_Solver.hh
 * \author Thomas M. Evans
 * \date   Mon Apr 13 18:14:03 2015
 * \brief  Keff_Solver class declaration.
 * \note   Copyright (c) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef MC_mc_Keff_Solver_hh
#define MC_mc_Keff_Solver_hh

#include "Fission_Matrix_Acceleration.hh"
#include "Keff_Tally.hh"
#include "Fission_Source.hh"
#include "Solver.hh"

namespace profugus
{

//===========================================================================//
/*!
 * \class Keff_Solver
 * \brief Base class for eigenvalue solvers.
 */
/*!
 * \example mc/test/tstKeff_Solver.cc
 *
 * Test of Keff_Solver.
 */
//===========================================================================//

template <class Geometry>
class Keff_Solver : public Solver<Geometry>
{
  public:
    //@{
    //! Typedefs.
    typedef Geometry                                Geometry_t;
    typedef Keff_Tally<Geometry_t>                  Keff_Tally_t;
    typedef std::shared_ptr<Keff_Tally_t>           SP_Keff_Tally;
    typedef Fission_Source<Geometry_t>              Fission_Source_t;
    typedef std::shared_ptr<Fission_Source_t>       SP_Fission_Source;
    typedef Fission_Matrix_Acceleration<Geometry>   FM_Acceleration_t;
    typedef std::shared_ptr<FM_Acceleration_t>      SP_FM_Acceleration;
    //@}

  protected:
    // >>> DATA

    // Eigenvalue tally
    SP_Keff_Tally b_keff_tally;

  public:
    //! Virtual Destructor.
    virtual ~Keff_Solver() = default;

    // >>> API

    //! Keff tally for tracking k-effective.
    SP_Keff_Tally keff_tally() const { return b_keff_tally; }

    //! Return acceleration.
    virtual SP_FM_Acceleration acceleration() const = 0;
};

} // end namespace profugus

#endif // MC_mc_Keff_Solver_hh

//---------------------------------------------------------------------------//
// end of MC/mc/Keff_Solver.hh
//---------------------------------------------------------------------------//
