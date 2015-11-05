//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/Solver.hh
 * \author Thomas M. Evans and Seth Johnson
 * \date   Tue May 13 14:39:00 2014
 * \brief  Solver class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef mc_Solver_hh
#define mc_Solver_hh

#include <memory>

#include "Tallier.hh"
#include "Physics.hh"

// Remove this once templated
#include "geometry/RTK_Geometry.hh"

namespace profugus
{

//===========================================================================//
/*!
 * \class Solver
 * \brief Base class for Monte Carlo top-level solvers.
 */
//===========================================================================//

class Solver
{
  public:
    //@{
    //! Typedefs.
    typedef Physics<Core>              Physics_t;
    typedef Physics_t::Geometry_t      Geometry_t;
    typedef Tallier                    Tallier_t;
    typedef std::shared_ptr<Tallier_t> SP_Tallier;
    //@}

  protected:
    // >>> DATA

    // Tally contoller.
    SP_Tallier b_tallier;

  public:
    // Virtual destructor for polymorphism.
    virtual ~Solver() = 0;

    //! Solve the problem.
    virtual void solve() = 0;

    //! Call to reset the solver and tallies for another calculation.
    virtual void reset() = 0;

    //! Get tallies.
    SP_Tallier tallier() const { return b_tallier; }
};

} // end namespace profugus

#endif // mc_Solver_hh

//---------------------------------------------------------------------------//
//                 end of Solver.hh
//---------------------------------------------------------------------------//
