//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   MC/mc/Solver.hh
 * \author Thomas M. Evans and Seth Johnson
 * \date   Tue May 13 14:39:00 2014
 * \brief  Solver class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef MC_mc_Solver_hh
#define MC_mc_Solver_hh

#include <memory>

#include "Tallier.hh"
#include "Physics.hh"

namespace profugus
{

//===========================================================================//
/*!
 * \class Solver
 * \brief Base class for Monte Carlo top-level solvers.
 */
//===========================================================================//

template <class Geometry>
class Solver
{
  public:
    //@{
    //! Typedefs.
    typedef Geometry                   Geometry_t;
    typedef Physics<Geometry_t>        Physics_t;
    typedef Tallier<Geometry_t>        Tallier_t;
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

#endif // MC_mc_Solver_hh

//---------------------------------------------------------------------------//
//                 end of Solver.hh
//---------------------------------------------------------------------------//
