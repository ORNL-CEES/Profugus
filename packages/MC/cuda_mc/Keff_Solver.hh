//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   CUDA_MC/cuda_mc/Keff_Solver.hh
 * \author Thomas M. Evans
 * \date   Mon Apr 13 18:14:03 2015
 * \brief  Keff_Solver class declaration.
 * \note   Copyright (c) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef CUDA_MC_cuda_mc_Keff_Solver_hh
#define CUDA_MC_cuda_mc_Keff_Solver_hh

#include "Keff_Tally.hh"
#include "Fission_Source.hh"
#include "Solver.hh"

namespace cuda_profugus
{

//===========================================================================//
/*!
 * \class Keff_Solver
 * \brief Base class for eigenvalue solvers.
 */
/*!
 * \example cuda_mc/test/tstKeff_Solver.cc
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
};

} // end namespace cuda_profugus

#endif // CUDA_MC_cuda_mc_Keff_Solver_hh

//---------------------------------------------------------------------------//
// end of CUDA_MC/cuda_mc/Keff_Solver.hh
//---------------------------------------------------------------------------//
