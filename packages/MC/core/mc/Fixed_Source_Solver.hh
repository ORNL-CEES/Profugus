//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/Fixed_Source_Solver.hh
 * \author Thomas M. Evans
 * \date   Tue May 13 14:40:06 2014
 * \brief  Fixed_Source_Solver class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef mc_Fixed_Source_Solver_hh
#define mc_Fixed_Source_Solver_hh

#include "Solver.hh"

#include "harness/DBC.hh"
#include "utils/Definitions.hh"
#include "Source_Transporter.hh"

namespace profugus
{

//===========================================================================//
/*!
 * \class Fixed_Source_Solver
 * \brief Solve a fixed-source problem using Monte Carlo.
 */
/*!
 * \example mc/test/tstFixed_Source_Solver.cc
 *
 * Test of Fixed_Source_Solver.
 */
//===========================================================================//

class Fixed_Source_Solver : public Solver
{
  public:
    //@{
    //! Typedefs.
    typedef Source_Transporter                    Source_Transporter_t;
    typedef std::shared_ptr<Source_Transporter_t> SP_Source_Transporter;
    typedef Source_Transporter_t::SP_Source       SP_Source;
    //@}

  private:
    // >>> DATA

    // Source transporter.
    SP_Source_Transporter d_transporter;

    // Source
    SP_Source d_source;

  public:
    // Constructor.
    Fixed_Source_Solver();

    // Set the underlying source transporter and source.
    void set(SP_Source_Transporter transporter, SP_Source source);

    // >>> INHERITED INTERFACE

    // Solve the fixed-source problem.
    void solve();

    // Reset the problem.
    void reset() { NOT_IMPLEMENTED("Resetting a fixed source calculation."); }

  private:
    // >>> IMPLEMENTATION

    // Total, global number of particles to run.
    def::size_type d_Np;

    // Processor/node info.
    int d_node, d_nodes;
};

} // end namespace profugus

#endif // mc_Fixed_Source_Solver_hh

//---------------------------------------------------------------------------//
//                 end of Fixed_Source_Solver.hh
//---------------------------------------------------------------------------//
