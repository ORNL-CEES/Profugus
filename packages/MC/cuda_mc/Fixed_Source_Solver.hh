//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/Fixed_Source_Solver.hh
 * \author Thomas M. Evans
 * \date   Tue May 13 14:40:06 2014
 * \brief  Fixed_Source_Solver class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_mc_Fixed_Source_Solver_hh
#define cuda_mc_Fixed_Source_Solver_hh

#include "Solver.hh"

#include "harness/DBC.hh"
#include "utils/Definitions.hh"
#include "Source_Transporter.hh"

namespace cuda_profugus
{

//===========================================================================//
/*!
 * \class Fixed_Source_Solver
 * \brief Solve a fixed-source problem using Monte Carlo.
 */
/*!
 * \example cuda_mc/test/tstFixed_Source_Solver.cc
 *
 * Test of Fixed_Source_Solver.
 */
//===========================================================================//

template <class Geometry>
class Fixed_Source_Solver : public Solver<Geometry>
{
    typedef Solver<Geometry> Base;
  public:
    //@{
    //! Typedefs.
    typedef Geometry                                    Geometry_t;
    typedef Source_Transporter<Geometry_t>              Source_Transporter_t;
    typedef std::shared_ptr<Source_Transporter_t>       SP_Source_Transporter;
    typedef typename Source_Transporter_t::SP_Source    SP_Source;
    //@}

  private:
    // >>> DATA

    using Base::b_tallier;

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
    void reset() { /* ... */ }

  private:
    // >>> IMPLEMENTATION

    // Total, global number of particles to run.
    def::size_type d_Np;

    // Processor/node info.
    int d_node, d_nodes;
};

} // end namespace cuda_profugus

#endif // cuda_mc_Fixed_Source_Solver_hh

//---------------------------------------------------------------------------//
//                 end of Fixed_Source_Solver.hh
//---------------------------------------------------------------------------//
