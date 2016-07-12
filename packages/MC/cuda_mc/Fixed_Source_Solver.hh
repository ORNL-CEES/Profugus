//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/Fixed_Source_Solver.hh
 * \author Steven Hamilton
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
#include "Source.cuh"
#include "Solver.hh"

namespace cuda_mc
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

template <class Geometry>
class Fixed_Source_Solver : public Solver<Geometry>
{
    typedef Solver<Geometry> Base;

  public:
    //@{
    //! Typedefs.
    typedef Geometry                                Geometry_t;
    typedef Source_Transporter<Geometry_t>          Source_Transporter_t;
    typedef std::shared_ptr<Source_Transporter_t>   SP_Source_Transporter;
    typedef std::shared_ptr<Source<Geometry>>       SP_Source;
    typedef Tallier<Geometry>                       Tallier_t;
    typedef std::shared_ptr<Tallier_t>              SP_Tallier;
    typedef Teuchos::RCP<Teuchos::ParameterList>    RCP_Std_DB;
    typedef def::size_type                          size_type;
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
    Fixed_Source_Solver(RCP_Std_DB db);

    // Set the underlying source transporter and source.
    void set(SP_Source_Transporter transporter,
             SP_Source             source,
             SP_Tallier            tallier);

    // >>> INHERITED INTERFACE

    // Solve the fixed-source problem.
    void solve();

    // Reset the problem.
    void reset() { INSIST(false,"Cannot reset a fixed source calculation."); }

  private:
    // >>> IMPLEMENTATION

    // Total, global number of particles to run.
    def::size_type d_Np;

    // Number of particles per batch for kernel launches
    def::size_type d_batch_size;

    // Processor/node info.
    int d_node, d_nodes;
};

} // end namespace cuda_mc

#endif // cuda_mc_Fixed_Source_Solver_hh

//---------------------------------------------------------------------------//
//                 end of Fixed_Source_Solver.hh
//---------------------------------------------------------------------------//
