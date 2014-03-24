//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   spn/Solver_Base.hh
 * \author Thomas M. Evans
 * \date   Mon Feb 17 20:36:54 2014
 * \brief  Solver_Base class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef spn_Solver_Base_hh
#define spn_Solver_Base_hh

#include "Teuchos_RCP.hpp"

#include "State.hh"
#include "Linear_System.hh"

namespace profugus
{

//===========================================================================//
/*!
 * \class Solver_Base
 * \brief Base class for SPN solvers (fixed-source and eigenvalue).
 */
/*!
 * \example spn/test/tstSolver_Base.cc
 *
 * Test of Solver_Base.
 */
//===========================================================================//

class Solver_Base
{
  public:
    //@{
    //! Typedefs.
    typedef Linear_System                      Linear_System_t;
    typedef Linear_System_t::Matrix_t          Matrix_t;
    typedef Linear_System_t::Vector_t          Vector_t;
    typedef Linear_System_t::RCP_ParameterList RCP_ParameterList;
    typedef Linear_System_t::RCP_Vector        RCP_Vector;
    typedef Linear_System_t::RCP_Dimensions    RCP_Dimensions;
    typedef Linear_System_t::RCP_Mat_DB        RCP_Mat_DB;
    typedef Linear_System_t::RCP_Mesh          RCP_Mesh;
    typedef Linear_System_t::RCP_Indexer       RCP_Indexer;
    typedef Linear_System_t::RCP_Global_Data   RCP_Global_Data;
    typedef State                              State_t;
    typedef Teuchos::RCP<Linear_System_t>      RCP_Linear_System;
    typedef typename State_t::View_Field       View_Field;
    //@}

  protected:
    // >>> DATA

    // Problem database.
    RCP_ParameterList b_db;

    // Linear system.
    RCP_Linear_System b_system;

  public:
    // Constructor.
    explicit Solver_Base(RCP_ParameterList db);

    // Virtual destructor.
    virtual ~Solver_Base() = 0;

    // Set up the solver.
    virtual void setup(RCP_Dimensions dim, RCP_Mat_DB mat, RCP_Mesh mesh,
                       RCP_Indexer indexer, RCP_Global_Data data) = 0;

    //! Write results of solve into state.
    virtual void write_state(State_t &state) = 0;

    // >>> ACCESSORS

    //! Get the linear system.
    const Linear_System_t& get_linear_system() const { return *b_system; }

    // >>> INHERITED METHODS

    // Write u-vector into the state.
    void write_u_into_state(const Vector_t &u, State_t &state);
};

} // end namespace profugus

#endif // spn_Solver_Base_hh

//---------------------------------------------------------------------------//
//                 end of Solver_Base.hh
//---------------------------------------------------------------------------//
