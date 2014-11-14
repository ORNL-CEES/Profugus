//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/Fission_Matrix_Acceleration.hh
 * \author Thomas M. Evans
 * \date   Wed Nov 12 14:54:34 2014
 * \brief  Fission_Matrix_Acceleration class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef mc_Fission_Matrix_Acceleration_hh
#define mc_Fission_Matrix_Acceleration_hh

#include "Teuchos_RCP.hpp"

#include "spn/State.hh"
#include "spn/Linear_System.hh"
#include "spn_driver/Problem_Builder.hh"

namespace profugus
{

//===========================================================================//
/*!
 * \class Fission_Matrix_Acceleration
 * \brief Do Fission Matrix acceleration of fission source.
 */
/*!
 * \example mc/test/tstFission_Matrix_Acceleration.cc
 *
 * Test of Fission_Matrix_Acceleration.
 */
//===========================================================================//

class Fission_Matrix_Acceleration
{
  public:
    // Typedefs.
    typedef spn::Problem_Builder                 Problem_Builder_t;
    typedef Problem_Builder_t::RCP_ParameterList RCP_ParameterList;
    typedef Problem_Builder_t::RCP_Mesh          RCP_Mesh;
    typedef Problem_Builder_t::RCP_Indexer       RCP_Indexer;
    typedef Problem_Builder_t::RCP_Global_Data   RCP_Global_Data;
    typedef Problem_Builder_t::RCP_Mat_DB        RCP_Mat_DB;
    typedef Teuchos::RCP<State>                  RCP_State;

  protected:
    // >>> DATA

    // Problem database.
    RCP_ParameterList b_db;

    // Mesh objects
    RCP_Mesh        b_mesh;
    RCP_Indexer     b_indexer;
    RCP_Global_Data b_gdata;

    // Material database.
    RCP_Mat_DB b_mat;

    // Forward and adjoint solutions to the SPN problem.
    RCP_State b_forward;
    RCP_State b_adjoint;

  public:
    // Constructor.
    Fission_Matrix_Acceleration() { /* ... */ }

    // Destructor.
    virtual ~Fission_Matrix_Acceleration() { /* ... */ }

    //! Build the SPN problem.
    virtual void build_problem(const Problem_Builder_t &builder) = 0;

    //! Initialize the acceleration.
    virtual void initialize() = 0;

    // >>> ACCESSORS

    //@{
    //! Get the forward and adjoint states.
    const State& forward() const { return *b_forward; }
    const State& adjoint() const { return *b_adjoint; }
    //@}
};

//===========================================================================//
/*!
 * \class Fission_Matrix_Acceleration_Impl
 * \brief Do implementation of fission matrix acceleration of fission source.
 */
//===========================================================================//

template<class T>
class Fission_Matrix_Acceleration_Impl : public Fission_Matrix_Acceleration
{
  public:
    // Typedefs.
    typedef profugus::Linear_System<T>               Linear_System_t;
    typedef typename Linear_System_t::RCP_Dimensions RCP_Dimensions;
    typedef Teuchos::RCP<Linear_System_t>            RCP_Linear_System;

  private:
    // >>> DATA

    // Problem dimensions.
    RCP_Dimensions d_dim;

    // Original SPN linear system.
    RCP_Linear_System d_system;

  public:
    // Constructor.
    Fission_Matrix_Acceleration_Impl();

    // Build the SPN problem.
    void build_problem(const Problem_Builder_t &builder);

    // Initialize the acceleration.
    void initialize();
};

} // end namespace profugus

#endif // mc_Fission_Matrix_Acceleration_hh

//---------------------------------------------------------------------------//
//                 end of Fission_Matrix_Acceleration.hh
//---------------------------------------------------------------------------//
