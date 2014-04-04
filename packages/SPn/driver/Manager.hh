//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   driver/Manager.hh
 * \author Thomas M. Evans
 * \date   Fri Mar 14 11:32:36 2014
 * \brief  Manager class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef driver_Manager_hh
#define driver_Manager_hh

#include <sstream>
#include <string>

#include "Teuchos_RCP.hpp"

#include "comm/P_Stream.hh"
#include "spn/Eigenvalue_Solver.hh"
#include "spn/Fixed_Source_Solver.hh"
#include "spn/Time_Dependent_Solver.hh"
#include "Problem_Builder.hh"

namespace profugus
{

//===========================================================================//
/*!
 * \class Manager
 * \brief Manager class that drives the SPN miniapp.
 */
//===========================================================================//

class Manager
{
  private:
    // Typedefs.
    typedef Problem_Builder::RCP_ParameterList   RCP_ParameterList;
    typedef Problem_Builder::RCP_Mesh            RCP_Mesh;
    typedef Problem_Builder::RCP_Indexer         RCP_Indexer;
    typedef Problem_Builder::RCP_Global_Data     RCP_Global_Data;
    typedef Problem_Builder::RCP_Mat_DB          RCP_Mat_DB;
    typedef Solver_Base::RCP_Dimensions          RCP_Dimensions;
    typedef Solver_Base::State_t                 State_t;
    typedef Teuchos::RCP<State_t>                RCP_State;
    typedef Teuchos::RCP<Solver_Base>            RCP_Solver_Base;
    typedef Teuchos::RCP<Fixed_Source_Solver>    RCP_Fixed_Source_Solver;
    typedef Teuchos::RCP<Eigenvalue_Solver>      RCP_Eigenvalue_Solver;
    typedef Teuchos::RCP<Time_Dependent_Solver>  RCP_Time_Dependent_Solver;
    typedef Fixed_Source_Solver::External_Source External_Source_t;
    typedef Teuchos::RCP<External_Source_t>      RCP_External_Source;

    // >>> DATA

    // Problem database.
    RCP_ParameterList d_db;

    // Mesh objects
    RCP_Mesh        d_mesh;
    RCP_Indexer     d_indexer;
    RCP_Global_Data d_gdata;

    // Material database.
    RCP_Mat_DB d_mat;

    // Problem state.
    RCP_State d_state;

    // Problem dimensions.
    RCP_Dimensions d_dim;

    // Solvers
    RCP_Solver_Base           d_solver_base;
    RCP_Fixed_Source_Solver   d_fixed_solver;
    RCP_Eigenvalue_Solver     d_eigen_solver;
    RCP_Time_Dependent_Solver d_time_dep_solver;

    // External source for fixed source problems.
    RCP_External_Source d_external_source;

  public:
    // Constructor.
    Manager();

    // Setup the problem.
    void setup(const std::string &xml_file);

    // Solve the problem.
    void solve();

    // Output.
    void output();

  private:
    // >>> IMPLEMENTATION

    // Processor.
    int d_node, d_nodes;

    // Problem name.
    std::string d_problem_name;

    //! Output messages in a common format.
#define SCREEN_MSG(stream)                            \
    {                                                 \
        std::ostringstream m;                         \
        m << ">>> " << stream;                        \
        profugus::pcout << m.str() << profugus::endl; \
    }
};

} // end namespace profugus

#endif // driver_Manager_hh

//---------------------------------------------------------------------------//
//                 end of Manager.hh
//---------------------------------------------------------------------------//
