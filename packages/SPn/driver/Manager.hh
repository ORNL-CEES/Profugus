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

#include "spn_tpetra/Eigenvalue_Solver.hh"
#include "spn_tpetra/Fixed_Source_Solver.hh"

namespace spn
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
    typedef Problem_Builder::RCP_ParameterList     RCP_ParameterList;
    typedef Problem_Builder::RCP_Mesh              RCP_Mesh;
    typedef Problem_Builder::RCP_Indexer           RCP_Indexer;
    typedef Problem_Builder::RCP_Global_Data       RCP_Global_Data;
    typedef Problem_Builder::RCP_Mat_DB            RCP_Mat_DB;
    typedef profugus::Solver_Base                  Solver_Base_t;
    typedef Solver_Base_t::RCP_Dimensions          RCP_Dimensions;
    typedef Teuchos::RCP<profugus::State>          RCP_State;
    typedef Teuchos::RCP<Solver_Base_t>            RCP_Solver_Base;

    typedef profugus::EpetraTypes                        EpetraTypes;
    typedef profugus::Fixed_Source_Solver<EpetraTypes>   Fixed_Source_Solver_t;
    typedef profugus::Eigenvalue_Solver<EpetraTypes>     Eigenvalue_Solver_t;
    typedef profugus::Time_Dependent_Solver<EpetraTypes> Time_Dependent_Solver_t;

    typedef Teuchos::RCP<Fixed_Source_Solver_t>    RCP_Fixed_Source_Solver;
    typedef Teuchos::RCP<Eigenvalue_Solver_t>      RCP_Eigenvalue_Solver;
    typedef Teuchos::RCP<Time_Dependent_Solver_t>  RCP_Time_Dependent_Solver;
    typedef typename Fixed_Source_Solver_t::External_Source External_Source_t;
    typedef Teuchos::RCP<External_Source_t>        RCP_External_Source;

    // Tpetra variants
    typedef profugus::tpetra::Solver_Base          Solver_Base_Tpetra_t;
    typedef Teuchos::RCP<Solver_Base_Tpetra_t>     RCP_Solver_Base_Tpetra;
    typedef profugus::tpetra::Fixed_Source_Solver  Fixed_Source_Solver_Tpetra_t;
    typedef profugus::tpetra::Eigenvalue_Solver    Eigenvalue_Solver_Tpetra_t;
    typedef Teuchos::RCP<Fixed_Source_Solver_Tpetra_t>
                RCP_Fixed_Source_Solver_Tpetra;
    typedef Teuchos::RCP<Eigenvalue_Solver_Tpetra_t>
                RCP_Eigenvalue_Solver_Tpetra;

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

    // Trilinos implementation (epetra or tpetra)
    std::string d_implementation;

    // Solvers
    RCP_Solver_Base           d_solver_base;
    RCP_Fixed_Source_Solver   d_fixed_solver;
    RCP_Eigenvalue_Solver     d_eigen_solver;
    RCP_Time_Dependent_Solver d_time_dep_solver;

    // Tpetra Solvers
    RCP_Solver_Base_Tpetra           d_solver_base_tpetra;
    RCP_Fixed_Source_Solver_Tpetra   d_fixed_solver_tpetra;
    RCP_Eigenvalue_Solver_Tpetra     d_eigen_solver_tpetra;

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

} // end namespace spn

#endif // driver_Manager_hh

//---------------------------------------------------------------------------//
//                 end of Manager.hh
//---------------------------------------------------------------------------//
