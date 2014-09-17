//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   spn/test/tstEnergy_Multigrid.cc
 * \author Thomas Evans, Steven Hamilton
 * \date   Monday March 10 12:38:40 2014
 * \brief  Energy Grid Transfer test.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <vector>
#include <string>
#include <iomanip>

#include "gtest/utils_gtest.hh"
#include <SPn/config.h>

#include "Teuchos_RCP.hpp"
#include "Teuchos_Array.hpp"
#include "AnasaziOperatorTraits.hpp"

#include "xs/Mat_DB.hh"
#include "mesh/Partitioner.hh"
#include "solvers/LinAlgTypedefs.hh"
#include "../Linear_System_FV.hh"
#include "../Dimensions.hh"
#include "../Energy_Multigrid.hh"
#include "../VectorTraits.hh"

#include "Test_XS.hh"

using Teuchos::RCP;
using Teuchos::rcp;
using profugus::EpetraTypes;

typedef Teuchos::ParameterList      ParameterList;
typedef Teuchos::RCP<ParameterList> RCP_ParameterList;
typedef profugus::Partitioner       Partitioner;

using namespace std;

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//

TEST(MultigridTest, Heuristic)
{
    typedef EpetraTypes T;
    typedef profugus::Energy_Multigrid<T> Energy_Multigrid;
    typedef typename T::MV MV;
    typedef typename T::OP OP;
    typedef Anasazi::OperatorTraits<double,MV,OP> OPT;

    int nodes = profugus::nodes();
    int node  = profugus::node();

    // Initialize database and set basic data
    RCP_ParameterList db = rcp(new ParameterList("Main"));
    db->set("delta_x", 1.0);
    db->set("delta_y", 1.0);
    db->set("delta_z", 10.0);
    db->set("num_cells_i", 4);
    db->set("num_cells_j", 4);
    db->set("num_cells_k", 4);

    if (nodes == 1)
    {
        db->set("num_blocks_i", 1);
        db->set("num_blocks_j", 1);
    }
    else if (nodes == 2)
    {
        db->set("num_blocks_i", 2);
        db->set("num_blocks_j", 1);
    }
    else if (nodes == 4)
    {
        db->set("num_blocks_i", 2);
        db->set("num_blocks_j", 2);
    }

    // Build mesh objects
    Partitioner p(db);
    p.build();
    Partitioner::RCP_Mesh mesh        = p.get_mesh();
    Partitioner::RCP_Indexer indexer  = p.get_indexer();
    Partitioner::RCP_Global_Data data = p.get_global_data();

    // Get mat_db
    int pn_order = 1;
    RCP<profugus::Mat_DB> mat =
        twelve_grp::make_mat(pn_order, mesh->num_cells());

    // Build SPN Dimensions
    int spn_order=3;
    RCP<profugus::Dimensions> dim = rcp( new profugus::Dimensions(spn_order) );

    // Set boundary conditions
    db->set("boundary", string("reflect"));
    Teuchos::Array<int> ref(6, 0);
    ref[0] = 1;
    ref[1] = 1;
    ref[2] = 1;
    ref[3] = 1;

    // Fine level linear system
    RCP<profugus::Linear_System<T> > system = rcp(
        new profugus::Linear_System_FV<T>(
            db, dim, mat, mesh, indexer, data) );
    system->build_Matrix();
    RCP<typename T::OP> matrix = system->get_Operator();

    // Create db for preconditioner
    RCP_ParameterList aztec_settings_db =
        rcp(new ParameterList("AztecOO Settings"));
    aztec_settings_db->set("Aztec Solver", string("GMRES"));
    aztec_settings_db->set("Aztec Preconditioner", string("Jacobi"));

    RCP_ParameterList forward_solve_db =
        rcp(new ParameterList("Forward Solve"));
    forward_solve_db->set("AztecOO Settings", *aztec_settings_db);

    RCP_ParameterList aztecoo_db =
        rcp(new ParameterList("AztecOO"));
    aztecoo_db->set("Forward Solve", *forward_solve_db);

    RCP_ParameterList solver_types_db =
        rcp(new ParameterList("Linear Solver Types"));
    solver_types_db->set("AztecOO", *aztecoo_db);

    RCP_ParameterList stratimikos_db =
        rcp(new ParameterList("Stratimikos"));
    stratimikos_db->set("Linear Solver Types", *solver_types_db);
    stratimikos_db->set("Linear Solver Type", string("AztecOO"));

    RCP_ParameterList smoother_db =
        rcp(new ParameterList("Smoother"));
    smoother_db->set("Stratimikos", *stratimikos_db);
    smoother_db->set("max_itr", 1);
    smoother_db->set("solver_type", string("stratimikos"));

    RCP_ParameterList prec_db =
        rcp(new ParameterList("Prec"));
    prec_db->set("Smoother", *smoother_db);

    // Create two vectors
    RCP<typename T::VECTOR> tmp_vec = system->get_RHS();
    RCP<typename T::MV> x =
        profugus::VectorTraits<T>::build_vector(system->get_Map());
    RCP<typename T::MV> y =
        profugus::VectorTraits<T>::build_vector(system->get_Map());

    // Create preconditioner
    Energy_Multigrid prec(db, prec_db, dim, mat, mesh, indexer, data, system);

    profugus::VectorTraits<T>::put_scalar(x,1.0);
    OPT::Apply(prec,*x,*y);

    vector<double> norm2(1);
    y->Norm2(&norm2[0]);
    {
        cout << "Norm of y after apply: " << setw(12) << scientific
             << setprecision(3) << norm2[0] << endl;
    }

    if (nodes == 1)
    {
        EXPECT_SOFTEQ(3.292e+02, norm2[0], 1.0e-3);
    }
    else if (nodes == 2)
    {
        EXPECT_SOFTEQ(2.999e+02, norm2[0], 1.0e-3);
    }
    else if (nodes == 4)
    {
        EXPECT_SOFTEQ(2.726e+02, norm2[0], 1.0e-3);
    }
}

//---------------------------------------------------------------------------//
//                 end of tstEnergy_Multigrid.cc
//---------------------------------------------------------------------------//
