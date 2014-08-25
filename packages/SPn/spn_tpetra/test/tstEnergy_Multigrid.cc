//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   spn_tpetra/test/tstEnergy_Multigrid.cc
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

#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Vector.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Array.hpp"

#include "xs/Mat_DB.hh"
#include "mesh/Partitioner.hh"
#include "spn/Dimensions.hh"
#include "../Linear_System_FV.hh"
#include "../Energy_Multigrid.hh"

#include "Test_XS.hh"

using Teuchos::RCP;
using Teuchos::rcp;

typedef profugus::tpetra::Energy_Multigrid  Energy_Multigrid;
typedef Energy_Multigrid::ParameterList     ParameterList;
typedef Energy_Multigrid::RCP_ParameterList RCP_ParameterList;
typedef profugus::Partitioner               Partitioner;
typedef KokkosClassic::SerialNode           Node;
typedef Tpetra::Operator<double,int,int,Node> Tpetra_Op;
typedef Tpetra::MultiVector<double,int,int,Node> Tpetra_MV;
typedef Tpetra::Vector<double,int,int,Node> Tpetra_Vector;

using namespace std;

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//

TEST(MultigridTest, Heuristic)
{
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
    RCP<profugus::tpetra::Linear_System> system = rcp(
        new profugus::tpetra::Linear_System_FV(db, dim, mat, mesh, indexer, data) );
    system->build_Matrix();
    RCP<Tpetra_Op> matrix = system->get_Operator();

    // Create db for preconditioner
    RCP_ParameterList block_gmres_db =
        rcp(new ParameterList("Block GMRES"));
    block_gmres_db->set("Convergence Tolerance",1e-4);
    block_gmres_db->set("Maximum Iterations",5);

    RCP_ParameterList belos_solver_types_db =
        rcp(new ParameterList("Solver Types"));
    belos_solver_types_db->set("Block GMRES",*block_gmres_db);

    RCP_ParameterList belos_db =
        rcp(new ParameterList("Belos"));
    belos_db->set("Solver Type",string("Block GMRES"));
    belos_db->set("Solver Types",*belos_solver_types_db);

    RCP_ParameterList solver_types_db =
        rcp(new ParameterList("Linear Solver Types"));
    solver_types_db->set("Belos", *belos_db);

    RCP_ParameterList stratimikos_db =
        rcp(new ParameterList("Stratimikos"));
    stratimikos_db->set("Linear Solver Types", *solver_types_db);
    stratimikos_db->set("Linear Solver Type", string("Belos"));
    stratimikos_db->set("Preconditioner Type",string("None"));

    RCP_ParameterList smoother_db =
        rcp(new ParameterList("Smoother"));
    smoother_db->set("Stratimikos", *stratimikos_db);
    smoother_db->set("max_itr", 5);
    smoother_db->set("solver_type", string("stratimikos"));

    RCP_ParameterList prec_db =
        rcp(new ParameterList("Prec"));
    prec_db->set("Smoother", *smoother_db);

    // Create two vectors
    RCP<Tpetra_Vector> tmp_vec = system->get_RHS();
    RCP<Tpetra_MV> x(
            Teuchos::rcp( new Tpetra_MV(*tmp_vec,Teuchos::Copy)));
    RCP<Tpetra_MV> y(
            Teuchos::rcp( new Tpetra_MV(*tmp_vec,Teuchos::Copy)));

    // Create preconditioner
    Energy_Multigrid prec(db, prec_db, dim, mat, mesh, indexer, data, system);

    x->putScalar(1.0);
    prec.apply(*x, *y);

    Teuchos::ArrayRCP<double> norm2(1);
    y->norm2(norm2());
    {
        cout << "Norm of y after apply: " << setw(12) << scientific
             << setprecision(3) << norm2[0] << endl;
    }

    // Heuristic test of output vector norm
    if (nodes == 1)
    {
        EXPECT_SOFTEQ(2.66046e2, norm2[0], 1.0e-3);
    }
    else if (nodes == 2)
    {
        EXPECT_SOFTEQ(2.899e+02, norm2[0], 1.0e-3);
    }
    else if (nodes == 4)
    {
        EXPECT_SOFTEQ(2.758546e+02, norm2[0], 1.0e-3);
    }
}

//---------------------------------------------------------------------------//
//                 end of tstEnergy_Multigrid.cc
//---------------------------------------------------------------------------//
