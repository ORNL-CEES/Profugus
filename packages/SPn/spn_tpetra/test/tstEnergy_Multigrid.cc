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
#include "Teuchos_XMLParameterListHelpers.hpp"

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

    /*
    const std::string plrefstr(
        );
        */

    const std::string plrefstr(
 "<ParameterList name='Multigrid Preconditioner'>                            \n\
   <ParameterList name='Smoother'>                                           \n\
    <Parameter name='Preconditioner' type='string' value='Ifpack2'/>         \n\
    <Parameter name='Ifpack2_Type' type='string' value='RILUK'/>             \n\
    <Parameter name='verbosity' type='string' value='None'/>                 \n\
    <Parameter name='max_itr' type='int' value='1'/>                         \n\
    <Parameter name='solver_type' type='string' value='stratimikos'/>        \n\
    <ParameterList name='Stratimikos'>                                       \n\
     <Parameter name='Linear Solver Type' type='string' value='Belos'/>      \n\
     <Parameter name='Preconditioner Type' type='string' value='None'/>      \n\
     <ParameterList name='Linear Solver Types'>                              \n\
      <ParameterList name='Belos'>                                           \n\
       <Parameter name='Solver Type' type='string' value='Block GMRES'/>     \n\
       <ParameterList name='Solver Types'>                                   \n\
        <ParameterList name='Block GMRES'>                                   \n\
         <Parameter name='Convergence Tolerance' type='double' value='1e-4'/>\n\
        </ParameterList>                                                     \n\
       </ParameterList>                                                      \n\
      </ParameterList>                                                       \n\
     </ParameterList>                                                        \n\
    </ParameterList>                                                         \n\
   </ParameterList>                                                          \n\
  </ParameterList>                                                           \n"
        );

    // Convert string to a Teuchos PL
    RCP_ParameterList prec_pl =
        Teuchos::getParametersFromXmlString(plrefstr);

    // Create two vectors
    RCP<Tpetra_Vector> tmp_vec = system->get_RHS();
    RCP<Tpetra_MV> x( new Tpetra_MV(tmp_vec->getMap(),1) );
    RCP<Tpetra_MV> y( new Tpetra_MV(tmp_vec->getMap(),1) );

    // Create preconditioner
    Energy_Multigrid prec(db, prec_pl, dim, mat, mesh, indexer, data, system);

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
        EXPECT_SOFTEQ(329.171, norm2[0], 1.0e-3);
    }
    else if (nodes == 2)
    {
        EXPECT_SOFTEQ(299.947, norm2[0], 1.0e-3);
    }
    else if (nodes == 4)
    {
        EXPECT_SOFTEQ(272.649, norm2[0], 1.0e-3);
    }
}

//---------------------------------------------------------------------------//
//                 end of tstEnergy_Multigrid.cc
//---------------------------------------------------------------------------//
