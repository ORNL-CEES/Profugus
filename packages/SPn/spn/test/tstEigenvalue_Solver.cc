//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   SPn/spn/test/tstEigenvalue_Solver.cc
 * \author Thomas M. Evans
 * \date   Mon Mar 10 23:18:53 2014
 * \brief  Unit-test for Eigenvalue_Solver.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "gtest/utils_gtest.hh"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_StandardCatchMacros.hpp"

#include "xs/Mat_DB.hh"
#include "mesh/Partitioner.hh"
#include "solvers/LinAlgTypedefs.hh"
#include "../Dimensions.hh"
#include "../Eigenvalue_Solver.hh"
#include "../Isotropic_Source.hh"

#include "Test_XS.hh"

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//

template <class T>
class Inf_Med_Eigenvalue_SolverTest : public testing::Test
{
  protected:
    // Typedefs usable inside the test fixture

    typedef profugus::Eigenvalue_Solver<T>            Solver;
    typedef Teuchos::RCP<Solver>                      RCP_Solver;
    typedef typename Solver::RCP_ParameterList        RCP_ParameterList;
    typedef typename Solver::Linear_System_t          Linear_System_t;
    typedef typename Linear_System_t::RCP_Mat_DB      RCP_Mat_DB;
    typedef typename Linear_System_t::RCP_Dimensions  RCP_Dimensions;
    typedef typename Linear_System_t::RCP_Mesh        RCP_Mesh;
    typedef typename Linear_System_t::RCP_Indexer     RCP_Indexer;
    typedef typename Linear_System_t::RCP_Global_Data RCP_Global_Data;
    typedef typename Linear_System_t::Matrix_t        Matrix_t;
    typedef profugus::Partitioner                     Partitioner;
    typedef Teuchos::RCP<profugus::State>             RCP_State;
    typedef profugus::Mat_DB                          Mat_DB_t;
    typedef Mat_DB_t::XS_t                            XS;
    typedef Mat_DB_t::RCP_XS                          RCP_XS;

  protected:
    // Initialization that are performed for each test
    void SetUp()
    {
        node  = profugus::node();
        nodes = profugus::nodes();

        // database
        db = Teuchos::rcp(new Teuchos::ParameterList("Test"));
    }

    void build(int order,
               int Ng)
    {
        num_groups = Ng;
        eqn_order  = order;

        // build 4x4x4 mesh

        db->set("delta_x", 1.0);
        db->set("delta_y", 1.0);
        db->set("delta_z", 1.0);

        db->set("num_cells_i", 3);
        db->set("num_cells_j", 3);
        db->set("num_cells_k", 3);

        if (nodes == 2)
        {
            db->set("num_blocks_i", 2);
        }
        if (nodes == 4)
        {
            db->set("num_blocks_i", 2);
            db->set("num_blocks_j", 2);
        }

        Partitioner p(db);
        p.build();

        mesh    = p.get_mesh();
        indexer = p.get_indexer();
        data    = p.get_global_data();

        solver  = Teuchos::rcp(new Solver(db));

        if (num_groups == 1)
            mat = one_grp::make_mat(3, mesh->num_cells());
        else
            mat = three_grp::make_mat(3, mesh->num_cells());

        EXPECT_FALSE(mat.is_null());
        EXPECT_FALSE(mesh.is_null());
        EXPECT_FALSE(indexer.is_null());
        EXPECT_FALSE(data.is_null());

        // add fission
        RCP_Mat_DB matf = Teuchos::rcp(new Mat_DB_t);
        {
            const XS &old = mat->xs();
            RCP_XS xs     = Teuchos::rcp(new XS);

            xs->set(old.pn_order(), old.num_groups());

            XS::OneDArray tot(num_groups), nusigf(num_groups), chi(num_groups);
            XS::TwoDArray P0(num_groups, num_groups),
                P1(num_groups, num_groups),
                P2(num_groups, num_groups),
                P3(num_groups, num_groups);

            const XS::Vector &sig = old.vector(0, XS::TOTAL);
            for (int g = 0; g < num_groups; ++g)
            {
                tot[g] = sig[g];
            }

            const XS::Matrix &sct0 = old.matrix(0, 0);
            const XS::Matrix &sct1 = old.matrix(0, 1);
            const XS::Matrix &sct2 = old.matrix(0, 2);
            const XS::Matrix &sct3 = old.matrix(0, 3);
            for (int g = 0; g < num_groups; ++g)
            {
                for (int gp = 0; gp < num_groups; ++gp)
                {
                    P0(g, gp) = sct0(g, gp);
                    P1(g, gp) = sct1(g, gp);
                    P2(g, gp) = sct2(g, gp);
                    P3(g, gp) = sct3(g, gp);
                }
            }

            if (num_groups == 1)
            {
                nusigf[0] = 0.5;
                chi[0]    = 1.0;
            }
            else if (num_groups == 3)
            {
                nusigf[0] = 0.3;
                nusigf[1] = 0.2;

                chi[0] = 0.2;
                chi[1] = 0.8;
            }

            xs->add(0, XS::TOTAL, tot);
            xs->add(0, XS::NU_SIG_F, nusigf);
            xs->add(0, XS::CHI, chi);

            xs->add(0, 0, P0);
            xs->add(0, 1, P1);
            xs->add(0, 2, P2);
            xs->add(0, 3, P3);

            xs->complete();
            matf->set(xs, mesh->num_cells());

            for (int n = 0; n < mesh->num_cells(); ++n)
            {
                matf->matid(n) = 0;
            }
        }

        bool success = true, verbose = true;
        try
        {
            dim = Teuchos::rcp(new profugus::Dimensions(eqn_order));
            solver->setup(dim, matf, mesh, indexer, data);
        }

        TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

        // make a state object
        state = Teuchos::rcp(new profugus::State(mesh, num_groups));
    }

  protected:
    RCP_ParameterList db;
    RCP_Mesh          mesh;
    RCP_Indexer       indexer;
    RCP_Global_Data   data;

    RCP_Mat_DB     mat;
    RCP_Dimensions dim;

    RCP_Solver solver;

    RCP_State state;

    int num_groups, eqn_order;

    int node, nodes;
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//
using profugus::EpetraTypes;
using profugus::TpetraTypes;
typedef ::testing::Types<EpetraTypes,TpetraTypes> MyTypes;
TYPED_TEST_CASE(Inf_Med_Eigenvalue_SolverTest, MyTypes);

TYPED_TEST(Inf_Med_Eigenvalue_SolverTest, 1Grp_SP1)
{
    typedef typename TestFixture::RCP_ParameterList RCP_ParameterList;

    // make an Anasazi database and set modified Gram-Schmidt
    // orthogonalization
    RCP_ParameterList edb = Teuchos::rcp(
        new Teuchos::ParameterList("eigenvalue_settings"));
    RCP_ParameterList adb = Teuchos::rcp(
        new Teuchos::ParameterList("anasazi_settings"));
    adb->set("eigensolver", std::string("Arnoldi"));
    edb->set("Anasazi", *adb);
    this->db->set("eigenvalue_db", *edb);

    this->build(1, 1);
    EXPECT_EQ(1, this->dim->num_equations());

    Teuchos::RCP<profugus::Isotropic_Source> q;

    this->solver->solve(q);

    EXPECT_SOFTEQ(0.8333333333333334, this->solver->get_eigenvalue(), 1.0e-6);
}

//---------------------------------------------------------------------------//

TYPED_TEST(Inf_Med_Eigenvalue_SolverTest, 1Grp_SP3)
{
    this->build(3, 1);
    EXPECT_EQ(2, this->dim->num_equations());

    Teuchos::RCP<profugus::Isotropic_Source> q;
    this->solver->solve(q);

    EXPECT_SOFTEQ(0.8333333333333334, this->solver->get_eigenvalue(), 1.0e-6);
}

//---------------------------------------------------------------------------//

TYPED_TEST(Inf_Med_Eigenvalue_SolverTest, 3Grp_SP1)
{
    typedef typename TestFixture::Linear_System_t Linear_System_t;
    typedef typename TypeParam::MV                MV;

    this->build(1, 3);
    EXPECT_EQ(1, this->dim->num_equations());

    Teuchos::RCP<profugus::Isotropic_Source> q;
    this->solver->solve(q);

    EXPECT_SOFTEQ(3.301149153942720, this->solver->get_eigenvalue(), 1.0e-6);

    Teuchos::RCP<const MV> ev = this->solver->get_eigenvector();
    Teuchos::ArrayRCP<const double> ev_data =
        profugus::VectorTraits<TypeParam>::get_data(ev);
    EXPECT_EQ(this->mesh->num_cells() * 3, ev_data.size());

    const Linear_System_t &system = this->solver->get_linear_system();

    double ref[] = {0.316914305060293, 0.867219274759711, 0.384052148455639};

    // normalization per cell
    for (int cell = 0; cell < this->mesh->num_cells(); ++cell)
    {
        double norm = 0.0;
        for (int g = 0; g < 3; ++g)
        {
            int index  = system.index(g, 0, cell);
            norm      += ev_data[index] * ev_data[index];
        }

        for (int g = 0; g < 3; ++g)
        {
            int index = system.index(g, 0, cell);
            // Take absolute value here because eigenvector can be negative
            // It would get normalized during write_state but here we're
            //  accessing the raw eigenvector.
            double v  = std::fabs( ev_data[index] / sqrt(norm) );
            EXPECT_SOFTEQ(ref[g], v, 1.0e-6);
        }
    }
}

//---------------------------------------------------------------------------//

TYPED_TEST(Inf_Med_Eigenvalue_SolverTest, 3Grp_SP3)
{
    typedef typename TestFixture::Linear_System_t Linear_System_t;
    typedef typename TypeParam::MV                MV;

    this->build(3, 3);
    EXPECT_EQ(2, this->dim->num_equations());

    Teuchos::RCP<profugus::Isotropic_Source> q;
    this->solver->solve(q);

    EXPECT_SOFTEQ(3.301149153942720, this->solver->get_eigenvalue(), 1.0e-6);

    Teuchos::RCP<const MV> ev = this->solver->get_eigenvector();
    EXPECT_EQ(this->mesh->num_cells() * 3 * 2,
              profugus::VectorTraits<TypeParam>::local_length(ev));

    const Linear_System_t &system = this->solver->get_linear_system();

    double ref[] = {0.316914305060293, 0.867219274759711, 0.384052148455639};

    // fill the state
    EXPECT_EQ(this->mesh->num_cells() * 3, this->state->size());

    this->solver->write_state(*this->state);
    profugus::State::View_Field phi = this->state->flux();
    EXPECT_EQ(this->mesh->num_cells() * 3, phi.size());

    // normalize the eigenvector on each cell
    double eps = 1.0e-4;
    for (int cell = 0; cell < this->mesh->num_cells(); ++cell)
    {
        int g0 = cell + 0 * this->mesh->num_cells();
        int g1 = cell + 1 * this->mesh->num_cells();
        int g2 = cell + 2 * this->mesh->num_cells();

        double norm = phi[g0]*phi[g0] + phi[g1]*phi[g1] + phi[g2]*phi[g2];
        norm        = 1.0 / sqrt(norm);

        EXPECT_SOFTEQ(ref[0], phi[g0]*norm, eps);
        EXPECT_SOFTEQ(ref[1], phi[g1]*norm, eps);
        EXPECT_SOFTEQ(ref[2], phi[g2]*norm, eps);
    }
}

//---------------------------------------------------------------------------//
//                 end of tstEigenvalue_Solver.cc
//---------------------------------------------------------------------------//
