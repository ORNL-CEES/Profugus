//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   MC/mc/test/tstFission_Matrix_Tally.cc
 * \author Thomas M. Evans
 * \date   Thu Sep 11 16:27:58 2014
 * \brief  Fission_Matrix_Tally unit-test.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <unordered_map>
#include <algorithm>
#include <memory>

#include "../Fission_Matrix_Tally.hh"

#include "gtest/utils_gtest.hh"

#include "utils/HDF5_Reader.hh"
#include "geometry/Mesh_Geometry.hh"
#include "TransporterTestBase.hh"

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//

class Fission_Matrix_TallyTest : public TransporterTestBase
{
  protected:
    // >>> TYPEDEFS
    typedef profugus::Core                              Geometry_t;
    typedef profugus::Fission_Matrix_Tally<Geometry_t>  Tally;
    typedef Tally::PL_Tally                             PL_Tally;
    typedef Tally::Src_Tally                            Src_Tally;
    typedef profugus::Mesh_Geometry                     Mesh_Geometry_t;
    typedef std::shared_ptr<Mesh_Geometry_t>            SP_Mesh_Geometry;
    typedef Mesh_Geometry_t::Vec_Dbl                    Vec_Dbl;
    typedef Particle_t::Geo_State_t                     State;
    typedef profugus::Fission_Matrix_Processor          Processor;
    typedef Processor::Idx                              Idx;
    typedef Processor::Ordered_Graph                    Graph;
    typedef Processor::Ordered_Matrix                   Matrix;
    typedef profugus::HDF5_Reader::Decomp               Decomp;

  protected:
    void SetUp()
    {
        TransporterTestBase::SetUp();

        // make the mesh geometry
        Vec_Dbl r = {0.0, 1.0, 2.0, 3.0, 4.0};
        mesh = std::make_shared<Mesh_Geometry_t>(r, r, r);

        // make fission matrix db
        auto fis_db = Teuchos::sublist(db, "fission_matrix_db");
        tally_db    = Teuchos::sublist(fis_db, "tally");

        db->set("num_cycles", 10);
    }

  protected:
    // >>> DATA

    SP_Mesh_Geometry mesh;
    RCP_Std_DB fis_db, tally_db;
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST_F(Fission_Matrix_TallyTest, construction)
{
    Tally fm(db, physics, mesh);

    EXPECT_TRUE(static_cast<bool>(fm.get_pl_tally()));
    EXPECT_TRUE(static_cast<bool>(fm.get_src_tally()));

    EXPECT_EQ(0, tally_db->get<int>("start_cycle"));
    EXPECT_EQ(10, tally_db->get<int>("output_cycle"));
}

//---------------------------------------------------------------------------//

TEST_F(Fission_Matrix_TallyTest, transport)
{
    Tally fm(db, physics, mesh);

    // get the pathlength and source tallies from the compound tally
    auto pl  = fm.get_pl_tally();
    auto src = fm.get_src_tally();

    Particle_t p;
    p.set_wt(0.5);
    p.set_group(0);
    p.set_matid(1);

    State &state = p.geo_state();
    state.d_dir = Space_Vector(1.0, 0.0, 0.0);
    state.d_r   = Space_Vector(1.5, 1.5, 1.5);

    // step in +x
    src->birth(p);
    pl->accumulate(2.0, p);

    // step in +y
    state.d_dir = Space_Vector(0.0, 1.0, 0.0);
    src->birth(p);
    pl->accumulate(2.0, p);

    // step in -z
    state.d_dir = Space_Vector(0.0, 0.0, -1.0);
    src->birth(p);
    pl->accumulate(1.0, p);

    // step in +x
    state.d_dir = Space_Vector(1.0, 0.0, 0.0);
    state.d_r   = Space_Vector(1.5, 1.5, 0.5);
    pl->accumulate(1.0, p);

    // build the matrix
    fm.build_matrix();

    // check it
    const auto &graph  = fm.processor().graph();
    const auto &matrix = fm.processor().matrix();

    EXPECT_EQ(7, graph.size());
    EXPECT_EQ(7, matrix.size());

    Graph gref = {Idx(5,  21),
                  Idx(6,  21),
                  Idx(21, 21),
                  Idx(22, 21),
                  Idx(23, 21),
                  Idx(25, 21),
                  Idx(29, 21)};

    EXPECT_EQ(gref, graph);

    double sigf = 2.4 * 4.0;
    double wt   = 1.5;

    Matrix mref = {sigf * 0.5 * (0.5 + 0.5) / wt,
                   sigf * 0.5 * 0.5 / wt,
                   sigf * 0.5 * (0.5 + 0.5 + 0.5) / wt,
                   sigf * 0.5 * 1.0 / wt,
                   sigf * 0.5 * 0.5 / wt,
                   sigf * 0.5 * 1.0 / wt,
                   sigf * 0.5 * 0.5 / wt};

    EXPECT_VEC_SOFTEQ(mref, matrix, 1.0e-14);
}

//---------------------------------------------------------------------------//

TEST_F(Fission_Matrix_TallyTest, end_cycle)
{
    tally_db->set("start_cycle", 2);
    Tally fm(db, physics, mesh);

    // get the pathlength and source tallies from the compound tally
    auto pl  = fm.get_pl_tally();
    auto src = fm.get_src_tally();

    Particle_t p;
    p.set_wt(0.5);
    p.set_group(0);
    p.set_matid(1);

    State &state = p.geo_state();
    state.d_dir = Space_Vector(1.0, 0.0, 0.0);
    state.d_r   = Space_Vector(1.5, 1.5, 1.5);

    EXPECT_FALSE(fm.tally_started());

    // step in +x
    src->birth(p);
    pl->accumulate(2.0, p);

    // we haven't built a matrix on the first 2 cycles
    fm.build_matrix();
    {
        const auto &graph  = fm.processor().graph();
        const auto &matrix = fm.processor().matrix();

        EXPECT_TRUE(graph.empty());
        EXPECT_TRUE(matrix.empty());
    }
    fm.end_cycle(1);

    EXPECT_FALSE(fm.tally_started());

    // step in +x
    src->birth(p);
    pl->accumulate(3.0, p);

    fm.build_matrix();
    {
        const auto &graph  = fm.processor().graph();
        const auto &matrix = fm.processor().matrix();

        EXPECT_TRUE(graph.empty());
        EXPECT_TRUE(matrix.empty());
    }
    fm.end_cycle(1);

    EXPECT_FALSE(fm.tally_started());

    // we haven't tallied anything yet!
    fm.build_matrix();
    {
        const auto &graph  = fm.processor().graph();
        const auto &matrix = fm.processor().matrix();

        EXPECT_TRUE(graph.empty());
        EXPECT_TRUE(matrix.empty());
    }

    // now it will tally

    // step in +y
    state.d_dir = Space_Vector(0.0, 1.0, 0.0);
    src->birth(p);
    pl->accumulate(2.0, p);

    EXPECT_TRUE(fm.tally_started());

    // build the matrix
    fm.build_matrix();
    {
        const auto &graph  = fm.processor().graph();
        const auto &matrix = fm.processor().matrix();

        EXPECT_EQ(3, graph.size());
        EXPECT_EQ(3, matrix.size());

        Graph gref = {Idx(21, 21),
                      Idx(25, 21),
                      Idx(29, 21)};

        EXPECT_EQ(gref, graph);

        double sigf = 2.4 * 4.0;

        Matrix mref = {sigf * 0.5 * 0.5 / 0.5,
                       sigf * 0.5 * 1.0 / 0.5,
                       sigf * 0.5 * 0.5 / 0.5};

        EXPECT_VEC_SOFTEQ(mref, matrix, 1.0e-14);
    }
}

//---------------------------------------------------------------------------//

TEST_F(Fission_Matrix_TallyTest, output)
{
    std::ostringstream m;
    m << "out" << nodes;

    tally_db->set("start_cycle", 1);
    tally_db->set("output_cycle", 2);
    db->set("problem_name", m.str());

    Tally fm(db, physics, mesh);

    // get the pathlength and source tallies from the compound tally
    auto pl  = fm.get_pl_tally();
    auto src = fm.get_src_tally();

    Particle_t p;
    p.set_wt(0.5);
    p.set_group(0);
    p.set_matid(1);

    State &state = p.geo_state();
    state.d_dir = Space_Vector(1.0, 0.0, 0.0);
    state.d_r   = Space_Vector(1.5, 1.5, 1.5);

    // step in +x
    src->birth(p);
    pl->accumulate(2.0, p);
    fm.end_cycle(1);

    // step in +x
    src->birth(p);
    pl->accumulate(2.0, p);
    fm.end_cycle(1);

    // step in +x
    src->birth(p);
    pl->accumulate(2.0, p);
    fm.end_cycle(1);

    // step in +y
    state.d_dir = Space_Vector(0.0, 1.0, 0.0);
    src->birth(p);
    pl->accumulate(2.0, p);
    fm.end_cycle(1);

    profugus::global_barrier();

    // check the file
    if (node == 0)
    {
        std::string f = m.str() + "_fm.h5";
        profugus::HDF5_Reader reader;
        reader.open(f);

        EXPECT_FALSE(reader.query("cycle_1"));
        EXPECT_TRUE(reader.query("cycle_2"));
        EXPECT_TRUE(reader.query("cycle_3"));
        EXPECT_FALSE(reader.query("cycle_4"));

        double sigf = 2.4 * 4.0;
        double wt   = 0.0;

        Matrix matrix;
        Graph  idx;
        int    n, nz;
        Decomp d;
        d.order = profugus::HDF5_Reader::ROW_MAJOR;

        reader.begin_group("cycle_2");

        reader.read("size", n);
        reader.read("non_zero", nz);
        EXPECT_EQ(64, n);
        EXPECT_EQ(3, nz);

        reader.get_decomposition("indices", d);
        EXPECT_EQ(2, d.ndims);
        EXPECT_EQ(3, d.global[0]);
        EXPECT_EQ(2, d.global[1]);

        d.local[0] = 3;
        d.local[1] = 2;

        idx.resize(3);

        reader.read("indices", d, &idx[0].first);

        EXPECT_EQ(Idx(21, 21), idx[0]);
        EXPECT_EQ(Idx(22, 21), idx[1]);
        EXPECT_EQ(Idx(23, 21), idx[2]);

        reader.read("matrix", matrix);

        wt = 1.0;

        Matrix m2ref = {sigf * 0.5 * 1.0 / wt,
                        sigf * 0.5 * 2.0 / wt,
                        sigf * 0.5 * 1.0 / wt};

        EXPECT_VEC_SOFTEQ(m2ref, matrix, 1.0e-14);

        reader.end_group();

        reader.begin_group("cycle_3");

        reader.read("size", n);
        reader.read("non_zero", nz);
        EXPECT_EQ(64, n);
        EXPECT_EQ(5, nz);

        reader.get_decomposition("indices", d);
        EXPECT_EQ(2, d.ndims);
        EXPECT_EQ(5, d.global[0]);
        EXPECT_EQ(2, d.global[1]);

        d.local[0] = 5;
        d.local[1] = 2;

        idx.resize(5);

        reader.read("indices", d, &idx[0].first);

        EXPECT_EQ(Idx(21, 21), idx[0]);
        EXPECT_EQ(Idx(22, 21), idx[1]);
        EXPECT_EQ(Idx(23, 21), idx[2]);
        EXPECT_EQ(Idx(25, 21), idx[3]);
        EXPECT_EQ(Idx(29, 21), idx[4]);

        reader.read("matrix", matrix);

        wt = 1.5;

        Matrix m3ref = {sigf * 0.5 * 1.5 / wt,
                        sigf * 0.5 * 2.0 / wt,
                        sigf * 0.5 * 1.0 / wt,
                        sigf * 0.5 * 1.0 / wt,
                        sigf * 0.5 * 0.5 / wt};

        EXPECT_VEC_SOFTEQ(m3ref, matrix, 1.0e-14);

        reader.end_group();

        reader.close();
    }
}

//---------------------------------------------------------------------------//
// PAIR HASH TEST
//---------------------------------------------------------------------------//

struct My_Hash
{
  public:

    std::hash<int> d_hash;
    int d_N;

    My_Hash(int N = 1000) : d_N(N) {/*...*/}

    size_t operator()(const std::pair<int, int> &x) const
    {
        return d_hash(x.first + d_N * x.second);
    }
};

TEST(Sparse_Matrix, proto)
{
    typedef std::pair<int, int> Key;
    typedef My_Hash             Hash;

    std::unordered_map<Key, double, Hash> m;

    Hash h(10);
    std::unordered_map<Key, double, Hash> matrix(m.bucket_count(), h);
    EXPECT_EQ(m.bucket_count(), matrix.bucket_count());

    // add some values
    matrix[Key(0, 1)] = 10.1;
    matrix[Key(1, 1)] = 11.1;
    matrix[Key(4, 2)] = 12.1;

    EXPECT_EQ(3, matrix.size());

    EXPECT_EQ(10.1, matrix.find(Key(0, 1))->second);
    EXPECT_EQ(11.1, matrix.find(Key(1, 1))->second);
    EXPECT_EQ(12.1, matrix.find(Key(4, 2))->second);

    std::vector<Key> keys(8);
    auto itr = keys.begin();

    for (const auto &e : matrix)
    {
        *itr = e.first;
        ++itr;
    }

    // push some more keys
    *itr = Key(9, 9); ++itr;
    *itr = Key(0, 3); ++itr;
    *itr = Key(4, 1); ++itr;
    *itr = Key(4, 5); ++itr;
    *itr = Key(1, 1); ++itr;
    EXPECT_EQ(keys.end(), itr);

    std::sort(keys.begin(), keys.end());
    EXPECT_EQ(8, keys.size());

    auto end = std::unique(keys.begin(), keys.end());
    EXPECT_EQ(7, std::distance(keys.begin(), end));

    std::vector<Key> ref = {Key(0, 1), Key(0, 3), Key(1, 1), Key(4, 1),
                            Key(4, 2), Key(4, 5), Key(9, 9)};

    for (auto k = keys.begin(), r = ref.begin(); k != end; ++k, ++r)
    {
        EXPECT_EQ(*r, *k);
    }
}

//---------------------------------------------------------------------------//
//                 end of tstFission_Matrix_Tally.cc
//---------------------------------------------------------------------------//
