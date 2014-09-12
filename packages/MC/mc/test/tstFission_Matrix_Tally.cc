//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/test/tstFission_Matrix_Tally.cc
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

#include "geometry/Mesh_Geometry.hh"
#include "TransporterTestBase.hh"

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//

class Fission_Matrix_TallyTest : public TransporterTestBase
{
  protected:
    // >>> TYPEDEFS
    typedef profugus::Fission_Matrix_Tally   Tally;
    typedef profugus::Mesh_Geometry          Mesh_Geometry_t;
    typedef std::shared_ptr<Mesh_Geometry_t> SP_Mesh_Geometry;
    typedef Mesh_Geometry_t::Vec_Dbl         Vec_Dbl;
    typedef Particle_t::Geo_State_t          State;

  protected:
    void SetUp()
    {
        TransporterTestBase::SetUp();

        // make the mesh geometry
        Vec_Dbl r = {0.0, 1.0, 2.0, 3.0, 4.0};
        mesh = std::make_shared<Mesh_Geometry_t>(r, r, r);

        // make fission matrix db
        fis_db = Teuchos::sublist(db, "fission_matrix_db");

        db->set("num_cycles", 10);
    }

  protected:
    // >>> DATA

    SP_Mesh_Geometry mesh;
    RCP_Std_DB fis_db;
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST_F(Fission_Matrix_TallyTest, construction)
{
    Tally fm(db, physics, mesh);
}

//---------------------------------------------------------------------------//

TEST_F(Fission_Matrix_TallyTest, transport)
{
    Tally fm(db, physics, mesh);

    Particle_t p;
    p.set_wt(0.5);
    p.set_group(0);
    p.set_matid(1);

    State &state = p.geo_state();
    state.d_dir = Space_Vector(1.0, 0.0, 0.0);
    state.d_r   = Space_Vector(1.5, 1.5, 1.5);

    // step in +x
    fm.birth(p);
    fm.accumulate(2.0, p);
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
