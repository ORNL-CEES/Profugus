//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   source/test/tstIsotropic_Source.cc
 * \author Thomas M. Evans
 * \date   Thu Oct 22 14:07:39 2009
 * \brief  Isotropic_Source test.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//


#include "gtest/utils_gtest.hh"
#include <vector>
#include <cmath>
#include <sstream>

#include "utils/Constants.hh"
#include "../Isotropic_Source.hh"

using namespace std;

using profugus::constants::inv_four_pi;
using profugus::Isotropic_Source;

typedef Isotropic_Source::ID_Field      ID_Field;
typedef Isotropic_Source::Shape         Shape;
typedef Isotropic_Source::Source_Shapes Source_Shapes;
typedef Isotropic_Source::Source_Field  Source_Field;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST(Source, all)
{
    // assume a 6 cell mesh with 4 groups and 3 source shapes
    Isotropic_Source q(6);

    EXPECT_EQ(0, q.num_groups());
    EXPECT_TRUE(!q.verify());

    // build the source
    {
        ID_Field ids(6, 0);
        Source_Shapes shape(3, Shape(4, 0.0));
        Source_Field source(6, 0.0);

        ids[1] = 1;
        ids[2] = 2;
        ids[3] = 1;
        ids[5] = 2;

        source[1] = 1.0;
        source[2] = 1.1;
        source[3] = 1.2;
        source[4] = 0.9;
        source[5] = 0.9;

        shape[0][0] = 0.0;
        shape[0][1] = 0.0;
        shape[0][2] = 0.0;
        shape[0][3] = 0.0;

        shape[1][0] = 1.5;
        shape[1][1] = 1.6;
        shape[1][2] = 1.7;
        shape[1][3] = 1.8;

        shape[2][0] = 2.5;
        shape[2][1] = 2.6;
        shape[2][2] = 2.7;
        shape[2][3] = 2.8;

        q.set(ids, shape, source);
    }

    EXPECT_EQ(4, q.num_groups());
    EXPECT_TRUE(q.verify());

    EXPECT_TRUE(soft_equiv(q.q_e(0, 0), 0.0));
    EXPECT_TRUE(soft_equiv(q.q_e(0, 1), 0.0));
    EXPECT_TRUE(soft_equiv(q.q_e(0, 2), 0.0));
    EXPECT_TRUE(soft_equiv(q.q_e(0, 3), 0.0));

    EXPECT_TRUE(soft_equiv(q.q_e(1, 0), 1.0 * 1.5 * inv_four_pi));
    EXPECT_TRUE(soft_equiv(q.q_e(1, 1), 1.0 * 1.6 * inv_four_pi));
    EXPECT_TRUE(soft_equiv(q.q_e(1, 2), 1.0 * 1.7 * inv_four_pi));
    EXPECT_TRUE(soft_equiv(q.q_e(1, 3), 1.0 * 1.8 * inv_four_pi));

    EXPECT_TRUE(soft_equiv(q.q_e(2, 0), 1.1 * 2.5 * inv_four_pi));
    EXPECT_TRUE(soft_equiv(q.q_e(2, 1), 1.1 * 2.6 * inv_four_pi));
    EXPECT_TRUE(soft_equiv(q.q_e(2, 2), 1.1 * 2.7 * inv_four_pi));
    EXPECT_TRUE(soft_equiv(q.q_e(2, 3), 1.1 * 2.8 * inv_four_pi));

    EXPECT_TRUE(soft_equiv(q.q_e(3, 0), 1.2 * 1.5 * inv_four_pi));
    EXPECT_TRUE(soft_equiv(q.q_e(3, 1), 1.2 * 1.6 * inv_four_pi));
    EXPECT_TRUE(soft_equiv(q.q_e(3, 2), 1.2 * 1.7 * inv_four_pi));
    EXPECT_TRUE(soft_equiv(q.q_e(3, 3), 1.2 * 1.8 * inv_four_pi));

    EXPECT_TRUE(soft_equiv(q.q_e(4, 0), 0.0));
    EXPECT_TRUE(soft_equiv(q.q_e(4, 1), 0.0));
    EXPECT_TRUE(soft_equiv(q.q_e(4, 2), 0.0));
    EXPECT_TRUE(soft_equiv(q.q_e(4, 3), 0.0));

    EXPECT_TRUE(soft_equiv(q.q_e(5, 0), 0.9 * 2.5 * inv_four_pi));
    EXPECT_TRUE(soft_equiv(q.q_e(5, 1), 0.9 * 2.6 * inv_four_pi));
    EXPECT_TRUE(soft_equiv(q.q_e(5, 2), 0.9 * 2.7 * inv_four_pi));
    EXPECT_TRUE(soft_equiv(q.q_e(5, 3), 0.9 * 2.8 * inv_four_pi));
}

//---------------------------------------------------------------------------//

TEST(Truncation, all)
{
    // assume a 6 cell mesh with 4 groups and 3 source shapes
    Isotropic_Source q(6);

    EXPECT_EQ(0, q.num_groups());
    EXPECT_TRUE(!q.verify());

    // build the source
    {
        ID_Field ids(6, 0);
        Source_Shapes shape(3, Shape(4, 0.0));
        Source_Field source(6, 0.0);

        ids[1] = 1;
        ids[2] = 2;
        ids[3] = 1;
        ids[5] = 2;

        source[1] = 1.0;
        source[2] = 1.1;
        source[3] = 1.2;
        source[4] = 0.9;
        source[5] = 0.9;

        shape[0][0] = 0.0;
        shape[0][1] = 0.0;
        shape[0][2] = 0.0;
        shape[0][3] = 0.0;

        shape[1][0] = 1.5;
        shape[1][1] = 1.6;
        shape[1][2] = 1.7;
        shape[1][3] = 1.8;

        shape[2][0] = 2.5;
        shape[2][1] = 2.6;
        shape[2][2] = 2.7;
        shape[2][3] = 2.8;

        q.set(ids, shape, source);
    }

    // full range truncation (no-op)
    q.truncate(0, 3);
    EXPECT_EQ(4, q.num_groups());

    // truncation to groups 1-2
    q.truncate(1, 2);

    EXPECT_EQ(2, q.num_groups());
    EXPECT_TRUE(q.verify());

    EXPECT_TRUE(soft_equiv(q.q_e(0, 0), 0.0));
    EXPECT_TRUE(soft_equiv(q.q_e(0, 1), 0.0));

    EXPECT_TRUE(soft_equiv(q.q_e(1, 0), 1.0 * 1.6 * inv_four_pi));
    EXPECT_TRUE(soft_equiv(q.q_e(1, 1), 1.0 * 1.7 * inv_four_pi));

    EXPECT_TRUE(soft_equiv(q.q_e(2, 0), 1.1 * 2.6 * inv_four_pi));
    EXPECT_TRUE(soft_equiv(q.q_e(2, 1), 1.1 * 2.7 * inv_four_pi));

    EXPECT_TRUE(soft_equiv(q.q_e(3, 0), 1.2 * 1.6 * inv_four_pi));
    EXPECT_TRUE(soft_equiv(q.q_e(3, 1), 1.2 * 1.7 * inv_four_pi));

    EXPECT_TRUE(soft_equiv(q.q_e(4, 0), 0.0));
    EXPECT_TRUE(soft_equiv(q.q_e(4, 1), 0.0));

    EXPECT_TRUE(soft_equiv(q.q_e(5, 0), 0.9 * 2.6 * inv_four_pi));
    EXPECT_TRUE(soft_equiv(q.q_e(5, 1), 0.9 * 2.7 * inv_four_pi));
}

//---------------------------------------------------------------------------//
//                        end of tstIsotropic_Source.cc
//---------------------------------------------------------------------------//
