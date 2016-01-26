//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   utils/test/tstMetaclass.cc
 * \author Seth R Johnson
 * \date   Tue Oct 15 13:43:00 2013
 * \brief  Test Metaclass
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../Metaclass.t.hh"

#include "gtest/utils_gtest.hh"
#include "utils/Vector_Lite.hh"

using profugus::Metaclass;

//---------------------------------------------------------------------------//
// Needed utility classes
//---------------------------------------------------------------------------//
// Empty classes for creating a unique dynamic struct class

struct Testing_Instantiator_A
{
};

struct Testing_Instantiator_B
{
};

struct Instance_Counter
{
    Instance_Counter()
    {
        ++num_instances;
    }

    ~Instance_Counter()
    {
        --num_instances;
    }

    // Whatever instance data
    double dummy;

    // >>> STATIC DATA
    static unsigned int num_instances;
};

unsigned int Instance_Counter::num_instances = 0;

//---------------------------------------------------------------------------//
// Force methods to instantiate
template class profugus::Metaclass<Testing_Instantiator_A>;
template class profugus::Metaclass<Testing_Instantiator_B>;

//---------------------------------------------------------------------------//
// Metaclass "instances", which are actually classes
typedef profugus::Metaclass<Testing_Instantiator_A> Metaclass_A;
typedef profugus::Metaclass<Testing_Instantiator_B> Metaclass_B;

//---------------------------------------------------------------------------//
// Basic metaclass test
//---------------------------------------------------------------------------//

class MetaclassTest : public ::testing::Test
{
  protected:

    void SetUp()
    {
    }

    void TearDown()
    {
        Metaclass_A::reset();
        Metaclass_B::reset();
    }
};

//---------------------------------------------------------------------------//

TEST_F(MetaclassTest, basic_members)
{
    typedef Metaclass_A Metaclass_t;

    // Add entries to our templated metaclass
    // (These would normally be stored as static variables of a class)
    const unsigned int a_idx = Metaclass_t::new_pod_member<int>("int_a");
    const unsigned int b_idx = Metaclass_t::new_pod_member<int>("int_b");
    const unsigned int c_idx = Metaclass_t::new_pod_member<float>("float_c");

    // Because we've called "reset", we expect a_idx, etc. to be 0, 1, 2.
    // However, in an actual code, the "new member" method may be called after
    // some prior method has added a new member, so the index that gets
    // returned should always be stored and used.
    EXPECT_EQ(0, a_idx);
    EXPECT_EQ(1, b_idx);
    EXPECT_EQ(2, c_idx);

    // Test static methods
    ASSERT_EQ(3, Metaclass_t::size());
    ASSERT_GE(Metaclass_t::storage_size(), 2 * sizeof(int) + sizeof(float));

    EXPECT_EQ("int_a"  , Metaclass_t::name(a_idx));
    EXPECT_EQ("int_b"  , Metaclass_t::name(b_idx));
    EXPECT_EQ("float_c", Metaclass_t::name(c_idx));

    // Create one instance of the metaclass
    Metaclass_t mystruct;

    // Set member data
    mystruct.access<int>(a_idx)   = 3;
    mystruct.access<int>(b_idx)   = 1023;
    mystruct.access<float>(c_idx) = 100.5f;

    // Check member data
    const Metaclass_t& mystruct_const = mystruct;
    EXPECT_EQ(3,      mystruct_const.access<int>(a_idx)  );
    EXPECT_EQ(1023,   mystruct_const.access<int>(b_idx)  );
    EXPECT_EQ(100.5f, mystruct_const.access<float>(c_idx));

    // Check copy construction
    {
        Metaclass_t duplicate(mystruct_const);

        EXPECT_EQ(3,      duplicate.access<int>(a_idx)  );
        EXPECT_EQ(1023,   duplicate.access<int>(b_idx)  );
        EXPECT_EQ(100.5f, duplicate.access<float>(c_idx));

        duplicate.access<int>(b_idx) = -127;
        EXPECT_EQ(-127,duplicate.access<int>(b_idx));
    }

    // Make sure that didn't change the original
    EXPECT_EQ(1023, mystruct_const.access<int>(b_idx));
}

//---------------------------------------------------------------------------//

TEST_F(MetaclassTest, assignment)
{
    typedef Metaclass_A Metaclass_t;
    const unsigned int a_idx = Metaclass_t::new_pod_member<int>("int_a");
    const unsigned int c_idx = Metaclass_t::new_pod_member<float>("float_c");

    // Create one instance of the metaclass
    Metaclass_t a, b;

    a.access<int>(a_idx)   = 3;
    a.access<float>(c_idx) = 100.5f;

    b = a;
    a.access<int>(a_idx) = 4; // change the original

    EXPECT_EQ(3,      b.access<int>(a_idx)  );
    EXPECT_EQ(100.5f, b.access<float>(c_idx));
}

//---------------------------------------------------------------------------//

TEST_F(MetaclassTest, complex_members)
{
    typedef Metaclass_A Metaclass_t;
    typedef profugus::Vector_Lite<float,3> Space_Vector;
    typedef profugus::Vector_Lite<int,3>   Dim_Vector;

    const unsigned int pos_idx
        = Metaclass_t::new_pod_member<Space_Vector>("pos");
    EXPECT_EQ(0, pos_idx);
    const unsigned int ijk_idx
        = Metaclass_t::new_pod_member<Dim_Vector>("ijk");
    EXPECT_EQ(1, ijk_idx);

    // Test static methods
    ASSERT_EQ(2, Metaclass_t::size());
    ASSERT_GE(Metaclass_t::storage_size(), 2 * sizeof(Space_Vector));

    // Create one instance of the metaclass
    Metaclass_t mystruct;

    // Set member data
    Space_Vector& pos = mystruct.access<Space_Vector>(pos_idx);
    Dim_Vector&   ijk = mystruct.access<Dim_Vector>(ijk_idx);

#ifdef CHECK_ON
    // Should be initialized to zero because of memset
    EXPECT_EQ(Space_Vector(0.,0.,0.), pos);
    EXPECT_EQ(Dim_Vector(0,0,0)     , ijk);
#else
    // Hide warnings if DBC disabled
    (void)sizeof(pos);
    (void)sizeof(ijk);
#endif
}

//---------------------------------------------------------------------------//

struct Equivalent_A
{
    int   a_a;
    float a_b;
    float a_c;
};

struct Equivalent_B
{
    char   b_a;
    double b_b;
};

//---------------------------------------------------------------------------//

TEST_F(MetaclassTest, multiple_classes)
{
    const unsigned int aa_idx = Metaclass_A::new_pod_member<int>("a_a");
    const unsigned int ab_idx = Metaclass_A::new_pod_member<float>("a_b");
    const unsigned int ac_idx = Metaclass_A::new_pod_member<float>("a_c");
    (void)ac_idx;

    const unsigned int ba_idx = Metaclass_B::new_pod_member<char>("b_a");
    const unsigned int bb_idx = Metaclass_B::new_pod_member<double>("b_b");

    ASSERT_EQ(3, Metaclass_A::size());
    ASSERT_EQ(2, Metaclass_B::size());

    ASSERT_GE(Metaclass_A::storage_size(), sizeof(int)  + 2 * sizeof(float) );
    ASSERT_GE(Metaclass_B::storage_size(), sizeof(char) +     sizeof(double));

    EXPECT_EQ(sizeof(Equivalent_A), Metaclass_A::storage_size());
    EXPECT_EQ(sizeof(Equivalent_B), Metaclass_B::storage_size());

    // Create one instance of the metaclass
    Metaclass_A a;

    a.access<int>(aa_idx)   = 3;
    a.access<float>(ab_idx) = 100.5f;

    Metaclass_B b;
    b.access<char>(ba_idx)   = '!';
    b.access<double>(bb_idx) = 4.;

    EXPECT_EQ(3 , a.access<int>(aa_idx));
    EXPECT_EQ(4., b.access<double>(bb_idx));
}


//---------------------------------------------------------------------------//

TEST_F(MetaclassTest, error_checking)
{
#ifndef REQUIRE_ON
    SKIP_TEST("DBC Require statements are disabled");
#endif

    const unsigned int aa_idx = Metaclass_A::new_pod_member<int>("a_a");
    const unsigned int ab_idx = Metaclass_A::new_pod_member<float>("a_b");
    (void)ab_idx;

    const unsigned int ba_idx = Metaclass_B::new_pod_member<char>("b_a");
    const unsigned int bb_idx = Metaclass_B::new_pod_member<double>("b_b");
    (void)ba_idx;
    (void)bb_idx;

    {
        // Create an instance
        Metaclass_A a;

        // we shouldn't be able to add a new member since there's an extant
        // instance
        ASSERT_THROW(Metaclass_A::new_pod_member<double>("oops"),
                     profugus::assertion);

        // nor should we be able to reset
        ASSERT_THROW(Metaclass_A::reset(), profugus::assertion);
    }

    ASSERT_NO_THROW(Metaclass_A::new_pod_member<double>("yay"));

    Metaclass_A a;
    a.access<int>(aa_idx)   = 3;

    // sizeof(int) != sizeof(char)
    ASSERT_THROW(a.access<char>(aa_idx), profugus::assertion);
    // out of bounds
    ASSERT_THROW(a.access<int>(5), profugus::assertion);
}

//---------------------------------------------------------------------------//

TEST_F(MetaclassTest, packing)
{
    const unsigned int a_idx = Metaclass_A::new_pod_member<char>("char");
    const unsigned int b_idx = Metaclass_A::new_pod_member<int>("int");
    const unsigned int c_idx = Metaclass_A::new_pod_member<double>("double");

    ASSERT_EQ(3, Metaclass_A::size());
    ASSERT_GE(Metaclass_A::storage_size(),
            sizeof(char) + sizeof(int) + sizeof(double));

    // Create one instance of the metaclass
    Metaclass_A a;

    a.access<char>(a_idx)   = 'S';
    a.access<int>(b_idx)    = 12345;
    a.access<double>(c_idx) = 3.14159;

    EXPECT_EQ(a.packed_size(), sizeof(char) + sizeof(int) + sizeof(double));

    std::vector<char> buffer(a.packed_size());
    a.pack(&buffer[0]);
    EXPECT_EQ('S', buffer[0]);

    // Test unpacking
    const Metaclass_A b(&buffer[0], buffer.size());
    EXPECT_EQ('S'    , b.access<char  >(a_idx));
    EXPECT_EQ(12345  , b.access<int   >(b_idx));
    EXPECT_EQ(3.14159, b.access<double>(c_idx));
}

//---------------------------------------------------------------------------//
// Vector tests
//---------------------------------------------------------------------------//

// Add vectors to each Metaclass_A instances so it looks like:
//
// struct Container
// {
//     int                       counter;
//     std::vector<float>        floats;
//     std::vector<Space_Vector> coords;
// };
class MetaclassVecTest : public ::testing::Test
{
  protected:
    typedef profugus::Vector_Lite<double,3> Space_Vector;

    typedef std::vector<float>        Vec_Float;
    typedef std::vector<Space_Vector> Vec_Space_Vector;

    typedef Metaclass_A Metaclass_t;

  protected:
    void SetUp()
    {
        counter_idx = Metaclass_t::new_pod_member<int>("counter");
        floats_idx  = Metaclass_t::new_vec_member<Vec_Float>("floats");
        coords_idx  = Metaclass_t::new_vec_member<Vec_Space_Vector>("coords");
    }

    void build(Metaclass_t& mc)
    {
        mc.access<int>(counter_idx) = 456;
        Vec_Float& flt_vec = mc.access<Vec_Float>(floats_idx);
        Vec_Space_Vector& coo_vec = mc.access<Vec_Space_Vector>(coords_idx);

        flt_vec.push_back(1.23f);
        flt_vec.push_back(4.f);

        coo_vec.push_back(Space_Vector(1.,0.,0.));
        coo_vec.push_back(Space_Vector(0.,2.,0.));
        coo_vec.push_back(Space_Vector(0.,0.,3.));
    }

    void check(const Metaclass_t& mc)
    {
        EXPECT_EQ(456, mc.access<int>(counter_idx));

        const Vec_Float&
            flt_vec = mc.access<Vec_Float>(floats_idx);
        ASSERT_EQ(2    , flt_vec.size());
        EXPECT_EQ(1.23f, flt_vec[0]);
        EXPECT_EQ(4.f  , flt_vec[1]);

        const Vec_Space_Vector&
            coo_vec = mc.access<Vec_Space_Vector>(coords_idx);
        ASSERT_EQ(3  , coo_vec.size());
        EXPECT_EQ(Space_Vector(1.,0.,0.), coo_vec[0]);
        EXPECT_EQ(Space_Vector(0.,2.,0.), coo_vec[1]);
        EXPECT_EQ(Space_Vector(0.,0.,3.), coo_vec[2]);
    }

    void TearDown()
    {
        Metaclass_t::reset();
    }

  protected:
    unsigned int counter_idx;
    unsigned int floats_idx;
    unsigned int coords_idx;
};

//---------------------------------------------------------------------------//

TEST_F(MetaclassVecTest, accessing)
{
    Metaclass_t a;
    build(a);

    check(a);
}

//---------------------------------------------------------------------------//

TEST_F(MetaclassVecTest, copying)
{
    Metaclass_t a;
    build(a);

    Metaclass_t b(a);

    // Delete stuff from first instance
    a.access<Vec_Float>(floats_idx).clear();
    a.access<Vec_Space_Vector>(coords_idx).clear();

    check(b);
}

//---------------------------------------------------------------------------//

TEST_F(MetaclassVecTest, assignment)
{
    Metaclass_t a;
    build(a);

    Metaclass_t b;
    EXPECT_TRUE(b.access<Vec_Float>(floats_idx).empty());

    b = a;
    check(b);
}

//---------------------------------------------------------------------------//

TEST_F(MetaclassVecTest, packing)
{
    Metaclass_t a;
    build(a);

    // Test packing
    EXPECT_EQ(a.packed_size(),
              sizeof(int)
            + sizeof(Metaclass_t::size_type) + 2 * sizeof(float)
            + sizeof(Metaclass_t::size_type) + 3 * sizeof(Space_Vector));

    std::vector<char> buffer(a.packed_size());
    a.pack(&buffer[0]);

    // Delete stuff from Metaclass a
    a.access<Vec_Float>(floats_idx).clear();
    a.access<Vec_Space_Vector>(coords_idx).clear();
    EXPECT_EQ(a.packed_size(),
              sizeof(int)
            + sizeof(Metaclass_t::size_type)
            + sizeof(Metaclass_t::size_type));

    // Test unpacking
    const Metaclass_t b(&buffer[0], buffer.size());
    check(b);
}

//---------------------------------------------------------------------------//
//                 end of tstMetaclass.cc
//---------------------------------------------------------------------------//
