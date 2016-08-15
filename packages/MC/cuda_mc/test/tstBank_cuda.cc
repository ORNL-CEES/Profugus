//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/test/tstBank.cc
 * \author Stuart Slattery
 * \date   Fri Apr 25 16:50:26 2014
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "rng/RNG_Control.hh"
#include "cuda_geometry/Mesh_Geometry.hh"
#include "../Bank.hh"
#include "Particle_Vector_Tester.hh"

#include "gtest/utils_gtest.hh"

//---------------------------------------------------------------------------//
// TESTS (simple test harness)
//---------------------------------------------------------------------------//

class BankTest : public ::testing::Test
{
  protected:
    typedef cuda_profugus::Mesh_Geometry     Geometry_t;
    typedef cuda_profugus::Bank<Geometry_t>  Bank_t;
    typedef Bank_t::Particle_t               Particle;
    typedef Bank_t::SP_Particle              SP_Particle;

    void SetUp()
    {
        // Initialization that are performed for each test
        // changes to these don't propagate between tests
        m_orig_p = std::make_shared<Particle>();
        m_orig_p->set_wt(1.23);
        m_orig_p->set_matid(1);
    }

    // data that get re-initialized between tests
    SP_Particle m_orig_p;
    Bank_t b;
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//
TEST_F( BankTest, cuda_pop )
{
    SP_Particle orig_p2 = std::make_shared<Particle>();
    orig_p2->set_wt(1.23);
    orig_p2->set_matid(2);

    // add two copies of two particles
    b.push(m_orig_p, 128);
    b.push(orig_p2, 128);

    EXPECT_TRUE(!b.empty());
    EXPECT_EQ(256, b.size());
    EXPECT_EQ(2, b.num_unique());

    // pop the particles into a particle vector
    int num_particle = 256;
    profugus::RNG_Control control( 3420239343 );
    Particle_Vector_Tester vector_tester( num_particle, control.rng() );
    cuda_utils::Shared_Device_Ptr<cuda_profugus::Particle_Vector<cuda_profugus::Mesh_Geometry> >
	particles = vector_tester.get_vector();
    b.pop( particles );

    // check the contents of the particle vector
    Teuchos::Array<double> weights = vector_tester.wt();
    Teuchos::Array<int> matids = vector_tester.matid();
    for ( int i = 0; i < num_particle; ++i )
    {
	EXPECT_EQ( weights[i], 1.23 );
	if ( i < 128 )
	{
	    EXPECT_EQ( matids[i], 2 );
	}
	else
	{
	    EXPECT_EQ( matids[i], 1 );
	}
    }
}

//---------------------------------------------------------------------------//
//                 end of tstBank.cc
//---------------------------------------------------------------------------//
