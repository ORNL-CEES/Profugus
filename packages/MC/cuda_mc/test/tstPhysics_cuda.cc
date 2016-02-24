//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/test/tstPhysics.cc
 * \author Thomas M. Evans
 * \date   Fri May 02 00:58:54 2014
 * \brief  Physics unit test.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../Physics.hh"

#include "gtest/utils_gtest.hh"

#include <vector>
#include <memory>
#include <iostream>
#include <iomanip>

#include "Teuchos_RCP.hpp"
#include "utils/Definitions.hh"
#include "rng/RNG_Control.hh"
#include "cuda_utils/Shared_Device_Ptr.hh"
#include "cuda_geometry/Mesh_Geometry.hh"

#include "Physics_Tester.hh"

using namespace std;

using def::X; using def::Y; using def::Z;

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//

template <class Geometry>
class PhysicsTest : public testing::Test
{
  protected:
    typedef Geometry                                    Geometry_t;
    typedef profugus::RNG_Control                       RNG_Control;
    typedef cuda_profugus::Physics<Geometry>            Physics_t;
    typedef typename Physics_t::ParameterList_t         ParameterList_t;
    typedef typename Physics_t::XS_t                    XS_t;
    typedef typename Physics_t::Fission_Site            Fission_Site;
    typedef typename Physics_t::Fission_Site_Container  Fission_Site_Container;
    typedef typename Geometry_t::Space_Vector           Space_Vector;

  protected:

    void SetUp()
    {
        build_geometry();

        build_xs();

        // make a rng
        int seed = 342412;
        profugus::RNG_Control control(seed);
        rng = control.rng();

        // make db
        db = Teuchos::rcp(new ParameterList_t("test"));
    }

    void build_geometry();

    void build_xs()
    {
        xs = Teuchos::rcp(new XS_t());
        xs->set(0, 5);

        vector<double> bnd(6, 0.0);
        bnd[0] = 100.0;
        bnd[1] = 10.0;
        bnd[2] = 1.0;
        bnd[3] = 0.1;
        bnd[4] = 0.01;
        bnd[5] = 0.001;
        xs->set_bounds(bnd);

        typename XS_t::OneDArray total(5);
        typename XS_t::TwoDArray scat(5, 5);

        double f[5] = {0.1, 0.4, 1.8, 5.7, 9.8};
        double c[5] = {0.3770, 0.4421, 0.1809, 0.0, 0.0};
        double n[5] = {2.4*f[0], 2.4*f[1], 2.4*f[2], 2.4*f[3], 2.4*f[4]};
        typename XS_t::OneDArray fission(begin(f), end(f));
        typename XS_t::OneDArray chi(begin(c), end(c));
        typename XS_t::OneDArray nus(begin(n), end(n));
        xs->add(1, XS_t::SIG_F, fission);
        xs->add(1, XS_t::NU_SIG_F, nus);
        xs->add(1, XS_t::CHI, chi);

        // mat 0
        total[0] = 5.2 ;
        total[1] = 11.4;
        total[2] = 18.2;
        total[3] = 29.9;
        total[4] = 27.3;
        xs->add(0, XS_t::TOTAL, total);

        // mat 1
        total[0] = 5.2  + f[0];
        total[1] = 11.4 + f[1];
        total[2] = 18.2 + f[2];
        total[3] = 29.9 + f[3];
        total[4] = 27.3 + f[4];
        xs->add(1, XS_t::TOTAL, total);

        scat(0, 0) = 1.2;
        scat(1, 0) = 0.9;
        scat(1, 1) = 3.2;
        scat(2, 0) = 0.4;
        scat(2, 1) = 2.8;
        scat(2, 2) = 6.9;
        scat(2, 3) = 1.5;
        scat(3, 0) = 0.1;
        scat(3, 1) = 2.1;
        scat(3, 2) = 5.5;
        scat(3, 3) = 9.7;
        scat(3, 4) = 2.1;
        scat(4, 1) = 0.2;
        scat(4, 2) = 1.3;
        scat(4, 3) = 6.6;
        scat(4, 4) = 9.9;
        xs->add(0, 0, scat);
        xs->add(1, 0, scat);

        xs->complete();
    }

  protected:
    def::Vec_Dbl edges;
    Teuchos::RCP<XS_t> xs;
    Teuchos::RCP<ParameterList_t> db;
    profugus::RNG_Control::RNG_t rng;
};

template <>
void PhysicsTest<cuda_profugus::Mesh_Geometry>::build_geometry()
{
    edges = {0.0, 100.0};
}

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

typedef ::testing::Types<cuda_profugus::Mesh_Geometry> MyTypes;
TYPED_TEST_CASE(PhysicsTest, MyTypes);

TYPED_TEST(PhysicsTest, Collisions_no_capture)
{
    typedef typename TestFixture::Physics_t     Physics_t;
    typedef typename TestFixture::Space_Vector  Space_Vector;

    // turn off implicit capture
    this->db->set("implicit_capture", false);

    // make the physics tester.
    int Np = 50000;
    Physics_Tester physics_tester( this->edges, this->edges, this->edges,
				   Np, this->rng, *(this->db), *(this->xs) );

    // initialize particles with the geometry and set to collide.
    Space_Vector r, d;
    r.x = 50.0;
    r.y = 50.0;
    r.z = 50.0;
    d.x = 1.0;
    d.y = 1.0;
    d.z = 1.0;
    physics_tester.geometry_initialize( r, d, 0 );

    // check distributions from analog collisions
    int scat[5]   = {0};
    int abs[5]    = {0};
    int ng[5]     = {0};
    int octant[8] = {0};

    vector< vector<int> > g2g(5, vector<int>(5, 0));

    // group CDF
    std::vector<double> cdf = {0.2, 0.4, 0.6, 0.8, 1.0};

    // Sample the particle groups
    physics_tester.sample_group( cdf );

    // Get the ingoing groups
    auto ing = physics_tester.particle_tester().group();

    // Count the ingoing groups
    for ( auto i : ing ) ng[i]++;

    // collide the particles
    physics_tester.physics().get_host_ptr()->collide( physics_tester.particles() );

    // Get the outgoing groups.
    auto outg = physics_tester.particle_tester().group();

    // Get the outgoing events
    auto events = physics_tester.particle_tester().event();

    // Get the outgoing states.
    auto geo_state = physics_tester.particle_tester().geo_state();

    // count the outgoing data
    for ( int i = 0; i < Np; ++i )
    {
	if (events[i] == cuda_profugus::events::ABSORPTION)
	{
	    abs[ing[i]]++;
	}
	else if (events[i] == cuda_profugus::events::SCATTER)
	{
	    scat[ing[i]]++;
	    g2g[ing[i]][outg[i]]++;

	    Space_Vector omega = geo_state[i].d_dir;

	    // check octant distribution (isotropic scattering)
	    if (omega.z > 0.0)
	    {
		if (omega.y > 0.0)
		{
		    if (omega.x > 0.0)
		    {
			octant[0]++;
		    }
		    else
		    {
			octant[1]++;
		    }
		}
		else
		{
		    if (omega.x > 0.0)
		    {
			octant[2]++;
		    }
		    else
		    {
			octant[3]++;
		    }
		}
	    }
	    else
	    {
		if (omega.y > 0.0)
		{
		    if (omega.x > 0.0)
		    {
			octant[4]++;
		    }
		    else
		    {
			octant[5]++;
		    }
		}
		else
		{
		    if (omega.x > 0.0)
		    {
			octant[6]++;
		    }
		    else
		    {
			octant[7]++;
		    }
		}
	    }
	}
    }

    int tots         = 0, tota = 0;
    double ref[]     = {0.5000, 0.7281, 0.7527, 0.5953, 0.4396};
    double g2gr[][5] = {{0.4615, 0.3462, 0.1538, 0.0385, 0.0000},
			{0.0000, 0.3855, 0.3373, 0.2530, 0.0241},
			{0.0000, 0.0000, 0.5036, 0.4015, 0.0949},
			{0.0000, 0.0000, 0.0843, 0.5449, 0.3708},
			{0.0000, 0.0000, 0.0000, 0.1750, 0.8250}};

    cout.precision(4);
    cout << setw(3) << "g" << setw(10) << "abs" << setw(10) << "scat"
	 << setw(10) << "c" << setw(10) << "diff" << endl;
    cout << "-------------------------------------------" << endl;
    for (int g = 0; g < 5; ++g)
    {
	double a = static_cast<double>(abs[g]) / ng[g];
	double s = static_cast<double>(scat[g]) / ng[g];

	cout << setw(3) << g
	     << setw(10) << fixed << a
	     << setw(10) << fixed << s
	     << setw(10) << fixed << ref[g]
	     << setw(10) << fixed << std::fabs(s - ref[g]) / ref[g]
	     << endl;

	EXPECT_SOFTEQ(s, ref[g], 0.04);

	tota += abs[g];
	tots += scat[g];
    }
    EXPECT_EQ(Np, tots + tota);

    cout << endl;
    cout << "Group-to-group Scattering" << endl;
    cout << "------------------------------------" << endl;
    for (int g = 0; g < 5; ++g)
    {
	for (int gp = 0; gp < 5; ++gp)
	{
	    double s = static_cast<double>(g2g[g][gp]) / scat[g];
	    double e = 0.0;
	    if (g2gr[g][gp] > 0.0)
		e = std::fabs(s - g2gr[g][gp]) / g2gr[g][gp];

	    cout << setw(1) << g << " -> " << setw(1) << gp
		 << setw(10) << fixed << s
		 << setw(10) << fixed << g2gr[g][gp]
		 << setw(10) << fixed << e
		 << endl;

	    EXPECT_SOFTEQ(s, g2gr[g][gp], 0.07);
	}
    }

    cout << endl;
    cout << "Isotropic scattering by octant" << endl;
    cout << "---------------------------------" << endl;
    int osum = 0;
    for (int i = 0; i < 8; ++i)
    {
	double o  = static_cast<double>(octant[i]) / tots;
	double r  = 1.0/8.0;
	osum     += octant[i];
	cout << setw(3) << i
	     << setw(10) << fixed << o
	     << setw(10) << fixed << r
	     << setw(10) << std::abs(o - r) / r << endl;

	EXPECT_SOFTEQ(o, r, 0.03);
    }
    EXPECT_EQ(tots, osum);
}

//---------------------------------------------------------------------------//
TYPED_TEST(PhysicsTest, Collisions_with_capture)
{
    typedef typename TestFixture::Physics_t     Physics_t;
    typedef typename TestFixture::Space_Vector  Space_Vector;

    // turn on implicit capture
    this->db->set("implicit_capture", true);

    // make the physics tester.
    int Np = 50000;
    Physics_Tester physics_tester( this->edges, this->edges, this->edges,
				   Np, this->rng, *(this->db), *(this->xs) );

    // initialize particles with the geometry and set to collide.
    Space_Vector r, d;
    r.x = 50.0;
    r.y = 50.0;
    r.z = 50.0;
    d.x = 1.0;
    d.y = 1.0;
    d.z = 1.0;
    physics_tester.geometry_initialize( r, d, 0 );
    
    // check distributions from implicit-capture collisions
    std::vector<double> c = {0.5000, 0.7281, 0.7527, 0.5953, 0.4396};

    // group CDF
    std::vector<double> cdf = {0.2, 0.4, 0.6, 0.8, 1.0};

    // Sample the particle groups
    physics_tester.sample_group( cdf );

    // Get the ingoing groups
    auto ing = physics_tester.particle_tester().group();

    // collide the particles
    physics_tester.physics().get_host_ptr()->collide( physics_tester.particles() );

    // Get the outgoing groups.
    auto outg = physics_tester.particle_tester().group();

    // Get the outgoing events
    auto events = physics_tester.particle_tester().event();

    // Get the outgoing states.
    auto geo_state = physics_tester.particle_tester().geo_state();

    // Get the weights.
    auto weights = physics_tester.particle_tester().wt();

    // check distributions from nonanalog collisions
    int octant[8] = {0};

    for (int n = 0; n < Np; ++n)
    {
	EXPECT_EQ( cuda_profugus::events::IMPLICIT_CAPTURE, events[n] );
	EXPECT_SOFTEQ( weights[n], c[ing[n]] * 0.9, 1.0e-4 );

	const Space_Vector &omega =geo_state[n].d_dir;

	// check octant distribution (isotropic scattering)
	if (omega.z > 0.0)
	{
	    if (omega.y > 0.0)
	    {
		if (omega.x > 0.0)
		{
		    octant[0]++;
		}
		else
		{
		    octant[1]++;
		}
	    }
	    else
	    {
		if (omega.x > 0.0)
		{
		    octant[2]++;
		}
		else
		{
		    octant[3]++;
		}
	    }
	}
	else
	{
	    if (omega.y > 0.0)
	    {
		if (omega.x > 0.0)
		{
		    octant[4]++;
		}
		else
		{
		    octant[5]++;
		}
	    }
	    else
	    {
		if (omega.x > 0.0)
		{
		    octant[6]++;
		}
		else
		{
		    octant[7]++;
		}
	    }
	}
    }

    cout << endl;
    cout << "Isotropic scattering by octant" << endl;
    cout << "---------------------------------" << endl;
    for (int i = 0; i < 8; ++i)
    {
	double o  = static_cast<double>(octant[i]) / Np;
	double r  = 1.0/8.0;
	cout << setw(3) << i
	     << setw(10) << fixed << o
	     << setw(10) << fixed << r
	     << setw(10) << std::abs(o - r) / r << endl;

	EXPECT_SOFTEQ(o, r, 0.04);
    }
}

//---------------------------------------------------------------------------//

TYPED_TEST(PhysicsTest, Access)
{
    typedef typename TestFixture::Physics_t     Physics_t;
    typedef typename TestFixture::Space_Vector  Space_Vector;

    using cuda_profugus::physics::TOTAL;
    using cuda_profugus::physics::SCATTERING;
    using cuda_profugus::physics::FISSION;

    // make the physics tester.
    int Np = 50;
    Physics_Tester physics_tester( this->edges, this->edges, this->edges,
				   Np, this->rng, *(this->db), *(this->xs) );

    EXPECT_FALSE(physics_tester.is_fissionable(0));
    EXPECT_TRUE(physics_tester.is_fissionable(1));

    // test data access without fission
    EXPECT_SOFTEQ(5.2, physics_tester.get_total(0,0,TOTAL), 1.e-12);
    EXPECT_SOFTEQ(2.6, physics_tester.get_total(0,0,SCATTERING), 1.e-12);
    EXPECT_SOFTEQ(0.0, physics_tester.get_total(0,0,FISSION), 1.e-12);

    EXPECT_SOFTEQ(11.4, physics_tester.get_total(0,1,TOTAL), 1.e-12);
    EXPECT_SOFTEQ(8.3, physics_tester.get_total(0,1,SCATTERING), 1.e-12);
    EXPECT_SOFTEQ(0.0, physics_tester.get_total(0,1,FISSION), 1.e-12);

    EXPECT_SOFTEQ(18.2, physics_tester.get_total(0,2,TOTAL), 1.e-12);
    EXPECT_SOFTEQ(13.7, physics_tester.get_total(0,2,SCATTERING), 1.e-12);
    EXPECT_SOFTEQ(0.0, physics_tester.get_total(0,2,FISSION), 1.e-12);

    EXPECT_SOFTEQ(29.9, physics_tester.get_total(0,3,TOTAL), 1.e-12);
    EXPECT_SOFTEQ(17.8, physics_tester.get_total(0,3,SCATTERING), 1.e-12);
    EXPECT_SOFTEQ(0.0, physics_tester.get_total(0,3,FISSION), 1.e-12);

    EXPECT_SOFTEQ(27.3, physics_tester.get_total(0,4,TOTAL), 1.e-12);
    EXPECT_SOFTEQ(12.0, physics_tester.get_total(0,4,SCATTERING), 1.e-12);
    EXPECT_SOFTEQ(0.0, physics_tester.get_total(0,4,FISSION), 1.e-12);

    // test data access with fission
    EXPECT_SOFTEQ(5.3, physics_tester.get_total(1,0,TOTAL), 1.e-12);
    EXPECT_SOFTEQ(2.6, physics_tester.get_total(1,0,SCATTERING), 1.e-12);
    EXPECT_SOFTEQ(0.1, physics_tester.get_total(1,0,FISSION), 1.e-12);

    EXPECT_SOFTEQ(11.8, physics_tester.get_total(1,1,TOTAL), 1.e-12);
    EXPECT_SOFTEQ(8.3, physics_tester.get_total(1,1,SCATTERING), 1.e-12);
    EXPECT_SOFTEQ(0.4, physics_tester.get_total(1,1,FISSION), 1.e-12);

    EXPECT_SOFTEQ(20.0, physics_tester.get_total(1,2,TOTAL), 1.e-12);
    EXPECT_SOFTEQ(13.7, physics_tester.get_total(1,2,SCATTERING), 1.e-12);
    EXPECT_SOFTEQ(1.8, physics_tester.get_total(1,2,FISSION), 1.e-12);

    EXPECT_SOFTEQ(35.6, physics_tester.get_total(1,3,TOTAL), 1.e-12);
    EXPECT_SOFTEQ(17.8, physics_tester.get_total(1,3,SCATTERING), 1.e-12);
    EXPECT_SOFTEQ(5.7, physics_tester.get_total(1,3,FISSION), 1.e-12);

    EXPECT_SOFTEQ(37.1, physics_tester.get_total(1,4,TOTAL), 1.e-12);
    EXPECT_SOFTEQ(12.0, physics_tester.get_total(1,4,SCATTERING), 1.e-12);
    EXPECT_SOFTEQ(9.8, physics_tester.get_total(1,4,FISSION), 1.e-12);
}

//---------------------------------------------------------------------------//

TYPED_TEST(PhysicsTest, initialization)
{
    typedef typename TestFixture::Physics_t     Physics_t;
    typedef typename TestFixture::Space_Vector  Space_Vector;

    // make the physics tester.
    int Np = 16;
    Physics_Tester physics_tester( this->edges, this->edges, this->edges,
				   Np, this->rng, *(this->db), *(this->xs) );

    // Create initializization energies.
    std::vector<double> energy =
	{ 100.0, 99.99, 
	  10.01, 10.0, 9.99, 
	  1.01, 1.0, 0.99,
	  0.101, 0.1, 0.099,
	  0.011, 0.01, 0.0099, 
	  0.0011, 0.001 };

    // Initialize the particles.
    physics_tester.physics().get_host_ptr()->initialize( energy, physics_tester.particles() );

    // Get the partice groups.
    auto groups = physics_tester.particle_tester().group();

    // check initializations
    EXPECT_EQ(0, groups[0]);
    EXPECT_EQ(0, groups[1]);
    EXPECT_EQ(0, groups[2]);
    EXPECT_EQ(0, groups[3]);
    EXPECT_EQ(1, groups[4]);
    EXPECT_EQ(1, groups[5]);
    EXPECT_EQ(1, groups[6]);
    EXPECT_EQ(2, groups[7]);
    EXPECT_EQ(2, groups[8]);
    EXPECT_EQ(2, groups[9]);
    EXPECT_EQ(3, groups[10]);
    EXPECT_EQ(3, groups[11]);
    EXPECT_EQ(3, groups[12]);
    EXPECT_EQ(4, groups[13]);
    EXPECT_EQ(4, groups[14]);
    EXPECT_EQ(4, groups[15]);

    // Check energy bounds
    double min, max;
    physics_tester.get_min_max_energy( min, max );
    EXPECT_EQ(0.001, min);
    EXPECT_EQ(100.0, max);
}

// //---------------------------------------------------------------------------//

// TYPED_TEST(PhysicsTest, fission_sampling)
// {
//     typedef typename TestFixture::Particle                  Particle;
//     typedef typename TestFixture::SP_Particle               SP_Particle;
//     typedef typename TestFixture::Physics_t                 Physics_t;
//     typedef typename TestFixture::Space_Vector              Space_Vector;
//     typedef typename TestFixture::Fission_Site_Container    FSC;

//     Physics_t physics(this->db, this->xs);
//     physics.set_geometry(this->geometry);

//     FSC fsites;

//     // make a particle
//     SP_Particle p(make_shared<Particle>());
//     this->geometry->initialize(Space_Vector(1.1, 0.5, 6.2),
//                                Space_Vector(1.0, 1.0, 1.0),
//                                p->geo_state());
//     p->set_matid(1);
//     p->set_wt(0.6);
//     p->set_rng(this->rng);
//     /*
//      * First 10 random numbers in sequence:
//      * 0.9709
//      * 0.3771
//      * 0.7536
//      * 0.1897
//      * 0.5297
//      * 0.8803
//      * 0.6286
//      * 0.3288
//      * 0.6362
//      * 0.8904
//      */

//     // put particle in group 3
//     p->set_group(3);

//     // sampling fission with these setting should result in 1 fission event
//     EXPECT_EQ(1, physics.sample_fission_site(*p, fsites, 1.04));
//     EXPECT_EQ(1, fsites.size());
//     EXPECT_EQ(1, fsites[0].m);
//     EXPECT_EQ(1.1, fsites[0].r[0]);
//     EXPECT_EQ(0.5, fsites[0].r[1]);
//     EXPECT_EQ(6.2, fsites[0].r[2]);

//     // this next one will fail
//     p->geo_state().d_r = Space_Vector(1.2, 0.3, 6.6);
//     EXPECT_EQ(0, physics.sample_fission_site(*p, fsites, 1.04));
//     EXPECT_EQ(1, fsites.size());
//     EXPECT_EQ(1, fsites[0].m);
//     EXPECT_EQ(1.1, fsites[0].r[0]);
//     EXPECT_EQ(0.5, fsites[0].r[1]);
//     EXPECT_EQ(6.2, fsites[0].r[2]);

//     // this one will pass
//     p->set_group(4);
//     p->set_wt(0.99);
//     EXPECT_EQ(3, physics.sample_fission_site(*p, fsites, 0.2));
//     EXPECT_EQ(4, fsites.size());
//     EXPECT_EQ(1, fsites[0].m);
//     EXPECT_EQ(1.1, fsites[0].r[0]);
//     EXPECT_EQ(0.5, fsites[0].r[1]);
//     EXPECT_EQ(6.2, fsites[0].r[2]);

//     // there are 3 fission sites at this location
//     EXPECT_EQ(1, fsites[1].m);
//     EXPECT_EQ(1.2, fsites[1].r[0]);
//     EXPECT_EQ(0.3, fsites[1].r[1]);
//     EXPECT_EQ(6.6, fsites[1].r[2]);
//     EXPECT_EQ(1, fsites[2].m);
//     EXPECT_EQ(1.2, fsites[2].r[0]);
//     EXPECT_EQ(0.3, fsites[2].r[1]);
//     EXPECT_EQ(6.6, fsites[2].r[2]);
//     EXPECT_EQ(1, fsites[2].m);
//     EXPECT_EQ(1.2, fsites[2].r[0]);
//     EXPECT_EQ(0.3, fsites[2].r[1]);
//     EXPECT_EQ(6.6, fsites[2].r[2]);

//     // test fission spectrum sample
//     {
//         EXPECT_TRUE(physics.initialize_fission(1, *p));
//         EXPECT_EQ(0, p->group());
//     }

//     // test initialization of particle from a fission site
//     EXPECT_TRUE(physics.initialize_fission(fsites[0], *p));
//     EXPECT_EQ(1, p->group());

//     EXPECT_EQ(4, fsites.size());

//     EXPECT_EQ(1.1, physics.fission_site(fsites[0])[0]);
//     EXPECT_EQ(0.5, physics.fission_site(fsites[0])[1]);
//     EXPECT_EQ(6.2, physics.fission_site(fsites[0])[2]);
//     EXPECT_EQ(1.2, physics.fission_site(fsites[1])[0]);
//     EXPECT_EQ(0.3, physics.fission_site(fsites[1])[1]);
//     EXPECT_EQ(6.6, physics.fission_site(fsites[1])[2]);

//     cout << "\n Size of fission-site = " << sizeof(fsites[0]) << " bytes"
//          << endl;

//     EXPECT_EQ(sizeof(fsites[0]), physics.fission_site_bytes());

//     // the size of the fission site container is 4 * 8 bytes (all of the
//     // elements of the struct are aligned along 64-bit boundaries because of
//     // the doubles in the space vector)
//     EXPECT_EQ(32, physics.fission_site_bytes());

//     // test packing of fission-container
//     vector<char> packed;
//     {
//         packed.resize(physics.fission_site_bytes() * 4);
//         memcpy(&packed[0], &fsites[0], packed.size());
//     }

//     {
//         FSC ufsc(4);
//         memcpy(&ufsc[0], &packed[0], packed.size());

//         EXPECT_EQ(4, ufsc.size());
//         EXPECT_EQ(1, ufsc[0].m);
//         EXPECT_EQ(1.1, ufsc[0].r[0]);
//         EXPECT_EQ(0.5, ufsc[0].r[1]);
//         EXPECT_EQ(6.2, ufsc[0].r[2]);
//         EXPECT_EQ(1, ufsc[1].m);
//         EXPECT_EQ(1.2, ufsc[1].r[0]);
//         EXPECT_EQ(0.3, ufsc[1].r[1]);
//         EXPECT_EQ(6.6, ufsc[1].r[2]);
//         EXPECT_EQ(1, ufsc[2].m);
//         EXPECT_EQ(1.2, ufsc[2].r[0]);
//         EXPECT_EQ(0.3, ufsc[2].r[1]);
//         EXPECT_EQ(6.6, ufsc[2].r[2]);
//         EXPECT_EQ(1, ufsc[3].m);
//         EXPECT_EQ(1.2, ufsc[3].r[0]);
//         EXPECT_EQ(0.3, ufsc[3].r[1]);
//         EXPECT_EQ(6.6, ufsc[3].r[2]);
//     }

//     // test null ops for no fission
//     p->set_matid(0);
//     EXPECT_FALSE(physics.initialize_fission(0, *p));
//     EXPECT_EQ(0, physics.sample_fission_site(*p, fsites, 0.1));
// }

//---------------------------------------------------------------------------//
//                 end of tstPhysics.cc
//---------------------------------------------------------------------------//
