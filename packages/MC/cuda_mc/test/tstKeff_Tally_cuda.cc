//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/test/tstKeff_Tally.cc
 * \author Thomas M. Evans
 * \date   Thu May 15 14:50:01 2014
 * \brief  Keff_Tally unit-test.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../Keff_Tally.hh"
#include "../Physics.hh"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "gtest/utils_gtest.hh"

#include "utils/Definitions.hh"
#include "rng/RNG_Control.hh"
#include "cuda_utils/Shared_Device_Ptr.hh"
#include "cuda_geometry/Mesh_Geometry.hh"

#include "Physics_Tester.hh"

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//

class Keff_TallyTest : public testing::Test
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
        node  = profugus::node();
        nodes = profugus::nodes();

        int seed = 342412;
        profugus::RNG_Control control(seed);
        rng = control.rng();

        db = Teuchos::rcp(new ParameterList_t("test"));
        build_physics();
    }

    void build_physics()
    {
        const int ng = 3;

        auto xs = Teuchos::rcp(new XS_t());
        xs->set(0, ng);

        // make group boundaries
        XS_t::OneDArray nbnd(4, 0.0);
        nbnd[0] = 100.0; nbnd[1] = 1.0; nbnd[2] = 0.01; nbnd[3] = 0.0001;
        xs->set_bounds(nbnd);

        double t1[] = {1.1, 1.6, 2.9};
        double t2[] = {10.0, 11.3, 16.2};

        XS_t::OneDArray tot1(std::begin(t1), std::end(t1));
        XS_t::OneDArray tot2(std::begin(t2), std::end(t2));

        xs->add(0, XS_t::TOTAL, tot1);
        xs->add(1, XS_t::TOTAL, tot2);

        double s1[][3] = {{0.7, 0.0, 0.0},
                          {0.2, 0.3, 0.0},
                          {0.1, 0.7, 1.9}};

        double s2[][3] = {{2.7, 0.0, 0.0},
                          {2.2, 2.3, 0.0},
                          {2.1, 2.7, 3.9}};

        XS_t::TwoDArray sct1(ng, ng, 0.0);
        XS_t::TwoDArray sct2(ng, ng, 0.0);

        for (int g = 0; g < 3; ++g)
        {
            for (int gp = 0; gp < 3; ++gp)
            {
                sct1(g, gp) = s1[g][gp];
                sct2(g, gp) = s2[g][gp];
            }
        }

        xs->add(0, 0, sct1);
        xs->add(1, 0, sct2);

        double c2[] = {0.4, 0.6, 0.0};
        double f2[] = {3.2, 4.2, 0.0};
        double n2[] = {2.4*3.2, 2.4*4.2, 0.0};

        XS_t::OneDArray chi2(std::begin(c2), std::end(c2));
        XS_t::OneDArray fis2(std::begin(f2), std::end(f2));
        XS_t::OneDArray nuf2(std::begin(n2), std::end(n2));

        xs->add(1, XS_t::CHI, chi2);
        xs->add(1, XS_t::SIG_F, fis2);
        xs->add(1, XS_t::NU_SIG_F, nuf2);

        xs->complete();

	int num_particle = 1000;
	std::vector<double> edges = {0.0,1.0};
	physics_tester = Teuchos::rcp( 
	    new Physics_Tester(edges,edges,edges,num_particle,rng,*db,*xs) );
    }

  protected:

    profugus::RNG_Control::RNG_t rng;
    RCP_Std_DB  db;
    Teuchos::RCP<Physics_Tester>  physics_tester;

    int node;
    int nodes;
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST_F(Keff_TallyTest, test_keff)
{
    // Make a keff tally
    Tally_t keff(1.0, physics_tester->physics());

    // Set reference k
    double ref_k = 0.0;

    keff.begin_cycle();
    EXPECT_EQ(0.0, keff.latest());

    // Tally some particles
    physics_tester->particle_tester()->set_wt( 0.65 );
    physics_tester->particle_tester()->set_matid( 1 );
    physics_tester->particle_tester()->set_group( 1 );
    physics_tester->particle_tester()->set_step( 1.4 );
    keff.accumulate( physics_tester->particles() );
    ref_k = 1.4 * 0.65 * 2.4 * 4.2;
    EXPECT_SOFTEQ(ref_k, keff.latest(), 1.0e-6);

    // Put particle in group 2
    physics_tester->particle_tester()->set_group( 2 ); // no fission in group 2
    physics_tester->particle_tester()->set_step( 1.8 );
    keff.accumulate( physics_tester->particles() );
    EXPECT_SOFTEQ(ref_k, keff.latest(), 1.0e-6);

    // Put particle in group 0
    physics_tester->particle_tester()->set_wt( 0.4 );
    physics_tester->particle_tester()->set_group( 0 );
    physics_tester->particle_tester()->set_step( 0.6 );
    keff.accumulate( physics_tester->particles() );
    ref_k += 0.4 * 0.6 * 2.4 * 3.2;
    EXPECT_SOFTEQ(ref_k, keff.latest(), 1.0e-6);

    // Put particle in material 0
    physics_tester->particle_tester()->set_wt( 0.45 );
    physics_tester->particle_tester()->set_matid( 0 );
    physics_tester->particle_tester()->set_step( 0.2 );
    keff.accumulate( physics_tester->particles() );
    EXPECT_SOFTEQ(ref_k, keff.latest(), 1.0e-6);

    // End the cycle
    int num_particle = 3.0 * physics_tester->particle_tester()->size();
    keff.end_cycle(num_particle);

    EXPECT_SOFTEQ(nodes * ref_k / num_particle, keff.latest(), 1.0e-6);

    // >>> TEST THE RESET METHOD
    keff.reset();

    // latest keff should stay the same
    EXPECT_SOFTEQ(nodes * ref_k / num_particle, keff.latest(), 1.0e-6);

    // path lengths and counter should be reset
    EXPECT_EQ(0., keff.keff_sum());
    EXPECT_EQ(0., keff.keff_sum_sq());
    EXPECT_EQ(0, keff.cycle_count());
}

//---------------------------------------------------------------------------//

TEST_F(Keff_TallyTest, test_keff_var)
{
    // Make a keff tally
    Tally_t keff(1.0, physics_tester->physics());
    int num_particle = 3.0 * physics_tester->particle_tester()->size();

    // Create reference variables
    double ref_k     = 0.0;
    double ref_kavg  = 0.0;
    double ref_kmom1 = 0.0;
    double ref_kmom2 = 0.0;
    double ref_kvar  = 0.0;

    /*** Begin INACTIVE CYCLE 1 ***/

    keff.begin_cycle();
    EXPECT_EQ(0.0, keff.latest());

    // Tally some particles
    physics_tester->particle_tester()->set_wt( 0.65 );
    physics_tester->particle_tester()->set_matid( 1 );
    physics_tester->particle_tester()->set_group( 1 );
    physics_tester->particle_tester()->set_step( 1.4 );
    keff.accumulate( physics_tester->particles() );
    ref_k = 1.4 * 0.65 * 2.4 * 4.2;
    EXPECT_SOFTEQ(ref_k, keff.latest(), 1.0e-6);

    // Put particle in group 2
    physics_tester->particle_tester()->set_group( 2 ); // no fission in group 2
    physics_tester->particle_tester()->set_step( 1.8 );
    keff.accumulate( physics_tester->particles() );
    EXPECT_SOFTEQ(ref_k, keff.latest(), 1.0e-6);

    // Put particle in group 0
    physics_tester->particle_tester()->set_wt( 0.4 );
    physics_tester->particle_tester()->set_group( 0 );
    physics_tester->particle_tester()->set_step( 0.6 );
    keff.accumulate( physics_tester->particles() );
    ref_k += 0.4 * 0.6 * 2.4 * 3.2;
    EXPECT_SOFTEQ(ref_k, keff.latest(), 1.0e-6);

    // Put particle in material 0
    physics_tester->particle_tester()->set_wt( 0.45 );
    physics_tester->particle_tester()->set_matid( 0 );
    physics_tester->particle_tester()->set_step( 0.2 );
    keff.accumulate( physics_tester->particles() );
    EXPECT_SOFTEQ(ref_k, keff.latest(), 1.0e-6);

    // End the cycle
    keff.end_cycle(num_particle);
    ref_k = nodes * ref_k / num_particle;
    EXPECT_SOFTEQ(ref_k, keff.latest(), 1.0e-6);

    /*** Begin INACTIVE CYCLE 2 ***/

    keff.begin_cycle();
    EXPECT_EQ(0.0, keff.latest());

    // Tally more particles
    physics_tester->particle_tester()->set_wt( 0.55 );
    physics_tester->particle_tester()->set_matid( 1 );
    physics_tester->particle_tester()->set_group( 1 );
    physics_tester->particle_tester()->set_step( 1.4 );
    keff.accumulate( physics_tester->particles() );
    ref_k = 1.4 * 0.55 * 2.4 * 4.2;
    EXPECT_SOFTEQ(ref_k, keff.latest(), 1.0e-6);

    // Put particle in group 2
    physics_tester->particle_tester()->set_group( 2 ); // no fission in group 2
    physics_tester->particle_tester()->set_step( 1.8 );
    keff.accumulate( physics_tester->particles() );
    EXPECT_SOFTEQ(ref_k, keff.latest(), 1.0e-6);

    // Put particle in group 0
    physics_tester->particle_tester()->set_wt( 0.2 );
    physics_tester->particle_tester()->set_group( 0 );
    physics_tester->particle_tester()->set_step( 0.6 );
    keff.accumulate( physics_tester->particles() );
    ref_k += 0.2 * 0.6 * 2.4 * 3.2;
    EXPECT_SOFTEQ(ref_k, keff.latest(), 1.0e-6);

    // Put particle in material 0
    physics_tester->particle_tester()->set_wt( 0.55 );
    physics_tester->particle_tester()->set_matid( 0 );
    physics_tester->particle_tester()->set_step( 0.2 );
    keff.accumulate( physics_tester->particles() );
    EXPECT_SOFTEQ(ref_k, keff.latest(), 1.0e-6);

    // End the cycle
    keff.end_cycle(num_particle);
    ref_k = nodes * ref_k / num_particle;
    EXPECT_SOFTEQ(ref_k, keff.latest(), 1.0e-6);

    /*** Begin ACTIVE CYCLE 1 ***/

    keff.begin_active_cycles();
    keff.begin_cycle();
    EXPECT_EQ(0.0, keff.latest());
    EXPECT_EQ(0.0, keff.keff_sum());
    EXPECT_EQ(0.0, keff.keff_sum_sq());

    // Tally some particles
    physics_tester->particle_tester()->set_wt( 0.65 );
    physics_tester->particle_tester()->set_matid( 1 );
    physics_tester->particle_tester()->set_group( 1 );
    physics_tester->particle_tester()->set_step( 1.4 );
    keff.accumulate( physics_tester->particles() );
    keff.accumulate(1.4, *p);
    ref_k = 1.4 * 0.65 * 2.4 * 4.2;
    EXPECT_SOFTEQ(ref_k, keff.latest(), 1.0e-6);

    // Put particle in group 2
    physics_tester->particle_tester()->set_group( 2 ); // no fission in group 2
    physics_tester->particle_tester()->set_step( 1.8 );
    keff.accumulate( physics_tester->particles() );
    EXPECT_SOFTEQ(ref_k, keff.latest(), 1.0e-6);

    // Put particle in group 0
    physics_tester->particle_tester()->set_wt( 0.4 );
    physics_tester->particle_tester()->set_group( 0 );
    physics_tester->particle_tester()->set_step( 0.6 );
    keff.accumulate( physics_tester->particles() );
    ref_k += 0.4 * 0.6 * 2.4 * 3.2;
    EXPECT_SOFTEQ(ref_k, keff.latest(), 1.0e-6);

    // Put particle in material 0
    physics_tester->particle_tester()->set_wt( 0.55 );
    physics_tester->particle_tester()->set_matid( 0 );
    physics_tester->particle_tester()->set_step( 0.2 );
    keff.accumulate( physics_tester->particles() );
    EXPECT_SOFTEQ(ref_k, keff.latest(), 1.0e-6);

    // End the cycle
    keff.end_cycle(num_particle);

    ref_k = nodes * ref_k / num_particle;
    ref_kmom1 += ref_k;
    ref_kmom2 += ref_k * ref_k;
    ref_kavg = ref_kmom1 / 1.0;

    EXPECT_SOFTEQ(ref_k, keff.latest(), 1.0e-6);
    EXPECT_SOFTEQ(ref_kavg, keff.mean(), 1.0e-6);

    /*** Begin ACTIVE CYCLE 2 ***/

    keff.begin_cycle();
    EXPECT_EQ(0.0, keff.latest());
    EXPECT_NE(0.0, keff.keff_sum());
    EXPECT_NE(0.0, keff.keff_sum_sq());

    // Tally some more particles
    physics_tester->particle_tester()->set_wt( 0.55 );
    physics_tester->particle_tester()->set_matid( 1 );
    physics_tester->particle_tester()->set_group( 1 );
    physics_tester->particle_tester()->set_step( 1.4 );
    keff.accumulate( physics_tester->particles() );
    ref_k = 1.4 * 0.55 * 2.4 * 4.2;
    EXPECT_SOFTEQ(ref_k, keff.latest(), 1.0e-6);

    // Put particle in group 2
    physics_tester->particle_tester()->set_group( 2 ); // no fission in group 2
    physics_tester->particle_tester()->set_step( 1.8 );
    keff.accumulate( physics_tester->particles() );
    EXPECT_SOFTEQ(ref_k, keff.latest(), 1.0e-6);

    // Put particle in group 0
    physics_tester->particle_tester()->set_wt( 0.2 );
    physics_tester->particle_tester()->set_group( 0 );
    physics_tester->particle_tester()->set_step( 0.6 );
    keff.accumulate( physics_tester->particles() );
    ref_k += 0.2 * 0.6 * 2.4 * 3.2;
    EXPECT_SOFTEQ(ref_k, keff.latest(), 1.0e-6);

    // Put particle in material 0
    physics_tester->particle_tester()->set_wt( 0.55 );
    physics_tester->particle_tester()->set_matid( 0 );
    physics_tester->particle_tester()->set_step( 0.2 );
    keff.accumulate( physics_tester->particles() );
    EXPECT_SOFTEQ(ref_k, keff.latest(), 1.0e-6);

    // End the cycle
    keff.end_cycle(num_particle);

    ref_k = nodes * ref_k / num_particle;
    ref_kmom1 += ref_k;
    ref_kmom2 += ref_k * ref_k;
    ref_kavg = ref_kmom1 / 2.0;
    ref_kvar = (ref_kmom2 - 2.0 * ref_kavg * ref_kavg) / 2.0 / 1.0;

    EXPECT_SOFTEQ(ref_k, keff.latest(), 1.0e-6);
    EXPECT_SOFTEQ(ref_kavg, keff.mean(), 1.0e-6);
    EXPECT_SOFTEQ(ref_kvar, keff.variance(), 1.0e-6);

    /*** Begin ACTIVE CYCLE 3 ***/

    keff.begin_cycle();
    EXPECT_EQ(0.0, keff.latest());
    EXPECT_NE(0.0, keff.keff_sum());
    EXPECT_NE(0.0, keff.keff_sum_sq());

    // Tally some more particles
    physics_tester->particle_tester()->set_wt( 0.8 );
    physics_tester->particle_tester()->set_matid( 1 );
    physics_tester->particle_tester()->set_group( 1 );
    physics_tester->particle_tester()->set_step( 1.4 );
    keff.accumulate( physics_tester->particles() );
    ref_k = 1.4 * 0.8 * 2.4 * 4.2;
    EXPECT_SOFTEQ(ref_k, keff.latest(), 1.0e-6);

    // Put particle in group 2
    physics_tester->particle_tester()->set_group( 2 ); // no fission in group 2
    physics_tester->particle_tester()->set_step( 1.8 );
    keff.accumulate( physics_tester->particles() );
    EXPECT_SOFTEQ(ref_k, keff.latest(), 1.0e-6);

    // Put particle in group 0
    physics_tester->particle_tester()->set_wt( 0.65 );
    physics_tester->particle_tester()->set_group( 0 );
    physics_tester->particle_tester()->set_step( 0.6 );
    keff.accumulate( physics_tester->particles() );
    ref_k += 0.65 * 0.6 * 2.4 * 3.2;
    EXPECT_SOFTEQ(ref_k, keff.latest(), 1.0e-6);

    // Put particle in material 0
    physics_tester->particle_tester()->set_wt( 0.55 );
    physics_tester->particle_tester()->set_matid( 0 );
    physics_tester->particle_tester()->set_step( 0.2 );
    keff.accumulate( physics_tester->particles() );
    EXPECT_SOFTEQ(ref_k, keff.latest(), 1.0e-6);

    // End the cycle
    keff.end_cycle(num_particle);

    ref_k = nodes * ref_k / num_particle;
    ref_kmom1 += ref_k;
    ref_kmom2 += ref_k * ref_k;
    ref_kavg = ref_kmom1 / num_particle;
    ref_kvar = (ref_kmom2 - num_particle * ref_kavg * ref_kavg) / num_particle / 2.0;

    EXPECT_SOFTEQ(ref_k, keff.latest(), 1.0e-6);
    EXPECT_SOFTEQ(ref_kavg, keff.mean(), 1.0e-6);
    EXPECT_SOFTEQ(ref_kvar, keff.variance(), 1.0e-6);
}

//---------------------------------------------------------------------------//
//                 end of tstKeff_Tally.cc
//---------------------------------------------------------------------------//
