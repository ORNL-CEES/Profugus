//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/test/tstDomain_Transporter_cuda.cc
 * \author Stuart Slattery
 * \brief  Domain tranpsporter test.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../Domain_Transporter.hh"
#include "../Cell_Tally.hh"
#include "../VR_Roulette.hh"

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
class DomainTransporterTest : public testing::Test
{
  protected:
    typedef Geometry                                    Geometry_t;
    typedef profugus::RNG_Control                       RNG_Control;
    typedef cuda_profugus::Domain_Transporter<Geometry> Domain_Transporter_t;
    typedef cuda_profugus::Physics<Geometry>            Physics_t;
    typedef cuda_profugus::Bank<Geometry>               Bank_t;
    typedef cuda_profugus::Tallier<Geometry>            Tallier_t;
    typedef cuda_profugus::Cell_Tally<Geometry>         Cell_Tally_t;
    typedef cuda_profugus::VR_Roulette<Geometry>        VR_Roulette_t;
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
        db->set("weight_cutoff", 0.001);
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
void DomainTransporterTest<cuda_profugus::Mesh_Geometry>::build_geometry()
{
    edges = {0.0, 20.0, 40.0, 60.0, 80.0, 100.0};
}

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

typedef ::testing::Types<cuda_profugus::Mesh_Geometry> MyTypes;
TYPED_TEST_CASE(DomainTransporterTest, MyTypes);

//---------------------------------------------------------------------------//

TYPED_TEST(DomainTransporterTest, take_step_roulette)
{
    typedef typename TestFixture::Tallier_t                 Tallier_t;
    typedef typename TestFixture::Cell_Tally_t              Cell_Tally_t;
    typedef typename TestFixture::VR_Roulette_t             VR_Roulette_t;
    typedef typename TestFixture::Physics_t                 Physics_t;
    typedef typename TestFixture::Bank_t                    Bank_t;
    typedef typename TestFixture::Domain_Transporter_t      Domain_Transporter_t;
    typedef typename TestFixture::Space_Vector              Space_Vector;
    typedef typename TestFixture::Fission_Site_Container    FSC;

    // turn off implicit capture
    this->db->set("implicit_capture", false);

    // make the physics tester.
    int Np = 15;
    Physics_Tester physics_tester( this->edges, this->edges, this->edges,
				   Np, this->rng, *(this->db), *(this->xs), 0 );
    
    // initialize the particles
    Space_Vector r, d;
    r.x = 50.0;
    r.y = 50.0;
    r.z = 50.0;
    d.x = 1.0;
    d.y = 1.0;
    d.z = 1.0;
    physics_tester.geometry_initialize( r, d, 0 );
    std::vector<double> cdf = {0.2, 0.4, 0.6, 0.8, 1.0};
    physics_tester.sample_group( cdf );
    physics_tester.particle_tester().live();
    physics_tester.particle_tester().set_batch( 0 );
    Teuchos::Array<cuda_profugus::events::Event> events( 
        Np, cuda_profugus::events::TAKE_STEP );
    physics_tester.particle_tester().set_event( events );
    physics_tester.particle_tester().sort_by_event();

    // make a tallier
    std::shared_ptr<Tallier_t> tallier = std::make_shared<Tallier_t>();    

    // add a cell tally
    int num_batch = 1;
    std::shared_ptr<Cell_Tally_t> cell_tally = 
        std::make_shared<Cell_Tally_t>( physics_tester.geometry(), num_batch );
    tallier->add_pathlength_tally( cell_tally );
    tallier->build();
    
    // make a roulette variance reduction
    std::shared_ptr<VR_Roulette_t> vr = std::make_shared<VR_Roulette_t>( *(this->db) );

    // create a domain transporter
    Domain_Transporter_t transporter;
    transporter.set( physics_tester.geometry(), physics_tester.physics() );
    transporter.set( vr );
    transporter.set( tallier );

    // Move the particles a step.
    cuda::Shared_Device_Ptr<Bank_t> bank;
    transporter.transport_step( physics_tester.particles(), bank );

    // Sort the particles.
    physics_tester.particles().get_host_ptr()->sort_by_event(
        physics_tester.particles().get_host_ptr()->size());

    // Check the we have only collisions and boundaries.
    events = physics_tester.particle_tester().event();
    for ( auto e : events )
    {
        EXPECT_TRUE( cuda_profugus::events::COLLISION == e ||
                     cuda_profugus::events::BOUNDARY == e );
    }

    // Post-process the step.
    transporter.process_step( physics_tester.particles(), bank );

    // Sort the particles.
    physics_tester.particles().get_host_ptr()->sort_by_event(
        physics_tester.particles().get_host_ptr()->size());

    // Check that we have only dead particles and particles that will take
    // another step.
    events = physics_tester.particle_tester().event();
    Teuchos::Array<int> live = physics_tester.particle_tester().alive();
    for ( int i = 0; i < Np; ++i )
    {
        if ( live[i] )
        {
            EXPECT_EQ( events[i], cuda_profugus::events::TAKE_STEP );
        }
        else
        {
            EXPECT_EQ( events[i], cuda_profugus::events::DEAD );
        }
    }
}

//---------------------------------------------------------------------------//
//                 end of tstDomain_Transporter_cuda.cc
//---------------------------------------------------------------------------//
