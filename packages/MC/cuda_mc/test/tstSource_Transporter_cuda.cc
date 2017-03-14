//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/test/tstSource_Transporter_cuda.cc
 * \author Stuart Slattery
 * \brief  Source tranpsporter test.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../Source_Transporter.hh"
#include "../Cell_Tally.hh"
#include "../Keff_Tally.hh"
#include "../VR_Roulette.hh"
#include "../Uniform_Source.hh"
#include "../Fission_Source.hh"

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
#include "mc/Global_RNG.hh"

#include "Physics_Tester.hh"

using namespace std;

using def::X; using def::Y; using def::Z;

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//

template <class Geometry>
class SourceTransporterTest : public testing::Test
{
  protected:
    typedef Geometry                                      Geometry_t;
    typedef profugus::RNG_Control                         RNG_Control;
    typedef cuda_profugus::Source_Transporter<Geometry>   Source_Transporter_t;
    typedef cuda_profugus::Physics<Geometry>              Physics_t;
    typedef cuda_profugus::Bank<Geometry>                 Bank_t;
    typedef cuda_profugus::Tallier<Geometry>              Tallier_t;
    typedef cuda_profugus::Box_Shape                      Shape;
    typedef cuda_profugus::Uniform_Source<Geometry,Shape> Uniform_Source_t;
    typedef cuda_profugus::Fission_Source<Geometry>       Fission_Source_t;
    typedef cuda_profugus::Cell_Tally<Geometry>           Cell_Tally_t;
    typedef cuda_profugus::Keff_Tally<Geometry>           Keff_Tally_t;
    typedef cuda_profugus::VR_Roulette<Geometry>          VR_Roulette_t;
    typedef typename Physics_t::ParameterList_t           ParameterList_t;
    typedef typename Physics_t::XS_t                      XS_t;
    typedef typename Physics_t::Fission_Site              Fission_Site;
    typedef typename Physics_t::Fission_Site_Container    Fission_Site_Container;
    typedef typename Geometry_t::Space_Vector             Space_Vector;

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
        db->set("particle_vector_size", 1000);
        db->set("Np", 10000);
        db->set("cuda_mc_diag_frac", 0.1);
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
void SourceTransporterTest<cuda_profugus::Mesh_Geometry>::build_geometry()
{
    edges = {0.0, 20.0, 40.0, 60.0, 80.0, 100.0};
}

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

typedef ::testing::Types<cuda_profugus::Mesh_Geometry> MyTypes;
TYPED_TEST_CASE(SourceTransporterTest, MyTypes);

//---------------------------------------------------------------------------//
TYPED_TEST(SourceTransporterTest, fixed_source)
{
    typedef typename TestFixture::Fission_Source_t          Fission_Source_t;
    typedef typename TestFixture::Uniform_Source_t          Uniform_Source_t;
    typedef typename TestFixture::Tallier_t                 Tallier_t;
    typedef typename TestFixture::Cell_Tally_t              Cell_Tally_t;
    typedef typename TestFixture::VR_Roulette_t             VR_Roulette_t;
    typedef typename TestFixture::Physics_t                 Physics_t;
    typedef typename TestFixture::Bank_t                    Bank_t;
    typedef typename TestFixture::Source_Transporter_t      Source_Transporter_t;
    typedef typename TestFixture::Space_Vector              Space_Vector;
    typedef typename TestFixture::Fission_Site_Container    FSC;

    // turn off implicit capture
    this->db->set("implicit_capture", false);

    // make the physics tester.
    int vector_size = 1000;
    Physics_Tester physics_tester( this->edges, this->edges, this->edges,
				   vector_size, this->rng, *(this->db), *(this->xs), 1 );

    // make a roulette variance reduction
    auto vr = std::make_shared<VR_Roulette_t>( *(this->db) );

    // make a tallier
    std::shared_ptr<Tallier_t> tallier = std::make_shared<Tallier_t>();

    // add a cell tally
    int num_batch = 2;
    auto cell_tally =
        std::make_shared<Cell_Tally_t>( this->db, physics_tester.geometry(),
                                        num_batch );
    tallier->add_pathlength_tally( cell_tally );
    tallier->build();

    // create a source.
    auto shape = physics_tester.source_shape();
    auto source = std::make_shared<Uniform_Source_t>( this->db,
                                                      physics_tester.geometry(),
                                                      this->xs->num_groups(),
                                                      num_batch );
    source->build_source( shape );

    // create a source transporter
    profugus::Global_RNG::d_rng = this->rng;
    Source_Transporter_t transporter(
        this->db, physics_tester.geometry(), physics_tester.physics() );
    transporter.set( vr );
    transporter.set( tallier );
    transporter.assign_source( source );

    // transport the source
    transporter.solve();

    // check tally results
    cell_tally->finalize( this->db->get<int>("Np") );
    Teuchos::Array<double> first_moment, second_moment;
    cell_tally->copy_moments_to_host( first_moment, second_moment );
    std::cout << first_moment << std::endl;
    std::cout << second_moment << std::endl;
}

//---------------------------------------------------------------------------//
TYPED_TEST(SourceTransporterTest, fission_source)
{
    typedef typename TestFixture::Fission_Source_t          Fission_Source_t;
    typedef typename TestFixture::Uniform_Source_t          Uniform_Source_t;
    typedef typename TestFixture::Tallier_t                 Tallier_t;
    typedef typename TestFixture::Keff_Tally_t              Keff_Tally_t;
    typedef typename TestFixture::VR_Roulette_t             VR_Roulette_t;
    typedef typename TestFixture::Physics_t                 Physics_t;
    typedef typename TestFixture::Bank_t                    Bank_t;
    typedef typename TestFixture::Source_Transporter_t      Source_Transporter_t;
    typedef typename TestFixture::Space_Vector              Space_Vector;
    typedef typename TestFixture::Fission_Site_Container    FSC;

    // turn off implicit capture
    this->db->set("implicit_capture", false);

    // make the physics tester.
    int vector_size = 1000;
    Physics_Tester physics_tester( this->edges, this->edges, this->edges,
				   vector_size, this->rng, *(this->db), *(this->xs), 1 );

    // make a roulette variance reduction
    auto vr = std::make_shared<VR_Roulette_t>( *(this->db) );

    // make a tallier
    std::shared_ptr<Tallier_t> tallier = std::make_shared<Tallier_t>();

    // add a keff tally
    double keff_init = 1.01;
    auto keff_tally =
        std::make_shared<Keff_Tally_t>( keff_init, physics_tester.physics(), vector_size);
    keff_tally->begin_active_cycles();
    tallier->add_pathlength_tally( keff_tally );
    tallier->build();

    // create a source.
    auto shape = physics_tester.source_shape();
    auto source = std::make_shared<Fission_Source_t>( this->db,
                                                      physics_tester.geometry(),
                                                      physics_tester.physics() );
    Teuchos::Array<double> density( physics_tester.geometry().get_host_ptr()->num_cells(),
                                    10*this->db->get<int>("Np") );
    Teuchos::ArrayView<const double> dens_view = density();
    source->build_initial_source( physics_tester.cart_mesh(), dens_view );
    auto fsites = source->create_fission_site_container();

    // create a source transporter
    profugus::Global_RNG::d_rng = this->rng;
    Source_Transporter_t transporter(
        this->db, physics_tester.geometry(), physics_tester.physics() );
    transporter.set( vr );
    transporter.set( tallier );
    transporter.assign_source( source );
    transporter.sample_fission_sites( fsites, keff_init );

    // transport the source
    transporter.solve();

    // check tally results
    keff_tally->end_cycle( this->db->get<int>("Np") );
    Teuchos::Array<double> first_moment, second_moment;
    std::cout << "Keff 1st Moment: " << keff_tally->keff_sum() << std::endl;
    std::cout << "Keff 2nd Moment: " << keff_tally->keff_sum_sq() << std::endl;

    // Look at fission sites.
    for ( auto fs : *fsites )
    {
        std::cout << "Fsite: " << fs.m << ", "
                  << fs.r.x << " "
                  << fs.r.y << " "
                  << fs.r.z << std::endl;
    }
}

//---------------------------------------------------------------------------//
//                 end of tstSource_Transporter_cuda.cc
//---------------------------------------------------------------------------//
