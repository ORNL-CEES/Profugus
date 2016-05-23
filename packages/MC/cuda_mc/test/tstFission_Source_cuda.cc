//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/test/tstFission_Source_cuda.cc
 * \author Stuart Slattery
 * \date   Tue May 06 11:54:26 2014
 * \brief  Fission_Source unit test.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <memory>

#include "Teuchos_Array.hpp"
#include "Teuchos_RCP.hpp"

#include "../Fission_Source.hh"
#include "cuda_geometry/Cartesian_Mesh.hh"
#include "utils/Definitions.hh"
#include "rng/RNG_Control.hh"
#include "mc/Global_RNG.hh"
#include "cuda_utils/Shared_Device_Ptr.hh"
#include "cuda_geometry/Mesh_Geometry.hh"

#include "gtest/utils_gtest.hh"
#include "Physics_Tester.hh"

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//

template <class Geometry>
class FissionSourceTest : public testing::Test
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
    typedef cuda_profugus::Fission_Source<Geometry>     Fission_Source;

  protected:

    void SetUp()
    {
        build_geometry();

        build_xs();

        // make a rng
        int seed = 342412;
        profugus::RNG_Control control(seed);
        rng = control.rng();
	profugus::Global_RNG::d_rng = rng;

        // make db
	vector_size = 100;
	np = 800;
        db = Teuchos::rcp(new ParameterList_t("test"));
	db->set( "Np", np );
        db->set( "particle_vector_size", vector_size );

	// build physics.
	physics_tester = std::make_shared<Physics_Tester>(
	    edges, edges, edges, vector_size, rng, *db, *xs, 1 );
    }

    void build_geometry();

    void build_xs()
    {
        xs = Teuchos::rcp(new XS_t());
        xs->set(0, 5);

	std::vector<double> bnd(6, 0.0);
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
        typename XS_t::OneDArray fission(std::begin(f), std::end(f));
        typename XS_t::OneDArray chi(std::begin(c), std::end(c));
        typename XS_t::OneDArray nus(std::begin(n), std::end(n));
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
    int np;
    int vector_size;
    def::Vec_Dbl edges;
    Teuchos::RCP<XS_t> xs;
    Teuchos::RCP<ParameterList_t> db;
    profugus::RNG_Control::RNG_t rng;
    std::shared_ptr<Physics_Tester> physics_tester;
};

template <>
void FissionSourceTest<cuda_profugus::Mesh_Geometry>::build_geometry()
{
    edges = {0.0, 100.0};
}

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

typedef ::testing::Types<cuda_profugus::Mesh_Geometry> MyTypes;
TYPED_TEST_CASE(FissionSourceTest, MyTypes);

//---------------------------------------------------------------------------//
// sample the initial source using the geometry
TYPED_TEST(FissionSourceTest, sample_geometry)
{
    typedef typename TestFixture::Fission_Source Fission_Source;

    // Create the fission source.
    Fission_Source source( this->db, 
			   this->physics_tester->geometry(),
			   this->physics_tester->physics() );

    EXPECT_EQ( this->np, source.Np() );
    EXPECT_EQ( this->edges[0], source.lower_coords().x );
    EXPECT_EQ( this->edges[0], source.lower_coords().y );
    EXPECT_EQ( this->edges[0], source.lower_coords().z );
    EXPECT_EQ( this->edges[1], source.width().x );
    EXPECT_EQ( this->edges[1], source.width().y );
    EXPECT_EQ( this->edges[1], source.width().z );

    // Initialize with no mesh - sample the geometry.
    source.build_initial_source();
    EXPECT_EQ( this->np / profugus::nodes(), source.num_to_transport() );
    EXPECT_EQ( this->np, source.total_num_to_transport() );
    EXPECT_EQ( this->np / profugus::nodes(), source.num_left() );
    EXPECT_EQ( 0, source.num_run() );
    EXPECT_FALSE( source.empty() );
    EXPECT_TRUE( source.is_initial_source() );

    // Set all the particles to dead.
    Teuchos::Array<cuda_profugus::events::Event> dead_event( 
	this->vector_size, cuda_profugus::events::DEAD );
    this->physics_tester->particle_tester().set_event( dead_event );

    // Get particles from the source.
    source.get_particles( this->physics_tester->particles() );
    EXPECT_EQ( source.num_to_transport() - this->vector_size, source.num_left() );
    EXPECT_EQ( this->vector_size, source.num_run() );
    
    // check particles
    Teuchos::Array<cuda_profugus::events::Event> events =
    	this->physics_tester->particle_tester().event();
    Teuchos::Array<int> matids = 
    	this->physics_tester->particle_tester().matid();
    Teuchos::Array<int> alive = 
    	this->physics_tester->particle_tester().alive();
    Teuchos::Array<double> wts = 
    	this->physics_tester->particle_tester().wt();
    for ( int i = 0; i < this->vector_size; ++i )
    {
    	EXPECT_EQ( cuda_profugus::events::TAKE_STEP, events[i] );
    	EXPECT_EQ( 1, matids[i] );
	EXPECT_TRUE( alive[i] );
    	EXPECT_EQ( 1.0, wts[i] );	
    }
}

//---------------------------------------------------------------------------//
// sample the initial source using the mesh
TYPED_TEST(FissionSourceTest, sample_mesh)
{
    typedef typename TestFixture::Fission_Source Fission_Source;

    // Create the fission source.
    Fission_Source source( this->db, 
			   this->physics_tester->geometry(),
			   this->physics_tester->physics() );

    EXPECT_EQ( this->np, source.Np() );
    EXPECT_EQ( this->edges[0], source.lower_coords().x );
    EXPECT_EQ( this->edges[0], source.lower_coords().y );
    EXPECT_EQ( this->edges[0], source.lower_coords().z );
    EXPECT_EQ( this->edges[1], source.width().x );
    EXPECT_EQ( this->edges[1], source.width().y );
    EXPECT_EQ( this->edges[1], source.width().z );

    // Initialize with mesh and fission density.
    Teuchos::Array<double> density( 1, this->np );
    Teuchos::ArrayView<const double> dens_view = density();
    source.build_initial_source( this->physics_tester->cart_mesh(), dens_view );
    EXPECT_EQ( this->np / profugus::nodes(), source.num_to_transport() );
    EXPECT_EQ( this->np, source.total_num_to_transport() );
    EXPECT_EQ( this->np / profugus::nodes(), source.num_left() );
    EXPECT_EQ( 0, source.num_run() );
    EXPECT_FALSE( source.empty() );
    EXPECT_TRUE( source.is_initial_source() );

    // Set all the particles to dead.
    Teuchos::Array<cuda_profugus::events::Event> dead_event( 
	this->vector_size, cuda_profugus::events::DEAD );
    this->physics_tester->particle_tester().set_event( dead_event );

    // Get particles from the source.
    source.get_particles( this->physics_tester->particles() );
    EXPECT_EQ( source.num_to_transport() - this->vector_size, source.num_left() );
    EXPECT_EQ( this->vector_size, source.num_run() );
    
    // check particles
    Teuchos::Array<cuda_profugus::events::Event> events =
    	this->physics_tester->particle_tester().event();
    Teuchos::Array<int> matids = 
    	this->physics_tester->particle_tester().matid();
    Teuchos::Array<int> alive = 
    	this->physics_tester->particle_tester().alive();
    Teuchos::Array<double> wts = 
    	this->physics_tester->particle_tester().wt();
    for ( int i = 0; i < this->vector_size; ++i )
    {
    	EXPECT_EQ( cuda_profugus::events::TAKE_STEP, events[i] );
    	EXPECT_EQ( 1, matids[i] );
	EXPECT_TRUE( alive[i] );
    	EXPECT_EQ( 1.0, wts[i] );	
    }
}

//---------------------------------------------------------------------------//
// sample the source with fission sites
TYPED_TEST(FissionSourceTest, sample_fission_sites)
{
    typedef typename TestFixture::Fission_Source Fission_Source;
    typedef typename TestFixture::Fission_Site Fission_Site;

    // Create the fission source.
    Fission_Source source( this->db, 
			   this->physics_tester->geometry(),
			   this->physics_tester->physics() );

    EXPECT_EQ( this->np, source.Np() );
    EXPECT_EQ( this->edges[0], source.lower_coords().x );
    EXPECT_EQ( this->edges[0], source.lower_coords().y );
    EXPECT_EQ( this->edges[0], source.lower_coords().z );
    EXPECT_EQ( this->edges[1], source.width().x );
    EXPECT_EQ( this->edges[1], source.width().y );
    EXPECT_EQ( this->edges[1], source.width().z );

    // Initialize with fission sites
    auto fission_sites = source.create_fission_site_container();
    Fission_Site site;
    site.m = 1;
    site.r.x = 50.0;
    site.r.y = 50.0;
    site.r.z = 50.0;
    for ( int n = 0; n < this->np; ++n )
    {
	fission_sites->push_back( site );
    }
    source.build_source( fission_sites );
    EXPECT_EQ( this->np, source.num_to_transport() );
    EXPECT_EQ( this->np * profugus::nodes(), source.total_num_to_transport() );
    EXPECT_EQ( this->np, source.num_left() );
    EXPECT_EQ( 0, source.num_run() );
    EXPECT_FALSE( source.empty() );
    EXPECT_TRUE( !source.is_initial_source() );

    // Set all the particles to dead.
    Teuchos::Array<cuda_profugus::events::Event> dead_event( 
	this->vector_size, cuda_profugus::events::DEAD );
    this->physics_tester->particle_tester().set_event( dead_event );

    // Get particles from the source.
    source.get_particles( this->physics_tester->particles() );
    EXPECT_EQ( source.num_to_transport() - this->vector_size, source.num_left() );
    EXPECT_EQ( this->vector_size, source.num_run() );
    
    // check particles
    Teuchos::Array<cuda_profugus::events::Event> events =
    	this->physics_tester->particle_tester().event();
    Teuchos::Array<int> matids = 
    	this->physics_tester->particle_tester().matid();
    Teuchos::Array<int> alive = 
    	this->physics_tester->particle_tester().alive();
    Teuchos::Array<double> wts = 
    	this->physics_tester->particle_tester().wt();
    for ( int i = 0; i < this->vector_size; ++i )
    {
    	EXPECT_EQ( cuda_profugus::events::TAKE_STEP, events[i] );
    	EXPECT_EQ( 1, matids[i] );
	EXPECT_TRUE( alive[i] );
    	EXPECT_EQ( 1.0/profugus::nodes(), wts[i] );	
    }
}

//---------------------------------------------------------------------------//
//                 end of tstFission_Source_cuda.cc
//---------------------------------------------------------------------------//
