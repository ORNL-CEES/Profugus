//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   MC/mc/test/tstKDE_Fission_Source.cc
 * \author Gregory G. Davidson
 * \date   Wed Dec 02 09:37:32 2015
 * \brief  KDE_Fission_Source class definitions.
 * \note   Copyright (c) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../KDE_Fission_Source.hh"

#include <Utils/config.h>
#include "Utils/gtest/utils_gtest.hh"
#include "SourceTestBase.hh"

using profugus::KDE_Fission_Source;

using def::X; using def::Y; using def::Z;
using def::Vec_Int;

inline void set_pos(double pos[3], double x, double y, double z)
{
    pos[X] = x;
    pos[Y] = y;
    pos[Z] = z;
}

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//
class KernelTest : public SourceTestBase
{
    typedef SourceTestBase Base;

  protected:
    typedef profugus::KDE_Fission_Source  KDE_Source;

    typedef KDE_Source::Fission_Site_Container Fission_Site_Container;
    typedef KDE_Source::Fission_Site           Fission_Site;
    typedef KDE_Source::SP_Physics             SP_Physics;
    typedef KDE_Source::SP_Fission_Sites       SP_Fission_Sites;
    typedef KDE_Source::SP_Particle            SP_Particle;

    virtual int get_seed() const
    {
        return 3421;
    }

    virtual void init_group_bounds()
    {
        Vec_Dbl n_bounds = {100.0, 0.001};

        b_group_bounds = std::make_shared<profugus::Group_Bounds>(n_bounds);
        ENSURE(b_group_bounds->num_groups() == 1);
    }

    /*
      Test Lattice:

      |-------|-------|
      |       |       |
      |  UO2  |  H2O  |
      |       |       |
      |-------|-------|
      |       |       |
      |  H2O  |  UO2  |
      |       |       |
      |-------|-------|

      PIN: FUEL = 1
           MOD  = 0

      LATTICE:
           UO2 - 1
           H2O - 2
     */
    virtual void init_geometry()
    {
        typedef Geometry_t::SP_Array SP_Core;
        typedef Geometry_t::Array_t  Core_t;
        typedef Core_t::SP_Object    SP_Lattice;
        typedef Core_t::Object_t     Lattice_t;
        typedef Lattice_t::SP_Object SP_Pin_Cell;
        typedef Lattice_t::Object_t  Pin_Cell_t;

        // make pin cells
        SP_Pin_Cell uo2(std::make_shared<Pin_Cell_t>(1, 0.54, 0, 1.26, 14.28));
        SP_Pin_Cell h2o(std::make_shared<Pin_Cell_t>(0, 1.26, 14.28));

        // make lattice
        SP_Lattice lat(std::make_shared<Lattice_t>(2, 2, 1, 3));
        lat->assign_object(uo2, 1);
        lat->assign_object(h2o, 2);
        lat->id(0, 0, 0) = 2; // H20
        lat->id(1, 0, 0) = 1; // UO2
        lat->id(0, 1, 0) = 1; // UO2
        lat->id(1, 1, 0) = 2; // H2O

        // complete lattice
        lat->complete(0.0, 0.0, 0.0);

        // make the core
        SP_Core core(std::make_shared<Core_t>(1, 1, 1, 1));
        core->assign_object(lat, 0);
        core->complete(0.0, 0.0, 0.0);

        // make the b_geometry
        b_geometry = std::make_shared<Geometry_t>(core);
    }

    /*
     - Mat 0 -> H20
     - Mat 1 -> FUEL
     */
    virtual void init_physics()
    {
        RCP_XS xs(Teuchos::rcp(new XS_t()));
        xs->set(0, 1);

        XS_t::OneDArray total0(1, 1.1),   total1(1, 10.0);
        XS_t::TwoDArray scat0(1, 1, 0.9), scat1(1, 1, 2.1);

        XS_t::OneDArray chi1(1, 1.0);
        XS_t::OneDArray sigf1(1, 4.2);
        XS_t::OneDArray nusigf1(1, 2.4*4.2);

        xs->add(0, XS_t::TOTAL, total0);
        xs->add(0, 0, scat0);

        xs->add(1, XS_t::TOTAL, total1);
        xs->add(1, XS_t::CHI, chi1);
        xs->add(1, XS_t::NU_SIG_F, nusigf1);
        xs->add(1, XS_t::SIG_F, sigf1);

        XS_t::OneDArray bounds(b_group_bounds->group_bounds());

        xs->set_bounds(bounds);

        xs->complete();

        b_physics = std::make_shared<Physics_t>(b_db, xs);

        EXPECT_FALSE(b_physics->is_fissionable(0));
        EXPECT_TRUE(b_physics->is_fissionable(1));
    }

    //------------------------------------------------------------------------//
    // Check that a given position is in one of our fuel pins
    //
    bool is_in_pin(const Space_Vector &pos)
    {
        // Pin origins are (1.89, 0.63) and (0.63, 1.89) with radius 0.54
        double dist_1 = sqrt(  (pos[0] - 1.89) * (pos[0] - 1.89)
                             + (pos[1] - 0.63) * (pos[1] - 0.63));
        double dist_2 = sqrt(  (pos[0] - 0.63) * (pos[0] - 0.63)
                             + (pos[1] - 1.89) * (pos[1] - 1.89));

        if (dist_1 > 0.54 && dist_2 > 0.54)
        {
            return false;
        }
        else
        {
            return true;
        }
    }
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST_F(KernelTest, simple_test)
{
    b_db->set<int>("Np", 12);
    b_db->set<std::string>("kernel_type", std::string("resample"));

    // make fission source
    KDE_Source source(b_db, b_geometry, b_physics, b_rcon);

    // get a fission source container
    SP_Fission_Sites fsrc = std::make_shared<Fission_Site_Container>();
    EXPECT_TRUE(fsrc->empty());

    // add 3 sites to the container (6 total fissions)
    Fission_Site site;
    site.m = 1;
    site.r = Space_Vector(0.63, 1.89, 12.3);
    fsrc->push_back(site);
    fsrc->push_back(site);
    fsrc->push_back(site);
    site.m = 1;
    site.r = Space_Vector(0.63, 1.89, 12.1);
    fsrc->push_back(site);
    site.m = 1;
    site.r = Space_Vector(0.63, 1.89, 11.1);
    fsrc->push_back(site);
    fsrc->push_back(site);
    EXPECT_EQ(6, fsrc->size());
    // Add 3 sites to the container for bottom-right pincell
    site.m = 1;
    site.r = Space_Vector(1.89, 0.63, 8.4);
    fsrc->push_back(site);
    site.r = Space_Vector(1.89, 0.63, 7.4);
    fsrc->push_back(site);
    site.r = Space_Vector(1.89, 0.63, 9.4);
    fsrc->push_back(site);

    // Before setting the fission source should have 9 entries
    EXPECT_EQ(9, fsrc->size());
    source.build_source(fsrc);
    // After setting, it should be empty (swapped with internal container)
    EXPECT_EQ(0, fsrc->size());

    // Calculate the bandwidth for cell 3
    {
        double sum   = 12.3*3.0 + 12.1 + 11.1*2.0;
        double sum_sq = (12.3*12.3)*3.0 + 12.1*12.1 + (11.1*11.1)*2;
        double variance = sum_sq/6.0 - (sum/6.0)*(sum/6.0);
        double ref_bandwidth = 1.06 * std::sqrt(variance) * std::pow(6,-0.20);
        EXPECT_SOFTEQ(ref_bandwidth, source.bandwidth(3), 1.0e-6);
    }
    // Calculate the bandwidth for cell 1
    {
        double sum   = 8.4 + 7.4 + 9.4;
        double sum_sq = (8.4*8.4) + (7.4*7.4) + (9.4*9.4);
        double variance = sum_sq/3.0 - (sum/3.0)*(sum/3.0);
        double ref_bandwidth = 1.06 * std::sqrt(variance) * std::pow(3,-0.20);
        EXPECT_SOFTEQ(ref_bandwidth, source.bandwidth(1), 1.0e-6);
    }

    EXPECT_EQ(12, source.Np());
    EXPECT_EQ(9, source.num_to_transport()); // local
    EXPECT_EQ(9 * nodes, source.total_num_to_transport());

    int ctr       = 0;
    double tot_wt = 0.0;
    while (!source.empty())
    {
        SP_Particle p = source.get_particle();
        ctr++;

        EXPECT_EQ(1, p->matid());
        EXPECT_SOFT_EQ(12.0 / (9 * nodes), p->wt());
        tot_wt += p->wt();
    }

    profugus::global_sum(tot_wt);

    EXPECT_SOFTEQ(12.0, tot_wt, 1.e-12);
    EXPECT_EQ(9, ctr);
    EXPECT_EQ(9, source.num_run());
    EXPECT_TRUE(source.empty());
}

//---------------------------------------------------------------------------//
#if 0
TEST_F(KernelTest, axial_kernel_test)
{
    b_db->set<double>("Np", 12);

    EXPECT_EQ(1, b_bmesh->num_blocks());
    EXPECT_EQ(nodes, b_bmesh->num_sets());

    // make fission source
    Initial_Source init_source(b_db, b_geometry, b_physics, b_rcon, b_bmesh);

    auto attr = std::make_shared<mc::Separable_Source_Attr>(
            std::make_shared<mc::Box_Shape>(0.0, 2.52, 0.0, 2.52, 0.0, 14.28),
            std::make_shared<mc::Line_Energy_Distribution>(100.0),
            std::make_shared<mc::Isotropic_Angular_Distribution>());

    init_source.add_source(std::make_shared<shift::Separable_Source>(
            attr, b_geometry, b_physics));
    init_source.complete();

    // Make the KDE Kernel
    auto kernel =
        std::make_shared<shift::Axial_Kernel>(b_geometry, b_physics);

    // Define the filename
    std::string filename = "kde_fiss_src-np" + std::to_string(nemesis::nodes())
                           + ".h5";
    std::vector<unsigned int> write_cycles = {0, 5, 10};

    // get a fission source container
    KDE_Source source(init_source, kernel, filename, write_cycles);

    Fission_Site_Container fsrc;
    EXPECT_TRUE(fsrc.empty());

    // add 3 sites to the container (6 total fissions) for upper-left
    // pincell
    Fission_Site site;
    site.physics_state.matid = 1;
    set_pos(site.geo_state.pos, 0.45, 2.1, 12.3);
    fsrc.push_back(site);
    fsrc.push_back(site);
    fsrc.push_back(site);
    site.physics_state.matid = 1;
    set_pos(site.geo_state.pos, 0.48, 1.9, 12.1);
    fsrc.push_back(site);
    site.physics_state.matid = 1;
    set_pos(site.geo_state.pos, 0.63, 1.8, 11.1);
    fsrc.push_back(site);
    fsrc.push_back(site);
    // Add 3 sites to the container for bottom-right pincell
    site.physics_state.matid = 1;
    set_pos(site.geo_state.pos, 1.89, 0.63, 8.4);
    fsrc.push_back(site);
    set_pos(site.geo_state.pos, 1.89, 0.63, 7.4);
    fsrc.push_back(site);
    set_pos(site.geo_state.pos, 1.89, 0.63, 9.4);
    fsrc.push_back(site);
    source.set_fission_sites(std::move(fsrc));
    EXPECT_EQ(9, fsrc.size());
    fsrc.clear();

    // Calculate the bandwidth for the upper-left pincell
    double sum_ul           = 12.3*3.0 + 12.1 + 11.1*2.0;
    double sum_sq_ul        = (12.3*12.3)*3.0 + 12.1*12.1 + (11.1*11.1)*2;
    double variance_ul      = sum_sq_ul/6.0 - (sum_ul/6.0)*(sum_ul/6.0);
    double ref_bandwidth_ul = 1.06*std::sqrt(variance_ul)*std::pow(6,-0.20);
    EXPECT_SOFTEQ(ref_bandwidth_ul, source.bandwidth(3), 1.0e-6);

    // Calculate the bandwidth for the bottom-right pincell
    double sum_br           = 8.4 + 7.4 + 9.4;
    double sum_sq_br        = (8.4*8.4) + (7.4*7.4) + (9.4*9.4);
    double variance_br      = sum_sq_br/3.0 - (sum_br/3.0)*(sum_br/3.0);
    double ref_bandwidth_br = 1.06*std::sqrt(variance_br)*std::pow(3,-0.20);
    EXPECT_SOFTEQ(ref_bandwidth_br, source.bandwidth(1), 1.0e-6);

    source.complete();

    EXPECT_EQ(12, source.num_requested());
    EXPECT_EQ(9, source.num_to_transport()); // local
    EXPECT_EQ(9, source.num_to_transport_in_set());
    EXPECT_EQ(9 * nodes, source.total_num_to_transport());

    int ctr       = 0;
    double tot_wt = 0.0;
    while (!source.empty())
    {
        Particle_t p;

        source.init_particle(p);
        ctr++;

        EXPECT_EQ(1, p.matid());
        EXPECT_SOFT_EQ(12.0 / (9 * nodes), p.wt());
        tot_wt += p.wt();
    }

    nemesis::global_sum(tot_wt);

    EXPECT_SOFTEQ(12.0, tot_wt, 1.e-12);
    EXPECT_EQ(9, ctr);
    EXPECT_EQ(9, source.num_run());
    EXPECT_TRUE(source.empty());
}
#endif

//---------------------------------------------------------------------------//
// end of MC/mc/test/tstKDE_Fission_Source.cc
//---------------------------------------------------------------------------//
