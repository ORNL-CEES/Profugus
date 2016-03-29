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

    //------------------------------------------------------------------------//
    // Check that a given position is in one of our fuel pins
    //
    void check_position(const Fission_Site_Container &sites,
                        const def::Space_Vector       pos,
                        double                        bandwidth)
    {
        // Find the fission site that has the same x and y position as the
        // sampled position
        bool found = false;
        int index = 0;
        while (index < sites.size() && !found)
        {
            if (soft_equiv(sites[index].r[0], pos[0], 1.0e-6) &&
                soft_equiv(sites[index].r[1], pos[1], 1.0e-6))
            {
                found = true;
            }
            else
            {
                ++index;
            }
        }
        EXPECT_FALSE(!found);

        std::cout << "Index: " << index << std::endl;
        std::cout << "Site: (" << sites[index].r[0] << "," << sites[index].r[1]
                  << "," << sites[index].r[2] << ")" << std::endl;
        std::cout << "Position: (" << pos[0] << "," << pos[1] << "," << pos[2]
                  << ")" << std::endl;
        std::cout << "Bandwidth: " << bandwidth << std::endl;
        std::cout << "Z-diff: " << std::fabs(pos[2]-sites[index].r[2])
                  << std::endl;
        std::cout << std::endl;

        // Verify that z is within bandwidth
        EXPECT_TRUE(std::fabs(pos[2]-sites[index].r[2]) <= bandwidth);
    }
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEST_F(KernelTest, simple_test)
{
    b_db->set<int>("Np", 12);

    // Make the KDE database
    ParameterList_t kde_db("kde_db");
    kde_db.set<std::string>("kernel_type", std::string("fission_rejection"));

    // Set the KDE database
    b_db->set("kde_db", kde_db);

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
    site.r = Space_Vector(0.62, 1.90, 12.1);
    fsrc->push_back(site);
    site.m = 1;
    site.r = Space_Vector(0.64, 1.88, 11.1);
    fsrc->push_back(site);
    fsrc->push_back(site);
    EXPECT_EQ(6, fsrc->size());
    // Add 3 sites to the container for bottom-right pincell
    site.m = 1;
    site.r = Space_Vector(1.89, 0.63, 8.4);
    fsrc->push_back(site);
    site.r = Space_Vector(1.88, 0.64, 7.4);
    fsrc->push_back(site);
    site.r = Space_Vector(1.90, 0.62, 9.4);
    fsrc->push_back(site);

    // Make a copy of the fission source sites.  Reverse them because we draw
    // particles from the end
    SP_Fission_Sites sites = std::make_shared<Fission_Site_Container>(*fsrc);
    std::reverse(sites->begin(), sites->end());

    // Before setting the fission source should have 9 entries
    EXPECT_EQ(9, fsrc->size());
    source.build_source(fsrc);
    // After setting, it should be empty (swapped with internal container)
    EXPECT_EQ(0, fsrc->size());

    // Calculate the bandwidth for cell 3
    {
        unsigned int nodes = profugus::nodes();
        double sum    = 0.0;
        double sum_sq = 0.0;
        for (unsigned int node = 0; node < nodes; ++node)
        {
            sum    += 12.3*3.0 + 12.1 + 11.1*2.0;
            sum_sq += (12.3*12.3)*3.0 + 12.1*12.1 + (11.1*11.1)*2;
        }
        unsigned int N = 6*nodes;
        double variance = sum_sq/N - (sum/N)*(sum/N);
        double ref_bandwidth = 1.06 * std::sqrt(variance) *
                               std::pow(N,-0.20);
        EXPECT_SOFTEQ(ref_bandwidth, source.bandwidth(3), 1.0e-6);
    }
    // Calculate the bandwidth for cell 1
    {
        unsigned int nodes = profugus::nodes();
        double sum = 0.0;
        double sum_sq = 0.0;
        for (unsigned int node = 0; node < nodes; ++node)
        {
            sum    += 8.4 + 7.4 + 9.4;
            sum_sq += (8.4*8.4) + (7.4*7.4) + (9.4*9.4);
        }
        unsigned int N = 3*nodes;
        double variance = sum_sq/N - (sum/N)*(sum/N);
        double ref_bandwidth = 1.06 * std::sqrt(variance) * std::pow(N,-0.20);
        EXPECT_SOFTEQ(ref_bandwidth, source.bandwidth(1), 1.0e-6);
    }

    EXPECT_EQ(12, source.Np());
    EXPECT_EQ(9, source.num_to_transport()); // local
    EXPECT_EQ(9 * nodes, source.total_num_to_transport());

    int ctr       = 0;
    double tot_wt = 0.0;

    for (const Fission_Site &fs : *sites)
    {
        EXPECT_FALSE(source.empty());

        SP_Particle p = source.get_particle();
        ctr++;

        // Get the bandwidth for this site
        double bw = source.bandwidth(b_geometry->cell(fs.r));

        // Verify that the particle is in the right location (x,y same, z
        // within bandwidth)
        check_position(*sites, p->geo_state().d_r, bw);

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
// end of MC/mc/test/tstKDE_Fission_Source.cc
//---------------------------------------------------------------------------//
