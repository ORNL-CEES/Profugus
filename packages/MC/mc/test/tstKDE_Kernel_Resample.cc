//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   MC/mc/test/tstKDE_Kernel_Resample.cc
 * \author Gregory G. Davidson
 * \date   Thu Jan 22 13:50:16 2015
 * \brief  Tests the KDE Kernels
 * \note   Copyright (C) 2008 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "../KDE_Kernel_Resample.hh"

#include "SourceTestBase.hh"
#include "Utils/gtest/utils_gtest.hh"

//---------------------------------------------------------------------------//

class KernelTest : public SourceTestBase
{
    typedef SourceTestBase Base;

  protected:
    // >>> TYPEDEFS
    typedef profugus::KDE_Kernel_Resample KDE_Kernel_Resample;

  protected:
    virtual int get_seed() const
    {
        return 3421;
    }

    virtual void init_group_bounds()
    {
        Vec_Dbl n_bounds = {100.0, 0.001};

        b_group_bounds =
            std::make_shared<profugus::Group_Bounds>(n_bounds);
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

    // 2 material definitions/1 group
    /*
     - Mat 0 -> Moderator
     - Mat 1 -> Fuel
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
TEST_F(KernelTest, axial_kernel_bandwidth_calc)
{
    typedef Physics_t::Fission_Site       Fission_Site;
    typedef profugus::geometry::cell_type cell_type;

    // Define the bandwidth exponent
    double coeff = 1.06;
    double exponent = -0.7;

    // Create a Axial KDE kernel
    KDE_Kernel_Resample kernel(b_geometry, b_physics, coeff, exponent);

    // Create a bunch of fission sites
    std::vector<Fission_Site> fis_sites;
    fis_sites.push_back(Fission_Site{0, Space_Vector(0.45, 2.1,  12.3)});
    fis_sites.push_back(Fission_Site{0, Space_Vector(0.45, 2.1,  12.3)});
    fis_sites.push_back(Fission_Site{0, Space_Vector(0.45, 2.1,  12.3)});
    fis_sites.push_back(Fission_Site{0, Space_Vector(0.48, 1.9,  12.1)});
    fis_sites.push_back(Fission_Site{0, Space_Vector(0.63, 1.8,  11.1)});
    fis_sites.push_back(Fission_Site{0, Space_Vector(0.63, 1.8,  11.1)});
    fis_sites.push_back(Fission_Site{0, Space_Vector(1.89, 0.63,  8.4)});
    fis_sites.push_back(Fission_Site{0, Space_Vector(1.89, 0.63,  7.4)});
    fis_sites.push_back(Fission_Site{0, Space_Vector(1.89, 0.63,  9.4)});

    // Add to these fission sites if in parallel
    std::vector<Fission_Site> total_fis_sites;
    for (unsigned int n = 0; n < profugus::nodes(); ++n)
    {
        // Make copy
        total_fis_sites.insert(total_fis_sites.end(), fis_sites.begin(),
                               fis_sites.end());
    }

    // Calculate the sums and sum squares
    std::map<cell_type, double> sum;
    std::map<cell_type, double> sum_sq;
    std::map<cell_type, unsigned int> num_sites;
    for (const Fission_Site &fs : total_fis_sites)
    {
        cell_type cell = b_geometry->cell(fs.r);

        sum[cell]    += fs.r[def::Z];
        sum_sq[cell] += fs.r[def::Z]*fs.r[def::Z];
        num_sites[cell] += 1;
    }

    std::map<cell_type, double> ref_bandwidths;
    auto sum_iter    = sum.cbegin();
    auto sum_sq_iter = sum_sq.cbegin();
    for (auto ns_iter = num_sites.cbegin(); ns_iter != num_sites.cend();
         ++ns_iter, ++sum_iter, ++sum_sq_iter)
    {
        CHECK(sum_iter != sum.cend());
        CHECK(sum_sq_iter != sum_sq.cend());

        cell_type cell = sum_iter->first;
        double sum     = sum_iter->second;
        double sum_sq  = sum_sq_iter->second;
        unsigned int N = ns_iter->second;

        double variance = sum_sq/N - (sum/N)*(sum/N);
        ref_bandwidths[cell] = 1.06*std::sqrt(variance)*std::pow(N, exponent);
    }

    // Calculate the bandwidth
    kernel.calc_bandwidths(fis_sites);

    for (const auto &elem : ref_bandwidths)
    {
        cell_type cell       = std::get<0>(elem);
        double ref_bandwidth = std::get<1>(elem);

        EXPECT_SOFTEQ(ref_bandwidth, kernel.bandwidth(cell), 1.0e-6);
    }
}

//---------------------------------------------------------------------------//
TEST_F(KernelTest, test_bounds)
{
    // Create a KDE kernel
    KDE_Kernel_Resample kernel(b_geometry, b_physics);

    // Set the bandwidth in the pin (cell 1)
    double bandwidth = 2.5;
    kernel.set_bandwidth(1, bandwidth);

    // Create random number generator
    profugus::RNG rng = b_rcon->rng();

    // Do 1000 samples
    for (unsigned int i = 0; i < 1000; ++i)
    {
        // Create a fission site at the middle of the top-right pin
        Space_Vector site_pos(0.63, 1.89, 7.14);

        // Sample the position
        Space_Vector sampled_pos = kernel.sample_position(site_pos, rng);

        // Verify that the x,y coordinates are the same, and the z-position is
        // within the bandwidth
        EXPECT_SOFTEQ(site_pos[0], sampled_pos[0], 1.0e-12);
        EXPECT_SOFTEQ(site_pos[1], sampled_pos[1], 1.0e-12);
        EXPECT_LE(site_pos[2]-bandwidth, sampled_pos[2]);
        EXPECT_GE(site_pos[2]+bandwidth, sampled_pos[2]);
    }
}

//---------------------------------------------------------------------------//
TEST_F(KernelTest, test_in_pin)
{
    // Create a KDE kernel
    KDE_Kernel_Resample kernel(b_geometry, b_physics);

    // Set the bandwidth
    double bandwidth = 2.5;
    kernel.set_bandwidth(1, bandwidth);

    // Create random number generator
    profugus::RNG rng = b_rcon->rng();

    // Do 1000 samples
    for (unsigned int i = 0; i < 1000; ++i)
    {
        // Create a fission site within the top-left pin
        // Sample the position
        Space_Vector site_pos(0.0, 0.0, 0.0);
        do
        {
            site_pos[0] = rng.ran() * 0.54 + 0.63;
            site_pos[1] = rng.ran() * 0.54 + 1.89;
            site_pos[2] = rng.ran() * 14.28;
        } while (!is_in_pin(site_pos));

        // Sample the position
        Space_Vector sampled_pos = kernel.sample_position(site_pos, rng);

        // Verify that the sampled position is in the pin
        EXPECT_TRUE(is_in_pin(sampled_pos));
    }
}

//---------------------------------------------------------------------------//
#if 0
TEST_F(KernelTest, heuristic_test_resample)
{
    // Create a KDE kernel
    KDE_Kernel_Resample kernel(b_geometry, b_physics);

    // Set the bandwidth
    double bandwidth = 2.5;
    kernel.set_bandwidth(bandwidth);

    // Create random number generator
    RNG_Control::RNG rng = b_rcon->rng();

    // Store the mean and variance of the sampled position
    std::vector<double> test_mean = {0.0, 0.0, 0.0};
    std::vector<double> test_var  = {0.0, 0.0, 0.0};
    std::vector<double> ref_mean  = {0.0, 0.0, 0.0};
    std::vector<double> ref_var   = {0.0, 0.0, 0.0};

    // Do 1000 test samples
    for (unsigned int i = 0; i < 1000; ++i)
    {
        // Create a fission site within the top-right pin
        // Sample the position
        Space_Vector site_pos(0.0, 0.0, 0.0);
        do
        {
            site_pos[0] = rng.ran()*0.59 + 1.555;
            site_pos[1] = rng.ran()*0.59 + 1.555;
            site_pos[2] = rng.ran()*14.28;
        } while (!is_in_pin(site_pos));

        // Sample the position
        Space_Vector sampled_pos = kernel.sample_position(site_pos, rng);

        // Add the test statistic
        test_mean[0] += sampled_pos[0];
        test_mean[1] += sampled_pos[1];
        test_mean[2] += sampled_pos[2];
        test_var[0] += sampled_pos[0]*sampled_pos[0];
        test_var[1] += sampled_pos[1]*sampled_pos[1];
        test_var[2] += sampled_pos[2]*sampled_pos[2];

        // Do the heuristic sample (using the RESAMPLE method)
        ref_mean[0] += site_pos[0];
        ref_mean[1] += site_pos[1];
        double new_z = 0.0;
        do
        {
            double eps = 0.0;
            double ran1 = 2.0 * rng.ran() - 1.0;
            double ran2 = 2.0 * rng.ran() - 1.0;
            double ran3 = 2.0 * rng.ran() - 1.0;
            if (std::fabs(ran3) >= std::fabs(ran2) &&
                std::fabs(ran3) >= std::fabs(ran1))
            {
                eps = ran2;
            }
            else
            {
                eps = ran3;
            }
            new_z = site_pos[2] + bandwidth*eps/2.0;
        } while (new_z < 0.0 && new_z > 14.28);
        ref_mean[2] += new_z;
        ref_var[0] += site_pos[0] * site_pos[0];
        ref_var[1] += site_pos[1] * site_pos[1];
        ref_var[2] += new_z*new_z;
    }

    // Normalize statistics
    for (int i = 0; i < 3; ++i)
    {
        test_mean[i] /= 1000.0;
        test_var[i]  /= 999.0;
        ref_mean[i]  /= 1000.0;
        ref_var[i]   /= 999.0;
    }

    EXPECT_VEC_SOFTEQ(ref_mean, test_mean, 1.0e-2);
    EXPECT_VEC_SOFTEQ(ref_var, test_var, 3.0e-2);
}

//---------------------------------------------------------------------------//

TEST_F(KernelTest, heuristic_test_fissite)
{
    // Create a KDE kernel
    KDE_Kernel_Resample kernel(b_geometry, b_physics);

    // Set the bandwidth
    double bandwidth = 2.5;
    kernel.set_bandwidth(bandwidth);

    // Create random number generator
    RNG_Control::RNG rng = b_rcon->rng();

    // Store the mean and variance of the sampled position
    std::vector<double> test_mean = {0.0, 0.0, 0.0};
    std::vector<double> test_var  = {0.0, 0.0, 0.0};
    std::vector<double> ref_mean  = {0.0, 0.0, 0.0};
    std::vector<double> ref_var   = {0.0, 0.0, 0.0};

    // Do 1000 test samples
    for (unsigned int i = 0; i < 1000; ++i)
    {
        // Create a fission site within the top-right pin
        // Sample the position
        Space_Vector site_pos(0.0, 0.0, 0.0);
        do
        {
            site_pos[0] = rng.ran()*0.59 + 1.555;
            site_pos[1] = rng.ran()*0.59 + 1.555;
            site_pos[2] = rng.ran()*14.28;
        } while (!is_in_pin(site_pos));

        // Sample the position
        Space_Vector sampled_pos = kernel.sample_position(site_pos, rng);

        // Add the test statistic
        test_mean[0] += sampled_pos[0];
        test_mean[1] += sampled_pos[1];
        test_mean[2] += sampled_pos[2];
        test_var[0] += sampled_pos[0]*sampled_pos[0];
        test_var[1] += sampled_pos[1]*sampled_pos[1];
        test_var[2] += sampled_pos[2]*sampled_pos[2];

        // Do the heuristic sample (using the FISSION_SITE method)
        ref_mean[0] += site_pos[0];
        ref_mean[1] += site_pos[1];
        double eps = 0.0;
        double ran1 = 2.0 * rng.ran() - 1.0;
        double ran2 = 2.0 * rng.ran() - 1.0;
        double ran3 = 2.0 * rng.ran() - 1.0;
        if (std::fabs(ran3) >= std::fabs(ran2) &&
            std::fabs(ran3) >= std::fabs(ran1))
        {
            eps = ran2;
        }
        else
        {
            eps = ran3;
        }
        double new_z = site_pos[2] + bandwidth*eps/2.0;
        // If outside boundary, use original site
        if (new_z < 0.0 || new_z > 14.28)
            new_z = site_pos[2];
        ref_mean[2] += new_z;
        ref_var[0]  += site_pos[0] * site_pos[0];
        ref_var[1]  += site_pos[1] * site_pos[1];
        ref_var[2]  += new_z*new_z;
    }

    // Normalize statistics
    for (int i = 0; i < 3; ++i)
    {
        test_mean[i] /= 1000.0;
        test_var[i]  /= 999.0;
        ref_mean[i]  /= 1000.0;
        ref_var[i]   /= 999.0;
    }

    EXPECT_VEC_SOFTEQ(ref_mean, test_mean, 1.0e-2);
    EXPECT_VEC_SOFTEQ(ref_var, test_var, 3.0e-2);
}

//---------------------------------------------------------------------------//

TEST_F(KernelTest, heuristic_test_bandwidth_rej)
{
    // For this test, we're simply going to try sampling from near the
    // boundary
    // Create a KDE kernel
    KDE_Kernel_Resample kernel(b_geometry, b_physics);

    // Set the bandwidth
    double bandwidth = 2.5;
    kernel.set_bandwidth(bandwidth);

    // Create random number generator
    RNG_Control::RNG rng = b_rcon->rng();

    {
        // Try sampling from near the bottom boundary
        Space_Vector site_pos(0.0, 0.0, 1.0);
        Space_Vector new_pos = kernel.sample_position(site_pos, rng);
        EXPECT_SOFTEQ(0.0, new_pos[0], 1.0e-6);
        EXPECT_SOFTEQ(0.0, new_pos[1], 1.0e-6);
        EXPECT_SOFTEQ(1.0, new_pos[2], 1.0e-6);
    }

    {
        // Try sampling from near the top boundary
        Space_Vector site_pos(0.0, 0.0, 14.0);
        Space_Vector new_pos = kernel.sample_position(site_pos, rng);
        EXPECT_SOFTEQ(0.0,  new_pos[0], 1.0e-6);
        EXPECT_SOFTEQ(0.0,  new_pos[1], 1.0e-6);
        EXPECT_SOFTEQ(14.0, new_pos[2], 1.0e-6);
    }
}
#endif

//---------------------------------------------------------------------------//
//                        end of tstKDE_Kernel_Resample.cc
//---------------------------------------------------------------------------//
