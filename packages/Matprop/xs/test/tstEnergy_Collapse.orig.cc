//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   material/test/tstEnergy_Collapse.cc
 * \author Steven Hamilton
 * \date   Thu Mar 28 13:26:47 2013
 * \brief  Energy_Collapse unit test.
 * \note   Copyright (C) 2008 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>

#include "harness/DBC.hh"
#include "harness/Soft_Equivalence.hh"
#include "comm/global.hh"
#include "comm/Parallel_Unit_Test.hh"
#include "harness/Soft_Equivalence.hh"
#include "release/Release.hh"

#include "utils/SP.hh"
#include "../Cross_Sections.hh"
#include "../Mat_DB.hh"
#include "../Energy_Collapse.hh"

using denovo::SP;
using nemesis::Parallel_Unit_Test;
using nemesis::soft_equiv;


#define ITFAILS ut.failure(__LINE__);
#define UNIT_TEST(a) if (!(a)) ut.failure(__LINE__);

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void multi_material_test(Parallel_Unit_Test &ut)
{
    // make database
    denovo::SP<database::Std_DB> db;
    {
        db = new database::Std_DB("mat");
        UNIT_TEST(db);

        db->new_key("Pn_order", 0);
        db->new_key("num_groups", 8);
        db->new_key("downscatter", false);
    }

    // make a mat db (for 10 cells)
    SP<denovo::Mat_DB> mat_db( new denovo::Mat_DB(10) );

    // set number of materials and groups
    int num_mat = 2;
    int num_grps = 8;
    mat_db->set_num(num_mat,num_grps);

    // build the cross sections
    std::vector< SP<denovo::Cross_Sections> > xs0(num_grps);
    std::vector< SP<denovo::Cross_Sections> > xs1(num_grps);

    for( int g=0; g<num_grps; ++g )
    {
        xs0[g] = new denovo::Cross_Sections(db,g);
        xs1[g] = new denovo::Cross_Sections(db,g);
        xs0[g]->sigma() = static_cast<double>(10*g+1);
        xs1[g]->sigma() = static_cast<double>(20*g+1);
        for( int gp=0; gp<=g; ++gp )
        {
            xs0[g]->sigma_s(gp,0) = static_cast<double>(g) + gp*0.01;
            xs1[g]->sigma_s(gp,0) = static_cast<double>(g) + gp*0.001;
        }
    }

    // Now a few upscattering entries
    xs0[3]->sigma_s(4,0) = 3.04;
    xs0[5]->sigma_s(7,0) = 5.07;
    xs0[6]->sigma_s(7,0) = 6.07;

    xs1[4]->sigma_s(5,0) = 4.005;
    xs1[5]->sigma_s(6,0) = 5.006;
    xs1[5]->sigma_s(7,0) = 5.007;
    xs1[6]->sigma_s(7,0) = 6.007;

    // add matids to cells
    std::set<int> mat0;
    std::set<int> mat1;

    mat0.insert(0);
    mat0.insert(3);
    mat0.insert(4);
    mat0.insert(5);
    mat0.insert(7);
    mat0.insert(8);

    mat1.insert(1);
    mat1.insert(2);
    mat1.insert(6);
    mat1.insert(9);

    mat_db->assign(0, mat0);
    mat_db->assign(1, mat1);

    for( int g=0; g<num_grps; ++g )
    {
        xs0[g]->complete();
        mat_db->assign(xs0[g], 0);

        xs1[g]->complete();
        mat_db->assign(xs1[g], 1);
    }

    // Build coarse mat_db
    std::vector<int> steer(3);
    steer[0] = 3;
    steer[1] = 2;
    steer[2] = 3;

    std::vector<double> weights(8);
    weights[0] = 1.0;
    weights[1] = 2.0;
    weights[2] = 3.0;
    weights[3] = 2.0;
    weights[4] = 1.0;
    weights[5] = 1.0;
    weights[6] = 2.0;
    weights[7] = 1.0;
    SP<denovo::Mat_DB> coarse_mat =
        denovo::Energy_Collapse::collapse_all_mats(mat_db,steer,weights);

    UNIT_TEST(coarse_mat->verify());

    // Now test collapse
    UNIT_TEST(soft_equiv(coarse_mat->get(0,0).sigma(), 86.0/6.0 ));
    UNIT_TEST(soft_equiv(coarse_mat->get(0,1).sigma(),103.0/3.0 ));
    UNIT_TEST(soft_equiv(coarse_mat->get(0,2).sigma(),244.0/4.0 ));
    UNIT_TEST(soft_equiv(coarse_mat->get(0,0).sigma_s(0,0),15.1/6.0 ));
    UNIT_TEST(!coarse_mat->get(0,0).has_upscatter(1) );
    UNIT_TEST(!coarse_mat->get(0,0).has_upscatter(2) );
    UNIT_TEST(soft_equiv(coarse_mat->get(0,1).sigma_s(0,0),42.16/6.0 ));
    UNIT_TEST(soft_equiv(coarse_mat->get(0,1).sigma_s(1,0),21.2/3.0 ));
    UNIT_TEST(!coarse_mat->get(0,1).has_upscatter(2) );
    UNIT_TEST(soft_equiv(coarse_mat->get(0,2).sigma_s(0,0),108.24/6.0 ));
    UNIT_TEST(soft_equiv(coarse_mat->get(0,2).sigma_s(1,0),54.3/3.0 ));
    UNIT_TEST(soft_equiv(coarse_mat->get(0,2).sigma_s(2,0),62.6/4.0 ));

    UNIT_TEST(soft_equiv(coarse_mat->get(1,0).sigma(),166.0/6.0 ));
    UNIT_TEST(soft_equiv(coarse_mat->get(1,1).sigma(),203.0/3.0 ));
    UNIT_TEST(soft_equiv(coarse_mat->get(1,2).sigma(),484.0/4.0 ));
    UNIT_TEST(soft_equiv(coarse_mat->get(1,0).sigma_s(0,0),15.01/6.0 ));
    UNIT_TEST(!coarse_mat->get(1,0).has_upscatter(1) );
    UNIT_TEST(!coarse_mat->get(1,0).has_upscatter(2) );
    UNIT_TEST(soft_equiv(coarse_mat->get(1,1).sigma_s(0,0),42.016/6.0 ));
    UNIT_TEST(soft_equiv(coarse_mat->get(1,1).sigma_s(1,0),18.016/3.0 ));
    UNIT_TEST(coarse_mat->get(1,1).has_upscatter(2) );
    UNIT_TEST(soft_equiv(coarse_mat->get(1,1).sigma_s(2,0),4.005/4.0 ));
    UNIT_TEST(soft_equiv(coarse_mat->get(1,2).sigma_s(0,0),108.024/6.0 ));
    UNIT_TEST(soft_equiv(coarse_mat->get(1,2).sigma_s(1,0),54.03/3.0 ));
    UNIT_TEST(soft_equiv(coarse_mat->get(1,2).sigma_s(2,0),72.072/4.0 ));

    // Make sure matids got assigned correctly
    UNIT_TEST(coarse_mat->matid(0) == 0);
    UNIT_TEST(coarse_mat->matid(3) == 0);
    UNIT_TEST(coarse_mat->matid(4) == 0);
    UNIT_TEST(coarse_mat->matid(5) == 0);
    UNIT_TEST(coarse_mat->matid(7) == 0);
    UNIT_TEST(coarse_mat->matid(8) == 0);

    UNIT_TEST(coarse_mat->matid(1) == 1);
    UNIT_TEST(coarse_mat->matid(2) == 1);
    UNIT_TEST(coarse_mat->matid(6) == 1);
    UNIT_TEST(coarse_mat->matid(9) == 1);

    if (ut.numFails == 0)
        ut.passes("Energy_Collapse works correctly for multiple materials.");
}



//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    Parallel_Unit_Test ut(argc, argv, denovo::release);

    try
    {
        // >>> UNIT TESTS
        if (nemesis::nodes() == 1)
        {
            multi_material_test(ut);
        }
        else
        {
            ut.passes("Test only designed to work on 1 processor.");
        }

    }
    catch (std::exception &err)
    {
        std::cout << "ERROR: While testing tstEnergy_Collapse, "
                  << err.what()
                  << std::endl;
        ut.numFails++;
    }
    catch( ... )
    {
        std::cout << "ERROR: While testing tstEnergy_Collapse, "
                  << "An unknown exception was thrown."
                  << std::endl;
        ut.numFails++;
    }
    return ut.numFails;
}

//---------------------------------------------------------------------------//
//                        end of tstEnergy_Collapse.cc
//---------------------------------------------------------------------------//
