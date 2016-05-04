//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   MC/cuda_mc/test/Test_XS.cc
 * \author Steven Hamilton
 * \date   Wed May 04 16:31:30 2016
 * \brief  Test_XS class definitions.
 * \note   Copyright (c) 2016 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "harness/DBC.hh"
#include "Test_XS.hh"

using profugus::XS;

Teuchos::RCP<profugus::XS> Test_XS::build_xs(int num_groups)
{
    VALIDATE(num_groups==1 || num_groups==2 || num_groups==3 || num_groups==5,
             "Invalid group request.");
    if (num_groups==1)
        return build_1grp_xs();
    else if (num_groups==2)
        return build_2grp_xs();
    else if (num_groups==3)
        return build_3grp_xs();
    else if (num_groups==5)
        return build_5grp_xs();
    return Teuchos::RCP<profugus::XS>();
}

Teuchos::RCP<XS> Test_XS::build_5grp_xs()
{
    Teuchos::RCP<XS> xs = Teuchos::rcp(new XS());
    xs->set(0, 5);

    std::vector<double> bnd(6, 0.0);
    bnd[0] = 100.0;
    bnd[1] = 10.0;
    bnd[2] = 1.0;
    bnd[3] = 0.1;
    bnd[4] = 0.01;
    bnd[5] = 0.001;
    xs->set_bounds(bnd);

    typename XS::OneDArray total(5);
    typename XS::TwoDArray scat(5, 5);

    double f[5] = {0.1, 0.4, 1.8, 5.7, 9.8};
    double c[5] = {0.3770, 0.4421, 0.1809, 0.0, 0.0};
    double n[5] = {2.4*f[0], 2.4*f[1], 2.4*f[2], 2.4*f[3], 2.4*f[4]};
    typename XS::OneDArray fission(std::begin(f), std::end(f));
    typename XS::OneDArray chi(std::begin(c), std::end(c));
    typename XS::OneDArray nus(std::begin(n), std::end(n));
    xs->add(1, XS::SIG_F, fission);
    xs->add(1, XS::NU_SIG_F, nus);
    xs->add(1, XS::CHI, chi);

    // mat 0
    total[0] = 5.2 ;
    total[1] = 11.4;
    total[2] = 18.2;
    total[3] = 29.9;
    total[4] = 27.3;
    xs->add(0, XS::TOTAL, total);

    // mat 1
    total[0] = 5.2  + f[0];
    total[1] = 11.4 + f[1];
    total[2] = 18.2 + f[2];
    total[3] = 29.9 + f[3];
    total[4] = 27.3 + f[4];
    xs->add(1, XS::TOTAL, total);

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
    return xs;
}

Teuchos::RCP<XS> Test_XS::build_3grp_xs()
{
    const int ng = 3;
    Teuchos::RCP<XS> xs = Teuchos::rcp(new XS());
    xs->set(0, ng);

    // make group boundaries
    XS::OneDArray nbnd(ng+1, 0.0);
    nbnd[0] = 100.0;
    nbnd[1] = 1.0;
    nbnd[2] = 0.01;
    nbnd[3] = 0.0001;
    xs->set_bounds(nbnd);

    double t1[ng] = {1.1, 1.6, 2.9};
    double t2[ng] = {10.0, 11.3, 16.2};

    XS::OneDArray tot1(std::begin(t1), std::end(t1));
    XS::OneDArray tot2(std::begin(t2), std::end(t2));

    xs->add(0, XS::TOTAL, tot1);
    xs->add(1, XS::TOTAL, tot2);

    double s1[][3] = {{0.7, 0.0, 0.0},
                      {0.2, 0.3, 0.0},
                      {0.1, 0.7, 1.9}};

    double s2[][3] = {{2.7, 0.0, 0.0},
                      {2.2, 2.3, 0.0},
                      {2.1, 2.7, 3.9}};

    XS::TwoDArray sct1(ng, ng, 0.0);
    XS::TwoDArray sct2(ng, ng, 0.0);

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

    XS::OneDArray chi2(std::begin(c2), std::end(c2));
    XS::OneDArray fis2(std::begin(f2), std::end(f2));
    XS::OneDArray nuf2(std::begin(n2), std::end(n2));

    xs->add(1, XS::CHI, chi2);
    xs->add(1, XS::SIG_F, fis2);
    xs->add(1, XS::NU_SIG_F, nuf2);

    xs->complete();
    return xs;
}

Teuchos::RCP<XS> Test_XS::build_2grp_xs()
{
    constexpr int ng = 2;
    Teuchos::RCP<XS> xs = Teuchos::rcp(new XS());
    xs->set(0, ng);

    // make group boundaries
    XS::OneDArray nbnd(ng+1, 0.0);
    nbnd[0] = 100.0;
    nbnd[1] = 1.0;
    nbnd[2] = 0.00001;
    xs->set_bounds(nbnd);

    double t1[ng] = {1.0, 0.1};

    XS::OneDArray tot1(std::begin(t1), std::end(t1));

    xs->add(0, XS::TOTAL, tot1);

    double s1[][ng] = {{0.5, 0.0},{0.1, 0.05}};

    XS::TwoDArray sct1(ng, ng, 0.0);

    for (int g = 0; g < ng; ++g)
    {
        for (int gp = 0; gp < ng; ++gp)
        {
            sct1(g, gp) = s1[g][gp];
        }
    }

    xs->add(0, 0, sct1);

    double c[] = {0.9, 0.10};
    double f[] = {0.2, 0.02};
    double n[] = {0.4, 0.05};

    XS::OneDArray chi(std::begin(c),std::end(c));
    XS::OneDArray fis(std::begin(f),std::end(f));
    XS::OneDArray nuf(std::begin(n),std::end(n));

    xs->add(0, XS::CHI,       chi);
    xs->add(0, XS::SIG_F,     fis);
    xs->add(0, XS::NU_SIG_F,  nuf);

    xs->complete();

    return xs;
}


Teuchos::RCP<XS> Test_XS::build_1grp_xs()
{
    const int ng = 1;
    Teuchos::RCP<XS> xs = Teuchos::rcp(new XS());
    xs->set(0, ng);

    // make group boundaries
    XS::OneDArray nbnd(ng+1, 0.0);
    nbnd[0] = 100.0;
    nbnd[1] = 0.00001;
    xs->set_bounds(nbnd);

    double t1[ng] = {1.0};

    XS::OneDArray tot1(std::begin(t1), std::end(t1));

    xs->add(0, XS::TOTAL, tot1);

    double s1[][3] = {{0.5}};

    XS::TwoDArray sct1(ng, ng, 0.0);

    for (int g = 0; g < ng; ++g)
    {
        for (int gp = 0; gp < ng; ++gp)
        {
            sct1(g, gp) = s1[g][gp];
        }
    }

    xs->add(0, 0, sct1);

    double c[] = {1.0};
    double f[] = {0.2};
    double n[] = {0.4};

    XS::OneDArray chi(std::begin(c),std::end(c));
    XS::OneDArray fis(std::begin(f),std::end(f));
    XS::OneDArray nuf(std::begin(n),std::end(n));

    xs->add(0, XS::CHI,       chi);
    xs->add(0, XS::SIG_F,     fis);
    xs->add(0, XS::NU_SIG_F,  nuf);



    xs->complete();
    return xs;
}

//---------------------------------------------------------------------------//
// end of MC/cuda_mc/test/Test_XS.cc
//---------------------------------------------------------------------------//
