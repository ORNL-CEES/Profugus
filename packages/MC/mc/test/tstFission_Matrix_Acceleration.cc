//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/test/tstFission_Matrix_Acceleration.cc
 * \author Thomas M. Evans
 * \date   Wed Nov 12 14:59:19 2014
 * \brief  Test for Fission_Matrix_Acceleration
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <memory>
#include <vector>
#include <cmath>

#include "solvers/LinAlgTypedefs.hh"
#include "spn/Moment_Coefficients.hh"
#include "../Fission_Matrix_Acceleration.hh"

#include "gtest/utils_gtest.hh"

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//

template<class T>
class FM_AccelerationTest : public ::testing::Test
{
  protected:
    // >>> TYPEDEFS
    typedef profugus::Fission_Matrix_Acceleration         Acceleration;
    typedef std::shared_ptr<Acceleration>                 SP_Acceleration;
    typedef Acceleration::Problem_Builder_t               SPN_Builder;
    typedef profugus::Fission_Matrix_Acceleration_Impl<T> Implementation;
    typedef std::shared_ptr<Implementation>               SP_Implementation;
    typedef typename Implementation::Linear_System_t      Linear_System_t;
    typedef std::vector<double>                           Vec_Dbl;

  protected:
    void SetUp()
    {
        // make the acceleration
        implementation =
            std::make_shared<profugus::Fission_Matrix_Acceleration_Impl<T>>();
        acceleration = implementation;
    }

    void get_adjoint(Vec_Dbl &phi)
    {
        auto u = implementation->adjoint();

        const Linear_System_t &system = implementation->spn_system();

        N  = system.get_dims()->num_equations();
        Ng = acceleration->mat().xs().num_groups();
        Nc = acceleration->mesh().num_cells();

        // SPN moments
        double u_m[4] = {0.0, 0.0, 0.0, 0.0};

        // resize phi
        phi.resize(Nc * Ng);

        // loop over groups
        for (int g = 0; g < Ng; ++g)
        {
            // loop over cells on this domain
            for (int cell = 0; cell < Nc; ++cell)
            {
                // loop over moments
                for (int n = 0; n < N; ++n)
                {
                    EXPECT_GT(u.size(), system.index(g, n, cell));

                    // get the moment from the solution vector
                    u_m[n] = u[system.index(g, n, cell)];
                }

                // assign the scalar flux (SPN Phi_0 moment)
                phi[cell + Nc * g] =  profugus::Moment_Coefficients::u_to_phi(
                    u_m[0], u_m[1], u_m[2], u_m[3]);
            }
        }
    }

  protected:
    // >>> DATA

    SP_Implementation implementation;
    SP_Acceleration   acceleration;
    SPN_Builder       builder;

    int N, Ng, Nc;
};

//---------------------------------------------------------------------------//
// SETUP TEMPLATED TESTS
//---------------------------------------------------------------------------//
// Currently only Epetra is supported (due to adjoints)

typedef testing::Types<profugus::EpetraTypes> Test_Types;

TYPED_TEST_CASE(FM_AccelerationTest, Test_Types);

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TYPED_TEST(FM_AccelerationTest, initialization)
{
    // setup the spn problem
    this->builder.setup("inf_med_k.xml");

    // build the spn problem
    this->acceleration->build_problem(this->builder);

    // initialize the spn problem
    this->acceleration->initialize();

    // check the reference
    typename TestFixture::Vec_Dbl adjoint;
    this->get_adjoint(adjoint);

    typename TestFixture::Vec_Dbl inf_med_vec(this->Ng, 0.0);
    double norm = 0.0;
    for (int g = 0; g < this->Ng; ++g)
    {
        for (int cell = 0; cell < this->Nc; ++cell)
        {
            inf_med_vec[g] += adjoint[cell + this->Nc * g];
        }
        inf_med_vec[g] /= static_cast<double>(this->Nc);
        norm           += inf_med_vec[g] * inf_med_vec[g];
    }
    norm = 1.0 / std::sqrt(norm);

    // normalize the vector
    for (auto &x : inf_med_vec)
    {
        x *= norm;
    }

    double ref[] = {0.33611574,  0.3155125 ,  0.31716512,  0.40852777,  
                    0.41532128, 0.4177149 ,  0.41594701};

    EXPECT_VEC_SOFTEQ(ref, inf_med_vec, 1.0e-5);

    EXPECT_SOFTEQ(1.1889633277246718, this->implementation->keff(), 1.0e-6);
}

//---------------------------------------------------------------------------//
//                 end of tstFission_Matrix_Acceleration.cc
//---------------------------------------------------------------------------//
