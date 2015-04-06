//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   solvers/test/tstAndersonSolver.cc
 * \author Steven Hamilton
 * \date   Mon Apr 06 08:50:55 2015
 * \brief  Test interface to NOX Anderson acceleration solver
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "gtest/utils_gtest.hh"

#include "../AndersonSolver.hh"
#include "spn/OperatorAdapter.hh"
#include "spn/MatrixTraits.hh"
#include "spn/VectorTraits.hh"
#include "LinAlgTraits.hh"

using profugus::AndersonSolver;

//---------------------------------------------------------------------------//
// Discretization of Chandrasekhar H-equation, cf. Kelley Eq. 5.22
//---------------------------------------------------------------------------//
template <class T>
class Chandrasekhar_H : public profugus::OperatorAdapter<T>
{
  public:

    typedef typename T::MV  MV;
    typedef typename T::MAP MAP;

    Chandrasekhar_H(double c, int N, Teuchos::RCP<const MAP> map)
        : profugus::OperatorAdapter<T>(map)
        , d_c(c)
        , d_N(N)
    {
        REQUIRE( d_c > 0.0 && d_c < 1.0 );
        REQUIRE( d_N > 0 );
    }

  private:

    void ApplyImpl(const MV &x, MV &y) const override
    {
        // Get all global values of x (all_gather)
        auto x_global = linalg_traits::get_global_copy<T>(
            Teuchos::rcpFromRef(x));

        // mu is discretized evaluation points
        std::vector<double> mu(d_N,0.0);
        for( int i=0; i<d_N; ++i )
        {
            mu[i] = (static_cast<double>(i+1)-0.5)/static_cast<double>(d_N);
        }

        // Evaluate H-equation
        std::vector<double> y_data(d_N,0.0);
        for( int i=0; i<d_N; ++i )
        {
            // Integral over mu
            double rsum=0.0;
            for( int j=0; j<d_N; ++j )
            {
                rsum += mu[i] * x_global[j] / (mu[i] + mu[j]);
            }
            y_data[i] = x_global[i] - 1.0 /
                (1.0 - (d_c/static_cast<double>(2*d_N)) * rsum);
        }

        linalg_traits::fill_vector<T>(Teuchos::rcpFromRef(y),y_data);
    }

    double d_c;
    int d_N;
};

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//

template <class T>
class AndersonSolverTest : public ::testing::Test
{
  protected:

    typedef typename T::OP  OP;
    typedef typename T::MV  MV;
    typedef typename T::MAP MAP;

    // Initialization that are performed for each test
    void SetUp()
    {
        // Build operator
        d_c = 0.5;
        d_N = 50;
        d_num_local = d_N / profugus::nodes();
        if( profugus::node() < (d_N % profugus::nodes()) )
            d_num_local++;
        d_map = profugus::MatrixTraits<T>::build_map(d_num_local,d_N);
        d_op = Teuchos::rcp( new Chandrasekhar_H<T>(d_c,d_N,d_map) );
    }

    double d_c;
    int d_N, d_num_local;
    Teuchos::RCP<const MAP> d_map;
    Teuchos::RCP<OP> d_op;

};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

using profugus::EpetraTypes;
using profugus::TpetraTypes;
typedef ::testing::Types<EpetraTypes,TpetraTypes> MyTypes;
TYPED_TEST_CASE(AndersonSolverTest, MyTypes);

TYPED_TEST(AndersonSolverTest, Heqn)
{
    typedef typename TypeParam::MV MV;

    Teuchos::RCP<Teuchos::ParameterList> pl( new Teuchos::ParameterList() );
    pl->set("tolerance",1e-6);
    pl->set("max_itr",50);
    Teuchos::RCP<Teuchos::ParameterList> anderson_pl =
        Teuchos::sublist(pl,"Anderson Parameters");
    anderson_pl->set("Storage Depth",5);

    // Create and initialize vectors
    Teuchos::RCP<MV> x = linalg_traits::build_vector<TypeParam>(this->d_N);
    std::vector<double> x_data(this->d_N,0.0);
    linalg_traits::fill_vector<TypeParam>(x,x_data);

    // Build solver
    profugus::AndersonSolver<TypeParam> solver(this->d_op,pl);

    // Solve
    solver.solve(x);

    // Check solution
    std::vector<double> expected =
       {1.012429,   1.030067,   1.044180,   1.056382,   1.067269,   1.077160,
        1.086252,   1.094680,   1.102543,   1.109917,   1.116860,   1.123420,
        1.129635,   1.135540,   1.141162,   1.146525,   1.151651,   1.156556,
        1.161258,   1.165771,   1.170108,   1.174280,   1.178298,   1.182171,
        1.185907,   1.189515,   1.193002,   1.196374,   1.199638,   1.202798,
        1.205861,   1.208831,   1.211712,   1.214509,   1.217225,   1.219865,
        1.222432,   1.224928,   1.227357,   1.229722,   1.232025,   1.234269,
        1.236456,   1.238589,   1.240669,   1.242699,   1.244680,   1.246614,
        1.248504,   1.250349};

    linalg_traits::test_vector<TypeParam>(x,expected);
}

//---------------------------------------------------------------------------//
//                 end of tstAndersonSolver.cc
//---------------------------------------------------------------------------//
