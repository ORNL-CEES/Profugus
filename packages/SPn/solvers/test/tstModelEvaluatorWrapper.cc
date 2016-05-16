//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   SPn/solvers/test/tstModelEvaluatorWrapper.cc
 * \author Steven Hamilton
 * \date   Thu Apr 02 09:51:48 2015
 * \brief  Test of ModelEvaluatorWrapper interface.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "gtest/utils_gtest.hh"
#include <vector>

#include "../ModelEvaluatorWrapper.hh"
#include "../LinAlgTypedefs.hh"
#include "../ThyraTraits.hh"
#include "LinAlgTraits.hh"

using profugus::ModelEvaluatorWrapper;

//---------------------------------------------------------------------------//
// Test fixture
//---------------------------------------------------------------------------//

// NOTE: the test class name must not contain underscores.
template <class T>
class ModelEvaluatorWrapperTest : public ::testing::Test
{
  protected:
    // Typedefs usable inside the test fixture

  protected:
    // Initialization that are performed for each test
    void SetUp()
    {
        d_N = 50;
        d_A = linalg_traits::build_matrix<T>("laplacian",d_N);
        d_x = linalg_traits::build_vector<T>(d_N);
        d_f = linalg_traits::build_vector<T>(d_N);
    }

  protected:
    // >>> Data that get re-initialized between tests

    int d_N;
    Teuchos::RCP<typename T::OP> d_A;
    Teuchos::RCP<typename T::MV> d_x;
    Teuchos::RCP<typename T::MV> d_f;
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

using profugus::EpetraTypes;
using profugus::TpetraTypes;
typedef ::testing::Types<EpetraTypes,TpetraTypes> MyTypes;
TYPED_TEST_CASE(ModelEvaluatorWrapperTest, MyTypes);

TYPED_TEST(ModelEvaluatorWrapperTest, Constant)
{
    typedef typename TypeParam::ST ST;
    // Initialize vectors
    std::vector<double> x_data(this->d_N,1.0);
    std::vector<double> f_data(this->d_N,0.0);
    linalg_traits::fill_vector<TypeParam>(this->d_x,x_data);
    linalg_traits::fill_vector<TypeParam>(this->d_f,f_data);

    // Create wrapper
    ModelEvaluatorWrapper<TypeParam> wrapper(this->d_A);

    // Wrap input and output vectors appropriately
    auto x_thyra = profugus::ThyraTraits<TypeParam>::buildThyraMV(
        this->d_x,wrapper.get_x_space());
    auto inputs = wrapper.createInArgs();
    inputs.set_x(x_thyra->col(0));

    auto f_thyra = profugus::ThyraTraits<TypeParam>::buildThyraMV(
        this->d_f,wrapper.get_f_space());
    Thyra::ModelEvaluatorBase::Evaluation<Thyra::VectorBase<ST> > f_eval(
        f_thyra->col(0));
    auto outputs = wrapper.createOutArgs();
    outputs.set_f(f_eval);

    // Apply model evaluator
    wrapper.evalModel(inputs,outputs);

    // Check output: x_i=0 except at boundaries
    f_data[0] = 1.0;
    f_data[this->d_N-1] = 1.0;
    for( int i=1; i<this->d_N-1; ++i )
        f_data[i] = 0.0;

    linalg_traits::test_vector<TypeParam>(this->d_f,f_data);
}

TYPED_TEST(ModelEvaluatorWrapperTest, Quadratic)
{
    typedef typename TypeParam::ST ST;

    // Initialize vectors: x_i = i^2
    std::vector<double> x_data(this->d_N,0.0);
    for( int i=0; i<this->d_N; ++i )
        x_data[i] = static_cast<double>((i+1)*(i+1));
    std::vector<double> f_data(this->d_N,0.0);
    linalg_traits::fill_vector<TypeParam>(this->d_x,x_data);
    linalg_traits::fill_vector<TypeParam>(this->d_f,f_data);

    // Create wrapper
    ModelEvaluatorWrapper<TypeParam> wrapper(this->d_A);

    // Wrap input and output vectors appropriately
    auto x_thyra = profugus::ThyraTraits<TypeParam>::buildThyraMV(
        this->d_x,wrapper.get_x_space());
    auto inputs = wrapper.createInArgs();
    inputs.set_x(x_thyra->col(0));

    auto f_thyra = profugus::ThyraTraits<TypeParam>::buildThyraMV(
        this->d_f,wrapper.get_f_space());
    Thyra::ModelEvaluatorBase::Evaluation<Thyra::VectorBase<ST> > f_eval(
        f_thyra->col(0));
    auto outputs = wrapper.createOutArgs();
    outputs.set_f(f_eval);

    // Apply model evaluator
    wrapper.evalModel(inputs,outputs);

    // Check output: f_i = -(i-1)^2 + 2*(i)^2 - (i+1)^2
    f_data[0] = -2.0;
    f_data[this->d_N-1] = static_cast<double>(this->d_N * this->d_N
            + 2 * this->d_N - 1 );
    for( int i=1; i<this->d_N-1; ++i )
        f_data[i] = -2.0;

    linalg_traits::test_vector<TypeParam>(this->d_f,f_data);
}

//---------------------------------------------------------------------------//
//                 end of tstModelEvaluatorWrapper.cc
//---------------------------------------------------------------------------//
