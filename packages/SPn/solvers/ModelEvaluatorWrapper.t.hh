//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   SPn/solvers/ModelEvaluatorWrapper.t.hh
 * \author hu4
 * \date   Wed Apr 01 12:39:40 2015
 * \brief  ModelEvaluatorWrapper template member definitions.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef SPn_solvers_ModelEvaluatorWrapper_t_hh
#define SPn_solvers_ModelEvaluatorWrapper_t_hh

#include "ModelEvaluatorWrapper.hh"

#include "comm/global.hh"
#include "spn/VectorTraits.hh"
#include "ThyraTraits.hh"

namespace profugus
{

//---------------------------------------------------------------------------//
//\brief Constructor
//---------------------------------------------------------------------------//
template <class T>
ModelEvaluatorWrapper<T>::ModelEvaluatorWrapper( Teuchos::RCP<const OP> op )
    : d_op(op)
{
    REQUIRE( d_op != Teuchos::null );

    // Wrap operator into Thyra LinearOpBase
    d_thyra_op = ThyraTraits<T>::buildThyraOP(d_op);
}

//---------------------------------------------------------------------------//
//\brief Return vector space corresponding to x
//---------------------------------------------------------------------------//
template <class T>
Teuchos::RCP<const Thyra::VectorSpaceBase<typename T::ST> >
ModelEvaluatorWrapper<T>::get_x_space() const
{
    return d_thyra_op->domain();
}

//---------------------------------------------------------------------------//
//\brief Return vector space corresponding to f
//---------------------------------------------------------------------------//
template <class T>
Teuchos::RCP<const Thyra::VectorSpaceBase<typename T::ST> >
ModelEvaluatorWrapper<T>::get_f_space() const
{
    return d_thyra_op->range();
}

//---------------------------------------------------------------------------//
//\brief Return InArgs supported by this evaluator
//---------------------------------------------------------------------------//
template <class T>
typename ModelEvaluatorWrapper<T>::InputArgs
ModelEvaluatorWrapper<T>::createInArgs() const
{
    // We only need x to evaluate function
    Thyra::ModelEvaluatorBase::InArgsSetup<double> inArgs;
    inArgs.setModelEvalDescription(this->description());
    inArgs.setSupports(Thyra::ModelEvaluatorBase::IN_ARG_x);
    return inArgs;
}

//---------------------------------------------------------------------------//
//\brief Return OutArgs supported by this evaluator
//---------------------------------------------------------------------------//
template <class T>
typename ModelEvaluatorWrapper<T>::OutputArgs
ModelEvaluatorWrapper<T>::createOutArgsImpl() const
{
    // We only evaluate f
    Thyra::ModelEvaluatorBase::OutArgsSetup<double> outArgs;
    outArgs.setModelEvalDescription(this->description());
    outArgs.setSupports(Thyra::ModelEvaluatorBase::OUT_ARG_f);
    return outArgs;
}

//---------------------------------------------------------------------------//
//\brief Evaluate model
//---------------------------------------------------------------------------//
template <class T>
void ModelEvaluatorWrapper<T>::evalModelImpl(const InputArgs  &inputs,
                                             const OutputArgs &outputs) const
{
    REQUIRE( inputs.supports( Thyra::ModelEvaluatorBase::IN_ARG_x) );
    REQUIRE( outputs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_f) );

    // Get x from inputs
    Teuchos::RCP<const Thyra::VectorBase<ST> > x = inputs.get_x();

    // Get f from outputs
    Teuchos::RCP<Thyra::VectorBase<ST> > f = outputs.get_f();

    // Apply operator to x
    Thyra::apply(*d_thyra_op,Thyra::NOTRANS,*x,f.ptr());
}

} // end namespace profugus

#endif // SPn_solvers_ModelEvaluatorWrapper_t_hh

//---------------------------------------------------------------------------//
//                 end of ModelEvaluatorWrapper.t.hh
//---------------------------------------------------------------------------//
