//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   solvers/ModelEvaluatorWrapper.hh
 * \author Steven Hamilton
 * \date   Wed Apr 01 12:39:40 2015
 * \brief  ModelEvaluatorWrapper class definition.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef solvers_ModelEvaluatorWrapper_hh
#define solvers_ModelEvaluatorWrapper_hh

#include "Thyra_StateFuncModelEvaluatorBase.hpp"
#include "Thyra_LinearOpBase.hpp"

namespace profugus
{

//===========================================================================//
/*!
 * \class ModelEvaluatorWrapper
 * \brief Wrap an Epetra/Tpetra operator into a Thyra ModelEvaluator.
 *
 * \sa ModelEvaluatorWrapper.t.hh for detailed descriptions.
 */
/*!
 * \example solvers/test/tstModelEvaluatorWrapper.cc
 *
 * Test of ModelEvaluatorWrapper.
 */
//===========================================================================//

template <class T>
class ModelEvaluatorWrapper :
    public Thyra::StateFuncModelEvaluatorBase<typename T::ST>
{
  public:

    typedef typename T::ST ST;
    typedef typename T::OP OP;
    typedef Thyra::ModelEvaluatorBase::InArgs<ST>  InputArgs;
    typedef Thyra::ModelEvaluatorBase::OutArgs<ST> OutputArgs;

    ModelEvaluatorWrapper( Teuchos::RCP<const OP> op );

    // Thyra interface functions
    Teuchos::RCP<const Thyra::VectorSpaceBase<ST> > get_x_space() const;
    Teuchos::RCP<const Thyra::VectorSpaceBase<ST> > get_f_space() const;
    InputArgs  createInArgs() const;
    OutputArgs createOutArgsImpl() const;
    void evalModelImpl(const InputArgs  &inputs,
                       const OutputArgs &outputs) const;

  private:

    Teuchos::RCP<const OP> d_op;
    Teuchos::RCP<const Thyra::LinearOpBase<ST> > d_thyra_op;
};

} // end namespace profugus

#endif // solvers_ModelEvaluatorWrapper_hh

//---------------------------------------------------------------------------//
//                 end of ModelEvaluatorWrapper.hh
//---------------------------------------------------------------------------//
