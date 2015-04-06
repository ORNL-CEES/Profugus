//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   solvers/ThyraTraits.hh
 * \author hu4
 * \date   Thu Apr 02 10:24:48 2015
 * \brief  ThyraTraits class definition.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef solvers_ThyraTraits_hh
#define solvers_ThyraTraits_hh

#include "Teuchos_RCP.hpp"
#include "Thyra_VectorSpaceBase.hpp"
#include "Thyra_VectorBase.hpp"
#include "Thyra_MultiVectorBase.hpp"
#include "Thyra_EpetraLinearOp.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_TpetraLinearOp.hpp"
#include "Thyra_TpetraThyraWrappers.hpp"

#include "harness/DBC.hh"
#include "LinAlgTypedefs.hh"

namespace profugus
{

template <class T>
class UndefinedThyraTraits
{
    void NotDefined(){T::this_class_is_missing_a_specialization();}
};

//===========================================================================//
/*!
 * \class ThyraTraits
 * \brief Traits class for interoperability of Epetra/Tpetra with Thyra.
 *
 * Considering a fundamental idea of Thyra is to provide interoperability
 * of different concrete operator/vector types, it is absurb that this class
 * needs to exist.
 *
 */
//===========================================================================//
template <class T>
class ThyraTraits
{
  public:

    typedef typename T::ST ST;
    typedef typename T::MV MV;
    typedef typename T::OP OP;
    typedef Thyra::MultiVectorBase<ST> ThyraMV;
    typedef Thyra::LinearOpBase<ST>    ThyraOP;

    //\brief Wrap MV into a Thyra MultiVector
    static Teuchos::RCP<ThyraMV> buildThyraMV( Teuchos::RCP<MV> x,
        Teuchos::RCP<const Thyra::VectorSpaceBase<ST> > space)
    {
        UndefinedThyraTraits<T>::NotDefined();
        return Teuchos::null;
    }

    //\brief Wrap const MV into a const Thyra MultiVector
    static Teuchos::RCP<const ThyraMV> buildThyraConstMV(
        Teuchos::RCP<const MV> x,
        Teuchos::RCP<const Thyra::VectorSpaceBase<ST> > space)
    {
        UndefinedThyraTraits<T>::NotDefined();
        return Teuchos::null;
    }

    //\brief
    static Teuchos::RCP<const ThyraOP> buildThyraOP(
        Teuchos::RCP<const OP> op)
    {
        UndefinedThyraTraits<T>::NotDefined();
        return Teuchos::null;
    }
};

// Specialization for Epetra
template <>
class ThyraTraits<EpetraTypes>
{
  public:

    typedef EpetraTypes::ST ST;
    typedef EpetraTypes::MV MV;
    typedef EpetraTypes::OP OP;
    typedef Thyra::MultiVectorBase<ST> ThyraMV;
    typedef Thyra::LinearOpBase<ST>    ThyraOP;

    //\brief Wrap MV into a Thyra MultiVector
    static Teuchos::RCP<ThyraMV> buildThyraMV( Teuchos::RCP<MV> x,
        Teuchos::RCP<const Thyra::VectorSpaceBase<ST> > space)
    {
        REQUIRE( x != Teuchos::null );
        REQUIRE( space != Teuchos::null );
        return Thyra::create_MultiVector(x,space);
    }

    //\brief Wrap const MV into a const Thyra MultiVector
    static Teuchos::RCP<const ThyraMV> buildThyraConstMV(
        Teuchos::RCP<const MV> x,
        Teuchos::RCP<const Thyra::VectorSpaceBase<ST> > space)
    {
        REQUIRE( x != Teuchos::null );
        REQUIRE( space != Teuchos::null );
        return Thyra::create_MultiVector(x,space);
    }

    //\brief Wrap operator into Thyra LinearOpBase
    static Teuchos::RCP<const ThyraOP> buildThyraOP(
        Teuchos::RCP<const OP> op)
    {
        REQUIRE( op != Teuchos::null );
        return Thyra::epetraLinearOp(op);
    }
};

// Specialization for Tpetra
template <>
class ThyraTraits<TpetraTypes>
{
  public:

    typedef TpetraTypes::ST   ST;
    typedef TpetraTypes::LO   LO;
    typedef TpetraTypes::GO   GO;
    typedef TpetraTypes::NODE NODE;
    typedef TpetraTypes::MV   MV;
    typedef TpetraTypes::OP   OP;
    typedef Thyra::MultiVectorBase<ST> ThyraMV;
    typedef Thyra::LinearOpBase<ST>    ThyraOP;

    //\brief Wrap MV into a Thyra MultiVector
    static Teuchos::RCP<ThyraMV> buildThyraMV( Teuchos::RCP<MV> x,
        Teuchos::RCP<const Thyra::VectorSpaceBase<ST> > space)
    {
        REQUIRE( x != Teuchos::null );
        REQUIRE( space != Teuchos::null );
        return Thyra::createMultiVector(x,space);
    }

    //\brief Wrap const MV into a const Thyra MultiVector
    static Teuchos::RCP<const ThyraMV> buildThyraConstMV(
        Teuchos::RCP<const MV> x,
        Teuchos::RCP<const Thyra::VectorSpaceBase<ST> > space)
    {
        REQUIRE( x != Teuchos::null );
        REQUIRE( space != Teuchos::null );
        return Thyra::createConstMultiVector(x,space);
    }

    //\brief Wrap operator into Thyra LinearOpBase
    static Teuchos::RCP<const ThyraOP> buildThyraOP(
        Teuchos::RCP<const OP> op)
    {
        REQUIRE( op != Teuchos::null );
        auto range  = Thyra::tpetraVectorSpace<ST>(op->getRangeMap());
        auto domain = Thyra::tpetraVectorSpace<ST>(op->getDomainMap());
        return Thyra::createConstLinearOp<ST,LO,GO,NODE>(op,range,domain);
    }
};

} // end namespace profugus

#endif // solvers_ThyraTraits_hh

//---------------------------------------------------------------------------//
//                 end of ThyraTraits.hh
//---------------------------------------------------------------------------//
