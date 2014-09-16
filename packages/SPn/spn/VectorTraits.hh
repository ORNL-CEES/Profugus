//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   spn/VectorTraits.hh
 * \author Steven Hamilton
 * \date   Mon Feb 17 13:10:27 2014
 * \brief  VectorTraits class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef spn_VectorTraits_hh
#define spn_VectorTraits_hh

#include "Teuchos_RCP.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_ArrayView.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_OrdinalTraits.hpp"

#include "harness/DBC.hh"
#include "solvers/LinAlgTypedefs.hh"


namespace profugus
{

//===========================================================================//
/*!
 * \class VectorTraits
 * \brief Traits class for Epetra/Tpetra matrices
 */
/*!
 * \example spn/test/tstVectorTraits.cc
 *
 * Test of VectorTraits.
 */
//===========================================================================//

template <class T>
class UndefinedVectorTraits
{
    void NotDefined(){T::this_class_is_missing_a_specialization();}
};

template <class T>
class VectorTraits
{
  public:
    //@{
    //! Typedefs.
    typedef typename T::MAP    Map_t;
    typedef typename T::VECTOR Vector_t;
    //@}

    static Teuchos::RCP<Vector_t> build_vector(Teuchos::RCP<const Map_t> map,
                                               int num_local)
    {
        UndefinedVectorTraits<T>::NotDefined();
        return Teuchos::null;
    }

    static int local_length(Teuchos::RCP<const Vector_t> vector)
    {
        UndefinedVectorTraits<T>::NotDefined();
        return 0;
    }

    static void put_scalar(Teuchos::RCP<Vector_t> vector, double val)
    {
        UndefinedVectorTraits<T>::NotDefined();
    }

    static Teuchos::ArrayView<double> get_data(
        Teuchos::RCP<Vector_t> vector)
    {
        UndefinedVectorTraits<T>::NotDefined();
        return Teuchos::ArrayView<double>();
    }

};

// Specialization on EpetraTypes
template <>
class VectorTraits<EpetraTypes>
{
  public:
    //@{
    //! Typedefs.
    typedef typename EpetraTypes::MAP    Map_t;
    typedef typename EpetraTypes::VECTOR Vector_t;
    //@}

    static Teuchos::RCP<Vector_t> build_vector(Teuchos::RCP<const Map_t> map,
                                               int num_local)
    {

        Teuchos::RCP<Vector_t> x(new Vector_t(*map));
        CHECK(x->MyLength() == num_local);
        return x;
    }

    static int local_length(Teuchos::RCP<const Vector_t> vector)
    {
        return vector->MyLength();
    }

    static void put_scalar(Teuchos::RCP<Vector_t> vector, double val)
    {
        vector->PutScalar(val);
    }

    static Teuchos::ArrayView<double> get_data(
        Teuchos::RCP<Vector_t> vector)
    {
        return Teuchos::arrayView<double>(vector->Values(),vector->MyLength());
    }

};

// Specialization on TpetraTypes
template <>
class VectorTraits<TpetraTypes>
{
  public:
    //@{
    //! Typedefs.
    typedef typename TpetraTypes::MAP    Map_t;
    typedef typename TpetraTypes::VECTOR Vector_t;
    //@}

    static Teuchos::RCP<Vector_t> build_vector(Teuchos::RCP<const Map_t> map,
                                               int num_local)
    {

        Teuchos::RCP<Vector_t> x(new Vector_t(map));
        CHECK(x->getLocalLength() == num_local);
        return x;
    }

    static int local_length(Teuchos::RCP<const Vector_t> vector)
    {
        return vector->getLocalLength();
    }

    static void put_scalar(Teuchos::RCP<Vector_t> vector, double val)
    {
        vector->putScalar(val);
    }

    static Teuchos::ArrayView<double> get_data(
        Teuchos::RCP<Vector_t> vector)
    {
        return vector->getDataNonConst()();
    }

};

} // end namespace profugus

#endif // spn_VectorTraits_hh

//---------------------------------------------------------------------------//
//                 end of VectorTraits.hh
//---------------------------------------------------------------------------//
