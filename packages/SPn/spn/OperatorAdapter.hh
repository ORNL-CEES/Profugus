//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   SPn/spn/OperatorAdapter.hh
 * \author Thomas M. Evans, Steven Hamilton
 * \date   Tue Feb 25 12:35:37 2014
 * \brief  OperatorAdapters class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef SPn_spn_OperatorAdapter_hh
#define SPn_spn_OperatorAdapter_hh

#include <vector>

#include "Teuchos_RCP.hpp"
#include "AnasaziMultiVecTraits.hpp"
#include "AnasaziEpetraAdapter.hpp"
#include "AnasaziTpetraAdapter.hpp"

#include "solvers/LinAlgTypedefs.hh"

namespace profugus
{

//===========================================================================//
/*!
 * \class OperatorAdapters
 * \brief Mix-In for implementing Epetra/Tpetra operators
 *
 * This class provides a unified interface for Epetra_Operator and
 * Tpetra::Operator by handling the various details of the operator
 * interfaces.  All the user needs to do is inherit from this class
 * (which is templated on either EpetraTypes or TpetraTypes) and implement
 * the ApplyImpl function, i.e.
 * \code
   template <class T>
   class ConcreteOperator : OperatorAdapter<T>
   {
      ConcreteOperator( RCP<T::MAP> map )
        : OperatorAdapter<T>(map)
      {}

      void ApplyImpl( const T::MV &x, T::MV &y ) const
      {
          // Do something with x and y...
      }
   };
   \endcode
 *
 * Note that the OperatorAdapter does not have a default constructor, a
 * map must be provided at construction time (and therefore the derived
 * class must have a map available at construction time).
 */
//===========================================================================//

template <class T>
class UndefinedOperatorAdapter
{
    void notDefined(){T::this_class_is_missing_a_specialization();}
};

template <class T>
class OperatorAdapter
{
    OperatorAdapter(){UndefinedOperatorAdapter<T>::notDefined();}
};

template <>
class OperatorAdapter<EpetraTypes> : public EpetraTypes::OP
{
  public:

    OperatorAdapter( Teuchos::RCP<const Epetra_Map>  domain_map,
                     Teuchos::RCP<const Epetra_Map>  range_map=Teuchos::null )
        : d_domain_map(domain_map)
        , d_range_map(range_map)
    {
        REQUIRE( d_domain_map != Teuchos::null );
    }

    virtual int Apply(const Epetra_MultiVector &x,
                            Epetra_MultiVector &y) const override
    {
        this->ApplyImpl(x,y);
        return 0;
    }

    // Required interface
    virtual int SetUseTranspose(bool use) override {return -1;}
    virtual bool UseTranspose() const override {return false;}
    int ApplyInverse(const Epetra_MultiVector &x,
                           Epetra_MultiVector &y ) const override
    {
        return -1;
    }
    bool HasNormInf()const override {return false;}
    double NormInf() const override {return 0.0;}
    const char * Label() const override {return "EpetraAdapter";}
    const Epetra_Comm & Comm() const override {return d_domain_map->Comm();}
    const Epetra_Map & OperatorDomainMap() const override
    {
        return *d_domain_map;
    }
    const Epetra_Map & OperatorRangeMap() const override
    {
        if( d_range_map != Teuchos::null )
            return *d_range_map;
        return *d_domain_map;
    }

  protected:

    virtual void ApplyImpl(const Epetra_MultiVector &x,
                                 Epetra_MultiVector &y) const = 0;

    Teuchos::RCP<const Epetra_Map> d_domain_map;
    Teuchos::RCP<const Epetra_Map> d_range_map;
};

template <>
class OperatorAdapter<TpetraTypes> : public TpetraTypes::OP
{
  public:

    typedef typename TpetraTypes::MAP                   Map_t;
    typedef typename TpetraTypes::MV                    MV;
    typedef typename Anasazi::MultiVecTraits<double,MV> MVT;

    OperatorAdapter( Teuchos::RCP<const Map_t>  domain_map,
                     Teuchos::RCP<const Map_t>  range_map=Teuchos::null )
        : d_domain_map(domain_map)
        , d_range_map(range_map)
    {
        REQUIRE( d_domain_map != Teuchos::null );
    }

    // Tpetra::Operator apply
    virtual void apply( const MV &x, MV &y,
        Teuchos::ETransp mode=Teuchos::NO_TRANS,
        double alpha=Teuchos::ScalarTraits<double>::one(),
        double beta=Teuchos::ScalarTraits<double>::zero()) const override
    {
        REQUIRE( mode == Teuchos::NO_TRANS );
        REQUIRE( MVT::GetNumberVecs(x) == MVT::GetNumberVecs(y) );

        // If alpha is zero, don't apply operator, just scale and return
        if( alpha==Teuchos::ScalarTraits<double>::zero() )
        {
            MVT::MvScale(y,beta);
            return;
        }

        // If beta is zero, do apply and scale
        if( beta==Teuchos::ScalarTraits<double>::zero() )
        {
            this->ApplyImpl(x,y);
            MVT::MvScale(y,alpha);
        }
        // For nonzero beta, need temporary vector
        else
        {
            Teuchos::RCP<MV> z = MVT::Clone(x,1);
            this->ApplyImpl(x,*z);
            MVT::MvAddMv(alpha,*z,beta,y,y);
        }

    }

    // Required Tpetra::Operator interface
    virtual bool hasTransposeApply() const override {return false;}
    Teuchos::RCP<const Map_t> getDomainMap() const override
    {
        REQUIRE( d_domain_map!= Teuchos::null );
        return d_domain_map;
    }
    Teuchos::RCP<const Map_t> getRangeMap() const override
    {
        if( d_range_map != Teuchos::null )
            return d_range_map;
        return d_domain_map;
    }

  protected:

    virtual void ApplyImpl(const MV &x, MV &y) const = 0;

    Teuchos::RCP<const Map_t> d_domain_map;
    Teuchos::RCP<const Map_t> d_range_map;
};

} // end namespace profugus

#endif // SPn_spn_OperatorAdapter_hh

//---------------------------------------------------------------------------//
// end of OperatorAdapter.hh
//---------------------------------------------------------------------------//
