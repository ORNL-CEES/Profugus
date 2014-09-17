//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   spn/OperatorAdapters.hh
 * \author Thomas M. Evans, Steven Hamilton
 * \date   Tue Feb 25 12:35:37 2014
 * \brief  OperatorAdapters class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef spn_OperatorAdapters_hh
#define spn_OperatorAdapters_hh

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
 * \brief Mix-Ins for implementing Epetra/Tpetra operators
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
                            Epetra_MultiVector &y) const = 0;

    // Required interface
    int SetUseTranspose(bool use){return -1;}
    int ApplyInverse(const Epetra_MultiVector &x,
                           Epetra_MultiVector &y ) const {return -1;}
    bool HasNormInf()const {return false;}
    double NormInf() const {return 0.0;}
    const char * Label() const {return "EpetraAdapter";}
    bool UseTranspose() const {return false;}
    const Epetra_Comm & Comm() const {return d_domain_map->Comm();}
    const Epetra_Map & OperatorDomainMap() const
    {
        return *d_domain_map;
    }
    const Epetra_Map & OperatorRangeMap() const
    {
        if( d_range_map != Teuchos::null )
            return *d_range_map;
        return *d_domain_map;
    }

  private:

    Teuchos::RCP<const Epetra_Map> d_domain_map;
    Teuchos::RCP<const Epetra_Map> d_range_map;
};

template <>
class OperatorAdapter<TpetraTypes> : public TpetraTypes::OP
{
  public:

    typedef typename TpetraTypes::MAP Map_t;
    typedef typename TpetraTypes::MV  MV;

    OperatorAdapter( Teuchos::RCP<const Map_t>  domain_map,
                     Teuchos::RCP<const Map_t>  range_map=Teuchos::null )
        : d_domain_map(domain_map)
        , d_range_map(range_map)
    {
        REQUIRE( d_domain_map != Teuchos::null );
    }

    virtual void apply( const MV &x, MV &y,
        Teuchos::ETransp mode=Teuchos::NO_TRANS,
        double alpha=Teuchos::ScalarTraits<double>::one(),
        double beta=Teuchos::ScalarTraits<double>::zero()) const = 0;

    // Required interface
    bool hasTransposeApply(){return false;}
    Teuchos::RCP<const Map_t> getDomainMap() const
    {
        REQUIRE( d_domain_map!= Teuchos::null );
        return d_domain_map;
    }
    Teuchos::RCP<const Map_t> getRangeMap() const
    {
        if( d_range_map != Teuchos::null )
            return d_range_map;
        return d_domain_map;
    }

  private:

    Teuchos::RCP<const Map_t> d_domain_map;
    Teuchos::RCP<const Map_t> d_range_map;
};

} // end namespace profugus

#endif // spn_OperatorAdapters_hh

//---------------------------------------------------------------------------//
//                 end of OperatorAdapters.hh
//---------------------------------------------------------------------------//
