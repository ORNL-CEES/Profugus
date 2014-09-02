//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   solvers/MueLuPreconditioner.hh
 * \author Thomas M. Evans, Steven Hamilton
 * \date   Fri Feb 21 13:41:13 2014
 * \brief  MueLuPreconditioner class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef solvers_MueLuPreconditioner_hh
#define solvers_MueLuPreconditioner_hh

#include <string>

#include <SPn/config.h>

#include "AnasaziMultiVecTraits.hpp"
#include "AnasaziOperatorTraits.hpp"
#include "AnasaziEpetraAdapter.hpp"
#include "AnasaziTpetraAdapter.hpp"
#include "Epetra_Operator.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Map.h"
#include "Teuchos_RCP.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Operator.hpp"

#include "Xpetra_TpetraCrsMatrix.hpp"
#include "Xpetra_TpetraMultiVector.hpp"
#include "Xpetra_EpetraCrsMatrix.hpp"
#include "Xpetra_EpetraMultiVector.hpp"

#include "TpetraTypedefs.hh"

#include "harness/DBC.hh"

namespace profugus
{

//===========================================================================//
/*
 *! \class MueLuPreconditionerBase
 *  \brief Implementation details of MueLuPreconditioner.
 */
//===========================================================================//
class MueLuPreconditionerBase
{
  public:

    MueLuPreconditionerBase();

  protected:

    void setup(Teuchos::RCP<Teuchos::ParameterList> pl);
    void ApplyImpl( Teuchos::RCP<const Xpetra_MultiVector> x,
                    Teuchos::RCP<Xpetra_MultiVector>       y ) const;

    Teuchos::RCP<Xpetra_CrsMatrix> d_matrix;
};

//===========================================================================//
/*!
 * \class MueLuPreconditioner
 * \brief Wrap MueLu into generic Epetra/Tpetra operator interface
 *
 * MueLu uses the Xpetra interface.  In order to have a standard operator
 * interface, we have to wrap the relevant MueLu objects and handle the
 * {E/T->X}petra conversion under the hood.
 */
/*!
 * \example solvers/test/tstMueLuPreconditioner.cc
 *
 * Test of MueLuPreconditioner.
 */
//===========================================================================//

// Dummy implementation for MV/OP combos lacking a specialization.
// Attempt to instantiate will cause a compile error.
template <class MV, class OP>
class MueLuPreconditioner
{
  public:

    MueLuPreconditioner();
};

// Implementation for Epetra_MultiVector/Operator
template <>
class MueLuPreconditioner<Epetra_MultiVector,Epetra_Operator>
    : public Epetra_Operator,
      public MueLuPreconditionerBase
{
  public:
    //@{
    //! Typedefs.
    typedef Epetra_MultiVector MV;
    typedef Epetra_Operator    OP;
    typedef Teuchos::RCP<OP>   RCP_Operator;
    //@}

  private:

    typedef MueLuPreconditionerBase Base;
    using Base::d_matrix;
    Teuchos::RCP<Epetra_CrsMatrix> d_A;

  public:

    // Constructor.
    explicit MueLuPreconditioner(Teuchos::RCP<Epetra_CrsMatrix> A,
                                 Teuchos::RCP<Teuchos::ParameterList> pl)
    {
        d_A = A;
        d_matrix = Teuchos::rcp( new Xpetra::EpetraCrsMatrix(A) );
        this->setup(pl);
    };

    // Apply (solve linear system)
    int Apply(const MV &x, MV &y) const
    {
        REQUIRE( x.MyLength() == y.MyLength() );

        Teuchos::RCP<Epetra_MultiVector> x_rcp =
            Teuchos::rcp_const_cast<Epetra_MultiVector>(
                Teuchos::rcpFromRef(x));
        Teuchos::RCP<const Xpetra_MultiVector> x_wrap =
            Xpetra::toXpetra(x_rcp);
        REQUIRE( x_wrap != Teuchos::null );
        Teuchos::RCP<Xpetra_MultiVector> y_wrap =
            Xpetra::toXpetra(Teuchos::rcpFromRef(y));
        REQUIRE( y_wrap != Teuchos::null );

        ApplyImpl(x_wrap,y_wrap);

        return 0;
    }

    // Required inherited interface.
    int SetUseTranspose(bool useTranspose){return -1;}
    int ApplyInverse(const MV &x, MV &b) const
    {
        return -1;
    }
    double NormInf() const { return 0.0; }
    const char *  Label() const { return "MueLuPreconditioner"; }
    bool UseTranspose() const { return false; }
    bool HasNormInf() const { return false; }
    const Epetra_Comm & Comm() const
    {
        REQUIRE(d_A != Teuchos::null);
        return d_A->Comm();
    }
    const Epetra_Map & OperatorDomainMap() const
    {
        REQUIRE(d_A != Teuchos::null);
        return d_A->OperatorDomainMap();
    }
    const Epetra_Map & OperatorRangeMap() const
    {
        REQUIRE(d_A != Teuchos::null);
        return d_A->OperatorRangeMap();
    }

};

// Implementation for Tpetra::MultiVector/Operator
template <>
class MueLuPreconditioner<Tpetra_MultiVector,Tpetra_Operator>
    : public Tpetra_Operator,
      public MueLuPreconditionerBase
{
  public:
    //@{
    //! Typedefs.
    typedef Tpetra_MultiVector MV;
    typedef Tpetra_Operator    OP;
    typedef Tpetra_Map         Map;
    typedef Teuchos::RCP<OP>   RCP_Operator;
    //@}

  private:

    typedef MueLuPreconditionerBase Base;
    using Base::d_matrix;
    Teuchos::RCP<Tpetra_CrsMatrix> d_A;

  public:

    // Constructor.
    explicit MueLuPreconditioner(Teuchos::RCP<Tpetra_CrsMatrix> A,
                                 Teuchos::RCP<Teuchos::ParameterList> pl)
    {
        d_A = A;
        d_matrix = Teuchos::rcp(
            new Xpetra::TpetraCrsMatrix<SCALAR,LO,GO,NODE>(A) );
        this->setup(pl);
    };

    // Apply (solve linear system)
    void apply(const MV &x, MV &y,
               Teuchos::ETransp mode=Teuchos::NO_TRANS,
               double alpha=1.0, double beta=0.0) const
    {
        REQUIRE( x.getLocalLength() == y.getLocalLength() );

        MV z(x,Teuchos::Copy);

        Teuchos::RCP<Tpetra_MultiVector> x_rcp =
            Teuchos::rcp_const_cast<Tpetra_MultiVector>(
                Teuchos::rcpFromRef(x));
        Teuchos::RCP<const Xpetra_MultiVector> x_wrap =
            Xpetra::toXpetra(x_rcp);
        REQUIRE( x_wrap != Teuchos::null );
        Teuchos::RCP<Xpetra_MultiVector> y_wrap =
            Xpetra::toXpetra(Teuchos::rcpFromRef(y));
        REQUIRE( y_wrap != Teuchos::null );

        ApplyImpl(x_wrap,y_wrap);

        if( beta == 0.0 )
        {
            y.scale(alpha);
        }
        else
        {
            y.update(beta,z,alpha);
        }
    }

    // Required inherited interface.
    bool hasTranposeApply() const {return false;}
    Teuchos::RCP<const Map> getDomainMap() const
    {
        REQUIRE(d_A != Teuchos::null);
        return d_A->getDomainMap();
    }
    Teuchos::RCP<const Map> getRangeMap() const
    {
        REQUIRE(d_A != Teuchos::null);
        return d_A->getRangeMap();
    }

};

} // end namespace profugus

#endif // solvers_MueLuPreconditioner_hh

//---------------------------------------------------------------------------//
//                 end of MueLuPreconditioner.hh
//---------------------------------------------------------------------------//
