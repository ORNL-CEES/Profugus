//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   SPn/solvers/MueLuPreconditioner.hh
 * \author Thomas M. Evans, Steven Hamilton
 * \date   Fri Feb 21 13:41:13 2014
 * \brief  MueLuPreconditioner class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef SPn_solvers_MueLuPreconditioner_hh
#define SPn_solvers_MueLuPreconditioner_hh

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

// Only include Xpetra headers if MueLu is enabled
#ifdef USE_MUELU
#include "Xpetra_TpetraCrsMatrix.hpp"
#include "Xpetra_TpetraMultiVector.hpp"
#include "Xpetra_EpetraCrsMatrix.hpp"
#include "Xpetra_EpetraMultiVector.hpp"

#include "MueLu.hpp"
#include "MueLu_ParameterListInterpreter.hpp"
#include "MueLu_EpetraOperator.hpp"
#include "MueLu_TpetraOperator.hpp"
#endif

#include "LinAlgTypedefs.hh"

#include "harness/DBC.hh"

namespace profugus
{

//===========================================================================//
/*
 *! \class MueLuPreconditionerBase
 *  \brief Implementation details of MueLuPreconditioner.
 *
 *  This class does nothing more than take an Xpetra-wrapped matrix and
 *  construct the MueLu Hierarchy.
 */
//===========================================================================//
class MueLuPreconditionerBase
{
  public:

    MueLuPreconditionerBase();

  protected:

    void setup(Teuchos::RCP<Teuchos::ParameterList> pl);

    typedef typename TpetraTypes::ST   ST;
    typedef typename TpetraTypes::LO   LO;
    typedef typename TpetraTypes::GO   GO;
    typedef typename TpetraTypes::NODE NODE;

#ifdef USE_MUELU
    Teuchos::RCP<XpetraTypes::MATRIX>   d_matrix;
    Teuchos::RCP<MueLu::Hierarchy<ST,LO,GO,NODE> > d_hierarchy;
#endif
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
template <class T>
class MueLuPreconditioner
{
  public:

    MueLuPreconditioner();
};

#ifdef USE_MUELU
// Implementation for Epetra_MultiVector/Operator
template <>
class MueLuPreconditioner<EpetraTypes>
    : public Epetra_Operator,
      public MueLuPreconditionerBase
{
  public:
    //@{
    //! Typedefs.
    typedef typename EpetraTypes::ST    ST;
    typedef typename EpetraTypes::LO    LO;
    typedef typename EpetraTypes::GO    GO;
    typedef typename EpetraTypes::NODE  NODE;
    typedef typename EpetraTypes::MV    MV;
    typedef typename EpetraTypes::OP    OP;
    typedef Teuchos::RCP<OP>            RCP_Operator;
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

        // Wrap Epetra_CrsMatrix into an Xpetra::EpetraCrsMatrix
        Teuchos::RCP<XpetraTypes::CRS_MATRIX> temp_matrix =
            Teuchos::rcp( new Xpetra::EpetraCrsMatrix(A) );

        // Now wrap Xpetra_CrsMatrix into an Xpetra::Matrix
        // Why wrap twice?  Because it's twice as awesome
        d_matrix = Teuchos::rcp( new Xpetra::CrsMatrixWrap<ST,LO,GO,NODE>(
            temp_matrix) );
        this->setup(pl);
    };

    // Apply (solve linear system)
    int Apply(const MV &x, MV &y) const
    {
        REQUIRE( x.MyLength() == y.MyLength() );

        MueLu::EpetraOperator op_wrap( d_hierarchy );
        int err = op_wrap.ApplyInverse(x,y);
        REQUIRE( 0 == err );

        return err;
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
class MueLuPreconditioner<TpetraTypes>
    : public TpetraTypes::OP,
      public MueLuPreconditionerBase
{
  public:
    //@{
    //! Typedefs.
    typedef TpetraTypes::ST         ST;
    typedef TpetraTypes::LO         LO;
    typedef TpetraTypes::GO         GO;
    typedef TpetraTypes::NODE       NODE;
    typedef TpetraTypes::MV         MV;
    typedef TpetraTypes::OP         OP;
    typedef TpetraTypes::MAP        MAP;
    typedef TpetraTypes::MATRIX     MATRIX;
    typedef XpetraTypes::MATRIX     X_MATRIX;
    typedef XpetraTypes::CRS_MATRIX X_CRS_MATRIX;
    //@}

  private:

    typedef MueLuPreconditionerBase Base;
    using Base::d_matrix;
    Teuchos::RCP<MATRIX> d_A;

  public:

    // Constructor.
    explicit MueLuPreconditioner(Teuchos::RCP<MATRIX> A,
                                 Teuchos::RCP<Teuchos::ParameterList> pl)
    {
        d_A = A;

        // Wrap Epetra_CrsMatrix into an Xpetra::EpetraCrsMatrix
        Teuchos::RCP<X_CRS_MATRIX> temp_matrix = Teuchos::rcp(
                new Xpetra::TpetraCrsMatrix<ST,LO,GO,NODE>(A) );

        // Now wrap Xpetra_CrsMatrix into an Xpetra::Matrix
        // Why wrap twice?  Because it's twice as awesome
        d_matrix = Teuchos::rcp( new Xpetra::CrsMatrixWrap<ST,LO,GO,NODE>(
            temp_matrix) );

        this->setup(pl);
    };

    // Apply (solve linear system)
    void apply(const MV &x, MV &y,
               Teuchos::ETransp mode=Teuchos::NO_TRANS,
               double alpha=1.0, double beta=0.0) const
    {
        REQUIRE( x.getLocalLength() == y.getLocalLength() );

        MueLu::TpetraOperator<ST,LO,GO,NODE> op_wrap(d_hierarchy);

        op_wrap.apply(x,y,mode,alpha,beta);
    }

    // Required inherited interface.
    bool hasTranposeApply() const {return false;}
    Teuchos::RCP<const MAP> getDomainMap() const
    {
        REQUIRE(d_A != Teuchos::null);
        return d_A->getDomainMap();
    }
    Teuchos::RCP<const MAP> getRangeMap() const
    {
        REQUIRE(d_A != Teuchos::null);
        return d_A->getRangeMap();
    }

};
#endif // USE_MUELU

} // end namespace profugus

#endif // SPn_solvers_MueLuPreconditioner_hh

//---------------------------------------------------------------------------//
//                 end of MueLuPreconditioner.hh
//---------------------------------------------------------------------------//
