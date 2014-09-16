//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   spn/Energy_Multigrid.hh
 * \author Thomas M. Evans, Steven Hamilton
 * \date   Tue Feb 25 13:04:00 2014
 * \brief  Energy_Multigrid class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef spn_Energy_Multigrid_hh
#define spn_Energy_Multigrid_hh

#include <vector>

#include "Teuchos_RCP.hpp"
#include "Epetra_Operator.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_RowMatrix.h"

#include "harness/DBC.hh"
#include "xs/Mat_DB.hh"
#include "mesh/Mesh.hh"
#include "mesh/LG_Indexer.hh"
#include "mesh/Global_Mesh_Data.hh"
#include "solvers/StratimikosSolver.hh"
#include "solvers/LinearSolver.hh"
#include "solvers/LinearSolverBuilder.hh"
#include "solvers/LinAlgTypedefs.hh"
#include "Dimensions.hh"
#include "Linear_System.hh"
#include "Energy_Restriction.hh"
#include "Energy_Prolongation.hh"

namespace profugus
{

//===========================================================================//
/*!
 * \class Energy_Multigrid
 * \brief Multigrid in energy preconditioner for SPN
 *
 * \sa Energy_Multigrid.cc for detailed descriptions.
 */
/*!
 * \example spn/test/tstEnergy_Multigrid.cc
 *
 * Test of Energy_Multigrid.
 */
//===========================================================================//

class Energy_Multigrid : public Epetra_Operator
{
  public:
    //@{
    //! Typedefs.
    typedef Epetra_Operator                   OP;
    typedef Epetra_MultiVector                MV;
    typedef LinearSolver<EpetraTypes>              LinearSolver_t;
    typedef LinearSolver_t::RCP_ParameterList RCP_ParameterList;
    typedef LinearSolver_t::ParameterList     ParameterList;
    //@}

  public:
    // Constructor.
    Energy_Multigrid( RCP_ParameterList              main_db,
                      RCP_ParameterList              prec_db,
                      Teuchos::RCP<Dimensions>       dim,
                      Teuchos::RCP<Mat_DB>           mat_db,
                      Teuchos::RCP<Mesh>             mesh,
                      Teuchos::RCP<LG_Indexer>       indexer,
                      Teuchos::RCP<Global_Mesh_Data> data,
                      Teuchos::RCP<Linear_System>    fine_system );

    int Apply( const Epetra_MultiVector &x,
                     Epetra_MultiVector &y ) const;

    // Required interface
    int SetUseTranspose(bool use){return -1;}
    int ApplyInverse(const Epetra_MultiVector &x,
                           Epetra_MultiVector &y ) const
    {
        INSIST(false,"Energy_Multigrid should use Apply, not ApplyInverse.");
        return -1;
    }
    bool HasNormInf()const {return false;}
    double NormInf() const {return 0.0;}
    const char * Label() const {return "Energy_Multigrid";}
    bool UseTranspose() const {return false;}
    const Epetra_Comm & Comm() const {return d_solutions[0]->Comm();}
    const Epetra_Map & OperatorDomainMap() const
    {
        REQUIRE( d_restrictions[0] != Teuchos::null );
        return d_operators[0]->OperatorDomainMap();
    }
    const Epetra_Map & OperatorRangeMap() const
    {
        REQUIRE( d_operators[0] != Teuchos::null );
        return d_operators[0]->OperatorRangeMap();
    }

  private:

    int d_num_levels;
    std::vector< Teuchos::RCP<Epetra_Operator> >  d_operators;
    std::vector< Teuchos::RCP<Epetra_Operator> >  d_restrictions;
    std::vector< Teuchos::RCP<Epetra_Operator> >  d_prolongations;
    std::vector< Teuchos::RCP<Epetra_Operator> >  d_preconditioners;
    std::vector< Teuchos::RCP<Epetra_Vector> >    d_solutions;
    std::vector< Teuchos::RCP<Epetra_Vector> >    d_residuals;
    std::vector< Teuchos::RCP<Epetra_Vector> >    d_rhss;
    std::vector< Teuchos::RCP<LinearSolver_t> >   d_smoothers;
};

} // end namespace profugus

#endif // spn_Energy_Multigrid_hh

//---------------------------------------------------------------------------//
//                 end of Energy_Multigrid.hh
//---------------------------------------------------------------------------//
