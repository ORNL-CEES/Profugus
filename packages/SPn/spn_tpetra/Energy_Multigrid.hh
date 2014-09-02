//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   spn_tpetra/Energy_Multigrid.hh
 * \author Thomas M. Evans, Steven Hamilton
 * \date   Tue Feb 25 13:04:00 2014
 * \brief  Energy_Multigrid class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef spn_tpetra_Energy_Multigrid_hh
#define spn_tpetra_Energy_Multigrid_hh

#include <vector>

#include "Teuchos_RCP.hpp"
#include "Tpetra_Operator.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_RowMatrix.hpp"

#include "harness/DBC.hh"
#include "xs/Mat_DB.hh"
#include "mesh/Mesh.hh"
#include "mesh/LG_Indexer.hh"
#include "mesh/Global_Mesh_Data.hh"
#include "solvers/StratimikosSolver.hh"
#include "solvers/LinearSolver.hh"
#include "solvers/LinearSolverBuilder.hh"
#include "spn/Dimensions.hh"
#include "Linear_System.hh"
#include "Energy_Restriction.hh"
#include "Energy_Prolongation.hh"

#include "solvers/TpetraTypedefs.hh"

namespace profugus
{
namespace tpetra
{

//===========================================================================//
/*!
 * \class Energy_Multigrid
 * \brief Multigrid in energy preconditioner for SPN
 *
 * \sa Energy_Multigrid.cc for detailed descriptions.
 */
/*!
 * \example spn_tpetra/test/tstEnergy_Multigrid.cc
 *
 * Test of Energy_Multigrid.
 */
//===========================================================================//

class Energy_Multigrid : public Tpetra_Operator
{
  public:
    //@{
    //! Typedefs.
    typedef Tpetra_Map                        Map_t;
    typedef Tpetra_Operator                   OP;
    typedef Tpetra_MultiVector                MV;
    typedef Tpetra_Vector                     Vector_t;
    typedef LinearSolver<MV,OP>               LinearSolver_t;
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

    void apply( const MV &x, MV &y, Teuchos::ETransp mode=Teuchos::NO_TRANS,
                double alpha=Teuchos::ScalarTraits<double>::one(),
                double beta=Teuchos::ScalarTraits<double>::zero()) const;

    // Required interface
    bool hasTransposeApply(){return false;}
    Teuchos::RCP<const Map_t> getDomainMap() const
    {
        REQUIRE( d_operators[0] != Teuchos::null );
        return d_operators[0]->getDomainMap();
    }
    Teuchos::RCP<const Map_t> getRangeMap() const
    {
        REQUIRE( d_operators[0] != Teuchos::null );
        return d_operators[0]->getRangeMap();
    }

  private:

    int d_num_levels;
    std::vector< Teuchos::RCP<OP> >  d_operators;
    std::vector< Teuchos::RCP<OP> >  d_restrictions;
    std::vector< Teuchos::RCP<OP> >  d_prolongations;
    std::vector< Teuchos::RCP<OP> >  d_preconditioners;
    std::vector< Teuchos::RCP<MV> >    d_solutions;
    std::vector< Teuchos::RCP<MV> >    d_residuals;
    std::vector< Teuchos::RCP<MV> >    d_rhss;
    std::vector< Teuchos::RCP<LinearSolver_t> >   d_smoothers;
};

} // end namespace profugus
} // end namespace tpetra

#endif // spn_tpetra_Energy_Multigrid_hh

//---------------------------------------------------------------------------//
//                 end of Energy_Multigrid.hh
//---------------------------------------------------------------------------//
