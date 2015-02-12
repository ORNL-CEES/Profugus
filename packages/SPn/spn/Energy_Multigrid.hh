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
#include "OperatorAdapter.hh"

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

template <class T>
class Energy_Multigrid : public OperatorAdapter<T>
{
  public:
    //@{
    //! Typedefs.
    typedef typename T::OP                             OP;
    typedef typename T::MV                             MV;
    typedef typename T::MAP                            MAP;
    typedef Anasazi::OperatorTraits<double,MV,OP>      OPT;
    typedef Anasazi::MultiVecTraits<double,MV>         MVT;
    typedef LinearSolver<T>                            LinearSolver_t;
    typedef typename LinearSolver_t::RCP_ParameterList RCP_ParameterList;
    typedef typename LinearSolver_t::ParameterList     ParameterList;
    //@}

  public:
    // Constructor.
    Energy_Multigrid( RCP_ParameterList               main_db,
                      RCP_ParameterList               prec_db,
                      Teuchos::RCP<Dimensions>        dim,
                      Teuchos::RCP<Mat_DB>            mat_db,
                      Teuchos::RCP<Mesh>              mesh,
                      Teuchos::RCP<LG_Indexer>        indexer,
                      Teuchos::RCP<Global_Mesh_Data>  data,
                      Teuchos::RCP<Linear_System<T> > fine_system );

    int Apply( const MV &x, MV &y ) const
    {
        ApplyImpl(x,y);
        return 0;
    }

    void apply( const MV &x, MV &y, Teuchos::ETransp mode=Teuchos::NO_TRANS,
                double alpha=Teuchos::ScalarTraits<double>::one(),
                double beta=Teuchos::ScalarTraits<double>::zero()) const
    {
        REQUIRE( alpha == 1.0 );
        REQUIRE( beta  == 0.0 );
        REQUIRE( mode == Teuchos::NO_TRANS );

        ApplyImpl(x,y);
    }


  private:

    void ApplyImpl(const MV &x, MV &y) const;

    int d_num_levels;
    std::vector< Teuchos::RCP<const MAP> >      d_maps;
    std::vector< Teuchos::RCP<OP> >             d_operators;
    std::vector< Teuchos::RCP<OP> >             d_restrictions;
    std::vector< Teuchos::RCP<OP> >             d_prolongations;
    std::vector< Teuchos::RCP<OP> >             d_preconditioners;
    std::vector< Teuchos::RCP<MV> >             d_solutions;
    std::vector< Teuchos::RCP<MV> >             d_residuals;
    std::vector< Teuchos::RCP<MV> >             d_rhss;
    std::vector< Teuchos::RCP<LinearSolver_t> > d_smoothers;
};

} // end namespace profugus

#endif // spn_Energy_Multigrid_hh

//---------------------------------------------------------------------------//
//                 end of Energy_Multigrid.hh
//---------------------------------------------------------------------------//
