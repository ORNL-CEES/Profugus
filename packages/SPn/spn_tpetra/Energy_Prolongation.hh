//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   spn_tpetra/Energy_Prolongation.hh
 * \author Thomas M. Evans, Steven Hamilton
 * \date   Tue Feb 25 12:35:44 2014
 * \brief  Energy_Prolongation class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef spn_tpetra_Energy_Prolongation_hh
#define spn_tpetra_Energy_Prolongation_hh

#include <vector>

#include "Teuchos_RCP.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Operator.hpp"
#include "Tpetra_Map.hpp"

#include "solvers/TpetraTypedefs.hh"

namespace profugus
{
namespace tpetra
{

//===========================================================================//
/*!
 * \class Energy_Prolongation
 * \brief Perform energy multigrid restriction
 *
 * This class performs a prolongation operation for the
 * SPN equations.  It is assumed that energy is the innermost variable.
 *
 * \sa Energy_Prolongation.cc for detailed descriptions.
 */
/*!
 * \example spn_tpetra/test/tstEnergy_Prolongation.cc
 *
 * Test of Energy_Prolongation.
 */
//===========================================================================//

class Energy_Prolongation : public Tpetra_Operator
{
  public:
    typedef Tpetra_Map         Map_t;
    typedef Tpetra_Operator    OP;
    typedef Tpetra_MultiVector MV;
    typedef Tpetra_Vector      Vector_t;

    Energy_Prolongation( Teuchos::RCP<const MV>  fine_vec,
                         Teuchos::RCP<const MV>  coarse_vec,
                         const std::vector<int> &steer_vec );

    void apply( const MV &x, MV &y, Teuchos::ETransp mode=Teuchos::NO_TRANS,
                double alpha=Teuchos::ScalarTraits<double>::one(),
                double beta=Teuchos::ScalarTraits<double>::zero()) const;

    // Required interface
    bool hasTransposeApply(){return false;}
    Teuchos::RCP<const Map_t> getDomainMap() const
    {
        REQUIRE( d_coarse_map != Teuchos::null );
        return d_coarse_map;
    }
    Teuchos::RCP<const Map_t> getRangeMap() const
    {
        REQUIRE( d_fine_map != Teuchos::null );
        return d_fine_map;
    }

  private:

    std::vector<int> d_steer_vec;

    int d_unks_per_grp;
    int d_fine_groups;
    int d_coarse_groups;

    Teuchos::RCP<const Map_t> d_coarse_map;
    Teuchos::RCP<const Map_t> d_fine_map;
};

} // end namespace profugus
} // end namespace tpetra

#endif // spn_tpetra_Energy_Prolongation_hh

//---------------------------------------------------------------------------//
//                 end of Energy_Prolongation.hh
//---------------------------------------------------------------------------//
