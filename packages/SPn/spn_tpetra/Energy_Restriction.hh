//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   spn_tpetra/Energy_Restriction.hh
 * \author Thomas M. Evans, Steven Hamilton
 * \date   Tue Feb 25 12:35:37 2014
 * \brief  Energy_Restriction class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef spn_tpetra_Energy_Restriction_hh
#define spn_tpetra_Energy_Restriction_hh

#include <vector>

#include "Teuchos_RCP.hpp"

#include "solvers/LinAlgTypedefs.hh"

namespace profugus
{
namespace tpetra
{

//===========================================================================//
/*!
 * \class Energy_Restriction
 * \brief Perform energy multigrid transfer functions.
 *
 * This class performs a restriction operation for the
 * SPN equations.  It is assumed that energy is the innermost variable.
 *
 * \sa Energy_Restriction.cc for detailed descriptions.
 */
/*!
 * \example spn_tpetra/test/tstEnergy_Restriction.cc
 *
 * Test of Energy_Restriction.
 */
//===========================================================================//

class Energy_Restriction : public TpetraTypes::OP
{
  public:

    typedef typename TpetraTypes::MAP    Map_t;
    typedef typename TpetraTypes::OP     OP;
    typedef typename TpetraTypes::MV     MV;
    typedef typename TpetraTypes::VECTOR Vector_t;

    Energy_Restriction( Teuchos::RCP<const MV>  fine_vec,
                        Teuchos::RCP<const MV>  coarse_vec,
                        const std::vector<int> &steer_vec );

    void apply( const MV &x, MV &y, Teuchos::ETransp mode=Teuchos::NO_TRANS,
                double alpha=Teuchos::ScalarTraits<double>::one(),
                double beta=Teuchos::ScalarTraits<double>::zero()) const;

    // Required interface
    bool hasTransposeApply(){return false;}
    Teuchos::RCP<const Map_t> getDomainMap() const
    {
        REQUIRE( d_fine_map != Teuchos::null );
        return d_fine_map;
    }
    Teuchos::RCP<const Map_t> getRangeMap() const
    {
        REQUIRE( d_coarse_map != Teuchos::null );
        return d_coarse_map;
    }



  private:

    std::vector<int> d_steer_vec;

    int d_unks_per_grp;
    int d_fine_groups;
    int d_coarse_groups;

    Teuchos::RCP<const Map_t> d_fine_map;
    Teuchos::RCP<const Map_t> d_coarse_map;
};

} // end namespace profugus
} // end namespace tpetra

#endif // spn_tpetra_Energy_Restriction_hh

//---------------------------------------------------------------------------//
//                 end of Energy_Restriction.hh
//---------------------------------------------------------------------------//
