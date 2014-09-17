//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   spn/Energy_Restriction.hh
 * \author Thomas M. Evans, Steven Hamilton
 * \date   Tue Feb 25 12:35:37 2014
 * \brief  Energy_Restriction class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef spn_Energy_Restriction_hh
#define spn_Energy_Restriction_hh

#include <vector>

#include "Teuchos_RCP.hpp"
#include "AnasaziMultiVecTraits.hpp"
#include "AnasaziEpetraAdapter.hpp"
#include "AnasaziTpetraAdapter.hpp"

#include "solvers/LinAlgTypedefs.hh"
#include "OperatorAdapter.hh"

namespace profugus
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
 * \example spn/test/tstEnergy_Restriction.cc
 *
 * Test of Energy_Restriction.
 */
//===========================================================================//

template <class T>
class Energy_Restriction : public OperatorAdapter<T>
{
  public:
    //@{
    //! Typedefs.
    typedef typename T::OP                        OP;
    typedef typename T::MV                        MV;
    typedef typename T::VECTOR                    VECTOR;
    typedef typename T::MAP                       MAP;
    typedef Anasazi::MultiVecTraits<double,MV>    MVT;
    //@}

    Energy_Restriction( Teuchos::RCP<const MAP>  fine_map,
                        Teuchos::RCP<const MAP>  coarse_map,
                        const std::vector<int>  &steer_vec );

    // Implement both the Epetra and Tpetra-style operator applies
    // Delegate actual computation to private function
    int Apply( const MV &x,
                     MV &y ) const
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

    void ApplyImpl( const MV &x, MV &y) const;

    std::vector<int> d_steer_vec;

    int d_unks_per_grp;
    int d_fine_groups;
    int d_coarse_groups;

    Teuchos::RCP<const MAP> d_fine_map;
    Teuchos::RCP<const MAP> d_coarse_map;
};

} // end namespace profugus

#endif // spn_Energy_Restriction_hh

//---------------------------------------------------------------------------//
//                 end of Energy_Restriction.hh
//---------------------------------------------------------------------------//
