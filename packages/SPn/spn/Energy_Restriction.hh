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

#include "solvers/LinAlgTypedefs.hh"

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

class Energy_Restriction : public Epetra_Operator
{
  public:
    //@{
    //! Typedefs.
    typedef EpetraTypes                           T;
    typedef typename T::OP                        OP;
    typedef typename T::MV                        MV;
    typedef typename T::VECTOR                    VECTOR;
    typedef typename T::MAP                       MAP;
    typedef Anasazi::MultiVecTraits<double,MV>    MVT;
    //@}

    Energy_Restriction( Teuchos::RCP<const MAP>  fine_map,
                        Teuchos::RCP<const MAP>  coarse_map,
                        const std::vector<int>  &steer_vec );

    int Apply( const Epetra_MultiVector &x,
                     Epetra_MultiVector &y ) const;

    // Required interface
    int SetUseTranspose(bool use){return -1;}
    int ApplyInverse(const Epetra_MultiVector &x,
                           Epetra_MultiVector &y ) const {return -1;}
    bool HasNormInf()const {return false;}
    double NormInf() const {return 0.0;}
    const char * Label() const {return "Energy_Restriction";}
    bool UseTranspose() const {return false;}
    const Epetra_Comm & Comm() const {return d_fine_map->Comm();}
    const Epetra_Map & OperatorDomainMap() const
    {
        return *d_fine_map;
    }
    const Epetra_Map & OperatorRangeMap() const
    {
        return *d_coarse_map;
    }

  private:

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
