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
#include "Epetra_MultiVector.h"
#include "Epetra_Operator.h"
#include "Epetra_Map.h"

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
    Energy_Restriction( const Epetra_MultiVector &fine_vec,
                        const Epetra_MultiVector &coarse_vec,
                        const std::vector<int>   &steer_vec );

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
    const Epetra_Comm & Comm() const {return d_fine_map.Comm();}
    const Epetra_Map & OperatorDomainMap() const
    {
        const Epetra_BlockMap *fine_blockmap = &d_fine_map;
        const Epetra_Map *fine_ptr =
            dynamic_cast<const Epetra_Map *>(fine_blockmap);
        return *fine_ptr;
    }
    const Epetra_Map & OperatorRangeMap() const
    {
        const Epetra_BlockMap *coarse_blockmap = &d_fine_map;
        const Epetra_Map *coarse_ptr =
            dynamic_cast<const Epetra_Map *>(coarse_blockmap);
        return *coarse_ptr;
    }

  private:

    std::vector<int> d_steer_vec;

    int d_unks_per_grp;
    int d_fine_groups;
    int d_coarse_groups;

    Epetra_BlockMap d_fine_map;
    Epetra_BlockMap d_coarse_map;
};

} // end namespace profugus

#endif // spn_Energy_Restriction_hh

//---------------------------------------------------------------------------//
//                 end of Energy_Restriction.hh
//---------------------------------------------------------------------------//
