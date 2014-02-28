//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   spn/Energy_Prolongation.hh
 * \author Thomas M. Evans, Steven Hamilton
 * \date   Tue Feb 25 12:35:44 2014
 * \brief  Energy_Prolongation class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef spn_Energy_Prolongation_hh
#define spn_Energy_Prolongation_hh

#include <vector>

#include "Teuchos_RCP.hpp"
#include "Epetra_MultiVector.h"
#include "Epetra_Operator.h"
#include "Epetra_Map.h"

namespace profugus
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
 * \example spn/test/tstEnergy_Prolongation.cc
 *
 * Test of Energy_Prolongation.
 */
//===========================================================================//

class Energy_Prolongation : public Epetra_Operator
{
  public:

    Energy_Prolongation( const Epetra_MultiVector &fine_vec,
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
    const char * Label() const {return "Energy_Prolongation";}
    bool UseTranspose() const {return false;}
    const Epetra_Comm & Comm() const {return d_fine_map.Comm();}
    const Epetra_Map & OperatorDomainMap() const
    {
        const Epetra_BlockMap *coarse_blockmap = &d_fine_map;
        const Epetra_Map *coarse_ptr =
            dynamic_cast<const Epetra_Map *>(coarse_blockmap);
        return *coarse_ptr;
    }
    const Epetra_Map & OperatorRangeMap() const
    {
        const Epetra_BlockMap *fine_blockmap = &d_fine_map;
        const Epetra_Map *fine_ptr =
            dynamic_cast<const Epetra_Map *>(fine_blockmap);
        return *fine_ptr;
    }


  private:

    std::vector<int> d_steer_vec;

    int d_unks_per_grp;
    int d_fine_groups;
    int d_coarse_groups;

    Epetra_BlockMap d_coarse_map;
    Epetra_BlockMap d_fine_map;
};

} // end namespace profugus

#endif // spn_Energy_Prolongation_hh

//---------------------------------------------------------------------------//
//                 end of Energy_Prolongation.hh
//---------------------------------------------------------------------------//
