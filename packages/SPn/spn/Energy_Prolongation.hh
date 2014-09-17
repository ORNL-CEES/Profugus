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
#include "AnasaziMultiVecTraits.hpp"
#include "AnasaziEpetraAdapter.hpp"

#include "solvers/LinAlgTypedefs.hh"

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

    //@{
    //! Typedefs.
    typedef EpetraTypes                           T;
    typedef typename T::OP                        OP;
    typedef typename T::MV                        MV;
    typedef typename T::VECTOR                    VECTOR;
    typedef typename T::MAP                       MAP;
    typedef Anasazi::MultiVecTraits<double,MV>    MVT;
    //@}

    Energy_Prolongation( Teuchos::RCP<const MAP> fine_vec,
                         Teuchos::RCP<const MAP> coarse_vec,
                         const std::vector<int> &steer_vec );

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
    const Epetra_Comm & Comm() const {return d_fine_map->Comm();}
    const Epetra_Map & OperatorDomainMap() const
    {
        return *d_coarse_map;
    }
    const Epetra_Map & OperatorRangeMap() const
    {
        return *d_fine_map;
    }


  private:

    std::vector<int> d_steer_vec;

    int d_unks_per_grp;
    int d_fine_groups;
    int d_coarse_groups;

    Teuchos::RCP<const MAP> d_coarse_map;
    Teuchos::RCP<const MAP> d_fine_map;
};

} // end namespace profugus

#endif // spn_Energy_Prolongation_hh

//---------------------------------------------------------------------------//
//                 end of Energy_Prolongation.hh
//---------------------------------------------------------------------------//
