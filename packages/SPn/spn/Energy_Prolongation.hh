//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   SPn/spn/Energy_Prolongation.hh
 * \author Thomas M. Evans, Steven Hamilton
 * \date   Tue Feb 25 12:35:44 2014
 * \brief  Energy_Prolongation class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef SPn_spn_Energy_Prolongation_hh
#define SPn_spn_Energy_Prolongation_hh

#include <vector>

#include "Teuchos_RCP.hpp"

#include "solvers/LinAlgTypedefs.hh"
#include "OperatorAdapter.hh"

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

template <class T>
class Energy_Prolongation : public OperatorAdapter<T>
{
  public:

    //@{
    //! Typedefs.
    typedef typename T::OP                        OP;
    typedef typename T::MV                        MV;
    typedef typename T::MAP                       MAP;
    typedef Anasazi::MultiVecTraits<double,MV>    MVT;
    //@}

    Energy_Prolongation( Teuchos::RCP<const MAP> coarse_map,
                         Teuchos::RCP<const MAP> fine_map,
                         const std::vector<int> &steer_vec );

  private:

    void ApplyImpl(const MV &x, MV &y) const;

    std::vector<int> d_steer_vec;

    int d_unks_per_grp;
    int d_fine_groups;
    int d_coarse_groups;

    Teuchos::RCP<const MAP> d_coarse_map;
    Teuchos::RCP<const MAP> d_fine_map;
};

} // end namespace profugus

#endif // SPn_spn_Energy_Prolongation_hh

//---------------------------------------------------------------------------//
//                 end of Energy_Prolongation.hh
//---------------------------------------------------------------------------//
