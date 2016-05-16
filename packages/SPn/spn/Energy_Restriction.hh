//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   SPn/spn/Energy_Restriction.hh
 * \author Thomas M. Evans, Steven Hamilton
 * \date   Tue Feb 25 12:35:37 2014
 * \brief  Energy_Restriction class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef SPn_spn_Energy_Restriction_hh
#define SPn_spn_Energy_Restriction_hh

#include <vector>

#include "Teuchos_RCP.hpp"

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
    typedef typename T::MAP                       MAP;
    typedef Anasazi::MultiVecTraits<double,MV>    MVT;
    //@}

    Energy_Restriction( Teuchos::RCP<const MAP>  fine_map,
                        Teuchos::RCP<const MAP>  coarse_map,
                        const std::vector<int>  &steer_vec );

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

#endif // SPn_spn_Energy_Restriction_hh

//---------------------------------------------------------------------------//
//                 end of Energy_Restriction.hh
//---------------------------------------------------------------------------//
