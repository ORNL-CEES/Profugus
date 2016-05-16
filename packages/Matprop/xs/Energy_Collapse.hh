//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Matprop/xs/Energy_Collapse.hh
 * \author Thomas M. Evans, Steven Hamilton
 * \date   Tue Feb 25 13:39:38 2014
 * \brief  Energy_Collapse class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef Matprop_xs_Energy_Collapse_hh
#define Matprop_xs_Energy_Collapse_hh

#include "Teuchos_RCP.hpp"

#include "utils/Definitions.hh"
#include "Mat_DB.hh"

namespace profugus
{

//===========================================================================//
/*!
 * \class Energy_Collapse
 * \brief Collapse Mat_DB to coarser energy group structure
 *
 * Create a Mat_DB from an existing Mat_DB with a collapsed energy group
 * structure.
 *
 * \sa Energy_Collapse.cc for detailed descriptions.
 */
/*!
 * \example xs/test/tstEnergy_Collapse.cc
 *
 * Test of Energy_Collapse.
 */
//===========================================================================//

class Energy_Collapse
{
  public:
    //@{
    //! Typedefs.
    typedef Mat_DB                 Mat_DB_t;
    typedef Teuchos::RCP<Mat_DB_t> RCP_Mat_DB;
    typedef def::Vec_Int           Vec_Int;
    typedef def::Vec_Dbl           Vec_Dbl;
    //@}

  private:
    // Prevent construction
    Energy_Collapse() { /* * */ }

  public:
    static RCP_Mat_DB collapse_all_mats(RCP_Mat_DB     fine_mat,
                                        const Vec_Int &collapse_vec,
                                        const Vec_Dbl &weights);
};

} // end namespace profugus

#endif // Matprop_xs_Energy_Collapse_hh

//---------------------------------------------------------------------------//
//                 end of Energy_Collapse.hh
//---------------------------------------------------------------------------//
