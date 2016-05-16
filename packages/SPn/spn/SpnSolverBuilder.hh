//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   SPn/spn/SpnSolverBuilder.hh
 * \author Steven Hamilton
 * \date   Mon Feb 17 13:10:27 2014
 * \brief  SpnSolverBuilder class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef SPn_spn_SpnSolverBuilder_hh
#define SPn_spn_SpnSolverBuilder_hh

#include "Teuchos_RCP.hpp"
#include "harness/DBC.hh"
#include "solvers/LinAlgTypedefs.hh"
#include "Solver_Base.hh"

namespace profugus
{

//===========================================================================//
/*!
 * \class SpnSolverBuilder
 * \brief Construct SPN solver of given type.
 */
//===========================================================================//

class SpnSolverBuilder
{
  public:

    typedef Teuchos::ParameterList      ParameterList;
    typedef Teuchos::RCP<ParameterList> RCP_ParameterList;

    static Teuchos::RCP<Solver_Base> build(std::string       problem,
                                           RCP_ParameterList db);
};

} // end namespace profugus

#endif // SPn_spn_SpnSolverBuilder_hh

//---------------------------------------------------------------------------//
//                 end of SpnSolverBuilder.hh
//---------------------------------------------------------------------------//
