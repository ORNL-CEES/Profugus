//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   solvers/PreconditionerBuilder.hh
 * \author Thomas M. Evans, Steven Hamilton
 * \date   Tue Feb 25 09:28:52 2014
 * \brief  PreconditionerBuilder class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef solvers_PreconditionerBuilder_hh
#define solvers_PreconditionerBuilder_hh

#include "harness/DBC.hh"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "LinAlgTypedefs.hh"

namespace profugus
{

//===========================================================================//
/*!
 * \class PreconditionerBuilder
 * \brief
 */
/*!
 * \example solvers/test/tstPreconditionerBuilder.cc
 *
 * Test of PreconditionerBuilder.
 */
//===========================================================================//

template <LinAlgType T>
class PreconditionerBuilder
{
  public:

    typedef typename LinAlgTypedefs<T>::OP       OP;
    typedef Teuchos::RCP<Teuchos::ParameterList> RCP_ParameterList;

    static Teuchos::RCP<OP> build_preconditioner(
        Teuchos::RCP<OP> op, RCP_ParameterList db );
};

} // end namespace profugus

#endif // solvers_PreconditionerBuilder_hh

//---------------------------------------------------------------------------//
//                 end of PreconditionerBuilder.hh
//---------------------------------------------------------------------------//
