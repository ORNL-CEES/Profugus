//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   SPn/spn/Energy_Restriction.pt.cc
 * \author Thomas M. Evans
 * \date   Fri Feb 21 14:41:20 2014
 * \brief  Energy_Restriction explicit instantiation.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Energy_Restriction.t.hh"
#include "solvers/LinAlgTypedefs.hh"

#include "AnasaziMultiVecTraits.hpp"
#include "AnasaziEpetraAdapter.hpp"
#include "AnasaziTpetraAdapter.hpp"

namespace profugus
{

template class Energy_Restriction<EpetraTypes>;
template class Energy_Restriction<TpetraTypes>;

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of Energy_Restriction.pt.cc
//---------------------------------------------------------------------------//
