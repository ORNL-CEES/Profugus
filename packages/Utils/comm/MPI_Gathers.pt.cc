//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   comm/MPI_Gathers.pt.cc
 * \author Thomas M. Evans
 * \date   Wed Jan  2 13:30:20 2008
 * \brief  MPI Explicit template instatiations.
 * \note   Copyright (C) 2008 Oak Ridge National Laboratory, UT-Battelle, LLC.
 *
 * Global reduction explicit instantiations.
 */
//---------------------------------------------------------------------------//

#include <comm/config.h>
#ifdef COMM_MPI

#include "MPI.t.hh"

namespace nemesis
{

//---------------------------------------------------------------------------//
// EXPLICIT INSTANTIATIONS OF GATHERS
//---------------------------------------------------------------------------//

template void all_gather(const short *, short *, int);
template void all_gather(const unsigned short *, unsigned short *, int);
template void all_gather(const int *, int *, int);
template void all_gather(const unsigned int *, unsigned int *, int);
template void all_gather(const long *, long *, int);
template void all_gather(const unsigned long *, unsigned long *, int);
template void all_gather(const float *, float *, int);
template void all_gather(const double *, double *, int);
template void all_gather(const long double *, long double *, int);

template void gather(const short *, short *, int, int);
template void gather(const unsigned short *, unsigned short *, int, int);
template void gather(const int *, int *, int, int);
template void gather(const unsigned int *, unsigned int *, int, int);
template void gather(const long *, long *, int, int);
template void gather(const unsigned long *, unsigned long *, int, int);
template void gather(const float *, float *, int, int);
template void gather(const double *, double *, int, int);
template void gather(const long double *, long double *, int, int);

} // end namespace nemesis

#endif // COMM_MPI

//---------------------------------------------------------------------------//
//                 end of MPI_Gathers.pt.cc
//---------------------------------------------------------------------------//
