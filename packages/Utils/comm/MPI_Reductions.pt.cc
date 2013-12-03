//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   comm/MPI_Reductions.pt.cc
 * \author Thomas M. Evans
 * \date   Wed Jan  2 13:30:20 2008
 * \brief  MPI Explicit template instatiations.
 * \note   Copyright (C) 2008 Oak Ridge National Laboratory, UT-Battelle, LLC.
 *
 * Global reduction explicit instantiations.
 */
//---------------------------------------------------------------------------//

#include <Utils/config.h>
#ifdef COMM_MPI

#include "MPI.t.hh"

namespace profugus
{

//---------------------------------------------------------------------------//
// EXPLICIT INSTANTIATIONS OF GLOBAL REDUCTIONS
//---------------------------------------------------------------------------//

template void global_sum(short &);
template void global_sum(unsigned short &);
template void global_sum(int &);
template void global_sum(unsigned int &);
template void global_sum(long &);
template void global_sum(unsigned long &);
template void global_sum(float &);
template void global_sum(double &);
template void global_sum(long double &);

template void global_prod(short &);
template void global_prod(unsigned short &);
template void global_prod(int &);
template void global_prod(unsigned int &);
template void global_prod(long &);
template void global_prod(unsigned long &);
template void global_prod(float &);
template void global_prod(double &);
template void global_prod(long double &);

template void global_max(short &);
template void global_max(unsigned short &);
template void global_max(int &);
template void global_max(unsigned int &);
template void global_max(long &);
template void global_max(unsigned long &);
template void global_max(float &);
template void global_max(double &);
template void global_max(long double &);

template void global_max(short &, const MPI_Comm&);
template void global_max(unsigned short &, const MPI_Comm&);
template void global_max(int &, const MPI_Comm&);
template void global_max(unsigned int &, const MPI_Comm&);
template void global_max(long &, const MPI_Comm&);
template void global_max(unsigned long &, const MPI_Comm&);
template void global_max(float &, const MPI_Comm&);
template void global_max(double &, const MPI_Comm&);
template void global_max(long double &, const MPI_Comm&);

template void global_min(short &);
template void global_min(unsigned short &);
template void global_min(int &);
template void global_min(unsigned int &);
template void global_min(long &);
template void global_min(unsigned long &);
template void global_min(float &);
template void global_min(double &);
template void global_min(long double &);

template void global_min(short &, const MPI_Comm&);
template void global_min(unsigned short &, const MPI_Comm&);
template void global_min(int &, const MPI_Comm&);
template void global_min(unsigned int &, const MPI_Comm&);
template void global_min(long &, const MPI_Comm&);
template void global_min(unsigned long &, const MPI_Comm&);
template void global_min(float &, const MPI_Comm&);
template void global_min(double &, const MPI_Comm&);
template void global_min(long double &, const MPI_Comm&);

template void global_sum(short *, int);
template void global_sum(unsigned short *, int);
template void global_sum(int *, int);
template void global_sum(unsigned int *, int);
template void global_sum(long *, int);
template void global_sum(unsigned long *, int);
template void global_sum(float *, int);
template void global_sum(double *, int);
template void global_sum(long double *, int);

template void global_prod(short *, int);
template void global_prod(unsigned short *, int);
template void global_prod(int *, int);
template void global_prod(unsigned int *, int);
template void global_prod(long *, int);
template void global_prod(unsigned long *, int);
template void global_prod(float *, int);
template void global_prod(double *, int);
template void global_prod(long double *, int);

template void global_max(short *, int);
template void global_max(unsigned short *, int);
template void global_max(int *, int);
template void global_max(unsigned int *, int);
template void global_max(long *, int);
template void global_max(unsigned long *, int);
template void global_max(float *, int);
template void global_max(double *, int);
template void global_max(long double *, int);

template void global_min(short *, int);
template void global_min(unsigned short *, int);
template void global_min(int *, int);
template void global_min(unsigned int *, int);
template void global_min(long *, int);
template void global_min(unsigned long *, int);
template void global_min(float *, int);
template void global_min(double *, int);
template void global_min(long double *, int);

//---------------------------------------------------------------------------//
// EXPLICIT INSTANTIATIONS OF REDUCTIONS
//---------------------------------------------------------------------------//

template void sum(short &, int);
template void sum(unsigned short &, int);
template void sum(int &, int);
template void sum(unsigned int &, int);
template void sum(long &, int);
template void sum(unsigned long &, int);
template void sum(float &, int);
template void sum(double &, int);
template void sum(long double &, int);

template void prod(short &, int);
template void prod(unsigned short &, int);
template void prod(int &, int);
template void prod(unsigned int &, int);
template void prod(long &, int);
template void prod(unsigned long &, int);
template void prod(float &, int);
template void prod(double &, int);
template void prod(long double &, int);

template void min(short &, int);
template void min(unsigned short &, int);
template void min(int &, int);
template void min(unsigned int &, int);
template void min(long &, int);
template void min(unsigned long &, int);
template void min(float &, int);
template void min(double &, int);
template void min(long double &, int);

template void max(short &, int);
template void max(unsigned short &, int);
template void max(int &, int);
template void max(unsigned int &, int);
template void max(long &, int);
template void max(unsigned long &, int);
template void max(float &, int);
template void max(double &, int);
template void max(long double &, int);

template void sum(short *, int, int);
template void sum(unsigned short *, int, int);
template void sum(int *, int, int);
template void sum(unsigned int *, int, int);
template void sum(long *, int, int);
template void sum(unsigned long *, int, int);
template void sum(float *, int, int);
template void sum(double *, int, int);
template void sum(long double *, int, int);

template void prod(short *, int, int);
template void prod(unsigned short *, int, int);
template void prod(int *, int, int);
template void prod(unsigned int *, int, int);
template void prod(long *, int, int);
template void prod(unsigned long *, int, int);
template void prod(float *, int, int);
template void prod(double *, int, int);
template void prod(long double *, int, int);

template void min(short *, int, int);
template void min(unsigned short *, int, int);
template void min(int *, int, int);
template void min(unsigned int *, int, int);
template void min(long *, int, int);
template void min(unsigned long *, int, int);
template void min(float *, int, int);
template void min(double *, int, int);
template void min(long double *, int, int);

template void max(short *, int, int);
template void max(unsigned short *, int, int);
template void max(int *, int, int);
template void max(unsigned int *, int, int);
template void max(long *, int, int);
template void max(unsigned long *, int, int);
template void max(float *, int, int);
template void max(double *, int, int);
template void max(long double *, int, int);

//---------------------------------------------------------------------------//
// EXPLICIT INSTANTIATIONS OF REDUCTION-SCATTERS
//---------------------------------------------------------------------------//

template void sum_scatter(const short *, short *, const int *);
template void sum_scatter(const unsigned short *, unsigned short *, const int *);
template void sum_scatter(const int *, int *, const int *);
template void sum_scatter(const unsigned int *, unsigned int *, const int *);
template void sum_scatter(const long *, long *, const int *);
template void sum_scatter(const unsigned long *, unsigned long *, const int *);
template void sum_scatter(const float *, float *, const int *);
template void sum_scatter(const double *, double *, const int *);
template void sum_scatter(const long double *, long double *, const int *);

template void prod_scatter(const short *, short *, const int *);
template void prod_scatter(const unsigned short *, unsigned short *, const int *);
template void prod_scatter(const int *, int *, const int *);
template void prod_scatter(const unsigned int *, unsigned int *, const int *);
template void prod_scatter(const long *, long *, const int *);
template void prod_scatter(const unsigned long *, unsigned long *, const int *);
template void prod_scatter(const float *, float *, const int *);
template void prod_scatter(const double *, double *, const int *);
template void prod_scatter(const long double *, long double *, const int *);

//---------------------------------------------------------------------------//
// EXPLICIT INSTANTIATIONS OF ALL-TO-ALLS
//---------------------------------------------------------------------------//

template void all_to_all(const short *, short *, int);
template void all_to_all(const unsigned short *, unsigned short *, int);
template void all_to_all(const int *, int *, int);
template void all_to_all(const unsigned int *, unsigned int *, int);
template void all_to_all(const long *, long *, int);
template void all_to_all(const unsigned long *, unsigned long *, int);
template void all_to_all(const float *, float *, int);
template void all_to_all(const double *, double *, int);
template void all_to_all(const long double *, long double *, int);


template void all_to_all(const short *, const int *, const int *,
                               short *, const int *, const int *);
template void all_to_all(const unsigned short *, const int *, const int *,
                               unsigned short *, const int *, const int *);
template void all_to_all(const int *, const int *, const int *,
                               int *, const int *, const int *);
template void all_to_all(const unsigned int *, const int *, const int *,
                               unsigned int *, const int *, const int *);
template void all_to_all(const long *, const int *, const int *,
                               long *, const int *, const int *);
template void all_to_all(const unsigned long *, const int *, const int *,
                               unsigned long *, const int *, const int *);
template void all_to_all(const float *, const int *, const int *,
                               float *, const int *, const int *);
template void all_to_all(const double *, const int *, const int *,
                               double *, const int *, const int *);
template void all_to_all(const long double *, const int *, const int *,
                               long double *, const int *, const int *);

} // end namespace profugus

#endif // COMM_MPI

//---------------------------------------------------------------------------//
//                 end of MPI_Reductions.pt.cc
//---------------------------------------------------------------------------//
