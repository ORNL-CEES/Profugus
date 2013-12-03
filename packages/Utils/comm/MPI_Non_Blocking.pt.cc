//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   comm/MPI_Non_Blocking.pt.cc
 * \author Thomas M. Evans
 * \date   Wed Jan  2 13:30:20 2008
 * \brief  MPI Explicit template instatiations.
 * \note   Copyright (C) 2008 Oak Ridge National Laboratory, UT-Battelle, LLC.
 *
 * Non-Blocking explicit instantiations.
 */
//---------------------------------------------------------------------------//

#include <comm/config.h>
#ifdef COMM_MPI

#include "MPI.t.hh"

namespace profugus
{

//---------------------------------------------------------------------------//
// EXPLICIT INSTANTIATIONS OF NON-BLOCKING SEND/RECEIVE
//---------------------------------------------------------------------------//

template Request send_async(const char *, int, int, int);
template Request send_async(const unsigned char *, int, int, int);
template Request send_async(const short *, int, int, int);
template Request send_async(const unsigned short *, int, int, int);
template Request send_async(const int *, int, int, int);
template Request send_async(const unsigned int *, int, int, int);
template Request send_async(const long *, int, int, int);
template Request send_async(const unsigned long *, int, int, int);
template Request send_async(const float *, int, int, int);
template Request send_async(const double *, int, int, int);
template Request send_async(const long double *, int, int, int);

template Request receive_async(char *, int, int, int);
template Request receive_async(unsigned char *, int, int, int);
template Request receive_async(short *, int, int, int);
template Request receive_async(unsigned short *, int, int, int);
template Request receive_async(int *, int, int, int);
template Request receive_async(unsigned int *, int, int, int);
template Request receive_async(long *, int, int, int);
template Request receive_async(unsigned long *, int, int, int);
template Request receive_async(float *, int, int, int);
template Request receive_async(double *, int, int, int);
template Request receive_async(long double *, int, int, int);

template void send_async(Request &, const char *, int, int, int);
template void send_async(Request &, const unsigned char *, int, int, int);
template void send_async(Request &, const short *, int, int, int);
template void send_async(Request &, const unsigned short *, int, int, int);
template void send_async(Request &, const int *, int, int, int);
template void send_async(Request &, const unsigned int *, int, int, int);
template void send_async(Request &, const long *, int, int, int);
template void send_async(Request &, const unsigned long *, int, int, int);
template void send_async(Request &, const float *, int, int, int);
template void send_async(Request &, const double *, int, int, int);
template void send_async(Request &, const long double *, int, int, int);

template void receive_async(Request &, char *, int, int, int);
template void receive_async(Request &, unsigned char *, int, int, int);
template void receive_async(Request &, short *, int, int, int);
template void receive_async(Request &, unsigned short *, int, int, int);
template void receive_async(Request &, int *, int, int, int);
template void receive_async(Request &, unsigned int *, int, int, int);
template void receive_async(Request &, long *, int, int, int);
template void receive_async(Request &, unsigned long *, int, int, int);
template void receive_async(Request &, float *, int, int, int);
template void receive_async(Request &, double *, int, int, int);
template void receive_async(Request &, long double *, int, int, int);

template void send_async_comm(Request &, const char *, const MPI_Comm&, int, int, int);
template void send_async_comm(Request &, const unsigned char *, const MPI_Comm&, int, int, int);
template void send_async_comm(Request &, const short *, const MPI_Comm&, int, int, int);
template void send_async_comm(Request &, const unsigned short *, const MPI_Comm&, int, int, int);
template void send_async_comm(Request &, const int *, const MPI_Comm&, int, int, int);
template void send_async_comm(Request &, const unsigned int *, const MPI_Comm&, int, int, int);
template void send_async_comm(Request &, const long *, const MPI_Comm&, int, int, int);
template void send_async_comm(Request &, const unsigned long *, const MPI_Comm&, int, int, int);
template void send_async_comm(Request &, const float *, const MPI_Comm&, int, int, int);
template void send_async_comm(Request &, const double *, const MPI_Comm&, int, int, int);
template void send_async_comm(Request &, const long double *, const MPI_Comm&, int, int, int);

template void receive_async_comm(Request &, char *, const MPI_Comm&, int, int, int);
template void receive_async_comm(Request &, unsigned char *, const MPI_Comm&, int, int, int);
template void receive_async_comm(Request &, short *, const MPI_Comm&, int, int, int);
template void receive_async_comm(Request &, unsigned short *, const MPI_Comm&, int, int, int);
template void receive_async_comm(Request &, int *, const MPI_Comm&, int, int, int);
template void receive_async_comm(Request &, unsigned int *, const MPI_Comm&, int, int, int);
template void receive_async_comm(Request &, long *, const MPI_Comm&, int, int, int);
template void receive_async_comm(Request &, unsigned long *, const MPI_Comm&, int, int, int);
template void receive_async_comm(Request &, float *, const MPI_Comm&, int, int, int);
template void receive_async_comm(Request &, double *, const MPI_Comm&, int, int, int);
template void receive_async_comm(Request &, long double *, const MPI_Comm&, int, int, int);

} // end namespace profugus

#endif // COMM_MPI

//---------------------------------------------------------------------------//
//                 end of MPI_Non_Blocking.pt.cc
//---------------------------------------------------------------------------//
