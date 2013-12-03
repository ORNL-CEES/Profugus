//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   comm/MPI_Blocking.pt.cc
 * \author Thomas M. Evans
 * \date   Wed Jan  2 13:30:20 2008
 * \brief  MPI Explicit template instatiations.
 * \note   Copyright (C) 2008 Oak Ridge National Laboratory, UT-Battelle, LLC.
 *
 * Blocking explicit instantiations.
 */
//---------------------------------------------------------------------------//

#include <Utils/config.h>
#ifdef COMM_MPI

#include "MPI.t.hh"

namespace profugus
{

//---------------------------------------------------------------------------//
// EXPLICIT INSTANTIATIONS OF BLOCKING SEND/RECEIVE AND BROADCAST
//---------------------------------------------------------------------------//

template int send(const char *, int, int, int);
template int send(const unsigned char *, int, int, int);
template int send(const short *, int, int, int);
template int send(const unsigned short *, int, int, int);
template int send(const int *, int, int, int);
template int send(const unsigned int *, int, int, int);
template int send(const long *, int, int, int);
template int send(const unsigned long *, int, int, int);
template int send(const float *, int, int, int);
template int send(const double *, int, int, int);
template int send(const long double *, int, int, int);

template int send_comm(const char *, const MPI_Comm&, int, int, int);
template int send_comm(const unsigned char *, const MPI_Comm&, int, int, int);
template int send_comm(const short *, const MPI_Comm&, int, int, int);
template int send_comm(const unsigned short *, const MPI_Comm&, int, int, int);
template int send_comm(const int *, const MPI_Comm&, int, int, int);
template int send_comm(const unsigned int *, const MPI_Comm&, int, int, int);
template int send_comm(const long *, const MPI_Comm&, int, int, int);
template int send_comm(const unsigned long *, const MPI_Comm&, int, int, int);
template int send_comm(const float *, const MPI_Comm&, int, int, int);
template int send_comm(const double *, const MPI_Comm&, int, int, int);
template int send_comm(const long double *, const MPI_Comm&, int, int, int);

template int receive(char *, int, int, int);
template int receive(unsigned char *, int, int, int);
template int receive(short *, int, int, int);
template int receive(unsigned short *, int, int, int);
template int receive(int *, int, int, int);
template int receive(unsigned int *, int, int, int);
template int receive(long *, int, int, int);
template int receive(unsigned long *, int, int, int);
template int receive(float *, int, int, int);
template int receive(double *, int, int, int);
template int receive(long double *, int, int, int);

template int receive_comm(char *, const MPI_Comm&, int, int, int);
template int receive_comm(unsigned char *, const MPI_Comm&, int, int, int);
template int receive_comm(short *, const MPI_Comm&, int, int, int);
template int receive_comm(unsigned short *, const MPI_Comm&, int, int, int);
template int receive_comm(int *, const MPI_Comm&, int, int, int);
template int receive_comm(unsigned int *, const MPI_Comm&, int, int, int);
template int receive_comm(long *, const MPI_Comm&, int, int, int);
template int receive_comm(unsigned long *, const MPI_Comm&, int, int, int);
template int receive_comm(float *, const MPI_Comm&, int, int, int);
template int receive_comm(double *, const MPI_Comm&, int, int, int);
template int receive_comm(long double *, const MPI_Comm&, int, int, int);

template int broadcast(char *, int, int);
template int broadcast(unsigned char *, int, int);
template int broadcast(short *, int, int);
template int broadcast(unsigned short *, int, int);
template int broadcast(int *, int, int);
template int broadcast(unsigned int *, int, int);
template int broadcast(long *, int, int);
template int broadcast(unsigned long *, int, int);
template int broadcast(float *, int, int);
template int broadcast(double *, int, int);
template int broadcast(long double *, int, int);

template int broadcast(char *, const MPI_Comm&, int, int);
template int broadcast(unsigned char *, const MPI_Comm&, int, int);
template int broadcast(short *, const MPI_Comm&, int, int);
template int broadcast(unsigned short *, const MPI_Comm&, int, int);
template int broadcast(int *, const MPI_Comm&, int, int);
template int broadcast(unsigned int *, const MPI_Comm&, int, int);
template int broadcast(long *, const MPI_Comm&, int, int);
template int broadcast(unsigned long *, const MPI_Comm&, int, int);
template int broadcast(float *, const MPI_Comm&, int, int);
template int broadcast(double *, const MPI_Comm&, int, int);
template int broadcast(long double *, const MPI_Comm&, int, int);

} // end namespace profugus

#endif // COMM_MPI

//---------------------------------------------------------------------------//
//                 end of MPI_Blocking.pt.cc
//---------------------------------------------------------------------------//
