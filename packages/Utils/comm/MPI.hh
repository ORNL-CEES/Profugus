//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   comm/MPI.hh
 * \author Thomas M. Evans
 * \date   Wed Jan  2 13:42:57 2008
 * \brief  MPI function declarations (defined in Functions.hh).
 * \note   Copyright (C) 2008 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef comm_MPI_hh
#define comm_MPI_hh

#include <comm/config.h>
#ifdef COMM_MPI

#include <algorithm>
#include <mpi.h>

#include "harness/DBC.hh"
#include "Functions.hh"

namespace profugus
{

//---------------------------------------------------------------------------//
// MPI Communicator
//---------------------------------------------------------------------------//

extern Communicator_t default_communicator;
extern Communicator_t communicator;

//---------------------------------------------------------------------------//
// BROADCASTS
//---------------------------------------------------------------------------//
/*!
 * Broadcast the range [first, last) from proc 0
 * into [result, ...) on all other processors.
 */
template<class ForwardIterator, class OutputIterator>
void broadcast(ForwardIterator first,
               ForwardIterator last,
               OutputIterator  result)
{
    typedef typename std::iterator_traits<ForwardIterator>::value_type
        value_type;
    typedef typename std::iterator_traits<ForwardIterator>::difference_type
        diff_type;

    // Proc 0 does not copy any data into the result iterator.

    if(node() == 0)
    {
        diff_type size = std::distance(first, last);

        value_type *buf = new value_type[size];
        std::copy(first, last, buf);

        for (int i=1; i < nodes(); ++i)
        {
            send(&size, 1, i);
            send(buf, size, i);
        }
        delete [] buf;
    }
    else
    {
        diff_type size;

        receive(&size, 1, 0);
        value_type *buf = new value_type[size];
        receive(buf,size,0);

        std::copy(buf, buf+size, result);

        delete [] buf;
    }
}

} // end namespace profugus

#endif // COMM_MPI

#endif // comm_MPI_hh

//---------------------------------------------------------------------------//
//              end of comm/MPI.hh
//---------------------------------------------------------------------------//
