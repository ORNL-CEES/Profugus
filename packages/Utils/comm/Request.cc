//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   comm/Request.cc
 * \author Thomas M. Evans
 * \date   Wed Jan  2 13:53:24 2008
 * \brief  Request class definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <iostream>

#include "harness/DBC.hh"
#include "Request.hh"

namespace profugus
{

//---------------------------------------------------------------------------//
// REQUEST MEMBER FUNCTIONS
//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 *
 * Register a new non blocking message request.
 */
Request::Request()
    : p(new RequestRef)
{
    ++p->n;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Copy constructor.
 *
 * Attach to an existing message request.
 */
Request::Request(const Request& req)
{
    if (req.inuse())
        p = req.p;
    else
        p = new RequestRef;
    ++p->n;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Destructor.
 *
 * Only delete if this is the last request holding the handle.
 */
Request::~Request()
{
    --p->n;
    if (p->n <= 0)
        delete p;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Assignment.
 *
 * Detach from our prior message request, waiting on it if necessary.  Then
 * attach to the new one.
 */
Request& Request::operator=(const Request& req)
{
    --p->n;
    if (p->n <= 0)
        delete p;

    if (req.inuse())
        p = req.p;
    else
        p = new RequestRef;

    ++p->n;

    return *this;
}

//---------------------------------------------------------------------------//
// REQUESTREF MEMBER FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 *
 * Register a new non blocking message request.
 */
Request::RequestRef::RequestRef()
    : n(0)
    , assigned(0)
{
    // empty
}

//---------------------------------------------------------------------------//
/*!
 * \brief Destructor.
 *
 * It is important that all existing requests are cleared before the
 * destructor is called.  We used to have a wait() in here; however, this
 * causes exception safety problems.  In any case, it is probably a bad idea
 * to clean up communication by going out of scope.
 */
Request::RequestRef::~RequestRef()
{
}

//---------------------------------------------------------------------------//
/*!
 * \brief Wait for an asynchronous message to complete.
 */
void Request::RequestRef::wait()
{
    if (assigned)
    {
#ifdef COMM_MPI
        MPI_Wait(&r, &s);
#endif
    }
    clear();
}

//---------------------------------------------------------------------------//
/*!
 * \brief Free request handle for a posted asynchronous receive.
 *
 * Note: once freed the handle must be reactivated to test for completeness or
 * to wait on it
 */
void Request::RequestRef::free()
{
#ifdef COMM_MPI
    if (assigned)
    {
        MPI_Cancel( &r );
        MPI_Request_free( &r );
    }
#endif
    clear();
}

//---------------------------------------------------------------------------//
/*!
 * \brief Tests for the completion of a non blocking operation.
 */
bool Request::RequestRef::complete()
{
#ifdef COMM_MPI
    int flag       = 0;
    bool indicator = false;
    if (assigned)
        MPI_Test( &r, &flag, &s );
    if (flag != 0)
    {
        clear();
        Check ( r == MPI_REQUEST_NULL);
        indicator = true;
    }
    return indicator;
#endif
#ifdef COMM_SCALAR
    throw profugus::assertion(
        "Send to self machinery has not been implemented in scalar mode.");
#endif
}

//---------------------------------------------------------------------------//
/*!
 * \brief Return the number of items returned on the last complete operation.
 */
unsigned Request::RequestRef::count()
{
#ifdef COMM_MPI
    int count;
    MPI_Get_count( &s, MPI_CHAR, &count );
    return count;
#else
    return 0;
#endif
}

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of Request.cc
//---------------------------------------------------------------------------//
