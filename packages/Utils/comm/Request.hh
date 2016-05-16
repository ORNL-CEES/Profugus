//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Utils/comm/Request.hh
 * \author Thomas M. Evans
 * \date   Wed Jan  2 13:53:24 2008
 * \brief  Request class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef Utils_comm_Request_hh
#define Utils_comm_Request_hh

#include <Utils/config.h>

#ifdef COMM_MPI
#include <mpi.h>
#endif

#include "Definitions.hh"

namespace profugus
{

//===========================================================================//
/*!
 * \class Request
 * \brief Non-blocking communication request class.
 *
 * This class provides an encapsulator for the message requests (MPI) which
 * are produced by non blocking calls.  Reference counting is used so that
 * these may be passed by value without accidentally triggering a program
 * stall.
 */
/*!
 * \example comm/test/tstRequest.cc
 */
//===========================================================================//

class Request
{
  private:
    // >>> MEMBER TYPES

    // Reference counting handler.
    struct RequestRef
    {
        // >>> DATA

        int n;
        int assigned;

#ifdef COMM_MPI
        MPI_Status  s;
        MPI_Request r;
#endif

        // >>> FUNCTIONAL INTERFACE

        RequestRef();
        ~RequestRef();

        void wait();
        void free();

        bool complete();

        unsigned count();

        int inuse() const { return assigned; }

        void set()   { assigned = 1; }
        void clear() { assigned = 0; }

      private:
        // >>> DISALLOWED METHODS

        RequestRef(const RequestRef& rep);
        RequestRef& operator=(const RequestRef& rep);
    };

  private:
    // >>> DATA

    //! Request handle.
    RequestRef *p;

  public:
    // Constructors.
    Request();
    Request(const Request& req);
    ~Request();
    Request& operator=(const Request& req);

    // >>> FUNCTIONAL INTERFACE

    //! Equivalence operator.
    bool operator==(const Request& right) { return (p == right.p); }

    //! Inequality operator.
    bool operator!=(const Request& right) { return (p != right.p ); }

    //! Wait for a request to complete.
    void wait() { p->wait(); }

    //! Free a request.
    void free() { p->free(); }

    //! Check to see if a request is complete.
    bool complete() { return p->complete(); }  // Should be const?

    //! Number of items returned on the last complete operation.
    unsigned count() { return p->count(); }

    //! Query to see if the request is active.
    int inuse() const { return p->inuse(); }

  private:
    // >>> IMPLEMENTATION

    // Set the handle.
    void set() { p->set(); }

#ifdef COMM_MPI
    // Get internal request state (for use by friends only).
    MPI_Request &r() { return p->r; }
#endif

    // FRIENDSHIP

    // Specific friend Comm functions that may need to manipulate the
    // RequestRef internals.

    template<class T>
    friend Request send_async(const T *buf, int nels, int dest, int tag);

    template<class T>
    friend Request receive_async(T *buf, int nels, int source, int tag);

    template<class T>
    friend void send_async(Request &r, const T *buf, int nels, int dest,
                           int tag);
    template<class T>
    friend void receive_async(Request &r, T *buf, int nels, int source,
                              int tag);

    template<class T>
    friend void send_async_comm(Request &r, const T *buf,
            const Communicator_t& comm, int nels, int dest, int tag);

    template<class T>
    friend void receive_async_comm(Request &r, T *buf,
            const Communicator_t& comm, int nels, int source, int tag);
};

} // end namespace profugus

#endif // Utils_comm_Request_hh

//---------------------------------------------------------------------------//
//              end of comm/Request.hh
//---------------------------------------------------------------------------//
