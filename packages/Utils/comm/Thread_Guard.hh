//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   utils/Thread_Guard.hh
 * \author Stuart R. Slattery
 * \date   Thu Sep 24 10:38:02 2015
 * \brief  Locally-scoped C++11 thread guard.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef comm_Thread_Guard_hh
#define comm_Thread_Guard_hh

#include <thread>

namespace profugus
{

//===========================================================================//
/*!
 * \class Thread_Guard
 * \brief Takes ownership of C++11 threads in the local scope and
 * automatically manages joining them for exception safety. 
 *
 * If, for example, some code in the local scope were to throw an execption:

 std::thread worker( some_func );
 call_code_that_throws_exception();
 worker.join();

 * Then join would never be called and all state managed by that thread would
 * not have its stack unwound. The local thread guard lets this be the case:

 Thread_Guard smart_worker( std::thread(some_func) );
 call_code_that_throws_execption();

 * When the exception is thrown, the smart_worker destructor is called,
 * joining the thread.
 *
 */
/*!
 * \example comm/test/tstThread_Guard.cc
 *
 * Local thread guard test.
 */
//===========================================================================//

class Thread_Guard
{
  private:
    // >>> DATA

    // The locally-scoped thread this thread guard is managing.
    std::thread d_thread;

  public:

    // Default constructor.
    Thread_Guard();
    
    // Thread constructor. This will take ownwership of the input thread.
    explicit Thread_Guard( std::thread local_thread );

    // Prevent compiler from creating copy constructor and assignment operator
    // to prevent copying of the move-only thread.
    Thread_Guard( const Thread_Guard& smart_thread ) = delete;
    Thread_Guard& operator=(
	const Thread_Guard& smart_thread ) = delete;

    // Allow for move assignment and move construction.
    Thread_Guard( Thread_Guard&& smart_thread ) = default;
    Thread_Guard& operator=( Thread_Guard&& smart_thread ) = default;

    // Destructor.
    ~Thread_Guard();

    // Get the id of the underlying thread.
    inline std::thread::id get_id() const;
};

//---------------------------------------------------------------------------//
// Local thread guard creation helper function.
//---------------------------------------------------------------------------//

// Create a local thread guard from a local thread. This will move the thread
// into the thread guard.
inline Thread_Guard thread_guard( std::thread local_thread );

//---------------------------------------------------------------------------//

} // end namespace profugus

//---------------------------------------------------------------------------//
// INLINE AND TEMPLATE DEFINITIONS
//---------------------------------------------------------------------------//

#include "Thread_Guard.i.hh"

//---------------------------------------------------------------------------//

#endif // comm_Thread_Guard_hh

//---------------------------------------------------------------------------//
//              end of comm/Thread_Guard.hh
//---------------------------------------------------------------------------//
