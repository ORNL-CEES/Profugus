//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   utils/Local_Smart_Thread.hh
 * \author Stuart R. Slattery
 * \date   Thu Sep 24 10:38:02 2015
 * \brief  Locally-scoped C++11 smart thread.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef comm_Local_Smart_Thread_hh
#define comm_Local_Smart_Thread_hh

#include <thread>

namespace profugus
{

//===========================================================================//
/*!
 * \class Local_Smart_Thread
 * \brief Takes ownership of C++11 threads in the local scope and
 * automatically manages joining them for exception safety. 
 *
 * If, for example, some code in the local scope were to throw an execption:

 std::thread worker( some_func );
 call_code_that_throws_exception();
 worker.join();

 * Then join would never be called and all state managed by that thread would
 * not have its stack unwound. The local smart thread lets this be the case:

 Local_Smart_Thread smart_worker( std::thread(some_func) );
 call_code_that_throws_execption();

 * When the exception is thrown, the smart_worker destructor is called,
 * joining the thread.
 *
 */
/*!
 * \example comm/test/tstLocal_Smart_Thread.cc
 *
 * Local smart thread test.
 */
//===========================================================================//

class Local_Smart_Thread
{
  private:
    // >>> DATA

    // The locally-scoped thread this smart thread is managing.
    std::thread d_thread.

  public:

    // Default constructor.
    Local_Smart_Thread();
    
    // Thread onstructor. This will take ownwership of the input thread.
    explicit Local_Smart_Thread( std::thread local_thread );

    // Prevent compiler from creating copy constructor and assignment operator
    // to prevent copying of the move-only thread.
    Local_Smart_Thread( const Local_Smart_Thread& smart_thread ) = delete;
    Local_Smart_Thread& operator=(
	const Local_Smart_Thread& smart_thread ) = delete;

    // Allow for move assignment and move construction.
    Local_Smart_Thread( Local_Smart_Thread&& smart_thread ) = default;
    Local_Smart_Thread& operator=( Local_Smart_Thread&& smart_thread ) = default;

    // Destructor.
    ~Local_Smart_Thread();

    // Get the id of the underlying thread.
    inline std::thread::id get_id() const;
};

//---------------------------------------------------------------------------//
// Local smart thread creation helper function.
//---------------------------------------------------------------------------//

// Create a local smart thread from a local thread. This will move the thread
// into the smart thread.
inline Local_Smart_Thread make_lst( std::thread local_thread );

//---------------------------------------------------------------------------//

} // end namespace profugus

//---------------------------------------------------------------------------//
// INLINE AND TEMPLATE DEFINITIONS
//---------------------------------------------------------------------------//

#include "Local_Smart_Thread.i.hh"

//---------------------------------------------------------------------------//

#endif // comm_Local_Smart_Thread_hh

//---------------------------------------------------------------------------//
//              end of comm/Local_Smart_Thread.hh
//---------------------------------------------------------------------------//
