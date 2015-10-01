//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   comm/Thread_Guard.i.hh
 * \author Stuart R. Slattery
 * \date   Thu Sep 24 10:38:02 2015
 * \brief  Locally-scoped C++11 thread guard.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef comm_Thread_Guard_i_hh
#define comm_Thread_Guard_i_hh

namespace profugus
{

//---------------------------------------------------------------------------//
/*!
 * \brief Get the id of the underlying local thread.
 */
std::thread::id Thread_Guard::get_id() const
{
    return d_thread.get_id();
}

//---------------------------------------------------------------------------//
/* 
 * \brief Create a local thread guard from a local thread. This will move the
 * thread into the thread guard.
 */
Thread_Guard thread_guard( std::thread local_thread )
{
    return Thread_Guard( std::move(local_thread) );
}

//---------------------------------------------------------------------------//

} // end namespace profugus

#endif // comm_Thread_Guard_i_hh

//---------------------------------------------------------------------------//
//              end of comm/Thread_Guard.i.hh
//---------------------------------------------------------------------------//
