//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   utils/Thread_Guard.cc
 * \author Stuart R. Slattery
 * \date   Thu Sep 24 10:38:02 2015
 * \brief  Locally-scoped thread guard.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Thread_Guard.hh"

#include "harness/DBC.hh"

namespace profugus
{

//---------------------------------------------------------------------------//
/*!
 * \brief Default constructor.
 */
Thread_Guard::Thread_Guard()
{ /* ... */ }
    
//---------------------------------------------------------------------------//
/* 
 * \brief Thread constructor. This will take ownwership of the input thread.
 */
Thread_Guard::Thread_Guard( std::thread local_thread )
    : d_thread( std::move(local_thread) )
{
    REQUIRE( d_thread.joinable() );
}

//---------------------------------------------------------------------------//
/*
 * \brief Destructor. If we have a valid thread we are managing, join it.
 */
Thread_Guard::~Thread_Guard()
{
    if ( d_thread.joinable() )
    {
	d_thread.join();
    }
}

//---------------------------------------------------------------------------//

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of Thread_Guard.cc
//---------------------------------------------------------------------------//
