//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   utils/Local_Smart_Thread.cc
 * \author Stuart R. Slattery
 * \date   Thu Sep 24 10:38:02 2015
 * \brief  Locally-scoped smart thread.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Local_Smart_Thread.hh"

#include "harness/DBC.hh"

namespace profugus
{

//---------------------------------------------------------------------------//
/*!
 * \brief Default constructor.
 */
Local_Smart_Thread::Local_Smart_Thread()
{ /* ... */ }
    
//---------------------------------------------------------------------------//
/* 
 * \brief Thread constructor. This will take ownwership of the input thread.
 */
Local_Smart_Thread::Local_Smart_Thread( std::thread local_thread )
    : d_thread( std::move(local_thread) )
{
    REQUIRE( d_thread.joinable() );
}

//---------------------------------------------------------------------------//
/*
 * \brief Destructor. If we have a valid thread we are managing, join it.
 */
Local_Smart_Thread::~Local_Smart_Thread()
{
    if ( d_thread.joinable() )
    {
	d_thread.join();
    }
}

//---------------------------------------------------------------------------//

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of Local_Smart_Thread.cc
//---------------------------------------------------------------------------//
