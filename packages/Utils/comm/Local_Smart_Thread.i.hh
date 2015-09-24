//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   comm/Local_Smart_Thread.i.hh
 * \author Stuart R. Slattery
 * \date   Thu Sep 24 10:38:02 2015
 * \brief  Locally-scoped C++11 smart thread.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef comm_Local_Smart_Thread_i_hh
#define comm_Local_Smart_Thread_i_hh

namespace profugus
{

//---------------------------------------------------------------------------//
/*!
 * \brief Get the id of the underlying local thread.
 */
std::thread::id Local_Smart_Thread::get_id() const
{
    return d_thread.get_id();
}

//---------------------------------------------------------------------------//
/* 
 * \brief Create a local smart thread from a local thread. This will move the
 * thread into the smart thread.
 */
Local_Smart_Thread make_lst( std::thread local_thread )
{
    return LocalSmartThread( std::move(local_thread) );
}

//---------------------------------------------------------------------------//

} // end namespace profugus

#endif // comm_Local_Smart_Thread_i_hh

//---------------------------------------------------------------------------//
//              end of comm/Local_Smart_Thread.i.hh
//---------------------------------------------------------------------------//
