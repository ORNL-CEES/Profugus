//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   comm/SpinLock.cc
 * \author Thomas M. Evans
 * \date   Wed Jan  2 16:14:33 2008
 * \brief  SpinLock class member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "harness/DBC.hh"
#include "global.hh"
#include "SpinLock.hh"

namespace profugus
{

//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 *
 * Waits for the preceeding processor to finish before continuing.
 */
SpinLock::SpinLock( int _lock /*=1*/ )
    : d_lock(_lock)
    , d_node(profugus::node())
    , d_nodes(profugus::nodes())
    , d_first(0)
    , d_last(d_nodes - 1)
{
    if (d_lock && d_node)
        receive( &d_trash, 0, d_node - 1, SL_Next );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Constructor with bounds.
 *
 * Sets up a spinlock over the processor range [begin, end).
 */
SpinLock::SpinLock(int begin,
                   int end)
    : d_lock(1)
    , d_node(profugus::node())
    , d_nodes(profugus::nodes())
    , d_first(begin)
    , d_last(end - 1)
{
    Require (d_last - d_first >= 0);
    Require (d_first >= 0 && d_first < d_nodes);
    Require (d_last  >= 0 && d_last  < d_nodes);

    if (d_node > d_first)
        receive( &d_trash, 0, d_node - 1, SL_Next );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Destructor
 *
 * Here we notify the next processor in the chain that he can proceed to
 * execute the block, and we go ahead about our business.
 */
SpinLock::~SpinLock()
{
    if (d_lock && d_node < d_last)
        send( &d_trash, 0, d_node + 1, SL_Next );
}

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of SpinLock.cc
//---------------------------------------------------------------------------//
