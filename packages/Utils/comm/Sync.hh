//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   comm/Sync.hh
 * \author Thomas M. Evans
 * \date   Wed Jan  2 16:08:10 2008
 * \brief  Sync class definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef comm_Sync_hh
#define comm_Sync_hh

#include "global.hh"
#include "Sync.hh"

namespace profugus
{

//===========================================================================//
/*!
 * \class HSync
 * \brief Head synchronizing class.
 *
 * Synchronizes processes at the head of a block by doing a global sync in the
 * ctor.
 */
//===========================================================================//

class HSync
{
    HSync(const HSync &);
    HSync& operator=(const HSync &);

  public:
    HSync(int s = 1) { if (s) global_barrier(); }
};

//===========================================================================//
/*!
 * \class TSync
 * \brief Tail synchronizing class.
 *
 * Synchronizes processes at the tail of a block by doing a global sync in the
 * dtor.
 */
//===========================================================================//

class TSync
{
    TSync(const TSync &);
    TSync& operator=(const TSync &);

    int sync;

  public:
    TSync(int s = 1) : sync(s) {}
    virtual ~TSync() { if (sync) global_barrier(); }
};

} // end namespace profugus

#endif // comm_Sync_hh

//---------------------------------------------------------------------------//
//              end of comm/Sync.hh
//---------------------------------------------------------------------------//
