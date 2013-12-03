//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   comm/SpinLock.hh
 * \author Thomas M. Evans
 * \date   Wed Jan  2 16:14:33 2008
 * \brief  SpinLock class definitions.
 * \note   Copyright (C) 2008 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef comm_SpinLock_hh
#define comm_SpinLock_hh

#include <iostream>
#include <cstdio>

#include "Sync.hh"

namespace nemesis
{

//===========================================================================//
/*!
 * \class SpinLock
 * \brief Serialize execution of a block.
 *
 * This class enables you to get a block of code to execute serially.  Each
 * processor begins to executes the block only after the one before it is
 * finished.  A processor sub-range, [begin, end), can be entered to sync
 * over.
 */
/*!
 * \example comm/test/tstSpinLock.cc
 */
//===========================================================================//

class SpinLock
{
    SpinLock( const SpinLock& );
    SpinLock& operator=( const SpinLock& );

    enum { SL_Next = 92874 };

    int d_trash;
    int d_lock;

    int d_node;
    int d_nodes;

    int d_first;
    int d_last;

  public:
    SpinLock(int _lock = 1);
    SpinLock(int begin, int end);
    virtual ~SpinLock();
};

//===========================================================================//
/*!
 * \class HSyncSpinLock
 * \brief Serialize a block, syncing at top.
 *
 * A spinlock that forces a global sync at the head of the block. A processor
 * sub-range, [begin, end), can be entered to sync over.
 */
//===========================================================================//

class HSyncSpinLock : public HSync, public SpinLock
{
    HSyncSpinLock(const HSyncSpinLock &);
    HSyncSpinLock& operator=(const HSyncSpinLock &);

  public:
    HSyncSpinLock(int l = 1) : HSync(l), SpinLock(l) {}
    HSyncSpinLock(int begin, int end) : HSync(1), SpinLock(begin, end) {}
};

//===========================================================================//
/*!
 * \class TSyncSpinLock
 * \brief Serialize a block, syncing at bottom.
 *
 * A spinlock that forces a global sync at the tail of the block. A processor
 * sub-range, [begin, end), can be entered to sync over.
 */
//===========================================================================//

class TSyncSpinLock : public TSync, public SpinLock
{
    TSyncSpinLock(const TSyncSpinLock &);
    TSyncSpinLock& operator=(const TSyncSpinLock &);

  public:
    TSyncSpinLock(int l = 1) : TSync(l), SpinLock(l) {}
    TSyncSpinLock(int begin, int end) : TSync(1), SpinLock(begin, end) {}
};

//===========================================================================//
/*!
 * \class HTSyncSpinLock
 * \brief Serialize a block, syncing at top and bottom.
 *
 * A spinlock that forces a global sync at the head and tail of the block. A
 * processor sub-range, [begin, end), can be entered to sync over.
 */
//===========================================================================//

class HTSyncSpinLock : public HSync, public TSync, public SpinLock
{
    HTSyncSpinLock(const HTSyncSpinLock &);
    HTSyncSpinLock& operator=(const HTSyncSpinLock &);

  public:
    HTSyncSpinLock(int l = 1) : HSync(l), TSync(l), SpinLock(l) {}
    HTSyncSpinLock(int b, int e) : HSync(1),  TSync(1), SpinLock(b, e) {}
};

} // end namespace nemesis

//! SPINLOCK macro.
#define SPINLOCK(a)            \
{                              \
    std::cout << std::flush;   \
    std::fflush(std::stdout);  \
    nemesis::HTSyncSpinLock h; \
    a;                         \
    std::cout << std::flush;   \
}

#endif // comm_SpinLock_hh

//---------------------------------------------------------------------------//
//              end of comm/SpinLock.hh
//---------------------------------------------------------------------------//
