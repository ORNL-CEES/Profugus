//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   rng/RNG.hh
 * \author Thomas M. Evans
 * \date   Fri Apr 25 14:57:29 2014
 * \brief  RNG class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef rng_RNG_hh
#define rng_RNG_hh

#include <vector>

#include <Utils/config.h>
#include "harness/DBC.hh"

//---------------------------------------------------------------------------//
// Declare SPRNG functions here instead of polluting namespace with all SPRNG
// definitions
//---------------------------------------------------------------------------//
extern "C" {
int    get_rn_int(int *igenptr);
float  get_rn_flt(int *igenptr);
double get_rn_dbl(int *igenptr);
int    free_rng(int *igenptr);
}
//---------------------------------------------------------------------------//

namespace profugus
{

//===========================================================================//
/*!
 * \class RNG
 * \brief Random-number generator class interface to SPRNG.
 *
 * The SPRNG random number class is a wrapper class for the <a
 * href="http://sprng.cs.fsu.edu/">SPRNG (Scalable Parallel Random Number
 * Generator)</a>) libraries.  The class is nominally controlled by the
 * Rnd_Control class; however, SPRNG random number objects can be created
 * directly from native SPRNG functions.  See the examples for more detail.
 *
 * The primary purpose of the SPRNG class is to manage the memory that is
 * controlled by the SPRNG library. SPRNG accomplishes this through reference
 * counting.  SPRNG objects can be copied an assigned freely without worrying
 * about memory.  When the last SPRNG object that accesses a particular random
 * number state is destroyed, the memory controlling that state is released.
 * Note that when a SPRNG is copied, both instantiations of the SPRNG object
 * access the same random number state.  See the examples for more info.
 *
 * The implementation of SPRNG is a modified version of Sprng-0.5 that has
 * been optimized for LFG generator contruction.  This code has been shipped
 * with Profugus in the rng/sprng sub-directory under the provisions of the
 * SPRNG opensource license.
 *
 * \sa <a href="http://sprng.cs.fsu.edu/">SPRNG</a>, RNG_Control
 */
/*!
 * \example rng/test/tstRNG.cc
 *
 * Test of RNG.
 */
//===========================================================================//

class RNG
{
  private:
    // >>> TYPES

    //! Reference Counting Class For SPRNG Random Numbers.
    struct RNGValue
    {
        // Counter and id to RNG library.
        int *id;
        int refcount;

        // Constructor.
        RNGValue(int *idnum) :  id(idnum), refcount(1) {}

        // Destructor.
        ~RNGValue() { ::free_rng(id); }
    };

  private:
    // >>> DATA

    // Pointer to memory in SPRNG library.
    RNGValue *d_streamid;

    // Number of this particular stream.
    int d_stream;

    // Size of the packed state
    static int d_packed_size;

  public:
    // Constructors
    inline RNG() : d_streamid(0), d_stream(0) {}
    inline RNG(int *, int);
    inline RNG(const RNG &);
    RNG(const std::vector<char> &);

    // Destructor, reclaim memory from SPRNG library.
    inline ~RNG();

    // Is RNG assigned?
    bool assigned() const { return d_streamid != 0; }

    // Assignment operator.
    RNG& operator=(const RNG &);

    // Pack RNG state.
    std::vector<char> pack() const;

    // >>> Services provided by RNG class.

    //! Get a uniform random number on [0,1] (closed on both sides)
    double ran() const { return uniform_impl(Type_Switch<double>()); }

    //! Get a uniform random number on [0,1]
    template<class RealType>
    inline RealType uniform() const
    {
        return uniform_impl(Type_Switch<RealType>());
    }

    //! Return the SPRNG ID pointer.
    int* get_id() const { Require (d_streamid); return d_streamid->id; }

    //! Return the stream number.
    int get_num() const { Require (d_streamid); return d_stream; }

    // Return the packed size
    int get_size() const;

    //! Do a diagnostic print.
    void print() const;

  private:
    template<class T>
    struct Type_Switch
    {
        typedef T value_type;
    };

    // Return a double-precision uniform value
    double uniform_impl(Type_Switch<double>) const
    {
        return ::get_rn_dbl(d_streamid->id);
    }

    // Return a single-precision uniform value
    float uniform_impl(Type_Switch<float>) const
    {
        return ::get_rn_flt(d_streamid->id);
    }
};

//---------------------------------------------------------------------------//
// INLINE RNG MEMBERS
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
RNG::RNG(int *idval, int number)
    : d_streamid(new RNGValue(idval))
    , d_stream(number)
{
}

//---------------------------------------------------------------------------//
/*!
 * \brief Copy constructor.
 */
RNG::RNG(const RNG &rhs)
    : d_streamid(rhs.d_streamid)
    , d_stream(rhs.d_stream)
{
    if (d_streamid)
        ++d_streamid->refcount;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Destructor
 */
RNG::~RNG()
{
    if (d_streamid && --d_streamid->refcount == 0)
        delete d_streamid;
}

} // end namespace profugus

#endif // rng_RNG_hh

//---------------------------------------------------------------------------//
//                 end of RNG.hh
//---------------------------------------------------------------------------//
