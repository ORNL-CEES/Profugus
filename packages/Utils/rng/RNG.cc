//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   rng/RNG.cc
 * \author Thomas M. Evans
 * \date   Fri Apr 25 14:57:29 2014
 * \brief  RNG member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "RNG.hh"

#include <cstdlib>

#include "utils/Packing_Utils.hh"
#include "sprng/sprng.h"

namespace profugus
{

//---------------------------------------------------------------------------//
// STATIC MEMBERS
//---------------------------------------------------------------------------//

int RNG::d_packed_size = 0;

//---------------------------------------------------------------------------//
// CONSTRUCTORS
//---------------------------------------------------------------------------//
/*!
 * \brief Unpacking constructor.
 */
RNG::RNG(const std::vector<char> &packed)
    : d_streamid(0)
    , d_stream(0)
{
    Require (packed.size() >= 2 * sizeof(int));

    // make an unpacker
    profugus::Unpacker u;

    // set the buffer
    u.set_buffer(packed.size(), &packed[0]);

    // unpack the stream num and size of state
    int rng_size = 0;
    u >> d_stream >> rng_size;
    Check (d_stream >= 0);
    Check (rng_size >= 0);

    // unpack the random stream state
    char *prng = new char[rng_size];
    for (int i = 0; i < rng_size; i++)
        u >> prng[i];

    // now rebuild the sprng object
    int *rnid = unpack_sprng(prng);

    // now make a new streamid
    d_streamid = new RNGValue(rnid);

    // reclaim memory
    delete [] prng;

    Ensure (u.get_ptr() == &packed[0] + packed.size());
    Ensure (d_streamid);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Pack a RNG object into a vector<char>.
 */
std::vector<char> RNG::pack() const
{
    Require (d_streamid);

    // make a packer
    profugus::Packer p;

    // first pack the random number object and determine the size
    char *prng   = 0;
    int rng_size = pack_sprng(d_streamid->id, &prng);
    int size     = rng_size + 2 * sizeof(int);
    Check (prng);

    // now set the buffer
    std::vector<char> packed(size);
    p.set_buffer(size, &packed[0]);

    // pack the stream number and rng size
    p << d_stream << rng_size;

    // pack the stream state
    for (int i = 0; i < rng_size; i++)
        p << prng[i];

    // free the prng buffer
    std::free(prng);

    Ensure (p.get_ptr() == &packed[0] + size);

    return packed;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Assignment operator.
 */
RNG& RNG::operator=(const RNG &rhs)
{
    // check to see if the values are the same
    if (d_streamid == rhs.d_streamid && d_stream == rhs.d_stream)
        return *this;

    // destroy this' value if it was the last random number in a particular
    // stream
    if (d_streamid && --d_streamid->refcount == 0)
        delete d_streamid;

    // do assignment
    d_streamid = rhs.d_streamid;
    d_stream   = rhs.d_stream;

    // increment count
    if (d_streamid)
        ++d_streamid->refcount;

    // return
    return *this;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Return the size of the packed character stream.
 *
 * \return the size
 */
int RNG::get_size() const
{
    Require(d_streamid);

    if (d_packed_size > 0)
        return d_packed_size;

    char *prng    = 0;
    int rng_size  = pack_sprng(d_streamid->id, &prng);
    d_packed_size = rng_size + 2 * sizeof(int);

    // clear memory
    std::free(prng);

    return d_packed_size;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Do a diagnostic print.
 */
void RNG::print() const
{
    Require (d_streamid);
    print_sprng(d_streamid->id);
}

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of RNG.cc
//---------------------------------------------------------------------------//
