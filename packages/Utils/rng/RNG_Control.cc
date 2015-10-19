//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   rng/RNG_Control.cc
 * \author Thomas M. Evans
 * \date   Fri Apr 25 14:56:51 2014
 * \brief  RNG_Control member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "RNG_Control.hh"

#include <cstdlib>
#include <vector>

#include "sprng/sprng.h"

namespace profugus
{
//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief RNG_Control class constructor.
 *
 * Constructor for the RNG_Control class only requires one argument with 3
 * additional options.
 *
 * \param seed seed value
 * \param number total number of independent streams; default = 1.0e9
 * \param stream stream index; default = 0
 * \param parameter SPRNG parameter that determines the size of the SPRNG
 *        state, see the SPRNG library documentation for details; default = 1
 */
RNG_Control::RNG_Control(int seed,
                         int number,
                         int parameter)
    : d_seed(seed)
    , d_number(number)
    , d_parameter(parameter)
{
    // make a spring object and pack it to determine the size
    int *id = init_sprng(0, d_number, d_seed, d_parameter);
    RNG temp( id, 0 );

    // pack it and set size
    std::vector<char> pack = temp.pack();
    d_size                 = pack.size();
    CHECK(d_size >= 0);
}

//---------------------------------------------------------------------------//
// MEMBER FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Create a SPRNG random number object at a given stream index.
 *
 * \param stream a user-requested random number index
 * \return SPRNG random number object
 *
 * Does not update the state of this object.
 */
RNG_Control::RNG_t RNG_Control::rng(int stream) const
{
    REQUIRE(stream <= d_number);

    // declare a stream
    int *id = init_sprng(stream, d_number, d_seed, d_parameter);

    // create a new RNG object
    RNG random( id, stream );

    // return the object
    return random;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Spawn a new random number object.
 *
 * \param random valid SPRNG random number object
 * \return SPRNG random number object
 *
 * This function spawns a new, independent random number object from an
 * existing SPRNG random number object.  The new SPRNG random number state
 * has as its index the index of the parent SPRNG object from which it has
 * been spawned.
 *
 * Memory control of the spawned stream is given to the returned SPRNG
 * object.  The spawned object is--from this point on--independent of the
 * stream from which it spawned.  This class is used to guarantee
 * reproducibility in random number codes.  By spawning only from existing
 * random number states, reproducibility in numeric operations can be
 * achieved.
*/
RNG_Control::RNG_t RNG_Control::spawn(const RNG_t &random) const
{
    // declare variables necessary to spawn a stream
    int **newstream;

    // spawn a new stream
    spawn_sprng(random.get_id(), 1, &newstream);

    // create a new SPRNG random number object with new stream
    RNG ran(newstream[0], random.get_num());

    // free the memory
    std::free(newstream);

    // return the new random object
    return ran;
}

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of RNG_Control.cc
//---------------------------------------------------------------------------//
