//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Utils/rng/RNG_Control.hh
 * \author Thomas M. Evans
 * \date   Fri Apr 25 14:56:51 2014
 * \brief  RNG_Control class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef Utils_rng_RNG_Control_hh
#define Utils_rng_RNG_Control_hh

#include "harness/DBC.hh"
#include "RNG.hh"

namespace profugus
{

//===========================================================================//
/*!
 * \class RNG_Control
 * \brief Initialize and control the creation of RNG random rumber objects.
 *
 * The RNG_Control class is used to manage the RNG random number classes.  The
 * class' main function is to create a new RNG random number object.  By using
 * the RNG_Control class, users are not required to "ever" access the SPRNG
 * random number library explicitly. RNG_Control calls all the appropriate
 * SPRNG library functions to set up a random number state.  Ownership of the
 * random number state is maintained by the created SPRNG class object. By
 * using the controller, all memory management \a vis-a-vis the SPRNG library
 * is taken care of automatically.
 *
 * The usage is simple, make a RNG_Control object and then query it for SPRNG
 * random number objects.  Once a SPRNG object is built, memory management of
 * the SPRNG library associated with that state is taken care of by the RNG
 * object.  The RNG_Control object keeps track of the number of independent *
 * streams created through the private member RNG_Control::streamnum (stream *
 * index).  This index is automatically incremented by one everytime a random
 * number is created. The linear progression of random number states can be
 * interupted by giving an optional stream index (integer) to the get_rn(int)
 * function or by reseting the random number stream index through the
 * set_num(int) function.
 *
 * \sa <a href="http://sprng.cs.fsu.edu/">SPRNG (Scalable
 * Parallel Random Number Generator Library)</a>, SPRNG
 */
/*!
 * \example rng/test/tstRNG_Control.cc
 *
 * Test of RNG_Control.
 */
//===========================================================================//

class RNG_Control
{
  public:
    //! Random number type.
    typedef RNG RNG_t;

    // Global RNG stream number counter.
    static int rn_stream;

  private:
    // >>> DATA

    // Seed for initialization of random number streams.
    int d_seed;

    // Total number of streams.
    int d_number;

    // Number of current stream.
    int d_stream;

    // Control parameter for stream inits.
    int d_parameter;

    // Size of packed stream state.
    int d_size;

  public:
    // Constructor.
    RNG_Control(int seed, int number = 1000000000, int stream = 0,
                int parameter = 1);

    // Create SPRNG objects.
    RNG_t rng();
    RNG_t rng(int stream);

    // Spawn a new random number object.
    RNG_t spawn(const RNG_t &) const;

    //! Query for the current random number stream index.
    int get_num() const { return d_stream; }

    //! Set (reset) the random number stream index.
    void set_num(int stream) { REQUIRE(stream < d_number); d_stream = stream; }

    //! Query size of a packed random number state.
    int get_size() const { return d_size; }

    //! Get the seed value used to initialize the SPRNG library.
    int get_seed() const { return d_seed; }

    //! Return the total number of current streams set.
    int get_number() const { return d_number; }

  private:
    // >>> IMPLEMENTATION

    // Make a generator.
    RNG_t make_random_number_generator() const;
};

} // end namespace profugus

#endif // Utils_rng_RNG_Control_hh

//---------------------------------------------------------------------------//
//                 end of RNG_Control.hh
//---------------------------------------------------------------------------//
