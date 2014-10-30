//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   acc/RNG.hh
 * \author Thomas M. Evans
 * \date   Thu Oct 30 13:18:58 2014
 * \brief  RNG class.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef acc_RNG_hh
#define acc_RNG_hh

namespace acc
{

//===========================================================================//
/*!
 * \class RNG
 * \brief A crappy random number generator for the gpu.
 */
//===========================================================================//

class RNG
{
  private:
    // >>> DATA

    // LCG constants.
    long d_a;
    long d_c;
    long d_m;

    // Seed
    long d_seed;

  public:
    explicit RNG(long seed);
    ~RNG();

    //! Get the random number.
    double ran()
    {
        d_seed = (d_a * d_seed + d_c) % d_m;
        return d_seed / static_cast<double>(d_m);
    }

    // Ran with seed.
    double ran(long &seed)
    {
        seed = (d_a * seed + d_c) % d_m;
        return seed / static_cast<double>(d_m);
    }
};

} // end namespace acc

#endif // acc_RNG_hh

//---------------------------------------------------------------------------//
//                 end of RNG.hh
//---------------------------------------------------------------------------//
