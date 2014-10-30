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
    const long d_a;
    const long d_c;
    const long d_m;

    // Seed
    long d_seed;

  public:
    RNG(long seed)
        : d_a(1664525)
        , d_c(1013904223)
        , d_m(4294967296)
        , d_seed(seed)
    {
#pragma acc enter data pcopyin(this)
    }

    ~RNG()
    {
#pragma acc exit data delete(this)
    }

    // Get the random number.
    double ran()
    {
        d_seed = (d_a * d_seed + d_c) % d_m;
        return d_seed / static_cast<double>(d_m);
    }
};

} // end namespace acc

#endif // acc_RNG_hh

//---------------------------------------------------------------------------//
//                 end of RNG.hh
//---------------------------------------------------------------------------//
