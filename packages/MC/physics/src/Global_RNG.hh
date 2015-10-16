//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/Global_RNG.hh
 * \author Thomas M. Evans
 * \date   Monday May 5 16:56:6 2014
 * \brief  Global_RNG class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef mc_Global_RNG_hh
#define mc_Global_RNG_hh

#include <memory>
#include "rng/RNG_Control.hh"

namespace profugus
{

//===========================================================================//
/*!
 * \class Global_RNG
 * \brief Static class that can be used to store globally-available RNGs.
 */
//===========================================================================//

class Global_RNG
{
  public:
    //@{
    //! Typedefs.
    typedef RNG_Control                    RNG_Control_t;
    typedef std::shared_ptr<RNG_Control_t> SP_RNG_Control;
    typedef RNG_Control_t::RNG_t           RNG_t;
    //@}

    //! RNG available to entire domain.
    static RNG_t d_rng;
};

} // end namespace profugus

#endif // mc_Global_RNG_hh

//---------------------------------------------------------------------------//
//              end of mc/Global_RNG.hh
//---------------------------------------------------------------------------//
