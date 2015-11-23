//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   MC/mc/kde/KDE_Kernel_Resample.cc
 * \author Gregory Davidson
 * \date   Mon Feb 16 14:21:15 2015
 * \brief  KDE_Kernel_Resample class definitions.
 * \note   Copyright (c) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef MC_mc_KDE_Kernel_Resample_hh
#define MC_mc_KDE_Kernel_Resample_hh

#include "KDE_Kernel.hh"

namespace profugus
{

class KDE_Kernel_Resample : public KDE_Kernel
{
    typedef KDE_Kernel  Base;

  public:
    //@{
    //! Useful typedefs
    typedef Base::SP_Geometry   SP_Geometry;
    typedef Base::SP_Physics    SP_Physics;
    typedef Base::Space_Vector  Space_Vector;
    //@}

  public:
    //! Constructor
    KDE_Kernel_Resample(SP_Geometry geometry,
                        SP_Physics  physics,
                        double      coefficient = 1.06,
                        double      exponent = -0.20);

    //! Sample a new position
    virtual Space_Vector sample_position(const Space_Vector &orig_position,
                                         RNG                &rng) const;
};

} // end namespace profugus

#endif // MC_mc_KDE_Kernel_Resample_hh

//---------------------------------------------------------------------------//
// end of MC/mc/kde/KDE_Kernel_Resample.hh
//---------------------------------------------------------------------------//
