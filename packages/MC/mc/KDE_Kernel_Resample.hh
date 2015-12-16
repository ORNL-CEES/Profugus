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


//===========================================================================//
/*!
 * \class KDE_Kernel_Resample
 * \brief Samples a KDE position.
 */
/*!
 * \example mc/test/tstKDE_Kernel.cc
 *
 * Test of KDE_Kernel.
 */
//===========================================================================//


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

    //! Enumeration describing the rejection method.  Fission rejection
    //! rejects if outside fissionable region.  Cell rejection rejects if
    //! outside original cell
    enum Reject_Method { FISSION_REJECTION, CELL_REJECTION };

  private:
    Reject_Method d_method;

  public:
    //! Constructor
    KDE_Kernel_Resample(SP_Geometry   geometry,
                        SP_Physics    physics,
                        Reject_Method method,
                        double        coefficient = 1.06,
                        double        exponent = -0.20);

    //! Return the rejection method
    Reject_Method reject_method() const { return d_method; }

    //! Sample a new position
    virtual Space_Vector sample_position(const Space_Vector &orig_position,
                                         RNG                &rng) const;

  private:
    // Sample using fission rejection
    Space_Vector sample_position_fiss_rej(const Space_Vector &orig_position,
                                          RNG                &rng) const;

    // Sample using cell rejection
    Space_Vector sample_position_cell_rej(const Space_Vector &orig_position,
                                          RNG                &rng) const;
};

} // end namespace profugus

#endif // MC_mc_KDE_Kernel_Resample_hh

//---------------------------------------------------------------------------//
// end of MC/mc/kde/KDE_Kernel_Resample.hh
//---------------------------------------------------------------------------//
