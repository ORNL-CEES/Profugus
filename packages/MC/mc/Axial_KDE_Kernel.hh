//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   MC/mc/kde/Axial_KDE_Kernel.cc
 * \author Gregory Davidson
 * \date   Mon Feb 16 14:21:15 2015
 * \brief  Axial_KDE_Kernel class definitions.
 * \note   Copyright (c) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef MC_mc_Axial_KDE_Kernel_hh
#define MC_mc_Axial_KDE_Kernel_hh

#include "KDE_Kernel.hh"

namespace profugus
{

//===========================================================================//
/*!
 * \class Axial_KDE_Kernel
 * \brief Samples a KDE position.
 */
/*!
 * \example mc/test/tstKDE_Kernel.cc
 *
 * Test of KDE_Kernel.
 */
//===========================================================================//

template<class Geometry>
class Axial_KDE_Kernel : public KDE_Kernel<Geometry>
{
    typedef KDE_Kernel<Geometry>  Base;

  public:
    //@{
    //! Useful typedefs
    typedef typename Base::SP_Geometry   SP_Geometry;
    typedef typename Base::SP_Physics    SP_Physics;
    typedef typename Base::Space_Vector  Space_Vector;
    typedef typename Base::size_type     size_type;
    typedef typename Base::cell_type     cell_type;
    //@}

    //! Enumeration describing the rejection method.  Fission rejection
    //! rejects if outside fissionable region.  Cell rejection rejects if
    //! outside original cell
    enum Reject_Method { FISSION_REJECTION, CELL_REJECTION };

  private:
    // Expose base class data members
    using Base::b_geometry;
    using Base::b_physics;
    using Base::b_coefficient;
    using Base::b_exponent;
    using Base::b_bndwidth_map;
    using Base::b_num_sampled;
    using Base::b_num_accepted;

    // Stores the rejection method
    Reject_Method d_method;

  public:
    //! Constructor
    Axial_KDE_Kernel(SP_Geometry   geometry,
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

#endif // MC_mc_Axial_KDE_Kernel_hh

//---------------------------------------------------------------------------//
// end of MC/mc/kde/Axial_KDE_Kernel.hh
//---------------------------------------------------------------------------//
