//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/VR_Analog.hh
 * \author Thomas M. Evans
 * \date   Fri May 09 13:09:17 2014
 * \brief  VR_Analog class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef mc_VR_Analog_hh
#define mc_VR_Analog_hh

#include "Variance_Reduction.hh"

namespace profugus
{

//===========================================================================//
/*!
 * \class VR_Analog
 * \brief Analog MC transport.
 */
//===========================================================================//

template <class Geometry>
class VR_Analog : public Variance_Reduction<Geometry>
{
    typedef Variance_Reduction<Geometry> Base;

  public:
    //@{
    //! Base-class typedefs.
    typedef typename Base::Particle_t Particle_t;
    typedef typename Base::Bank_t     Bank_t;
    //@}

  public:
    // Constructor.
    explicit VR_Analog()
        : Base()
    {
        Base::b_splitting = false;
    }

    //! Do nothing at surfaces
    void post_surface(Particle_t& particle, Bank_t& bank) const { /* * */ }

    //! Do nothing at collisions
    void post_collision(Particle_t& particle, Bank_t& bank) const { /* * */ }
};

} // end namespace profugus

#endif // mc_VR_Analog_hh

//---------------------------------------------------------------------------//
//                 end of VR_Analog.hh
//---------------------------------------------------------------------------//
