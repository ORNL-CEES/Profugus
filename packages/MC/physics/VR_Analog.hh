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

class VR_Analog : public Variance_Reduction
{
    typedef Variance_Reduction Base;

  public:
    //@{
    //! Base-class typedefs.
    typedef Base::Particle_t Particle_t;
    typedef Base::Bank_t     Bank_t;
    //@}

  public:
    // Constructor.
    explicit VR_Analog()
        : Base()
    {
        Base::b_splitting = false;
    }

    //! Do nothing at surfaces
    void post_surface(Particle_t& particle, events::Event& event, Bank_t& bank) 
	const { /* * */ }

    //! Do nothing at collisions
    void post_collision(Particle_t& particle, events::Event& event, Bank_t& bank) 
	const { /* * */ }
};

} // end namespace profugus

#endif // mc_VR_Analog_hh

//---------------------------------------------------------------------------//
//                 end of VR_Analog.hh
//---------------------------------------------------------------------------//
