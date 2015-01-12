//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   acc/Domain_Transporter.hh
 * \author Thomas M. Evans
 * \date   Sat Nov 01 10:56:35 2014
 * \brief  Domain_Transporter class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef acc_Domain_Transporter_hh
#define acc_Domain_Transporter_hh

#include <memory>

#include "core/geometry/Geometry.hh"
#include "Geometry.hh"
#include "Particle.hh"

namespace acc
{

//===========================================================================//
/*!
 * \class Domain_Transporter
 * \brief ACC-wrapper for the Domain_Transporter.
 */
//===========================================================================//

class Domain_Transporter
{
  public:
    //@{
    //! ACC typedefs.
    typedef std::shared_ptr<Geometry> SP_ACC_Geometry;
    typedef std::vector<Particle>     Vec_ACC_Particles;
    //@}

    //@{
    //! Standard typedefs.
    typedef std::shared_ptr<profugus::Core> SP_Geometry;
    //@}

  private:
    // >>> CPU DATA

    // ACC Geometry.
    SP_ACC_Geometry d_geometry;

    // >>> DEVICE DATA

  public:
    // Domain transporter.
    Domain_Transporter();

    // Destructor.
    ~Domain_Transporter();

    // >>> SET FUNCTIONS

    // Set the geometry (builds the ACC device geometry).
    void set(SP_Geometry geometry);

    // >>> TRANSPORT

    // Transport a vector of ACC-device particles.
    void transport(Vec_ACC_Particles &particles);
};

} // end namespace acc

#endif // acc_Domain_Transporter_hh

//---------------------------------------------------------------------------//
//                 end of Domain_Transporter.hh
//---------------------------------------------------------------------------//
