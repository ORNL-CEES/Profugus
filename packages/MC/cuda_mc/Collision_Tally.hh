//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/Collision_Tally.hh
 * \author Stuart Slattery
 * \brief  Collision class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_mc_Collision_Tally_hh
#define cuda_mc_Collision_Tally_hh

#include <Teuchos_Array.hpp>

#include "utils/Definitions.hh"

#include "cuda_utils/Shared_Device_Ptr.hh"

#include "mc/Definitions.hh"

#include "Particle_Vector.hh"

namespace cuda_profugus
{
//===========================================================================//
/*!
 * \class Collision_Tally
 * \brief Collision_Tally class for MC transport with batch statistics.
 */
/*!
 * \example mc/test/tstCollision_Tally.cc
 *
 * Test of Particle.
 */
//===========================================================================//
template <class Geometry>
class Collision_Tally
{
  public:
    //@{
    //! Typedefs.
    typedef typename Geometry::Geo_State_t Geo_State_t;
    typedef profugus::events::Event Event_t;
    //@}

  private:
    // >>> DATA

    // Geometry over which we are tallying.
    cuda::Shared_Device_Ptr<Geometry> d_geometry;

    // Number of statistical batches in the tally.
    int d_num_batch;

    // Number of cells in the tally.
    int d_num_cells;

    // Cell tallies. Indexed as d_tally[batch][cell]. On-device.
    double* d_tally;

  public:

    // >>> HOST API

    // Constructor
    Collision_Tally( const cuda::Shared_Device_Ptr<Geometry>& geometry, 
		     const int num_batch );
    
    // Destructor.
    ~Collision_Tally();

    // Get the number of cells in the tally.
    int num_cells() const { return d_num_cells; }
    
    // Get the number of batches in the tally.
    int num_batch() const { return d_num_batch; }

    // Tally the particles in a vector.
    void tally( const cuda::Shared_Device_Ptr<Particle_Vector>& particles );

    // Finalize the tally.
    void finalize( const std::size_t total_num_particle );

    // Copy the first and second tally moments from the device to the
    // host. The moments are lazy-evaluated in this function and indexed by
    // cell.
    void copy_moments_to_host( Teuchos::Array<double>& first_moment,
			       Teuchos::Array<double>& second_moment ) const;
};

//---------------------------------------------------------------------------//

} // end namespace cuda_profugus

#endif // cuda_mc_Collision_Tally_hh

//---------------------------------------------------------------------------//
//                 end of Collision_Tally.hh
//---------------------------------------------------------------------------//
