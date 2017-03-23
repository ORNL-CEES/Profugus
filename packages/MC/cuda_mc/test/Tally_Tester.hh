//---------------------------------*-CUDA-*----------------------------------//
/*!
 * \file   cuda_mc/test/Tally_Tester.hh
 * \author Stuart Slattery
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_mc_test_Tally_Tester_hh
#define cuda_mc_test_Tally_Tester_hh

#include "cuda_utils/Shared_Device_Ptr.hh"

#include "rng/RNG.hh"

#include "cuda_geometry/Mesh_Geometry.hh"

#include "../Particle_Vector.hh"

#include <vector>

//---------------------------------------------------------------------------//
class Tally_Tester
{
  public:

    typedef cuda_profugus::Mesh_Geometry_DMM Geometry_DMM;
    typedef cuda_profugus::Mesh_Geometry Geometry;
    typedef cuda_profugus::Particle_Vector<Geometry> Particle_Vector;
    typedef typename Particle_Vector::Event_t Event_t;
    typedef typename Geometry::Geo_State_t Geo_State_t;

    // Constructor. Will build mesh geometry and particle vector.
    Tally_Tester( const std::vector<double>& x_edges,
			    const std::vector<double>& y_edges,
			    const std::vector<double>& z_edges,
			    const int num_particle,
			    const profugus::RNG& rng,
			    const int num_batch );

    // Get the geometry.
    cuda_utils::Shared_Device_Ptr<Geometry> geometry() const
    { return d_geometry; }

    // Get the particles
    cuda_utils::Shared_Device_Ptr<Particle_Vector> particles() const
    { return d_particles; }

  private:

    cuda_utils::Shared_Device_Ptr<Geometry> d_geometry;
    cuda_utils::Shared_Device_Ptr<Particle_Vector> d_particles;
};

//---------------------------------------------------------------------------//

#endif // cuda_mc_test_Tally_Tester_hh

//---------------------------------------------------------------------------//
//                 end of cuda_utils/Tally_Tester.hh
//---------------------------------------------------------------------------//
