//---------------------------------*-CUDA-*----------------------------------//
/*!
 * \file   cuda_mc/test/Physics_Tester.hh
 * \author Stuart Slattery
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_mc_test_Physics_Tester_hh
#define cuda_mc_test_Physics_Tester_hh

#include "cuda_utils/Shared_Device_Ptr.hh"

#include "rng/RNG.hh"

#include "cuda_geometry/Mesh_Geometry.hh"

#include "xs/XS.hh"

#include "../Physics.hh"

#include "Particle_Vector_Tester.hh"

#include <Teuchos_Array.hpp>
#include <Teuchos_ParameterList.hpp>

#include <vector>

//---------------------------------------------------------------------------//
class Physics_Tester
{
  public:

    typedef cuda_profugus::Mesh_Geometry Geometry;
    typedef cuda_profugus::Particle_Vector<Geometry> Particle_Vector;
    typedef typename Particle_Vector::Event_t Event_t;
    typedef cuda_profugus::Physics<Geometry> Physics;
    typedef typename Geometry::Space_Vector Space_Vector;

    // Constructor. Will build mesh geometry and particle vector.
    Physics_Tester( const std::vector<double>& x_edges,
		    const std::vector<double>& y_edges,
		    const std::vector<double>& z_edges,
		    const int vector_size,
		    const profugus::RNG& rng,
		    Teuchos::ParameterList& db,
		    const profugus::XS& xs );

    // Get the physics.
    cuda::Shared_Device_Ptr<Physics>& physics()
    { return d_physics; }

    // Get the particles
    cuda::Shared_Device_Ptr<Particle_Vector>& particles()
    { return d_particle_tester.get_vector(); }

    // Get the particle vector tester.
    Particle_Vector_Tester& particle_tester()
    { return d_particle_tester; }

    // Initialize particles with the geometry.
    void geometry_initialize( const Space_Vector r, const Space_Vector d,
			      const int matid );

    // Sample a cdf and set the particle group.
    void sample_group( const std::vector<double>& cdf );

    // Check if a matid is fissionable.
    bool is_fissionable( const int matid ) const;

    // Get a total cross section.
    double get_total( const int matid,
		      const int group,
		      const cuda_profugus::physics::Reaction_Type type ) const;

  private:
    
    cuda::Shared_Device_Ptr<Physics> d_physics;
    cuda::Shared_Device_Ptr<Geometry> d_geometry;
    Particle_Vector_Tester d_particle_tester;
    int d_size;    
};

//---------------------------------------------------------------------------//

#endif // cuda_mc_test_Physics_Tester_hh

//---------------------------------------------------------------------------//
//                 end of cuda_utils/Physics_Tester.hh
//---------------------------------------------------------------------------//
