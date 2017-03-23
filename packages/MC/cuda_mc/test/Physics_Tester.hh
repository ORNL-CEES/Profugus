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
#include "../Box_Shape.hh"

#include "Particle_Vector_Tester.hh"

#include <Teuchos_Array.hpp>
#include <Teuchos_ParameterList.hpp>

#include <vector>

//---------------------------------------------------------------------------//
class Physics_Tester
{
  public:

    typedef cuda_profugus::Cartesian_Mesh Cartesian_Mesh;
    typedef cuda_profugus::Mesh_Geometry_DMM Geometry_DMM;
    typedef cuda_profugus::Mesh_Geometry Geometry;
    typedef cuda_profugus::Particle_Vector<Geometry> Particle_Vector;
    typedef typename Particle_Vector::Event_t Event_t;
    typedef cuda_profugus::Physics<Geometry> Physics;
    typedef Physics::Fission_Site Fission_Site;
    typedef typename Geometry::Space_Vector Space_Vector;
    typedef cuda_profugus::Box_Shape Shape;

    // Constructor. Will build mesh geometry and particle vector.
    Physics_Tester( const std::vector<double>& x_edges,
		    const std::vector<double>& y_edges,
		    const std::vector<double>& z_edges,
		    const int vector_size,
		    const profugus::RNG& rng,
		    Teuchos::ParameterList& db,
		    const profugus::XS& xs,
		    const int matid );

    // Get the physics.
    cuda_utils::Shared_Device_Ptr<Physics>& physics()
    { return d_physics; }

    // Get the geometry.
    cuda_utils::Shared_Device_Ptr<Geometry>& geometry()
    { return d_geometry; }

    // Get the cartesian mesh under the geometry.
    cuda_utils::Shared_Device_Ptr<Cartesian_Mesh>& cart_mesh()
    { return d_cart_mesh; }

    // Get the particles
    cuda_utils::Shared_Device_Ptr<Particle_Vector>& particles()
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

    // Get the min and max particle energies.
    void get_min_max_energy( double& min, double& max ) const;

    // Initialize a particle from a fission spectrum.
    void initialize_fission_from_spectrum( const int matid,
					   const double ran,
					   int& group,
					   bool& sampled ) const;

    // Initialize a particle from a fission site.
    void initialize_fission_from_site( const Fission_Site &fs,
				       const double ran,
				       int& group,
				       bool& sampled ) const;

    // Get the source shape.
    cuda_utils::Shared_Device_Ptr<Shape> source_shape() const
    { return d_shape; }

  private:

    cuda_utils::Shared_Device_Ptr<Physics> d_physics;
    cuda_utils::Shared_Device_Ptr<Geometry> d_geometry;
    cuda_utils::Shared_Device_Ptr<Cartesian_Mesh> d_cart_mesh;
    cuda_utils::Shared_Device_Ptr<Shape> d_shape;
    Particle_Vector_Tester d_particle_tester;
    int d_size;
};

//---------------------------------------------------------------------------//

#endif // cuda_mc_test_Physics_Tester_hh

//---------------------------------------------------------------------------//
//                 end of cuda_utils/Physics_Tester.hh
//---------------------------------------------------------------------------//
