//---------------------------------*-CUDA-*----------------------------------//
/*!
 * \file   cuda_mc/test/Uniform_Source_Tester.hh
 * \author Stuart Slattery
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_mc_test_Uniform_Source_Tester_hh
#define cuda_mc_test_Uniform_Source_Tester_hh

#include "cuda_utils/Shared_Device_Ptr.hh"

#include "rng/RNG.hh"

#include "cuda_geometry/Mesh_Geometry.hh"

#include "../Particle_Vector.hh"
#include "../Box_Shape.hh"

#include <Teuchos_Array.hpp>

#include <vector>

//---------------------------------------------------------------------------//
class Uniform_Source_Tester
{
  public:

    typedef cuda_profugus::Mesh_Geometry_DMM Geometry_DMM;
    typedef cuda_profugus::Mesh_Geometry Geometry;
    typedef cuda_profugus::Particle_Vector<Geometry> Particle_Vector;
    typedef typename Particle_Vector::Event_t Event_t;
    typedef cuda_profugus::Box_Shape Shape;

    // Constructor. Will build mesh geometry and particle vector.
    Uniform_Source_Tester( const std::vector<double>& x_edges,
			   const std::vector<double>& y_edges,
			   const std::vector<double>& z_edges,
			   const int matid,
			   const int vector_size,
			   const profugus::RNG& rng,
			   const int num_group );

    // Get the geometry.
    cuda_utils::Shared_Device_Ptr<Geometry> geometry() const
    { return d_geometry; }

    // Get the source shape.
    cuda_utils::Shared_Device_Ptr<Shape> shape() const
    { return d_shape; }

    // Get the particles
    cuda_utils::Shared_Device_Ptr<Particle_Vector> particles()
    { return d_particles; }

    // get a vector of matids.
    Teuchos::Array<int> matid();

    // get a vector of weights.
    Teuchos::Array<double> wt();

    // get a vector of alive status.
    Teuchos::Array<int> alive();

    // get a vector of groups.
    Teuchos::Array<int> group();

    // get a vector of events.
    Teuchos::Array<Event_t> event();

    // get a vector of batches.
    Teuchos::Array<int> batch();

    // kill all the particles so we can make a new batch.
    void kill_particles();

  private:

    int d_size;
    std::shared_ptr<Geometry_DMM> d_host_geom;
    cuda_utils::Shared_Device_Ptr<Geometry> d_geometry;
    cuda_utils::Shared_Device_Ptr<Shape> d_shape;
    cuda_utils::Shared_Device_Ptr<Particle_Vector> d_particles;
};

//---------------------------------------------------------------------------//

#endif // cuda_mc_test_Uniform_Source_Tester_hh

//---------------------------------------------------------------------------//
//                 end of cuda_utils/Uniform_Source_Tester.hh
//---------------------------------------------------------------------------//
