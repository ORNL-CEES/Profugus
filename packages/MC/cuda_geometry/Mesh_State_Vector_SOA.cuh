//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   MC/cuda_geometry/Mesh_State_Vector_SOA.cuh
 * \author Steven Hamilton
 * \date   Tue Apr 11 11:22:34 2017
 * \brief  Mesh_State_Vector_SOA class declaration.
 * \note   Copyright (c) 2017 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef MC_cuda_geometry_Mesh_State_Vector_SOA_cuh
#define MC_cuda_geometry_Mesh_State_Vector_SOA_cuh

#include "CudaUtils/cuda_utils/Definitions.hh"
#include "CudaUtils/cuda_utils/Device_View_Field.hh"

namespace cuda_profugus
{

//===========================================================================//
/*!
 * \class Mesh_State_Vector_SOA
 * \brief Class for managing a collection of mesh states.
 */
//===========================================================================//

class Mesh_State_Vector_SOA
{
  public:
    //@{
    //! Typedefs
    typedef cuda_utils::Coordinates                 Coordinates;
    typedef cuda_utils::Space_Vector                Space_Vector;
    typedef cuda::Device_View_Field<int>            Int_View;
    typedef cuda::Device_View_Field<double>         Double_View;
    typedef cuda::Device_View_Field<Coordinates>    Coordinate_View;
    typedef cuda::Device_View_Field<Coordinates>    Space_Vector_View;
    //@}

  private:

    // >>> DATA
    Coordinate_View   d_ijk;
    Space_Vector_View d_r;
    Space_Vector_View d_dir;
    Coordinate_View   d_next_ijk;
    Double_View       d_next_dist;
    Int_View          d_exiting_face;
    Int_View          d_reflecting_face;

  public:

    // Constructor
    Mesh_State_Vector_SOA(
        Coordinate_View   ijk,
        Space_Vector_View r,
        Space_Vector_View dir,
        Coordinate_View   next_ijk,
        Double_View       next_dist,
        Int_View          exiting_face,
        Int_View          reflecting_face)
      : d_ijk(ijk)
      , d_r(r)
      , d_dir(dir)
      , d_next_ijk(ijk)
      , d_next_dist(next_dist)
      , d_exiting_face(exiting_face)
      , d_reflecting_face(reflecting_face)
    {
    }

    // Number of states allocated
    __device__ int size() const {return d_ijk.size();}

    // Coordinates
    __device__ Coordinates& ijk(int pid)       
    {
        DEVICE_REQUIRE(pid<size());
        return d_ijk[pid];
    }
    __device__ const Coordinates& ijk(int pid) const 
    {
        DEVICE_REQUIRE(pid<size());
        return d_ijk[pid];
    }

    // Position
    __device__ Space_Vector& pos(int pid)       
    {
        DEVICE_REQUIRE(pid<size());
        return d_r[pid];
    }
    __device__ const Space_Vector& pos(int pid) const 
    {
        DEVICE_REQUIRE(pid<size());
        return d_r[pid];
    }

    // Direction
    __device__ Space_Vector& dir(int pid)       
    {
        DEVICE_REQUIRE(pid<size());
        return d_dir[pid];
    }
    __device__ const Space_Vector& dir(int pid) const 
    {
        DEVICE_REQUIRE(pid<size());
        return d_dir[pid];
    }

    // Next coordinates
    __device__ Coordinates& next_ijk(int pid)       
    {
        DEVICE_REQUIRE(pid<size());
        return d_next_ijk[pid];
    }
    __device__ const Coordinates& next_ijk(int pid) const 
    {
        DEVICE_REQUIRE(pid<size());
        return d_next_ijk[pid];
    }

    // Next distance
    __device__ double& next_dist(int pid)             
    {
        DEVICE_REQUIRE(pid<size());
        return d_next_dist[pid];
    }
    __device__ const double& next_dist(int pid) const 
    {
        DEVICE_REQUIRE(pid<size());
        return d_next_dist[pid];
    }

    // Exiting face
    __device__ int& exiting_face(int pid)       
    {
        DEVICE_REQUIRE(pid<size());
        return d_exiting_face[pid];
    }
    __device__ const int& exiting_face(int pid) const 
    {
        DEVICE_REQUIRE(pid<size());
        return d_exiting_face[pid];
    }

    // Reflecting face
    __device__ int& reflecting_face(int pid)       
    {
        DEVICE_REQUIRE(pid<size());
        return d_reflecting_face[pid];
    }
    __device__ const int& reflecting_face(int pid) const 
    {
        DEVICE_REQUIRE(pid<size());
        return d_reflecting_face[pid];
    }

};

//===========================================================================//
/*!
 * \class Mesh_State_Vector_SOA_DMM
 * \brief Device memory manager for Mesh_State_Vector_SOA.
 */
//===========================================================================//

class Mesh_State_Vector_SOA_DMM :
    public cuda::Device_Memory_Manager<Mesh_State_Vector_SOA>
{
  public:

    typedef Mesh_State_Vector_SOA Geo_State_Vector_t;

    //! Typedefs
    typedef cuda_utils::Coordinates             Coordinates;
    typedef cuda_utils::Space_Vector            Space_Vector;
    typedef thrust::device_vector<int>          Int_Vector;
    typedef thrust::device_vector<double>       Double_Vector;
    typedef thrust::device_vector<Coordinates>  Coordinate_Vector;
    typedef thrust::device_vector<Coordinates>  Space_Vector_Vector;

    // Constructor
    Mesh_State_Vector_SOA_DMM(){}

    // Number of states currently allocated
    int size() const {return d_ijk.size();}

    // Memory manager interface
    Geo_State_Vector_t device_instance()
    {
        return Geo_State_Vector_t(cuda::make_view(d_ijk),
                                  cuda::make_view(d_r),
                                  cuda::make_view(d_dir),
                                  cuda::make_view(d_next_ijk),
                                  cuda::make_view(d_next_dist),
                                  cuda::make_view(d_exiting_face),
                                  cuda::make_view(d_reflecting_face));
    }

    // Initialize vectors for specified number of states
    void initialize(int num_states)
    {
        REQUIRE(num_states > 0);
        d_ijk.resize(num_states);
        d_r.resize(num_states);
        d_dir.resize(num_states);
        d_next_ijk.resize(num_states);
        d_next_dist.resize(num_states);
        d_exiting_face.resize(num_states);
        d_reflecting_face.resize(num_states);
    }

  private:

    // >>> DATA
    Coordinate_Vector   d_ijk;
    Space_Vector_Vector d_r;
    Space_Vector_Vector d_dir;
    Coordinate_Vector   d_next_ijk;
    Double_Vector       d_next_dist;
    Int_Vector          d_exiting_face;
    Int_Vector          d_reflecting_face;
};

//---------------------------------------------------------------------------//
} // end namespace cuda_profugus

//---------------------------------------------------------------------------//
#endif // MC_cuda_geometry_Mesh_State_Vector_SOA_cuh

//---------------------------------------------------------------------------//
// end of MC/cuda_geometry/Mesh_State_Vector_SOA.cuh
//---------------------------------------------------------------------------//
