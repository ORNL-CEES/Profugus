//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   MC/cuda_geometry/Mesh_State_Vector_AOS.cuh
 * \author Steven Hamilton
 * \date   Tue Apr 11 11:22:34 2017
 * \brief  Mesh_State_Vector_AOS class declaration.
 * \note   Copyright (c) 2017 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef MC_cuda_geometry_Mesh_State_Vector_AOS_cuh
#define MC_cuda_geometry_Mesh_State_Vector_AOS_cuh

#include "CudaUtils/cuda_utils/Definitions.hh"
#include "CudaUtils/cuda_utils/Device_View_Field.hh"
#include "Mesh_State.hh"

namespace cuda_profugus
{

//===========================================================================//
/*!
 * \class Mesh_State_Vector_AOS
 * \brief Class for managing a collection of mesh states.
 */
//===========================================================================//

class Mesh_State_Vector_AOS
{
  public:
    //@{
    //! Typedefs
    typedef cuda::Device_View_Field<Mesh_State> State_View;
    typedef cuda_utils::Coordinates             Coordinates;
    typedef cuda_utils::Space_Vector            Space_Vector;
    //@}

  private:

    // >>> DATA
    State_View d_states;

  public:

    // Constructor
    Mesh_State_Vector_AOS(State_View states)
      : d_states(states)
    {
    }

    // Number of states allocated
    __device__ int size() const {return d_states.size();}

    // Coordinates
    __device__ Coordinates& ijk(int pid)       
    {
        DEVICE_REQUIRE(pid<size());
        return d_states[pid].ijk;
    }
    __device__ const Coordinates& ijk(int pid) const 
    {
        DEVICE_REQUIRE(pid<size());
        return d_states[pid].ijk;
    }

    // Position
    __device__ Space_Vector& pos(int pid)       
    {
        DEVICE_REQUIRE(pid<size());
        return d_states[pid].d_r;
    }
    __device__ const Space_Vector& pos(int pid) const 
    {
        DEVICE_REQUIRE(pid<size());
        return d_states[pid].d_r;
    }

    // Direction
    __device__ Space_Vector& dir(int pid)       
    {
        DEVICE_REQUIRE(pid<size());
        return d_states[pid].d_dir;
    }
    __device__ const Space_Vector& dir(int pid) const 
    {
        DEVICE_REQUIRE(pid<size());
        return d_states[pid].d_dir;
    }

    // Next coordinates
    __device__ Coordinates& next_ijk(int pid)       
    {
        DEVICE_REQUIRE(pid<size());
        return d_states[pid].next_ijk;
    }
    __device__ const Coordinates& next_ijk(int pid) const 
    {
        DEVICE_REQUIRE(pid<size());
        return d_states[pid].next_ijk;
    }

    // Next distance
    __device__ double& next_dist(int pid)             
    {
        DEVICE_REQUIRE(pid<size());
        return d_states[pid].next_dist;
    }
    __device__ const double& next_dist(int pid) const 
    {
        DEVICE_REQUIRE(pid<size());
        return d_states[pid].next_dist;
    }

    // Exiting face
    __device__ int& exiting_face(int pid)       
    {
        DEVICE_REQUIRE(pid<size());
        return d_states[pid].exiting_face;
    }
    __device__ const int& exiting_face(int pid) const 
    {
        DEVICE_REQUIRE(pid<size());
        return d_states[pid].exiting_face;
    }

    // Reflecting face
    __device__ int& reflecting_face(int pid)       
    {
        DEVICE_REQUIRE(pid<size());
        return d_states[pid].reflecting_face;
    }
    __device__ const int& reflecting_face(int pid) const 
    {
        DEVICE_REQUIRE(pid<size());
        return d_states[pid].reflecting_face;
    }

};

//===========================================================================//
/*!
 * \class Mesh_State_Vector_AOS_DMM
 * \brief Device memory manager for Mesh_State_Vector_AOS.
 */
//===========================================================================//

class Mesh_State_Vector_AOS_DMM :
    public cuda::Device_Memory_Manager<Mesh_State_Vector_AOS>
{
  public:

    typedef Mesh_State_Vector_AOS Geo_State_Vector_t;

    // Constructor
    Mesh_State_Vector_AOS_DMM(){}

    // Number of states currently allocated
    int size() const {return d_states.size();}

    // Memory manager interface
    Geo_State_Vector_t device_instance()
    {
        return Geo_State_Vector_t(cuda::make_view(d_states));
    }

    // Initialize vectors for specified number of states
    void initialize(int num_states)
    {
        REQUIRE(num_states > 0);
        d_states.resize(num_states);
    }

  private:

    using State_Vector = thrust::device_vector<Mesh_State>;

    // >>> DATA
    State_Vector d_states;

};


//---------------------------------------------------------------------------//
} // end namespace cuda_profugus

//---------------------------------------------------------------------------//
#endif // MC_cuda_geometry_Mesh_State_Vector_AOS_cuh

//---------------------------------------------------------------------------//
// end of MC/cuda_geometry/Mesh_State_Vector_AOS.cuh
//---------------------------------------------------------------------------//
