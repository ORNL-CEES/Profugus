//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   MC/cuda_rtk/RTK_State_Vector_AOS.cuh
 * \author Steven Hamilton
 * \date   Tue Apr 11 11:22:34 2017
 * \brief  RTK_State_Vector_AOS class declaration.
 * \note   Copyright (c) 2017 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef MC_cuda_rtk_RTK_State_Vector_AOS_cuh
#define MC_cuda_rtk_RTK_State_Vector_AOS_cuh

#include "CudaUtils/cuda_utils/Definitions.hh"
#include "CudaUtils/cuda_utils/Device_View_Field.hh"

namespace cuda_profugus
{

//===========================================================================//
/*!
 * \class RTK_State_Vector_AOS
 * \brief Class for managing a collection of mesh states.
 */
//===========================================================================//

class RTK_State_Vector_AOS
{
  public:
    //@{
    //! Typedefs
    typedef cuda_utils::Coordinates                 Coordinates;
    typedef cuda_utils::Space_Vector                Space_Vector;
    typedef cuda::Device_View_Field<RTK_State>      State_View;
    //@}

  private:

    // >>> DATA
    State_View d_states;

  public:

    // Constructor
    RTK_State_Vector_AOS(State_View states)
      : d_states(states)
    {
    }

    // Number of states allocated
    __device__ int size() const {return d_states.size();}

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

    // Region
    __device__ int & region(int pid)       
    {
        DEVICE_REQUIRE(pid<size());
        return d_states[pid].region;
    }
    __device__ const int & region(int pid) const 
    {
        DEVICE_REQUIRE(pid<size());
        return d_states[pid].region;
    }

    // Segment
    __device__ int & segment(int pid)       
    {
        DEVICE_REQUIRE(pid<size());
        return d_states[pid].segment;
    }
    __device__ const int & segment(int pid) const 
    {
        DEVICE_REQUIRE(pid<size());
        return d_states[pid].segment;
    }

    // Face
    __device__ int & face(int pid)       
    {
        DEVICE_REQUIRE(pid<size());
        return d_states[pid].face;
    }
    __device__ const int & face(int pid) const 
    {
        DEVICE_REQUIRE(pid<size());
        return d_states[pid].face;
    }

    // Next egion
    __device__ int & next_region(int pid)       
    {
        DEVICE_REQUIRE(pid<size());
        return d_states[pid].next_region;
    }
    __device__ const int & next_region(int pid) const 
    {
        DEVICE_REQUIRE(pid<size());
        return d_states[pid].next_region;
    }

    // Next segment
    __device__ int & next_segment(int pid)       
    {
        DEVICE_REQUIRE(pid<size());
        return d_states[pid].next_segment;
    }
    __device__ const int & next_segment(int pid) const 
    {
        DEVICE_REQUIRE(pid<size());
        return d_states[pid].next_segment;
    }

    // Next face
    __device__ int & next_face(int pid)       
    {
        DEVICE_REQUIRE(pid<size());
        return d_states[pid].next_face;
    }
    __device__ const int & next_face(int pid) const 
    {
        DEVICE_REQUIRE(pid<size());
        return d_states[pid].next_face;
    }

    // Dist to next region
    __device__ double& dist_to_next_region(int pid)             
    {
        DEVICE_REQUIRE(pid<size());
        return d_states[pid].dist_to_next_region;
    }
    __device__ const double& dist_to_next_region(int pid) const 
    {
        DEVICE_REQUIRE(pid<size());
        return d_states[pid].dist_to_next_region;
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

    // Escaping face
    __device__ int& escaping_face(int pid)       
    {
        DEVICE_REQUIRE(pid<size());
        return d_states[pid].escaping_face;
    }
    __device__ const int& escaping_face(int pid) const 
    {
        DEVICE_REQUIRE(pid<size());
        return d_states[pid].escaping_face;
    }

    // Exiting level
    __device__ int& exiting_level(int pid, int level)       
    {
        DEVICE_REQUIRE(pid<size());
        DEVICE_REQUIRE(level<3);
        return d_states[pid].exiting_level[level];
    }
    __device__ const int& exiting_level(int pid, int level) const 
    {
        DEVICE_REQUIRE(pid<size());
        DEVICE_REQUIRE(level<3);
        return d_states[pid].exiting_level[level];
    }

    // Level coordinates
    __device__ Coordinates& level_coord(int pid, int level)       
    {
        DEVICE_REQUIRE(pid<size());
        DEVICE_REQUIRE(level<3);
        return d_states[pid].level_coord[level];
    }
    __device__ const Coordinates& level_coord(int pid, int level) const 
    {
        DEVICE_REQUIRE(pid<size());
        DEVICE_REQUIRE(level<3);
        return d_states[pid].level_coord[level];
    }

};

//===========================================================================//
/*!
 * \class RTK_State_Vector_AOS_DMM
 * \brief Device memory manager for RTK_State_Vector_AOS.
 */
//===========================================================================//

class RTK_State_Vector_AOS_DMM :
    public cuda::Device_Memory_Manager<RTK_State_Vector_AOS>
{
  public:

    typedef RTK_State_Vector_AOS Geo_State_Vector_t;

    //! Typedefs
    typedef thrust::device_vector<RTK_State>  State_Vector;

    // Constructor
    RTK_State_Vector_AOS_DMM(){}

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

    // >>> DATA
    State_Vector d_states;
};

//---------------------------------------------------------------------------//
} // end namespace cuda_profugus

//---------------------------------------------------------------------------//
#endif // MC_cuda_rtk_RTK_State_Vector_AOS_cuh

//---------------------------------------------------------------------------//
// end of MC/cuda_rtk/RTK_State_Vector_AOS.cuh
//---------------------------------------------------------------------------//
