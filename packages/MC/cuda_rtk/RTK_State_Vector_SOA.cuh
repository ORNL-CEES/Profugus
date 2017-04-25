//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   MC/cuda_rtk/RTK_State_Vector_SOA.cuh
 * \author Steven Hamilton
 * \date   Tue Apr 11 11:22:34 2017
 * \brief  RTK_State_Vector_SOA class declaration.
 * \note   Copyright (c) 2017 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef MC_cuda_rtk_RTK_State_Vector_SOA_cuh
#define MC_cuda_rtk_RTK_State_Vector_SOA_cuh

#include "CudaUtils/cuda_utils/Definitions.hh"
#include "CudaUtils/cuda_utils/Device_View_Field.hh"

namespace cuda_profugus
{

//===========================================================================//
/*!
 * \class RTK_State_Vector_SOA
 * \brief Class for managing a collection of mesh states.
 */
//===========================================================================//

class RTK_State_Vector_SOA
{
  public:
    //@{
    //! Typedefs
    typedef cuda_utils::Coordinates                 Coordinates;
    typedef cuda_utils::Space_Vector                Space_Vector;
    typedef cuda::Device_View_Field<int>            Int_View;
    typedef cuda::Device_View_Field<double>         Double_View;
    typedef cuda::Device_View_Field<Coordinates>    Coordinate_View;
    typedef cuda::Device_View_Field<Space_Vector>   Space_Vector_View;
    //@}

  private:

    // >>> DATA
    Space_Vector_View d_r;
    Space_Vector_View d_dir;
    Int_View          d_exiting_face;
    Int_View          d_escaping_face;
    Int_View          d_reflecting_face;
    Coordinate_View   d_level_coord;   // Sized 3x state vector length
    Int_View          d_exiting_level; // Sized 3x state vector length

    // Pincell stuff
    Int_View    d_region;
    Int_View    d_segment;
    Int_View    d_face;
    Int_View    d_next_region;
    Int_View    d_next_segment;
    Int_View    d_next_face;
    Double_View d_dist_to_next_region;

  public:

    // Constructor
    RTK_State_Vector_SOA(
        Space_Vector_View   r,
        Space_Vector_View   dir,
        Int_View            exiting_face,
        Int_View            escaping_face,
        Int_View            reflecting_face,
        Coordinate_View     level_coord,
        Int_View            exiting_level,
        Int_View            region,
        Int_View            segment,
        Int_View            face,
        Int_View            next_region,
        Int_View            next_segment,
        Int_View            next_face,
        Double_View         dist_to_next_region)
      : d_r(r)
      , d_dir(dir)
      , d_exiting_face(exiting_face)
      , d_escaping_face(escaping_face)
      , d_reflecting_face(reflecting_face)
      , d_level_coord(level_coord)
      , d_exiting_level(exiting_level)
      , d_region(region)
      , d_segment(segment)
      , d_face(face)
      , d_next_region(next_region)
      , d_next_segment(next_segment)
      , d_next_face(next_face)
      , d_dist_to_next_region(dist_to_next_region)
    {
    }

    // Number of states allocated
    __device__ int size() const {return d_r.size();}

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

    // Region
    __device__ int & region(int pid)       
    {
        DEVICE_REQUIRE(pid<size());
        return d_region[pid];
    }
    __device__ const int & region(int pid) const 
    {
        DEVICE_REQUIRE(pid<size());
        return d_region[pid];
    }

    // Segment
    __device__ int & segment(int pid)       
    {
        DEVICE_REQUIRE(pid<size());
        return d_segment[pid];
    }
    __device__ const int & segment(int pid) const 
    {
        DEVICE_REQUIRE(pid<size());
        return d_segment[pid];
    }

    // Face
    __device__ int & face(int pid)       
    {
        DEVICE_REQUIRE(pid<size());
        return d_face[pid];
    }
    __device__ const int & face(int pid) const 
    {
        DEVICE_REQUIRE(pid<size());
        return d_face[pid];
    }

    // Next egion
    __device__ int & next_region(int pid)       
    {
        DEVICE_REQUIRE(pid<size());
        return d_next_region[pid];
    }
    __device__ const int & next_region(int pid) const 
    {
        DEVICE_REQUIRE(pid<size());
        return d_next_region[pid];
    }

    // Next segment
    __device__ int & next_segment(int pid)       
    {
        DEVICE_REQUIRE(pid<size());
        return d_next_segment[pid];
    }
    __device__ const int & next_segment(int pid) const 
    {
        DEVICE_REQUIRE(pid<size());
        return d_next_segment[pid];
    }

    // Next face
    __device__ int & next_face(int pid)       
    {
        DEVICE_REQUIRE(pid<size());
        return d_next_face[pid];
    }
    __device__ const int & next_face(int pid) const 
    {
        DEVICE_REQUIRE(pid<size());
        return d_next_face[pid];
    }

    // Dist to next region
    __device__ double& dist_to_next_region(int pid)             
    {
        DEVICE_REQUIRE(pid<size());
        return d_dist_to_next_region[pid];
    }
    __device__ const double& dist_to_next_region(int pid) const 
    {
        DEVICE_REQUIRE(pid<size());
        return d_dist_to_next_region[pid];
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

    // Escaping face
    __device__ int& escaping_face(int pid)       
    {
        DEVICE_REQUIRE(pid<size());
        return d_escaping_face[pid];
    }
    __device__ const int& escaping_face(int pid) const 
    {
        DEVICE_REQUIRE(pid<size());
        return d_escaping_face[pid];
    }

    // Exiting level
    __device__ int& exiting_level(int pid, int level)       
    {
        DEVICE_REQUIRE(pid<size());
        DEVICE_REQUIRE(level<3);
        return d_exiting_level[pid + level*size()];
    }
    __device__ const int& exiting_level(int pid, int level) const 
    {
        DEVICE_REQUIRE(pid<size());
        DEVICE_REQUIRE(level<3);
        return d_exiting_level[pid + level*size()];
    }

    // Level coordinates
    __device__ Coordinates& level_coord(int pid, int level)       
    {
        DEVICE_REQUIRE(pid<size());
        DEVICE_REQUIRE(level<3);
        return d_level_coord[level*size() + pid];
    }
    __device__ const Coordinates& level_coord(int pid, int level) const 
    {
        DEVICE_REQUIRE(pid<size());
        DEVICE_REQUIRE(level<3);
        return d_level_coord[level*size() + pid];
    }

};

//===========================================================================//
/*!
 * \class RTK_State_Vector_SOA_DMM
 * \brief Device memory manager for RTK_State_Vector_SOA.
 */
//===========================================================================//

class RTK_State_Vector_SOA_DMM :
    public cuda::Device_Memory_Manager<RTK_State_Vector_SOA>
{
  public:

    typedef RTK_State_Vector_SOA Geo_State_Vector_t;

    //! Typedefs
    typedef cuda_utils::Coordinates             Coordinates;
    typedef cuda_utils::Space_Vector            Space_Vector;
    typedef thrust::device_vector<int>          Int_Vector;
    typedef thrust::device_vector<double>       Double_Vector;
    typedef thrust::device_vector<Coordinates>  Coordinate_Vector;
    typedef thrust::device_vector<Space_Vector> Space_Vector_Vector;

    // Constructor
    RTK_State_Vector_SOA_DMM(){}

    // Number of states currently allocated
    int size() const {return d_r.size();}

    // Memory manager interface
    Geo_State_Vector_t device_instance()
    {
        return Geo_State_Vector_t(cuda::make_view(d_r),
                                  cuda::make_view(d_dir),
                                  cuda::make_view(d_exiting_face),
                                  cuda::make_view(d_escaping_face),
                                  cuda::make_view(d_reflecting_face),
                                  cuda::make_view(d_level_coord),
                                  cuda::make_view(d_exiting_level),
                                  cuda::make_view(d_region),
                                  cuda::make_view(d_segment),
                                  cuda::make_view(d_face),
                                  cuda::make_view(d_next_region),
                                  cuda::make_view(d_next_segment),
                                  cuda::make_view(d_next_face),
                                  cuda::make_view(d_dist_to_next_region));
    }

    // Initialize vectors for specified number of states
    void initialize(int num_states)
    {
        REQUIRE(num_states > 0);
        d_r.resize(num_states);
        d_dir.resize(num_states);
        d_exiting_face.resize(num_states);
        d_escaping_face.resize(num_states);
        d_reflecting_face.resize(num_states);
        d_level_coord.resize(3*num_states);
        d_exiting_level.resize(3*num_states);
        d_region.resize(num_states);
        d_segment.resize(num_states);
        d_face.resize(num_states);
        d_next_region.resize(num_states);
        d_next_segment.resize(num_states);
        d_next_face.resize(num_states);
        d_dist_to_next_region.resize(num_states);
    }

  private:

    // >>> DATA
    Space_Vector_Vector d_r;
    Space_Vector_Vector d_dir;
    Int_Vector          d_exiting_face;
    Int_Vector          d_escaping_face;
    Int_Vector          d_reflecting_face;
    Coordinate_Vector   d_level_coord;   // Sized 3x state vector length
    Int_Vector          d_exiting_level; // Sized 3x state vector length

    // Pincell stuff
    Int_Vector    d_region;
    Int_Vector    d_segment;
    Int_Vector    d_face;
    Int_Vector    d_next_region;
    Int_Vector    d_next_segment;
    Int_Vector    d_next_face;
    Double_Vector d_dist_to_next_region;
};

//---------------------------------------------------------------------------//
} // end namespace cuda_profugus

//---------------------------------------------------------------------------//
#endif // MC_cuda_rtk_RTK_State_Vector_SOA_cuh

//---------------------------------------------------------------------------//
// end of MC/cuda_rtk/RTK_State_Vector_SOA.cuh
//---------------------------------------------------------------------------//
