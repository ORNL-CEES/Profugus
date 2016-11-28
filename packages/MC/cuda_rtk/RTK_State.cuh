//---------------------------------*-CUDA-*----------------------------------//
/*!
 * \file   MC/cuda_rtk/RTK_State.cuh
 * \author Thomas Evans
 * \date   Fri Nov 18 15:06:45 2016
 * \brief  RTK_State kernel declarations.
 * \note   Copyright (c) 2016 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef MC_cuda_rtk_RTK_State_cuh
#define MC_cuda_rtk_RTK_State_cuh

#include "CudaUtils/cuda_utils/Device_Vector_Lite.hh"

namespace cuda_profugus
{

//===========================================================================//
/*!
 * \struct RTK_State
 * \brief CUDA implementation of handle to basic Reactor ToolKit pin-cell,
 * core geometry package
 */
//===========================================================================//

struct RTK_State
{
    //@{
    //! Types.
    using Space_Vector = Device_Vector_Lite<double, 3>;
    using Coordinates  = Device_Vector_Lite<int, 3>;
    //@}

    //! Position.
    Space_Vector d_r;

    //! Direction.
    Space_Vector d_dir;

    // >>> GEOMETRIC STATE

    //! Faces in pin-cell and vessel.
    enum Faces {MODERATOR = 0,
                NONE      = 900,
                INTERNAL  = 901,
                MINUS_X   = 1000,
                PLUS_X    = 1001,
                MINUS_Y   = 1002,
                PLUS_Y    = 1003,
                MINUS_Z   = 1004,
                PLUS_Z    = 1005,
                R0_VESSEL = 2000,
                R1_VESSEL = 2001,
                VESSEL    = 2002};

    //@{
    //! Plus/minus faces.
    __device__
    static int plus_face(int d)
    {
        if (d == 0)
            return PLUS_X;
        else if (d == 1)
            return PLUS_Y;
        else if (d == 2)
            return PLUS_Z;

        return -1;
    }

    __device__
    static int minus_face(int d)
    {
    if (d == 0)
        return MINUS_X;
    else if (d == 1)
        return MINUS_Y;
    else if (d == 2)
        return MINUS_Z;

    return -1;
    }
    //@}

    //@{
    //! Pin-cell semantics.
    int    region;
    int    segment;
    int    face;
    int    next_region;
    int    next_segment;
    int    next_face;
    double dist_to_next_region;
    //@}

    //! Exiting face indicator.
    int exiting_face;

    // >>> Max levels supported = 3

    //! Coordinates in array at each level.
    Device_Vector_Lite<Coordinates, 3> level_coord;

    //! Crossing boundary indicator by level.
    Device_Vector_Lite<int, 3> exiting_level;

    //! Escaping face in geometry.
    int escaping_face;

    //! Reflecting face in geometry.
    int reflecting_face;
};

} // end namespace cuda_profugus

#endif // MC_cuda_rtk_RTK_State_cuh

//---------------------------------------------------------------------------//
// end of MC/cuda_rtk/RTK_State.cuh
//---------------------------------------------------------------------------//
