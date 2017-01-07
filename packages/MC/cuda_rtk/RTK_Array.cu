//---------------------------------*-CUDA-*----------------------------------//
/*!
 * \file   MC/cuda_rtk/RTK_Array.cu
 * \author Tom Evans
 * \date   Wed Jan 04 15:43:43 2017
 * \brief  RTK_Array member and kernel definitions.
 * \note   Copyright (c) 2017 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <vector>

#include "RTK_Array.cuh"

#include "Utils/harness/DBC.hh"

namespace cuda_profugus
{

//---------------------------------------------------------------------------//
// RTK_ARRAY MEMBERS
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
template<class T>
RTK_Array<T>::RTK_Array(
    int          level,
    Dim_Vector   N,
    View_Int     layout,
    Object_View  objects,
    View_Dbl     x,
    View_Dbl     y,
    View_Dbl     z,
    Space_Vector corner,
    Space_Vector length)
    : d_N(std::move(N))
    , d_layout(std::move(layout))
    , d_objects(std::move(objects))
    , d_x(std::move(x))
    , d_y(std::move(y))
    , d_z(std::move(z))
    , d_corner(std::move(corner))
    , d_length(std::move(length))
    , d_level(level)
{
}

//---------------------------------------------------------------------------//
// RTK_ARRAY EXPLICIT INSTANTIATIONS
//---------------------------------------------------------------------------//

template class RTK_Array<RTK_Cell>;
template class RTK_Array< RTK_Array<RTK_Cell> >;

//---------------------------------------------------------------------------//
// RTK_LATTICE_ARRAY_DMM MEMBERS
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
RTK_Lattice_Array_DMM::RTK_Lattice_Array_DMM(
    const Host_Lattice_Array &host_array)
{
    using def::X; using def::Y; using def::Z;
    using def::I; using def::J; using def::K;

    // Get the number of pincells in this array
    auto num_pins = host_array.num_objects();

    // Make managers for the pincells
    DV::Managers managers(num_pins);

    // Build the managers
    for (int n = 0; n < num_pins; ++n)
    {
        // Get the pin on the host
        const auto &pin = host_array.object(n);

        // Make the device manager
        managers[n] = std::make_shared<RTK_Cell_DMM>(pin);
    }

    // Build the device view to the pins in this array
    d_objects = std::make_shared<DV>(managers);

    // Build the array edges
    const auto &x = host_array.edges(X);
    const auto &y = host_array.edges(Y);
    const auto &z = host_array.edges(Z);
    d_x = thrust::device_vector<double>(x.begin(), x.end());
    d_y = thrust::device_vector<double>(y.begin(), y.end());
    d_z = thrust::device_vector<double>(z.begin(), z.end());

    // Build the number of cells
    d_N = {host_array.size(X), host_array.size(Y), host_array.size(Z)};
    CHECK(num_pins == d_N[X] * d_N[Y] * d_N[Z]);

    // Get the corner and length
    Host_Lattice_Array::Space_Vector lo, hi;
    host_array.get_extents(lo, hi);

    d_corner = Space_Vector::from_host(lo);
    d_length = Space_Vector::from_host(hi - lo);

    // Build the layout
    std::vector<int> layout(num_pins, -1);
    for (int k = 0; k < d_N[K]; ++k)
    {
        for (int j = 0; j < d_N[J]; ++j)
        {
            for (int i = 0; i < d_N[I]; ++i)
            {
                CHECK(host_array.index(i, j, k) < num_pins);
                layout[host_array.index(i, j, k)] = host_array.id(i, j, k);
            }
        }
    }
    d_layout = thrust::device_vector<int>(layout.begin(), layout.end());
}

//---------------------------------------------------------------------------//
/*!
 * \brief Create a device instance.
 */
Lattice_Array RTK_Lattice_Array_DMM::device_instance()
{
    return Lattice_Array(
        0,
        d_N,
        cuda::make_view(d_layout),
        d_objects->get_view(),
        cuda::make_view(d_x),
        cuda::make_view(d_y),
        cuda::make_view(d_z),
        d_corner,
        d_length);
}

//---------------------------------------------------------------------------//
// RTK_CORE_ARRAY_DMM MEMBERS
//---------------------------------------------------------------------------//
#if 0
/*!
 * \brief Constructor.
 */
RTK_Core_Array_DMM::RTK_Core_Array_DMM(
    const Host_Core_Array &host_array)
{
}

//---------------------------------------------------------------------------//
/*!
 * \brief Create a device instance.
 */
Core_Array RTK_Core_Array_DMM::device_instance()
{
    return Core_Array();
}
#endif
} // end namespace cuda_profugus

//---------------------------------------------------------------------------//
// end of MC/cuda_rtk/RTK_Array.cu
//---------------------------------------------------------------------------//
