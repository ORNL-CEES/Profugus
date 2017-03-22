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
// TRAITS FOR RTK DMM MANAGERS
//---------------------------------------------------------------------------//

namespace
{

//---------------------------------------------------------------------------//

template<class Device_T>
class RTK_DMM_Traits
{
};

//---------------------------------------------------------------------------//
// Specialization on RTK_Cell

template<>
class RTK_DMM_Traits<RTK_Cell>
{
  public:
    using DMM = RTK_Cell_DMM;
};

//---------------------------------------------------------------------------//
// Specialization on RTK_Lattice

template<>
class RTK_DMM_Traits< RTK_Array<RTK_Cell> >
{
  public:
    using DMM = RTK_Array_DMM<
      profugus::RTK_Array<profugus::RTK_Cell>,
      RTK_Array<RTK_Cell> >;
};

} // end unnamed namespace

//---------------------------------------------------------------------------//
// RTK_ARRAY MEMBERS
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
template<class T>
RTK_Array<T>::RTK_Array(
    int              level,
    Dim_Vector       N,
    View_Int         layout,
    Object_View      objects,
    View_Dbl         x,
    View_Dbl         y,
    View_Dbl         z,
    View_Int         num_cells,
    View_Int         cell_offsets,
    Space_Vector     corner,
    Space_Vector     length,
    Reflecting_Faces reflect)
    : d_N(std::move(N))
    , d_layout(std::move(layout))
    , d_objects(std::move(objects))
    , d_x(std::move(x))
    , d_y(std::move(y))
    , d_z(std::move(z))
    , d_num_cells(std::move(num_cells))
    , d_cell_offsets(std::move(cell_offsets))
    , d_corner(std::move(corner))
    , d_length(std::move(length))
    , d_level(level)
    , d_reflect(reflect)
{
}

//---------------------------------------------------------------------------//
// RTK_ARRAY EXPLICIT INSTANTIATIONS
//---------------------------------------------------------------------------//

template class RTK_Array<RTK_Cell>;
template class RTK_Array< RTK_Array<RTK_Cell> >;

//---------------------------------------------------------------------------//
// RTK_ARRAY_DMM_MEMBERS
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
template<class Host_Array_T, class Device_Array_T>
RTK_Array_DMM<Host_Array_T,Device_Array_T>::RTK_Array_DMM(
    const Host_Array_T &host_array)
{
    using Object_DMM = typename RTK_DMM_Traits<Object_t>::DMM;
    using def::X; using def::Y; using def::Z;
    using def::I; using def::J; using def::K;

    // Get the number of objects in this array
    auto num_objects = host_array.num_objects();

    // Make map from host object index to device object index
    std::vector<int> h2d(num_objects, -1);

    // Device object index
    int devidx = 0;

    // Make managers for the objects
    typename DV::Managers managers;

    // Get the host-to-device layout index mapping and build the managers
    for (int n = 0; n < num_objects; ++n)
    {
        if (host_array.has_object(n))
        {
            // Make the mapping
            h2d[n] = devidx;

            // Get the object
            const auto &host_object = host_array.object(n);

            // Make the device manager
            managers.push_back(std::make_shared<Object_DMM>(host_object));

            // Advance the device index
            ++devidx;
        }
    }
    CHECK(devidx == managers.size());

    // Build the device view to the pins in this array
    d_objects = std::make_shared<DV>(managers);

    // Build the array edges
    const auto &x = host_array.edges(X);
    const auto &y = host_array.edges(Y);
    const auto &z = host_array.edges(Z);
    d_x = thrust::device_vector<double>(x.begin(), x.end());
    d_y = thrust::device_vector<double>(y.begin(), y.end());
    d_z = thrust::device_vector<double>(z.begin(), z.end());

    // Build the number of cells in each axis
    d_N = {host_array.size(X), host_array.size(Y), host_array.size(Z)};
    CHECK(host_array.size() == d_N[X] * d_N[Y] * d_N[Z]);

    // Get the corner and length
    typename Host_Array_T::Space_Vector lo, hi;
    host_array.get_extents(lo, hi);

    d_corner = Space_Vector::from_host(lo);
    d_length = Space_Vector::from_host(hi - lo);

    // Build the layout
    std::vector<int> layout(host_array.size(), -1);
    for (int k = 0; k < d_N[K]; ++k)
    {
        for (int j = 0; j < d_N[J]; ++j)
        {
            for (int i = 0; i < d_N[I]; ++i)
            {
                CHECK(host_array.id(i, j, k) < num_objects);

                // Get the device index of the object
                devidx = h2d[host_array.id(i,j,k)];
                CHECK(devidx >= 0);
                CHECK(devidx < d_objects->size());

                // Assign the layout
                layout[host_array.index(i, j, k)] = devidx;
            }
        }
    }
    d_layout = thrust::device_vector<int>(layout.begin(), layout.end());

    // Set the level
    d_level = host_array.level();

    // Store the reflecting faces.
    const auto &r = host_array.reflecting_faces();
    d_reflect     = {r[0], r[1], r[2], r[3], r[4], r[5]};

    // Compute cell counts and offset
    std::vector<int> host_cells;
    for (int iz = 0; iz < host_array.size(K); ++iz)
    {
        for (int iy = 0; iy < host_array.size(J); ++iy)
        {
            for (int ix = 0; ix < host_array.size(I); ++ix)
            {
                const auto& obj = host_array.object(ix,iy,iz);
                host_cells.push_back(obj.num_cells());
            }
        }
    }

    // Create cell offsets from cell counts
    std::vector<int> host_offsets(host_cells.size(),0);
    std::partial_sum(host_cells.begin(),host_cells.end()-1,
                     host_offsets.begin()+1);

    // Copy to device
    d_num_cells    = host_cells;
    d_cell_offsets = host_offsets;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Create a device instance.
 */
template<class Host_Array_T, class Device_Array_T>
Device_Array_T RTK_Array_DMM<Host_Array_T,Device_Array_T>::device_instance()
{
    return Device_Array_T(
        d_level,
        d_N,
        cuda::make_view(d_layout),
        d_objects->get_view(),
        cuda::make_view(d_x),
        cuda::make_view(d_y),
        cuda::make_view(d_z),
        cuda::make_view(d_num_cells),
        cuda::make_view(d_cell_offsets),
        d_corner,
        d_length,
        d_reflect);
}

//---------------------------------------------------------------------------//
// RTK_ARRAY_DMM EXPLICIT INSTANTIATIONS
//---------------------------------------------------------------------------//

template class RTK_Array_DMM<
    profugus::RTK_Array<profugus::RTK_Cell>,
    RTK_Array<RTK_Cell> >;
template class RTK_Array_DMM<
    profugus::RTK_Array< profugus::RTK_Array<profugus::RTK_Cell> >,
    RTK_Array< RTK_Array<RTK_Cell> > >;

} // end namespace cuda_profugus

//---------------------------------------------------------------------------//
// end of MC/cuda_rtk/RTK_Array.cu
//---------------------------------------------------------------------------//
