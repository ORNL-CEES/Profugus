//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/Group_Bounds.hh
 * \author Stuart Slattery
 * \brief  Group_Bounds class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_mc_Group_Bounds_hh
#define cuda_mc_Group_Bounds_hh

#include <vector>

#include <cuda_utils/CudaDBC.hh>
#include <cuda_utils/CudaMacros.hh>
#include <cuda_utils/Utility_Functions.hh>

namespace cuda_profugus
{

//===========================================================================//
/*!
 * \class Group_Bounds
 * \brief Group bounds container for tallies and cross section sampling.
 *
 * The group bounds are defined as follows
 *
 * \verbatim

  bound        |10.  5.   2.   1.|
  group #      |   0    1    2   |

  \endverbatim
 */
/*!
 * \example mc/test/tstGroup_Bounds.cc
 *
 * Test of Group_Bounds.
 */
//===========================================================================//

class Group_Bounds
{
  public:
    //@{
    //! Typedefs.
    typedef std::vector<double>           Vec_Dbl;
    //@}

  private:
    // >>> DATA

    // Bounds. On device.
    double* d_bounds;

    // Size of bounds.
    int d_size;

  public:

    // Construct with just neutron and photon bounds for now.
    Group_Bounds(const Vec_Dbl &bounds);

    // Destructor.
    ~Group_Bounds();

    // >>> ACCESSORS

    //! Access group bounds for a single particle.
    PROFUGUS_DEVICE_FUNCTION
    const double* group_bounds() const { return d_bounds; }

    //! Get the total number of groups.
    PROFUGUS_HOST_DEVICE_FUNCTION
    unsigned int num_groups() const { return d_size - 1; }

    // Access upper and lower energy bound for a flattened group.
    PROFUGUS_DEVICE_FUNCTION
    void get_energy(int group_index, double &lower, double &upper) const
    {
	DEVICE_REQUIRE(group_index + 1 < d_size );
	upper = d_bounds[group_index];
	lower = d_bounds[group_index + 1];
	DEVICE_ENSURE(lower < upper);
    }

    // Get the linearized group index for a particle at a given energy
    PROFUGUS_DEVICE_FUNCTION
    bool find(const double energy, int &group_index) const
    {
	DEVICE_REQUIRE(energy >= 0.);

	if ((energy > d_bounds[0]) || (energy < d_bounds[d_size-1]))
	    return false;

	// Find the group index; use std::greater because it's in descending
	// order
	group_index =
	    cuda_utils::utility::lower_bound( d_bounds, (d_bounds + d_size), energy,
					cuda_utils::utility::greater_than<double>() )
	    - d_bounds - 1;

	if (group_index == -1)
	    ++group_index;

	DEVICE_CHECK(group_index >= 0);
	DEVICE_CHECK(group_index < d_size - 1);
	return true;
    }
};

//---------------------------------------------------------------------------//

} // end namespace cuda_profugus

#endif // cuda_mc_Group_Bounds_hh

//---------------------------------------------------------------------------//
//                 end of Group_Bounds.hh
//---------------------------------------------------------------------------//
