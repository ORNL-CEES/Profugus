//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   MC/mc/Group_Bounds.hh
 * \author Thomas M. Evans
 * \date   Wed Apr 30 14:05:49 2014
 * \brief  Group_Bounds class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef MC_mc_Group_Bounds_hh
#define MC_mc_Group_Bounds_hh

#include <vector>
#include <memory>

namespace profugus
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
    typedef std::shared_ptr<Group_Bounds> SP_Group_Bounds;
    //@}

  private:
    // >>> DATA

    // Bounds.
    Vec_Dbl d_bounds;

  public:

    // Factory method for building a logarithmic group bounds structure.
    static SP_Group_Bounds build_logarithmic(double lower, double higher,
                                             int num_bins);

  public:
    // Construct with just neutron and photon bounds for now.
    Group_Bounds(const Vec_Dbl &bounds);

    // >>> ACCESSORS

    //! Access group bounds for a single particle.
    const Vec_Dbl& group_bounds() const { return d_bounds; }

    //! Get the total number of groups.
    unsigned int num_groups() const { return d_bounds.size() - 1; }

    // Access upper and lower energy bound for a flattened group.
    void get_energy(int group_index, double &lower, double &upper) const;

    // Get the linearized group index for a particle at a given energy
    bool find(const double energy, int &group_index) const;
};

} // end namespace profugus

#endif // MC_mc_Group_Bounds_hh

//---------------------------------------------------------------------------//
//                 end of Group_Bounds.hh
//---------------------------------------------------------------------------//
