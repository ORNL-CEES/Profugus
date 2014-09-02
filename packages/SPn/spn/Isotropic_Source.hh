//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   spn/Isotropic_Source.hh
 * \author Thomas M. Evans
 * \date   Fri Feb 14 18:14:05 2014
 * \brief  Isotropic_Source class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef spn_Isotropic_Source_hh
#define spn_Isotropic_Source_hh

#include <string>
#include <vector>
#include "harness/DBC.hh"
#include "utils/Definitions.hh"

namespace profugus
{

//===========================================================================//
/*!
 * \class Isotropic_Source
 * \brief Defines an isotropic, external source in each cell.
 */
/*!
 * \example spn/test/tstIsotropic_Source.cc
 *
 * Test of Isotropic_Source.
 */
//===========================================================================//

class Isotropic_Source
{
  public:
    //@{
    //! Typedefs.
    typedef def::Vec_Dbl       Source_Field;
    typedef def::Vec_Int       ID_Field;
    typedef def::Vec_Dbl       Shape;
    typedef std::vector<Shape> Source_Shapes;
    //@}

  private:
    // >>> DATA

    // Source shape ids in each cell.
    ID_Field d_ids;

    // Source shapes.
    Source_Shapes d_shapes;

    // Source strength in particles/cm^3 in each cell.
    Source_Field d_source;

  public:
    // Constructor.
    explicit Isotropic_Source(int num_cells);

    // >>> DERIVED INTERFACE

    //! Return the number of groups.
    int num_groups() const { return d_shapes.empty() ? 0 : d_shapes[0].size(); }

    //! Get source in particles/cm^3/str.
    double q_e(int cell, int g) const
    {
        REQUIRE(!d_shapes.empty());

        REQUIRE(cell >= 0 && cell < d_ids.size());
        REQUIRE(d_ids[cell] >= 0 && d_ids[cell] < d_shapes.size());
        REQUIRE(g >= 0 && g < d_shapes[d_ids[cell]].size());

        return d_source[cell] * d_shapes[d_ids[cell]][g];
    }

    // Truncate the source into a subset of groups.
    void truncate(int g_first, int g_last);

    // Verify that the sources have been assigned.
    bool verify() const;

    //! Source Label
    std::string label() const { return "isotropic"; }

    // Whether we contain only isotropic sources
    bool is_isotropic() const { return true; }

    // >>> CLASS INTERFACE

    // Set the source
    void set(const ID_Field &ids, const Source_Shapes &shapes,
             const Source_Field &source);

    //! Number of source spectral shapes.
    int num_shapes() const { return d_shapes.size(); }
};

} // end namespace profugus

#endif // spn_Isotropic_Source_hh

//---------------------------------------------------------------------------//
//                 end of Isotropic_Source.hh
//---------------------------------------------------------------------------//
