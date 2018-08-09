//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   MC/mc/Assembly_Model.hh
 * \author Steven Hamilton
 * \date   Wed Aug 08 16:13:28 2018
 * \brief  Assembly_Model class declaration.
 * \note   Copyright (c) 2018 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef MC_mc_Assembly_Model_hh
#define MC_mc_Assembly_Model_hh

namespace mc
{

//===========================================================================//
/*!
 * \class Assembly_Model
 * \brief High-level description of fuel assembly for surrogate models
 */
/*!
 * \example mc/test/tstAssembly_Model.cc
 *
 * Test of Assembly_Model.
 */
//===========================================================================//

class Assembly_Model
{
  public:
    //@{
    //! Typedefs
    <++>
    //@}

    enum PIN_TYPE {FUEL, GUIDE};

  private:
    // >>> DATA
    int d_Nx;
    int d_Ny;
    std::vector<PIN_TYPE> d_pin_map;
    std::vector<double> d_x_edges;
    std::vector<double> d_y_edges;
    double d_height;

    // Cylinder radii
    double d_fuel_radius;
    double d_clad_radius;
    double d_guide_radius;

  public:

    // Constructor
    Assembly_Model(const std::vector<PIN_TYPE>& pin_map,
                   const std::vector<double>&   x_edges,
                   const std::vector<double>&   y_edges,
                   double                       height);

    // Set fuel pin radius
    void set_fuel_radius(double fr)
    {
        REQUIRE(fr > 0);
        d_fuel_radius = fr;
    }

    // Set clad outer radius
    void set_clad_radius(double cr)
    {
        REQUIRE(cr > 0);
        d_clad_radius = cr;
    }

    // Set guide tube outer radius
    void set_guide_radius(double gr)
    {
        REQUIRE(gr > 0);
        d_guidel_radius = gr;
    }

    // Accessors
    double fuel_radius() const {return d_fuel_radius;}
    double clad_radius() const {return d_clad_radius;}
    double guide_radius() const {return d_guide_radius;}
    const std::vector<double>& x_edges() const {return d_x_edges;}
    const std::vector<double>& y_edges() const {return d_y_edges;}

    // Convert (i,j) to cardinal pin index
    int pin_id(int i, int j) const
    {
        REQUIRE(i < d_Nx);
        REQUIRE(j < d_Ny);
        return i + d_Nx * j;
    }

    // Type of pin at (i,j) location
    PIN_TYPE pin_type(int i, int j) const
    {
        REQUIRE(i < d_Nx);
        REQUIRE(j < d_Ny);
        return d_pin_map[pin_id(i,j)];
    }

    // Assembly height
    double height() const {return d_height;}

    // Flow area (cm^2)
    double flow_area(int i, int j) const;
};

//---------------------------------------------------------------------------//
} // end namespace mc

//---------------------------------------------------------------------------//
// INLINE DEFINITIONS
//---------------------------------------------------------------------------//
// #include "Assembly_Model.i.hh"
//---------------------------------------------------------------------------//
#endif // MC_mc_Assembly_Model_hh

//---------------------------------------------------------------------------//
// end of MC/mc/Assembly_Model.hh
//---------------------------------------------------------------------------//
