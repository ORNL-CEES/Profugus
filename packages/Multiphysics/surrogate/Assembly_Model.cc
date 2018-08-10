//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   MC/mc/Assembly_Model.cc
 * \author Steven Hamilton
 * \date   Wed Aug 08 16:13:28 2018
 * \brief  Assembly_Model class definitions.
 * \note   Copyright (c) 2018 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Assembly_Model.hh"

#include "Utils/utils/Constants.hh"

namespace mc
{
//---------------------------------------------------------------------------//
// Constructor
//---------------------------------------------------------------------------//
Assembly_Model::Assembly_Model(const std::vector<PIN_TYPE>& pin_map,
                               const std::vector<double>&   x_edges,
                               const std::vector<double>&   y_edges,
                               double                       height)
    : d_pin_map(pin_map)
    , d_x_edges(x_edges)
    , d_y_edges(y_edges)
    , d_height(height)
{
    d_Nx = d_x_edges.size() - 1;
    d_Ny = d_y_edges.size() - 1;
    REQUIRE(d_pin_map.size() == d_Nx * d_Ny);
}

//---------------------------------------------------------------------------//
// Constructor
//---------------------------------------------------------------------------//
double Assembly_Model::flow_area(int i, int j) const
{
    REQUIRE(i < d_Nx);
    REQUIRE(j < d_Ny);

    double dx = d_x_edges[i+1] - d_x_edges[i];
    double dy = d_y_edges[j+1] - d_y_edges[j];

    auto type = pin_type(i,j);
    double r = 0.0;
    if (type == FUEL)
        r = d_clad_radius;
    else if (type == GUIDE)
        r = d_guide_radius;
    CHECK(r > 0.0);

    double area = dx * dy - profugus::constants::pi * r * r;
    ENSURE(area > 0);
    return area;
}

//---------------------------------------------------------------------------//
} // end namespace mc

//---------------------------------------------------------------------------//
// end of MC/mc/Assembly_Model.cc
//---------------------------------------------------------------------------//
