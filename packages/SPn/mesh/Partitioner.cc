//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mesh/Partitioner.cc
 * \author Thomas M. Evans
 * \date   Wed Feb 12 09:54:53 2014
 * \brief  Partitioner member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "harness/Soft_Equivalence.hh"
#include "utils/Container_Functions.hh"
#include "Partitioner.hh"

namespace profugus
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 *
 * \param pl
 */
Partitioner::Partitioner(RCP_ParameterList pl)
{
    using def::I; using def::J; using def::K;

    Require (!pl.is_null());

    // initialize the parameter list
    init_pl(pl);

    // set dimension
    d_dimension = pl->get<int>("dimension");
    Insist(d_dimension == 2 || d_dimension == 3, "Dimension must be 2 or 3");

    // initialize blocks
    d_Nb[I]      = pl->get<int>("num_blocks_i");
    d_Nb[J]      = pl->get<int>("num_blocks_j");
    d_k_blocks   = pl->get<int>("num_z_blocks");
    d_num_blocks = d_Nb[I] * d_Nb[J];

    // build global edges of mesh
    Check (d_edges[I].empty() && d_edges[J].empty() && d_edges[K].empty());

    // determine whether the mesh grid is set with uniform spacings,
    // by supplying the cell-edges explicitly for each direction,
    // or nothing (allow for derivative classes to set d_edges)

    int num_cells = 0;
    double delta  = 0.0;

    if (pl->isParameter("num_cells_i"))
    {
        num_cells = pl->get<int>("num_cells_i");
        delta     = pl->get<double>("delta_x");
        Validate(num_cells > 0, "num_cells_i must be postitive");
        Validate(delta > 0., "delta_x must be postitive");

        build_uniform_edges(num_cells, delta, d_edges[I]);
    }
    else
    {
        d_edges[I] = pl->get<def::Vec_Dbl>("x_edges");
        Validate(profugus::is_sorted(d_edges[I].begin(), d_edges[I].end()),
                 "Mesh edges along X axis are not monotonically increasing.");
    }

    if (pl->isParameter("num_cells_j"))
    {
        num_cells = pl->get<int>("num_cells_j");
        delta     = pl->get<double>("delta_y");
        Validate(num_cells > 0, "num_cells_j must be postitive");
        Validate(delta > 0., "delta_y must be postitive");

        build_uniform_edges(num_cells, delta, d_edges[J]);
    }
    else
    {
        d_edges[J] = pl->get<def::Vec_Dbl>("y_edges");
        Validate(profugus::is_sorted(d_edges[J].begin(), d_edges[J].end()),
                 "Mesh edges along Y axis are not monotonically increasing.");
    }

    // process K if 3-dimensional
    if (d_dimension == 3)
    {
        if (pl->isParameter("num_cells_k"))
        {
            num_cells = pl->get<int>("num_cells_k");
            delta     = pl->get<double>("delta_z");
            Insist(num_cells > 0, "num_cells_k must be postitive");
            Insist(delta > 0., "delta_z must be postitive");

            // this->build_uniform_edges(num_cells, delta, d_edges[K]);
        }
        else
        {
            d_edges[K] = pl->get<def::Vec_Dbl>("z_edges");
            Validate(profugus::is_sorted(d_edges[K].begin(), d_edges[K].end()),
                     "Mesh edges along Z axis are not monotonically "
                     "increasing.");
        }
    }


    Ensure (d_Nb[I] > 0);
    Ensure (d_Nb[J] > 0);
    Ensure (d_dimension > 0);
    Ensure(d_edges[I].size() > 1);
    Ensure(d_edges[J].size() > 1);
    Ensure(dimension() == 3 ? d_edges[K].size() > 1 : d_edges[K].empty());
}

//---------------------------------------------------------------------------//
// PRIVATE FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Initialize the parameter list with defaults.
 */
void Partitioner::init_pl(RCP_ParameterList pl)
{
    // make a parameter list of defaults
    ParameterList defaults;

    // number of blocks
    defaults.set("num_blocks_i", 1);
    defaults.set("num_blocks_j", 1);
    defaults.set("num_z_blocks", 1);

    // dimension
    defaults.set("dimension", 3);

    // update the parameter list
    pl->setParametersNotAlreadySet(defaults);

    // determine dimension by examing the existing list
    bool has_k = (pl->isParameter("num_cells_k") ||
                  pl->isParameter("z_edges"));
    if (!has_k)
    {
        // set dimension to 2
        pl->set("dimension", 2);
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Build global uniform cell edges.
 */
void Partitioner::build_uniform_edges(int      num_cells,
                                      double   delta,
                                      Vec_Dbl& edges)
{
    Require (num_cells > 0);
    Require (delta > 0.);

    edges.resize(num_cells + 1);
    for (size_t i = 0; i < edges.size(); ++i)
        edges[i] = delta * i;

    Ensure (edges.size() == num_cells + 1);
    Ensure (profugus::soft_equiv(edges.front(), 0.));
    Ensure (profugus::soft_equiv(edges.back(), delta * num_cells));
}

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of Partitioner.cc
//---------------------------------------------------------------------------//
