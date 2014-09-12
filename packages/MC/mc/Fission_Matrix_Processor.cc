//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/Fission_Matrix_Processor.cc
 * \author Thomas M. Evans
 * \date   Fri Sep 12 11:08:38 2014
 * \brief  Fission_Matrix_Processor member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Fission_Matrix_Processor.hh"

#include "comm/global.hh"

namespace profugus
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
Fission_Matrix_Processor::Fission_Matrix_Processor()
    : d_type(INTERNAL)
    , d_parent(-1)
    , d_children(2, -1)
    , d_node(profugus::node())
    , d_nodes(profugus::nodes())
{
    // calculate the children and parents
    d_parent      = (d_node - 1) / 2;
    d_children[1] = (d_node + 1) * 2;
    d_children[0] = d_children[1] - 1;

    // node 0 is the root
    if (d_node == 0)
        d_parent = -1;

    // check children
    if (d_children[0] > d_nodes - 1)
    {
        d_children = -1;
        d_type     = EXTERNAL;
    }
    else if (d_children[1] > d_nodes - 1)
    {
        d_children[1] = -1;
    }
}

//---------------------------------------------------------------------------//
// PUBLIC FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Globally reduce and build the fission matrix.
 *
 * This builds a flattened fission matrix.  The order is stored in the graph.
 */
void Fission_Matrix_Processor::build_matrix(
    const Sparse_Matrix &local_matrix,
    const Denominator   &local_denominator)
{
}

//---------------------------------------------------------------------------//
/*!
 * \brief Reset internal fission matrix memory.
 */
void Fission_Matrix_Processor::reset()
{
}

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of Fission_Matrix_Processor.cc
//---------------------------------------------------------------------------//
