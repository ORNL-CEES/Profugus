//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mesh/Partitioner.hh
 * \author Thomas M. Evans
 * \date   Wed Feb 12 09:54:53 2014
 * \brief  Partitioner class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef mesh_Partitioner_hh
#define mesh_Partitioner_hh

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

#include "harness/DBC.hh"
#include "utils/Vector_Lite.hh"
#include "utils/Definitions.hh"

#include "Global_Mesh_Data.hh"
#include "LG_Indexer.hh"
#include "Mesh.hh"

namespace profugus
{

//===========================================================================//
/*!
 * \class Partitioner
 * \brief Partition a Mesh object.
 *
 * This partitioner does a very simple partitioning where the client specifies
 * the number of partitions in \e (i,j).  The partitioner tries to make each
 * block (mesh on a processor) the same size.  If the number of cells in each
 * direction does not divide evenly, then cells are added to each direction
 * starting at \e (i=0,j=0). The mesh is decomposed into \e B blocks.
 */
/*!
 * \example mesh/test/tstPartitioner.cc
 *
 * Test of Partitioner.
 */
//===========================================================================//

class Partitioner
{
  public:
    //@{
    //! Typedefs.
    typedef Teuchos::ParameterList              ParameterList;
    typedef Teuchos::RCP<ParameterList>         RCP_ParameterList;
    typedef Mesh                                Mesh_t;
    typedef short unsigned int                  Dimension;
    typedef Mesh_t::size_type                   size_type;
    typedef Mesh_t::dim_type                    dim_type;
    typedef Mesh_t::Dim_Vector                  Dim_Vector;
    typedef Mesh_t::Space_Vector                Space_Vector;
    typedef def::Vec_Int                        Vec_Int;
    typedef def::Vec_Dbl                        Vec_Dbl;
    typedef profugus::Vector_Lite<size_type, 2> IJ_Set;
    typedef profugus::Vector_Lite<Vec_Dbl, 3>   IJK_Vec_Dbl;
    typedef profugus::Vector_Lite<Vec_Int, 2>   IJ_Vec_Int;
    //@}

    //@{
    //! Object typedefs.
    typedef Teuchos::RCP<Mesh_t>           RCP_Mesh;
    typedef Teuchos::RCP<LG_Indexer>       RCP_Indexer;
    typedef Teuchos::RCP<Global_Mesh_Data> RCP_Global_Data;
    //@}

  private:
    // >>> DATA

    //! Mesh, created during \c build().
    RCP_Mesh d_mesh;

    //! LG_Indexer, created during \c build().
    RCP_Indexer d_indexer;

    //! Global mesh data, created during \c build().
    RCP_Global_Data d_data;

    //! Global cell edges in each direction.
    IJK_Vec_Dbl d_edges;

  public:
    // Constructor.
    Partitioner(RCP_ParameterList pl);

    // Partition the mesh.
    void build();

    // >>> ACCESSORS

    //! Get the mesh.
    RCP_Mesh get_mesh() const { return d_mesh; }

    //! Get the local-global indexer.
    RCP_Indexer get_indexer() const { return d_indexer; }

    //! Get global mesh data.
    RCP_Global_Data get_global_data() const { return d_data; }

    // Dimensionality of this partitioner
    Dimension dimension() const { return d_dimension; }

    //! Get the number of blocks per set.
    size_type num_blocks() const { return d_num_blocks; }

    //! Get the number of processor blocks in the I/J direction.
    size_type num_blocks(int dir) const { return d_Nb[dir]; }

  private:
    // >>> IMPLEMENTATION

    //  Initialize parameterlist.
    void init_pl(RCP_ParameterList pl);

    // Build global uniform edges of the mesh.
    void build_uniform_edges(int num_cells, double delta, Vec_Dbl& edges);

    //! Local domain (node)
    size_type d_domain;

    //! Local spatial block, initialized during \c build().
    size_type d_block;

    //! Number of processors (domains).
    size_type d_nodes;

    //! Number of blocks
    size_type d_num_blocks;

    //! Number of k-blocks.
    size_type d_k_blocks;

    //! Number of blocks in I/J directions.
    IJ_Set d_Nb;

    //! Dimension: 3 for 3-D, 2 for 2-D
    Dimension d_dimension;
};

} // end namespace profugus

#endif // mesh_Partitioner_hh

//---------------------------------------------------------------------------//
//                 end of Partitioner.hh
//---------------------------------------------------------------------------//
