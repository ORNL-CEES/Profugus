//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/Fission_Matrix_Processor.hh
 * \author Thomas M. Evans
 * \date   Fri Sep 12 11:08:38 2014
 * \brief  Fission_Matrix_Processor class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef mc_Fission_Matrix_Processor_hh
#define mc_Fission_Matrix_Processor_hh

#include <utility>
#include <unordered_map>
#include <vector>

#include "harness/DBC.hh"
#include "utils/Vector_Lite.hh"

namespace profugus
{

//===========================================================================//
/*!
 * \class Fission_Matrix_Processor
 * \brief Process the fission matrix.
 */
/*!
 * \example mc/test/tstFission_Matrix_Processor.cc
 *
 * Test of Fission_Matrix_Processor.
 */
//===========================================================================//

class Fission_Matrix_Processor
{
  public:
    //@{
    //! Typedefs.
    typedef profugus::Vector_Lite<int, 2> Children;

    //! Hash table for pair of ints.
    struct Idx_Hash
    {
      public:
        std::hash<int> d_hash;
        int            d_N;

        //! Constructor.
        Idx_Hash(int N = 0) : d_N(N) {/*...*/}

        size_t operator()(const std::pair<int, int> &x) const
        {
            REQUIRE(d_N > 0);
            return d_hash(x.first + d_N * x.second);
        }
    };

    //@{
    //! Sparse matrix storage for FM Tally.
    typedef std::pair<int, int>                       Idx;
    typedef std::unordered_map<Idx, double, Idx_Hash> Sparse_Matrix;
    typedef std::vector<double>                       Denominator;
    //@}

    //@{
    //! Flattened, ordered containers for fission matrix.
    typedef std::vector<Idx>    Ordered_Graph;
    typedef std::vector<double> Ordered_Matrix;
    //@}

    //! Node types.
    enum Node_Types
    {
        NONE     = -1,
        INTERNAL = 0,
        EXTERNAL = 1,
    };

  private:
    // >>> DATA

    // Graph of full-fission matrix (ordered).
    Ordered_Graph d_graph;

    // Fission matrix (ordered).
    Ordered_Matrix d_matrix;

  public:
    // Constructor.
    Fission_Matrix_Processor();

    // Build the fission matrix from local, on-processor contributions.
    void build_matrix(const Sparse_Matrix &local_matrix,
                      const Denominator &local_denominator);

    // Reset internal fission matrix memory.
    void reset();

    // >>> ACCESSORS

    //! Get the flattened, ordered, globally-reduced fission matrix.
    const Ordered_Matrix& matrix() const { return d_matrix; }

    //! Get the flattened, ordered, graph of of the global fission matrix.
    const Ordered_Graph& graph() const { return d_graph; }

    //! Get parent node.
    int parent() const { return d_parent; }

    //! Get children nodes.
    const Children& children() const { return d_children; }

    //! Node type.
    const Node_Types& node_type() const { return d_type; }

  private:
    // >>> IMPLEMENTATION

    // Node type.
    Node_Types d_type;

    // Parent node.
    int d_parent;

    // Children.
    Children d_children;

    // Number of domains and domain id.
    int d_node, d_nodes;
};

} // end namespace profugus

#endif // mc_Fission_Matrix_Processor_hh

//---------------------------------------------------------------------------//
//                 end of Fission_Matrix_Processor.hh
//---------------------------------------------------------------------------//
