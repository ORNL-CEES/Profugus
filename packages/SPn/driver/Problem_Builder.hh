//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   driver/Problem_Builder.hh
 * \author Thomas M. Evans
 * \date   Wed Mar 12 22:25:22 2014
 * \brief  Problem_Builder class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef driver_Problem_Builder_hh
#define driver_Problem_Builder_hh

#include <string>
#include <vector>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_TwoDArray.hpp"

#include "xs/Mat_DB.hh"
#include "mesh/Partitioner.hh"
#include "SPn/Isotropic_Source.hh"

namespace spn
{

//===========================================================================//
/*!
 * \class Problem_Builder
 * \brief Read and initialize an SPn problem,
 */
//===========================================================================//

class Problem_Builder
{
  public:
    //@{
    //! Typedefs.
    typedef Teuchos::ParameterList                   ParameterList;
    typedef Teuchos::RCP<ParameterList>              RCP_ParameterList;
    typedef Teuchos::RCP<profugus::Mat_DB>           RCP_Mat_DB;
    typedef profugus::Partitioner::RCP_Mesh          RCP_Mesh;
    typedef profugus::Partitioner::RCP_Indexer       RCP_Indexer;
    typedef profugus::Partitioner::RCP_Global_Data   RCP_Global_Data;
    typedef Teuchos::RCP<profugus::Isotropic_Source> RCP_Source;
    //@}

  private:
    // >>> DATA

    // Problem-parameterlist (talks to solver components).
    RCP_ParameterList d_db;

    // Mesh objects.
    RCP_Mesh        d_mesh;
    RCP_Indexer     d_indexer;
    RCP_Global_Data d_gdata;

    // Material database.
    RCP_Mat_DB d_mat;

    // External source.
    RCP_Source d_source;

  public:
    // Constructor.
    Problem_Builder();

    // Setup the problem.
    void setup(const std::string &xml_file);

    // >>> ACCESSORS

    //! Get problem database.
    RCP_ParameterList problem_db() const { return d_db; }

    //@{
    //! Get the mesh objects.
    RCP_Mesh mesh() const { return d_mesh; }
    RCP_Indexer indexer() const { return d_indexer; }
    RCP_Global_Data global_data() const { return d_gdata; }
    //@}

    //! Get the material database.
    RCP_Mat_DB mat_db() const { return d_mat; }

    //! Get the external source (could be null).
    RCP_Source source() const { return d_source; }

  private:
    // >>> IMPLEMENTATION

    // Typedefs.
    typedef Teuchos::Comm<int>          Comm;
    typedef Teuchos::RCP<const Comm>    RCP_Comm;
    typedef Teuchos::Array<int>         OneDArray_int;
    typedef Teuchos::Array<double>      OneDArray_dbl;
    typedef Teuchos::Array<std::string> OneDArray_str;
    typedef Teuchos::TwoDArray<int>     TwoDArray_int;
    typedef Teuchos::TwoDArray<double>  TwoDArray_dbl;

    // Build implementation.
    void build_mesh();
    void build_matids();
    void calc_axial_matids(int level, TwoDArray_int &matids);
    void build_matdb();
    void build_source(const ParameterList &source_db);

    // Number of assemblies and pins per assembly.
    int d_Na[2];
    int d_Np[2];

    // Local material ids.
    std::vector<int> d_matids;

    // Teuchos communicator.
    RCP_Comm d_comm;

    // Problem-setup parameterlists.
    RCP_ParameterList d_coredb;
    RCP_ParameterList d_assblydb;
    RCP_ParameterList d_matdb;
};

} // end namespace spn

#endif // driver_Problem_Builder_hh

//---------------------------------------------------------------------------//
//                 end of Problem_Builder.hh
//---------------------------------------------------------------------------//
