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

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_TwoDArray.hpp"

#include "xs/Mat_DB.hh"
#include "mesh/Partitioner.hh"

namespace profugus
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
    typedef Teuchos::ParameterList       ParameterList;
    typedef Teuchos::RCP<ParameterList>  RCP_ParameterList;
    typedef Teuchos::RCP<Mat_DB>         RCP_Mat_DB;
    typedef Partitioner::RCP_Mesh        RCP_Mesh;
    typedef Partitioner::RCP_Indexer     RCP_Indexer;
    typedef Partitioner::RCP_Global_Data RCP_Global_Data;
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

  private:
    // >>> IMPLEMENTATION

    // Typedefs.
    typedef Teuchos::Comm<int>          Comm;
    typedef Teuchos::RCP<const Comm>    RCP_Comm;
    typedef Teuchos::Array<int>         OneDArray_int;
    typedef Teuchos::Array<double>      OneDArray_dbl;
    typedef Teuchos::Array<std::string> OneDArray_str;
    typedef Teuchos::TwoDArray<int>     TwoDArray_int;

    // Build implementation.
    void build_mesh();
    void build_matdb();

    // Teuchos communicator.
    RCP_Comm d_comm;

    // Validating Parameterlist.
    RCP_ParameterList d_validator;

    // Problem-setup parameterlists.
    RCP_ParameterList d_coredb;
    RCP_ParameterList d_assblydb;
    RCP_ParameterList d_matdb;
    RCP_ParameterList d_meshdb;
};

} // end namespace profugus

#endif // driver_Problem_Builder_hh

//---------------------------------------------------------------------------//
//                 end of Problem_Builder.hh
//---------------------------------------------------------------------------//
