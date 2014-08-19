//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   spn_tpetra/Linear_System.hh
 * \author Thomas M. Evans
 * \date   Sun Oct 28 18:37:01 2012
 * \brief  Linear_System class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef spn_tpetra_Linear_System_hh
#define spn_tpetra_Linear_System_hh

#include <string>
#include <SPn/config.h>

// Trilinos Includes
#include "Teuchos_RCP.hpp"
#include "Teuchos_Array.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_RowMatrix.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_RowGraph.hpp"

#include "mesh/Mesh.hh"
#include "mesh/LG_Indexer.hh"
#include "mesh/Global_Mesh_Data.hh"
#include "spn/Isotropic_Source.hh"
#include "spn/Moment_Coefficients.hh"

namespace profugus
{
namespace tpetra
{

//===========================================================================//
/*!
 * \class Linear_System
 * \brief Base class that defines the Linear System for an SPN problem.
 *
 * This class builds the SPN linear system for both fixed source:
 * \f[
   \mathbf{A}\mathbf{u} = \mathbf{Q}\:,
 * \f]
 * and eigenvalue problems:
 * \f[
   \mathbf{A}\mathbf{u} = \frac{1}{k}\mathbf{B}\mathbf{u}\:.
 * \f]
 * The operators \b A and \b B along with vectors \b u and \b Q are defined in
 * the SPN technical note and Denovo methods manual.
 *
 * This class is equivalent to spn_tpetra/Linear_System but uses Tpetra instead.
 */
//===========================================================================//

class Linear_System
{
  public:
    //@{
    //! Typedefs.
    typedef Moment_Coefficients::Serial_Matrix     Serial_Matrix;
    typedef Moment_Coefficients::RCP_Mat_DB        RCP_Mat_DB;
    typedef Moment_Coefficients::RCP_Dimensions    RCP_Dimensions;
    typedef Moment_Coefficients::RCP_ParameterList RCP_ParameterList;
    typedef Moment_Coefficients::RCP_Timestep      RCP_Timestep;
    typedef Teuchos::RCP<Moment_Coefficients>      RCP_Moment_Coefficients;
    typedef KokkosClassic::SerialNode              Node;
    typedef Tpetra::Map<int,int,Node>              Map_t;
    typedef Tpetra::RowGraph<int,int,Node>         Graph_t;
    typedef Tpetra::RowMatrix<double,int,int,Node> Matrix_t;
    typedef Tpetra::CrsMatrix<double,int,int,Node> CrsMatrix_t;
    typedef Tpetra::Operator<double,int,int,Node>  Operator_t;
    typedef Tpetra::Vector<double,int,int,Node>    Vector_t;
    typedef Teuchos::RCP<Map_t>                    RCP_Map;
    typedef Teuchos::RCP<const Graph_t>            RCP_Graph;
    typedef Teuchos::RCP<Matrix_t>                 RCP_Matrix;
    typedef Teuchos::RCP<Operator_t>               RCP_Operator;
    typedef Teuchos::RCP<Vector_t>                 RCP_Vector;
    typedef Isotropic_Source                       External_Source;
    typedef Teuchos::RCP<Mesh>                     RCP_Mesh;
    typedef Teuchos::RCP<LG_Indexer>               RCP_Indexer;
    typedef Teuchos::RCP<Global_Mesh_Data>         RCP_Global_Data;
    typedef Teuchos::Array<int>                    Array_Int;
    typedef Teuchos::Array<double>                 Array_Dbl;
    //@}

    //@{
    //! Tpetra communicator type.
    typedef Teuchos::Comm<int>    Comm;
    //@}

  protected:
    // >>> SHARED DATA

    // Problem database.
    RCP_ParameterList b_db;

    // SPN dimensions.
    RCP_Dimensions b_dim;

    // Material database.
    RCP_Mat_DB b_mat;

    // Timestep.
    RCP_Timestep b_dt;

    // Moment-coefficient generator.
    RCP_Moment_Coefficients b_mom_coeff;

    // Epetra objects.
    RCP_Map      b_map;
    RCP_Operator b_operator; // SPN matrix
    RCP_Operator b_fission;  // Fission matrix
    RCP_Vector   b_rhs;

    // Isotropic source coefficients.
    double b_src_coefficients[4];

    // Isotropic boundary source coefficients.
    double b_bnd_coefficients[4];

    // Processor nodes.
    int b_node, b_nodes;

  public:
    // Constructor.
    Linear_System(RCP_ParameterList db, RCP_Dimensions dim, RCP_Mat_DB mat,
                  RCP_Timestep dt);

    // Destructor.
    virtual ~Linear_System() = 0;

    //! Make the matrix.
    virtual void build_Matrix() = 0;

    //! Build the right-hand-side from an external, isotropic source.
    virtual void build_RHS(const External_Source &q) = 0;

    //! Build the right-hand-side fission matrix.
    virtual void build_fission_matrix() = 0;

    // >>> ACCESSORS

    //! Get an RCP to the communication map.
    RCP_Map get_Map() const { return b_map; }

    //! Get an RCP to the full LHS Operator.
    RCP_Operator get_Operator() const { return b_operator; }

    //! Get an RCP to the LHS matrix (may not be full matrix)
    virtual RCP_Matrix get_Matrix() const { return Teuchos::null; }

    //! Get an RCP to the Right-Hand-Side vector.
    RCP_Vector get_RHS() const { return b_rhs; }

    //! Get an RCP to the Right-Hand-Side vector.
    RCP_Operator get_fission_matrix() const { return b_fission; }

    //! Get problem dimensions.
    RCP_Dimensions get_dims() const { return b_dim; }

    //! Get the local/global index for (group, eqn, spatial_unknown).
    virtual int index(int g, int eqn, int spatial_unknown) const = 0;
};

} // end namespace profugus
} // end namespace tpetra

#endif // spn_tpetra_Linear_System_hh

//---------------------------------------------------------------------------//
//              end of spn_tpetra/Linear_System.hh
//---------------------------------------------------------------------------//
