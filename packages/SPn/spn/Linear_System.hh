//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   spn/Linear_System.hh
 * \author Thomas M. Evans
 * \date   Sun Oct 28 18:37:01 2012
 * \brief  Linear_System class definition.
 * \note   Copyright (C) 2012 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef spn_Linear_System_hh
#define spn_Linear_System_hh

#include <string>
#include <SPn/config.h>

// Trilinos Includes
#include "Teuchos_RCP.hpp"
#include "Teuchos_Array.hpp"
#include "Epetra_Map.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_Vector.h"

#include "mesh/Mesh.hh"
#include "mesh/LG_Indexer.hh"
#include "mesh/Global_Mesh_Data.hh"
#include "Isotropic_Source.hh"
#include "Moment_Coefficients.hh"

namespace profugus
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
 * For now this uses Epetra, but it will probably get converted to Thyra at a
 * later date.
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
    typedef Teuchos::RCP<Moment_Coefficients>      RCP_Moment_Coefficients;
    typedef Epetra_Map                             Map_t;
    typedef Epetra_RowMatrix                       Matrix_t;
    typedef Epetra_Vector                          Vector_t;
    typedef Teuchos::RCP<Epetra_Map>               RCP_Map;
    typedef Teuchos::RCP<Epetra_RowMatrix>         RCP_Matrix;
    typedef Teuchos::RCP<Epetra_Operator>          RCP_Operator;
    typedef Teuchos::RCP<Epetra_Vector>            RCP_Vector;
    typedef Isotropic_Source                       External_Source;
    typedef Teuchos::RCP<Mesh>                     RCP_Mesh;
    typedef Teuchos::RCP<LG_Indexer>               RCP_Indexer;
    typedef Teuchos::RCP<Global_Mesh_Data>         RCP_Global_Data;
    typedef Teuchos::Array<int>                    Array_Int;
    typedef Teuchos::Array<double>                 Array_Dbl;
    //@}

    //@{
    //! Epetra communicator type.
#ifdef COMM_MPI
    typedef Epetra_MpiComm    Comm;
#else
    typedef Epetra_SerialComm Comm;
#endif
    //@}

  protected:
    // >>> SHARED DATA

    // Problem database.
    RCP_ParameterList b_db;

    // SPN dimensions.
    RCP_Dimensions b_dim;

    // Material database.
    RCP_Mat_DB b_mat;

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
    Linear_System(RCP_ParameterList db, RCP_Dimensions dim, RCP_Mat_DB mat);

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

#endif // spn_Linear_System_hh

//---------------------------------------------------------------------------//
//              end of spn/Linear_System.hh
//---------------------------------------------------------------------------//
