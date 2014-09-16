//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   spn_tpetra/Linear_System_FV.hh
 * \author Thomas M. Evans
 * \date   Fri Feb 14 19:58:19 2014
 * \brief  Linear_System_FV class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef spn_tpetra_Linear_System_FV_hh
#define spn_tpetra_Linear_System_FV_hh

#include <vector>

#include "harness/DBC.hh"
#include "utils/Vector_Lite.hh"
#include "utils/Definitions.hh"
#include "solvers/LinAlgTypedefs.hh"
#include "spn/FV_Bnd_Indexer.hh"
#include "spn/FV_Gather.hh"
#include "Linear_System.hh"

// Additional Trilinos pieces
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_LAPACK.hpp"

namespace profugus
{
namespace tpetra
{

//===========================================================================//
/*!
 * \class Linear_System_FV
 * \brief Build a linear SPN system based on a Cartesian Finite-Volume
 * discretization.
 */
/*!
 * \example spn_tpetra/test/tstLinear_System_FV.cc
 *
 * Test of Linear_System_FV.
 */
//===========================================================================//

class Linear_System_FV : public Linear_System
{
  private:
    //! Base class typedef.
    typedef Linear_System Base;

  public:
    //@{
    //! Typedefs.
    typedef def::Vec_Int                 Vec_Int;
    typedef def::Vec_Dbl                 Vec_Dbl;
    typedef Teuchos::RCP<FV_Bnd_Indexer> RCP_Bnd_Indexer;
    typedef typename TpetraTypes::MATRIX CrsMatrix_t;
    //@}

  private:
    // >>> DATA

    // Mesh.
    RCP_Mesh d_mesh;

    // L-G indexer for Cartesian mesh.
    RCP_Indexer d_indexer;

    // Element SPN and fission matrices.
    Teuchos::RCP<CrsMatrix_t> d_matrix;
    Teuchos::RCP<CrsMatrix_t> d_fission;

  public:
    // Constructor.
    Linear_System_FV(RCP_ParameterList db, RCP_Dimensions dim, RCP_Mat_DB mat,
                     RCP_Mesh mesh, RCP_Indexer indexer, RCP_Global_Data data,
                     RCP_Timestep dt = Teuchos::null);

    // Make the matrix.
    void build_Matrix();

    //! Build the right-hand-side fission matrix.
    void build_fission_matrix();

    // Build the right-hand-side from an external, isotropic source.
    void build_RHS(const External_Source &q);

    // Local/Global index in vector space corresponding to (group, eqn, cell).
    inline int index(int g, int eqn, int spatial_unknown) const;

    // >>> ACCESSORS

    //! Get an RCP to the LHS matrix (may not be full matrix)
    RCP_Matrix get_Matrix() const { return d_matrix; }

    // Get the boundary indexer for a face.
    inline RCP_Bnd_Indexer bnd_indexer(int face) const;

    //! Local number of volume unknowns (blocks) in the matrix.
    int vol_unknowns() const { return d_Nv_local; }

    //! Local number of boundary unknowns (blocks) in the matrix.
    int bnd_unknowns() const { return d_Nb_local; }

  private:
    // >>> IMPLEMENTATION

    // Types.
    typedef profugus::Vector_Lite<int, 3> Size_Vec;
    typedef FV_Gather::RCP_Face_Field     RCP_Face_Field;

    // Calculate boundary information.
    void calc_bnd_sizes();

    // Add boundary sources to the RHS.
    void add_boundary_sources(int face_id, int &face, const char *phi_b_str,
                              int num_face_cells);

    // Map local-to-global boundary elements on a face.
    void map_bnd_l2g(RCP_Bnd_Indexer indexer, int N_abscissa, int N_ordinate,
                     Vec_Int &l2g);

    // Insert spatially coupled elements.
    void spatial_coupled_element(int n, int i, int j, int k, int g_i, int g_j,
                                 RCP_Face_Field Dx_low, RCP_Face_Field Dx_high,
                                 RCP_Face_Field Dy_low, RCP_Face_Field Dy_high);

    // Add spatial element to the matrix.
    void add_spatial_element(int eqn, int row_cell, int col_cell,
                             double delta_l, double delta_r, double delta_c,
                             const Serial_Matrix &Dl, const Serial_Matrix &Dr);

    // Add boundary element to matrix.
    void build_bnd_element(int eqn, int global_row, int global_col,
                           double delta_c);

    // Add boundary equations to the matrix.
    void add_boundary_equations(int face, int global_cell, int local_cell,
                                double delta_c);

    // Insert a block matrix (GXG) into the matrix.
    void insert_block_matrix(int row_n, int row_cell, int row_off,
                             int col_m, int col_cell, int col_off,
                             const Serial_Matrix &M, CrsMatrix_t &matrix);

    // Gather object.
    FV_Gather d_gather;

    // Boundary indexers for each face (-x,+x,-y,+y,-z,+z).
    std::vector<RCP_Bnd_Indexer> d_bnd_index;

    // Number of groups.
    int d_Ng;

    // Number of equations.
    int d_Ne;

    // Dimensions of local mesh.
    Size_Vec d_N;
    int      d_Nc;

    // Dimensions of global mesh.
    Size_Vec d_G;
    int      d_Gc;

    // First and last global cell indices by dimension.
    int d_first, d_last_I, d_last_J, d_last_K;

    // Global and local block-size of volume unknowns.
    int d_Nv_global, d_Nv_local;

    // Global and local block-size of boundary equations.
    int d_Nb_global,  d_Nb_local;

    // Unknowns per computational cell and matrix-block size.
    int d_unknowns_per_cell, d_block_size;

    // Non-reflecting faces (-x,+x,-y,+y,-z,+z) -> 0 if reflecting, number of
    // global/local faces if vacuum/source
    int d_bc_global[6], d_bc_local[6];

    // Serial work matrices.
    Serial_Matrix d_D_c; // diffusion matrix in local cell for moment n
    Serial_Matrix d_C_c; // local cell matrix coefficient
    Serial_Matrix d_D;   // diffusion matrix in a neighbor cell for moment n
    Serial_Matrix d_C;   // matrix coefficient from neighbor cell
    Serial_Matrix d_W;   // work matrix

    // Global cell widths by dimension.
    profugus::Vector_Lite<def::Vec_Dbl, 3> d_widths;

    // LAPACK object.
    Teuchos::LAPACK<int, double> d_lapack;

    // LAPACK work arrays.
    def::Vec_Dbl d_work;
    def::Vec_Int d_ipiv;
    int          d_info;

    // Indices used to make matrix rows.
    Vec_Int d_indices; // the max will always be (6 (coupled cells) + 4
                       // (number of equations)) X number of groups; however
                       // we only insert one G x G block at a time, so this is
                       // set to G

    // Values for a matrix row.
    Vec_Dbl d_values;
};

//---------------------------------------------------------------------------//
// INLINE FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Local or global index in vector space corresponding to a (group,
 * eqn, cell).
 *
 * The cardinal index to an unknown location in the SPN vector space is:
 * \f[
   n = g + N_g(e + N_e(c))
 * \f]
 *
 * The index is global if the cell is global; it is local otherwise.
 *
 * \param g group index
 * \param eqn equation index in range \c [0,4)
 * \param cell local or global cell; if it local the index is local
 */
int Linear_System_FV::index(int g,
                            int eqn,
                            int cell) const
{
    REQUIRE(g < d_Ng);
    REQUIRE(eqn < d_Ne);
    return g + d_Ng * (eqn + d_Ne * (cell));
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get the boundary indexer for a face.
 *
 * \param face face index in the range \c [0,5] where the index corresponds to
 * \c -x,+x,-y,+y,-z,+z
 *
 * \return null if no boundary unknowns on the requested face
 */
Linear_System_FV::RCP_Bnd_Indexer Linear_System_FV::bnd_indexer(int face) const
{
    REQUIRE(face >= 0 && face < 6);
    return d_bnd_index[face];
}

} // end namespace profugus
} // end namespace tpetra

#endif // spn_tpetra_Linear_System_FV_hh

//---------------------------------------------------------------------------//
//                 end of Linear_System_FV.hh
//---------------------------------------------------------------------------//
