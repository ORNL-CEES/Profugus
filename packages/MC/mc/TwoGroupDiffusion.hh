//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   MC/mc/TwoGroupDiffusion.hh
 * \author Steven Hamilton
 * \date   Tue Aug 07 10:31:04 2018
 * \brief  TwoGroupDiffusion class declaration.
 * \note   Copyright (c) 2018 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef MC_mc_TwoGroupDiffusion_hh
#define MC_mc_TwoGroupDiffusion_hh

#include <vector>

#include "Utils/harness/DBC.hh"
#include "Utils/comm/global.hh"
#include "TwoGroupCrossSections.hh"

// Trilinos includes
#include "Teuchos_RCP.hpp"
#include "Tpetra_Operator.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_RowMatrix.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_CrsGraph.hpp"
#include "Kokkos_DefaultNode.hpp"

namespace mc
{

//===========================================================================//
/*!
 * \class TwoGroupDiffusion
 * \brief Two-group, 3D neutron diffusion solver.
 */
/*!
 * \example mc/test/tstTwoGroupDiffusion.cc
 *
 * Test of TwoGroupDiffusion.
 */
//===========================================================================//

class TwoGroupDiffusion
{
  public:

    enum BC_TYPE {VACUUM, REFLECT};
    enum FACE    {LO_X, HI_X, LO_Y, HI_Y, LO_Z, HI_Z};

    //@{
    //! Typedefs
    using XS      = mc::TwoGroupCrossSections;
    using XS_Data = XS::XS_Data;
    using Vec_Dbl = std::vector<double>;
    using Vec_Int = std::vector<int>;
    using Vec_BC  = std::vector<BC_TYPE>;
    using Pin_Map = std::vector<XS::XS_Type>;
    using ST = double;
    using LO = int;
    using GO = int;
    using NODE   = KokkosClassic::DefaultNode::DefaultNodeType;
    using MV     = Tpetra::MultiVector<ST,LO,GO,NODE>;
    using VECTOR = Tpetra::Vector<ST,LO,GO,NODE>;
    using OP     = Tpetra::Operator<ST,LO,GO,NODE>;
    using MAP    = Tpetra::Map<LO,GO,NODE>;
    using MATRIX = Tpetra::CrsMatrix<ST,LO,GO,NODE>;
    using GRAPH  = Tpetra::CrsGraph<LO,GO,NODE>;
    //@}

  private:
    // >>> DATA
    Vec_Dbl d_dx;
    Vec_Dbl d_dy;
    Vec_Dbl d_dz;
    Pin_Map d_pin_map;
    XS      d_xs;
    size_t d_Nx, d_Ny, d_Nz, d_num_cells;

    // Solver parameters
    double d_tol;
    int    d_max_iters;

    // Boundary conditions
    std::vector<BC_TYPE> d_bcs;

    // Data stored over single group
    Teuchos::RCP<MAP>    d_map;
    Teuchos::RCP<GRAPH>  d_graph;
    Teuchos::RCP<MATRIX> d_A_therm, d_A_fast;
    Teuchos::RCP<VECTOR> d_x_therm, d_x_fast;
    Teuchos::RCP<VECTOR> d_b_therm, d_b_fast;

  public:

    // Constructor
    TwoGroupDiffusion(const Vec_Dbl& dx,
                      const Vec_Dbl& dy,
                      const Vec_Dbl& dz,
                      const Pin_Map& map,
                      const Vec_BC&  bcs);

    // Solve
    void solve(const Vec_Dbl& temperatures,
               const Vec_Dbl& densities,
                     Vec_Dbl& power);

  private:

    int cellid(int ix, int iy, int iz)
    {
        REQUIRE(ix >= 0);
        REQUIRE(ix < d_Nx);
        REQUIRE(iy >= 0);
        REQUIRE(iy < d_Ny);
        REQUIRE(iz >= 0);
        REQUIRE(iz < d_Nz);
        int cell = ix + d_Nx * (iy + d_Ny * iz);
        ENSURE(cell >= 0);
        ENSURE(cell <  d_num_cells);
        return cell;
    }

    void build_matrices(const std::vector<XS_Data>& xs_data);

    // Compute 2-norm of vector
    double vec_norm(const Teuchos::ArrayRCP<double>& vec)
    {
        double nrm = 0.0;
        for (auto val : vec)
            nrm += val * val;
        return std::sqrt(nrm);
    }

    // Apply scale factor to vector data
    void scale_vec(Teuchos::ArrayRCP<double>& vec, double factor)
    {
        for (auto& val : vec)
            val *= factor;
    }
};

//---------------------------------------------------------------------------//
} // end namespace mc

//---------------------------------------------------------------------------//
#endif // MC_mc_TwoGroupDiffusion_hh

//---------------------------------------------------------------------------//
// end of MC/mc/TwoGroupDiffusion.hh
//---------------------------------------------------------------------------//
