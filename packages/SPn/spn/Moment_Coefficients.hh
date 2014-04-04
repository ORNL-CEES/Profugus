//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   spn/Moment_Coefficients.hh
 * \author Thomas M. Evans
 * \date   Mon Feb 10 13:07:08 2014
 * \brief  Moment_Coefficients class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef spn_Moment_Coefficients_hh
#define spn_Moment_Coefficients_hh

#include <vector>

#include <spn/config.h>
#include "Teuchos_RCP.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_LAPACK.hpp"
#include "Teuchos_ParameterList.hpp"

#include "utils/Definitions.hh"
#include "utils/Vector_Lite.hh"
#include "utils/Static_Map.hh"
#include "xs/Mat_DB.hh"
#include "SPN_Constants.hh"
#include "Dimensions.hh"
#include "Timestep.hh"

namespace profugus
{

//===========================================================================//
/*!
 * \class Moment_Coefficients
 * \brief Coefficients used to couple moments in the SPN equations.
 */
/*!
 * \example spn/test/tstMoment_Coefficients.cc
 *
 * Test of Moment_Coefficients.
 */
//===========================================================================//

class Moment_Coefficients
{
  public:
    //@{
    //! Typedefs.
    typedef Mat_DB                                        Mat_DB_t;
    typedef Teuchos::RCP<Mat_DB_t>                        RCP_Mat_DB;
    typedef Mat_DB_t::XS_t                                XS_t;
    typedef Mat_DB_t::RCP_XS                              RCP_XS;
    typedef Teuchos::RCP<Teuchos::ParameterList>          RCP_ParameterList;
    typedef Teuchos::RCP<Dimensions>                      RCP_Dimensions;
    typedef Teuchos::RCP<Timestep>                        RCP_Timestep;
    typedef Teuchos::SerialDenseMatrix<int, double>       Serial_Matrix;
    typedef Teuchos::RCP<Serial_Matrix>                   RCP_Serial_Matrix;
    typedef Static_Map<def::size_type, RCP_Serial_Matrix> Hash_Table;
    typedef Teuchos::RCP<Hash_Table>                      RCP_Hash_Table;
    //@}

  private:
    // >>> DATA

    // SPN problem dimensions.
    RCP_Dimensions d_dim;

    // Material database.
    RCP_Mat_DB d_mat;

    // Timestep.
    RCP_Timestep d_dt;

  public:
    // Constructor.
    Moment_Coefficients(RCP_ParameterList db, RCP_Dimensions dim,
                        RCP_Mat_DB mat, RCP_Timestep dt = Teuchos::null);

    // Make cross section (Sigma_n) matrix.
    void make_Sigma(int n, int matid, Serial_Matrix &S);

    // Make diffusion matrices.
    void make_D(int n, int cell, Serial_Matrix &D);

    // Make A-matrix block entries.
    void make_A(int n, int m, int cell, Serial_Matrix &A);

    // Get B-matrix diagonal block entries.
    void make_B(int n, int m, Serial_Matrix &B);

    // Get F fission matrix block entries.
    void make_F(int n, int m, int cell, Serial_Matrix &F);

    // >>> ACCESSORS

    //! Number of groups.
    int num_groups() const { return d_Ng; }

    //! Number of equations.
    int num_equations() const { return d_dim->num_equations(); }

    //! Minimum number of moments in cross section data across all materials.
    int min_scattering_moments() const { return d_min_moments; }

    // >>> STATIC METHODS

    // Convert u->phi.
    static double u_to_phi(double u0, double u1, double u2, double u3)
    {
        return u0 - fractions::rat_2_3 * u1 + fractions::rat_8_15 * u2 -
            fractions::rat_16_35 * u3;
    }

  private:
    // >>> TYPES

    typedef def::Vec_Int           Vec_Int;
    typedef def::Vec_Dbl           Vec_Dbl;
    typedef Vector_Lite<double, 4> Coefficients;
    typedef XS_t::Vector           Vector;
    typedef XS_t::Matrix           Matrix;

    // >>> DATA

    // Number of groups.
    int d_Ng;

    // Work block matrix.
    Serial_Matrix d_W;

    // Mapping of equation to diffusion-moment order.
    int d_d[4];

    // Diffusion coefficients.
    double d_alpha[4];

    // Mapping of equation to SPN order for A-matrix.
    int d_a[4];

    // A-matrix cofficients.
    std::vector< std::vector<Coefficients> > d_c;

    // B-matrix coefficients.
    double d_b[4][4];

    // F-matrix coefficients.
    double d_f[4][4];

    // Storage for all Sigma matrices
    RCP_Hash_Table d_Sigma;

    // Minimum scattering moments in cross section data across all materials.
    int d_min_moments;

    // Flag to turn on Pn outscatter correction
    bool d_outscatter_correction;

    // LAPACK object.
    Teuchos::LAPACK<int, double> d_lapack;

    // Work areas for LAPACK.
    Vec_Dbl d_work;
    Vec_Int d_ipiv;
    int     d_info;
};

} // end namespace profugus

#endif // spn_Moment_Coefficients_hh

//---------------------------------------------------------------------------//
//                 end of Moment_Coefficients.hh
//---------------------------------------------------------------------------//
