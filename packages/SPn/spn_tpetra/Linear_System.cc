//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   spn_tpetra/Linear_System.cc
 * \author Thomas M. Evans
 * \date   Sun Oct 28 18:37:01 2012
 * \brief  Linear_System member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Teuchos_ParameterList.hpp"

#include "harness/DBC.hh"
#include "comm/global.hh"
#include "spn/SPN_Constants.hh"
#include "Linear_System.hh"

namespace profugus
{
namespace tpetra
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR AND DESTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor
 *
 * \param dim SPN dimensions object
 * \param mat material database
 */
Linear_System::Linear_System(RCP_ParameterList db,
                             RCP_Dimensions    dim,
                             RCP_Mat_DB        mat,
                             RCP_Timestep      dt)
    : b_db(db)
    , b_dim(dim)
    , b_mat(mat)
    , b_dt(dt)
    , b_mom_coeff(Teuchos::rcp(new Moment_Coefficients(db, dim, mat, dt)))
    , b_node(profugus::node())
    , b_nodes(profugus::nodes())
{
    REQUIRE(!db.is_null());
    REQUIRE(!dim.is_null());
    REQUIRE(!mat.is_null());

    // source coefficients
    b_src_coefficients[0] =  1.0;
    b_src_coefficients[1] = -fractions::rat_2_3;
    b_src_coefficients[2] =  fractions::rat_8_15;
    b_src_coefficients[3] = -fractions::rat_16_35;

    // bnd source coefficients
    b_bnd_coefficients[0] =  fractions::rat_1_2;
    b_bnd_coefficients[1] = -fractions::rat_1_8;
    b_bnd_coefficients[2] =  fractions::rat_1_16;
    b_bnd_coefficients[3] = -fractions::rat_5_128;

    if (!b_db->isParameter("boundary"))
    {
        b_db->set("boundary", "reflect");
    }

    // add the default boundary treatment for reflecting
    if (b_db->get<std::string>("boundary") == "reflect")
    {
        if (!b_db->isSublist("boundary_db"))
        {
            // set to all reflecting b.c. if undefined
            Teuchos::ParameterList boundary("boundary");
            Array_Int reflect(6, 1);
            boundary.set("reflect", reflect);

            b_db->set("boundary_db", boundary);
        }
    }

    ENSURE(!b_mom_coeff.is_null());
}

//---------------------------------------------------------------------------//
/*!
 * \brief Virtual destructor.
 */
Linear_System::~Linear_System()
{
}

} // end namespace profugus
} // end namespace tpetra

//---------------------------------------------------------------------------//
//                 end of Linear_System.cc
//---------------------------------------------------------------------------//
