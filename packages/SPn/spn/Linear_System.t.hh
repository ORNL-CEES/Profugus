//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   SPn/spn/Linear_System.t.hh
 * \author Thomas M. Evans
 * \date   Sun Oct 28 18:37:01 2012
 * \brief  Linear_System template member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef SPn_spn_Linear_System_t_hh
#define SPn_spn_Linear_System_t_hh

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_Ptr.hpp"
#include "Epetra_CrsMatrix.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_RowMatrixTransposer.h"

#include "harness/DBC.hh"
#include "harness/Warnings.hh"
#include "comm/global.hh"
#include "SPN_Constants.hh"
#include "Linear_System.hh"

namespace profugus
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
template <class T>
Linear_System<T>::Linear_System(RCP_ParameterList db,
                                RCP_Dimensions    dim,
                                RCP_Mat_DB        mat,
                                RCP_Timestep      dt)
    : b_db(db)
    , b_dim(dim)
    , b_mat(mat)
    , b_dt(dt)
    , b_mom_coeff(Teuchos::rcp(new Moment_Coefficients(db, dim, mat, dt)))
    , b_adjoint(false)
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
template <class T>
Linear_System<T>::~Linear_System()
{
}

//---------------------------------------------------------------------------//
// PUBLIC FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Setup system to build adjoints.
 *
 * If this is set to true, then build adjoint operators from the forward.
 * Note that the forward operators are never changed.  The forward operators
 * can be retained by subsequent calls to this function with adjoint set to
 * false.
 *
 * \pre only works with Epetra
 */
template <class T>
void Linear_System<T>::set_adjoint(bool adjoint)
{
    // no op for anything but Trilinos
    ADD_WARNING("Adjoint only supported for Epetra, turning off adjoint");
}

//---------------------------------------------------------------------------//
/*!
 * \brief Epecialization on adjoints for Epetra.
 */
template<>
void Linear_System<EpetraTypes>::set_adjoint(bool adjoint)
{
    REQUIRE(Teuchos::nonnull(b_operator));
    REQUIRE(Teuchos::nonnull(b_fission));

    // set adjoint flag
    b_adjoint = adjoint;

    // if this is adjoint build adjoint operators
    if (b_adjoint)
    {
        // cast the operator to a Epetra_RowMatrix
        Teuchos::RCP<Epetra_RowMatrix> a =
            Teuchos::rcp_dynamic_cast<Epetra_RowMatrix>(b_operator);
        Teuchos::RCP<Epetra_RowMatrix> b =
            Teuchos::rcp_dynamic_cast<Epetra_RowMatrix>(b_fission);

        // if this is a Row_Matrix then we can use the row-matrix transposer
        if (Teuchos::nonnull(a) && Teuchos::nonnull(b))
        {
            Epetra_RowMatrixTransposer a_trans(a.getRawPtr());
            Epetra_RowMatrixTransposer b_trans(b.getRawPtr());

            // make raw pointers for the new matrices
            Epetra_CrsMatrix *aT = nullptr, *bT = nullptr;

            // make the transposes
            a_trans.CreateTranspose(true, aT);
            b_trans.CreateTranspose(true, bT);
            CHECK(aT != nullptr);
            CHECK(bT != nullptr);

            // assign to the adjoint operators
            b_adjoint_operator = Teuchos::rcp(aT);
            b_adjoint_fission  = Teuchos::rcp(bT);
        }

        // otherwise set use transpose and hope for the best
        else
        {
            b_adjoint_operator = b_operator;
            b_adjoint_fission  = b_fission;
        }

        ENSURE(Teuchos::nonnull(b_adjoint_operator));
        ENSURE(Teuchos::nonnull(b_adjoint_fission));
    }

    // set the regular operators transpose flag
    b_operator->SetUseTranspose(b_adjoint);
    b_fission->SetUseTranspose(b_adjoint);
}

} // end namespace profugus

#endif // SPn_spn_Linear_System_t_hh

//---------------------------------------------------------------------------//
//                 end of Linear_System.t.hh
//---------------------------------------------------------------------------//
