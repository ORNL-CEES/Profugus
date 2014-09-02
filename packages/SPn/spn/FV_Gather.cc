//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   spn/FV_Gather.cc
 * \author Thomas M. Evans
 * \date   Tue Oct 30 11:04:52 2012
 * \brief  FV_Gather member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Dimensions.hh"
#include "FV_Gather.hh"

namespace profugus
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
FV_Gather::FV_Gather(RCP_Mesh                 mesh,
                     RCP_Moment_Coefficients  coefficients,
                     const Indexer_t         &indexer)
    : d_mesh(mesh)
    , d_coefficients(coefficients)
    , d_Nb(indexer.num_blocks(def::I), indexer.num_blocks(def::J))
    , d_ij(mesh->block(def::I), mesh->block(def::J))
    , d_domain(profugus::node())
    , d_domains(profugus::nodes())
{
    using def::I; using def::J; using def::K;
    using def::X; using def::Y; using def::Z;
    using def::PROBLEM_BOUNDARY;

    // only allow 1-set problems
    Insist (indexer.num_sets() == 1,
            "SPN solver does not support multiple sets.");

    // for now, only allow 1-block in Z
    Insist (mesh->block(K) == 1, "Only allow 1 Z-block in SPN.");

    REQUIRE(!mesh.is_null());
    REQUIRE(!coefficients.is_null());
    REQUIRE(d_domain == d_ij[I] + d_ij[J] * d_Nb[I]);
    REQUIRE(d_domains == d_Nb[I] * d_Nb[J]);
    REQUIRE(d_domains == 1 ? d_ij[I] == 0 && d_ij[J] == 0 : true);

    // number of groups
    int Ng = d_coefficients->num_groups();

    // initialize all neighbors to problem boundary
    d_neighbor_I[LO] = PROBLEM_BOUNDARY;
    d_neighbor_I[HI] = PROBLEM_BOUNDARY;
    d_neighbor_J[LO] = PROBLEM_BOUNDARY;
    d_neighbor_J[HI] = PROBLEM_BOUNDARY;

    // first/last blocks in (i,j)
    int first  = 0;
    int last_I = d_Nb[I] - 1;
    int last_J = d_Nb[J] - 1;

    // make face fields and define neighbor blocks
    if (d_ij[I] > first)
    {
        d_neighbor_I[LO] = convert(d_ij[I] - 1, d_ij[J]);
        d_incoming_I[LO] = Teuchos::rcp(new SDM_Face_Field(*d_mesh, X, Ng));
        d_outgoing_I[LO] = Teuchos::rcp(new SDM_Face_Field(*d_mesh, X, Ng));

        CHECK(d_neighbor_I[LO] < d_domains);
        CHECK(d_neighbor_I[LO] >= 0);
    }
    if (d_ij[I] < last_I)
    {
        d_neighbor_I[HI] = convert(d_ij[I] + 1, d_ij[J]);
        d_incoming_I[HI] = Teuchos::rcp(new SDM_Face_Field(*d_mesh, X, Ng));
        d_outgoing_I[HI] = Teuchos::rcp(new SDM_Face_Field(*d_mesh, X, Ng));

        CHECK(d_neighbor_I[HI] < d_domains);
        CHECK(d_neighbor_I[HI] >= 0);
    }
    if (d_ij[J] > first)
    {
        d_neighbor_J[LO] = convert(d_ij[I], d_ij[J] - 1);
        d_incoming_J[LO] = Teuchos::rcp(new SDM_Face_Field(*d_mesh, Y, Ng));
        d_outgoing_J[LO] = Teuchos::rcp(new SDM_Face_Field(*d_mesh, Y, Ng));

        CHECK(d_neighbor_J[LO] < d_domains);
        CHECK(d_neighbor_J[LO] >= 0);
    }
    if (d_ij[J] < last_J)
    {
        d_neighbor_J[HI] = convert(d_ij[I], d_ij[J] + 1);
        d_incoming_J[HI] = Teuchos::rcp(new SDM_Face_Field(*d_mesh, Y, Ng));
        d_outgoing_J[HI] = Teuchos::rcp(new SDM_Face_Field(*d_mesh, Y, Ng));

        CHECK(d_neighbor_J[HI] < d_domains);
        CHECK(d_neighbor_J[HI] >= 0);
    }

    ENSURE(d_domains == 1 ? d_incoming_I[LO].is_null() : true);
    ENSURE(d_domains == 1 ? d_incoming_J[LO].is_null() : true);
    ENSURE(d_domains == 1 ? d_incoming_I[HI].is_null() : true);
    ENSURE(d_domains == 1 ? d_incoming_J[HI].is_null() : true);
    ENSURE(d_domains == 1 ? d_neighbor_I[LO] == PROBLEM_BOUNDARY : true);
    ENSURE(d_domains == 1 ? d_neighbor_J[LO] == PROBLEM_BOUNDARY : true);
    ENSURE(d_domains == 1 ? d_neighbor_I[HI] == PROBLEM_BOUNDARY : true);
    ENSURE(d_domains == 1 ? d_neighbor_J[HI] == PROBLEM_BOUNDARY : true);
}

//---------------------------------------------------------------------------//
// PUBLIC FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Gather data on all blocks.
 *
 * \param eqn equation order of diffusion coefficients in range [0,4)
 */
void FV_Gather::gather(int eqn)
{
    using def::I; using def::J; using def::PROBLEM_BOUNDARY;

    REQUIRE(eqn >= 0 && eqn < Dimensions::max_num_equations());

    // return immediately if only running on 1 domain (although it works
    // without this, we use if for efficiency only)
    if (d_domains == 1) return;

    // post receives
    post_receives();

    // fill low faces with local diffusion coefficients
    fill_I_face(eqn, d_outgoing_I[LO], 0);
    fill_J_face(eqn, d_outgoing_J[LO], 0);
    fill_I_face(eqn, d_outgoing_I[HI], d_mesh->num_cells_dim(I) - 1);
    fill_J_face(eqn, d_outgoing_J[HI], d_mesh->num_cells_dim(J) - 1);

    // send out data on low sides
    if (!d_outgoing_I[LO].is_null())
    {
        CHECK(!d_incoming_I[LO].is_null());
        CHECK(d_neighbor_I[LO] != PROBLEM_BOUNDARY);
        profugus::send(
            d_outgoing_I[LO]->data_pointer(), d_outgoing_I[LO]->data_size(),
            d_neighbor_I[LO], 452);
    }
    if (!d_outgoing_J[LO].is_null())
    {
        CHECK(!d_incoming_J[LO].is_null());
        CHECK(d_neighbor_J[LO] != PROBLEM_BOUNDARY);
        profugus::send(
            d_outgoing_J[LO]->data_pointer(), d_outgoing_J[LO]->data_size(),
            d_neighbor_J[LO], 453);
    }

    // send out data on high sides
    if (!d_outgoing_I[HI].is_null())
    {
        CHECK(!d_incoming_I[HI].is_null());
        CHECK(d_neighbor_I[HI] != PROBLEM_BOUNDARY);
        profugus::send(
            d_outgoing_I[HI]->data_pointer(), d_outgoing_I[HI]->data_size(),
            d_neighbor_I[HI], 450);
    }
    if (!d_outgoing_J[HI].is_null())
    {
        CHECK(!d_incoming_J[HI].is_null());
        CHECK(d_neighbor_J[HI] != PROBLEM_BOUNDARY);
        profugus::send(
            d_outgoing_J[HI]->data_pointer(), d_outgoing_J[HI]->data_size(),
            d_neighbor_J[HI], 451);
    }

    // wait on all of the posted receives
    d_request_I[LO].wait();
    d_request_I[HI].wait();
    d_request_J[LO].wait();
    d_request_J[HI].wait();

    ENSURE(!d_request_I[LO].inuse());
    ENSURE(!d_request_I[HI].inuse());
    ENSURE(!d_request_J[LO].inuse());
    ENSURE(!d_request_J[HI].inuse());
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get diffusion matrices from the low side neighbor.
 *
 * \param face I or J enumeration indicating face direction
 *
 * \return face field of diffusion coefficients from the low-side neighbor; it
 * could be unassigned if the face is adjacent to a problem boundary
 */
FV_Gather::RCP_Face_Field FV_Gather::low_side_D(int face) const
{
    using def::I; using def::J; using def::K;

    REQUIRE(face < K);

    // return the appropriate field
    if (face == I)
    {
        return d_incoming_I[LO];
    }
    return d_incoming_J[LO];
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get diffusion matrices from the low side neighbor.
 *
 * \param face I or J enumeration indicating face direction
 *
 * \return face field of diffusion coefficients from the low-side neighbor; it
 * could be unassigned if the face is adjacent to a problem boundary
 */
FV_Gather::RCP_Face_Field FV_Gather::high_side_D(int face) const
{
    using def::I; using def::J; using def::K;

    REQUIRE(face < K);

    // return the appropriate field
    if (face == I)
    {
        return d_incoming_I[HI];
    }
    return d_incoming_J[HI];
}

//---------------------------------------------------------------------------//
// PRIVATE IMPLEMENTATION
//---------------------------------------------------------------------------//
/*!
 * \brief Post receives.
 */
void FV_Gather::post_receives()
{
    REQUIRE(!d_request_I[LO].inuse());
    REQUIRE(!d_request_I[HI].inuse());
    REQUIRE(!d_request_J[LO].inuse());
    REQUIRE(!d_request_J[HI].inuse());

    // post receives on this block

    // low sides
    if (!d_incoming_I[LO].is_null())
    {
        CHECK(!d_outgoing_I[LO].is_null());
        profugus::receive_async(
            d_request_I[LO], d_incoming_I[LO]->data_pointer(),
            d_incoming_I[LO]->data_size(), d_neighbor_I[LO], 450);
    }
    if (!d_incoming_J[LO].is_null())
    {
        CHECK(!d_outgoing_J[LO].is_null());
        profugus::receive_async(
            d_request_J[LO], d_incoming_J[LO]->data_pointer(),
            d_incoming_J[LO]->data_size(), d_neighbor_J[LO], 451);
    }

    // high sides
    if (!d_incoming_I[HI].is_null())
    {
        CHECK(!d_outgoing_I[HI].is_null());
        profugus::receive_async(
            d_request_I[HI], d_incoming_I[HI]->data_pointer(),
            d_incoming_I[HI]->data_size(), d_neighbor_I[HI], 452);
    }
    if (!d_incoming_J[HI].is_null())
    {
        CHECK(!d_outgoing_J[HI].is_null());
        profugus::receive_async(
            d_request_J[HI], d_incoming_J[HI]->data_pointer(),
            d_incoming_J[HI]->data_size(), d_neighbor_J[HI], 453);
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Fill data on an I-face.
 */
void FV_Gather::fill_I_face(int            eqn,
                            RCP_Face_Field field,
                            int            i)
{
    using def::J; using def::K;

    // return if the field doesn't exist (meaning there is no data to transfer
    // here
    if (field.is_null()) return;

    REQUIRE(field->abscissa() == d_mesh->num_cells_dim(J));
    REQUIRE(field->ordinate() == d_mesh->num_cells_dim(K));

    // loop through cells and build the diffusion coefficients and assign them
    // to the face field
    for (int k = 0, Nk = field->ordinate(); k < Nk; ++k)
    {
        for (int j = 0, Nj = field->abscissa(); j < Nj; ++j)
        {
            // get a view of the block matrix at this location
            Serial_Matrix D = field->view(j, k);

            // calculate the diffusion coefficient in this cell
            d_coefficients->make_D(eqn, d_mesh->convert(i, j, k), D);
        }
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Fill data on an J-face.
 */
void FV_Gather::fill_J_face(int            eqn,
                            RCP_Face_Field field,
                            int            j)
{
    using def::I; using def::K;

    // return if the field doesn't exist (meaning there is no data to transfer
    // here
    if (field.is_null()) return;

    REQUIRE(field->abscissa() == d_mesh->num_cells_dim(I));
    REQUIRE(field->ordinate() == d_mesh->num_cells_dim(K));

    // loop through cells and build the diffusion coefficients and assign them
    // to the face field
    for (int k = 0, Nk = field->ordinate(); k < Nk; ++k)
    {
        for (int i = 0, Ni = field->abscissa(); i < Ni; ++i)
        {
            // get a view of the block matrix at this location
            Serial_Matrix D = field->view(i, k);

            // calculate the diffusion coefficient in this cell
            d_coefficients->make_D(eqn, d_mesh->convert(i, j, k), D);
        }
    }
}

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of FV_Gather.cc
//---------------------------------------------------------------------------//
