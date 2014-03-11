//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   spn/Solver_Base.cc
 * \author Thoams M. Evans
 * \date   Mon Feb 17 20:57:20 2014
 * \brief  Solver_Base member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "harness/DBC.hh"
#include "Moment_Coefficients.hh"
#include "Solver_Base.hh"

namespace profugus
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR/DESTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
Solver_Base::Solver_Base(RCP_ParameterList db)
    : b_db(db)
{
    Require (!b_db.is_null());
}

//---------------------------------------------------------------------------//
/*!
 * \brief Null virtual destructor.
 */
Solver_Base::~Solver_Base()
{
}

//---------------------------------------------------------------------------//
// PROTECTED METHODS
//---------------------------------------------------------------------------//
/*!
 * \brief Write a u-vector into the state.
 */
void Solver_Base::write_u_into_state(const Vector_t &u,
                                     State_t        &state)
{
    Require (state.num_cells() * b_system->get_dims()->num_equations()
             * state.num_groups() <= u.MyLength());
    Require (state.num_cells() == state.mesh().num_cells());

    // SPN moments
    double u_m[4] = {0.0, 0.0, 0.0, 0.0};

    // number of equations (moments) in the SPN solution
    int N = b_system->get_dims()->num_equations();

    // number of groups
    int Ng = state.num_groups();

    // number of local cells
    int Nc = state.num_cells();

    // the state field is ordered group->cell whereas the u vector is
    // ordered cell->moments->group

    // loop over groups
    for (int g = 0; g < Ng; ++g)
    {
        // get all the fluxes for this group
        View_Field flux = state.flux(g, g);
        Check (flux.size() == Nc);

        // loop over cells on this domain
        for (int cell = 0; cell < Nc; ++cell)
        {
            // loop over moments
            for (int n = 0; n < N; ++n)
            {
                Check (b_system->index(g, n, cell) < u.MyLength());

                // get the moment from the solution vector
                u_m[n] = u[b_system->index(g, n, cell)];
            }

            // assign the scalar flux (SPN Phi_0 moment) to the state ->
            // CURRENTLY WE ONLY SUPPORT 1 UNKNOWN PER CELL
            flux[cell] =  Moment_Coefficients::u_to_phi(
                u_m[0], u_m[1], u_m[2], u_m[3]);
        }
    }
}

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of Solver_Base.cc
//---------------------------------------------------------------------------//
