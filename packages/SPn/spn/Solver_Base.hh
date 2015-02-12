//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   spn/Solver_Base.hh
 * \author Thomas M. Evans
 * \date   Mon Feb 17 20:36:54 2014
 * \brief  Solver_Base class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef spn_Solver_Base_hh
#define spn_Solver_Base_hh

#include "Teuchos_RCP.hpp"

#include "State.hh"
#include "Linear_System.hh"
#include "VectorTraits.hh"
#include "Isotropic_Source.hh"
#include "solvers/LinAlgTypedefs.hh"

namespace profugus
{

//===========================================================================//
/*!
 * \class Solver_Base
 * \brief Base class for SPN solvers (fixed-source and eigenvalue).
 */
/*!
 * \example spn/test/tstSolver_Base.cc
 *
 * Test of Solver_Base.
 */
//===========================================================================//

class Solver_Base
{
  public:
    //@{
    //! Typedefs.
    typedef Teuchos::RCP<Dimensions>           RCP_Dimensions;
    typedef Teuchos::RCP<Mat_DB>               RCP_Mat_DB;
    typedef Teuchos::RCP<Mesh>                 RCP_Mesh;
    typedef Teuchos::RCP<LG_Indexer>           RCP_Indexer;
    typedef Teuchos::RCP<Global_Mesh_Data>     RCP_Global_Data;
    typedef Isotropic_Source                   External_Source;
    typedef State::View_Field                  View_Field;
    //@}

  public:
    // Constructor.
    Solver_Base(){};

    // Virtual destructor.
    virtual ~Solver_Base(){};

    // Set up the solver.
    virtual void setup(RCP_Dimensions dim, RCP_Mat_DB mat, RCP_Mesh mesh,
                       RCP_Indexer indexer, RCP_Global_Data data) = 0;

    // Solve problem
    virtual void solve(Teuchos::RCP<const External_Source> q) = 0;

    //! Write results of solve into state.
    virtual void write_state(State &state) = 0;

    //! Write problem matrix (or matrices to file)
    virtual void write_problem_to_file() const = 0;

};

//===========================================================================//
/*!
 * \class Solver_Base_Tmpl
 * \brief Templated base class for SPN solvers (fixed-source and eigenvalue).
 */
//===========================================================================//
template <class T>
class Solver_Base_Tmpl : public Solver_Base
{
  protected:

    using Solver_Base::RCP_Dimensions;
    using Solver_Base::RCP_Mat_DB;
    using Solver_Base::RCP_Mesh;
    using Solver_Base::RCP_Indexer;
    using Solver_Base::RCP_Global_Data;
    using Solver_Base::View_Field;

    typedef typename T::MV                MV;
    typedef Teuchos::ParameterList        ParameterList;
    typedef Teuchos::RCP<ParameterList>   RCP_ParameterList;
    typedef Linear_System<T>              Linear_System_t;
    typedef Teuchos::RCP<Linear_System_t> RCP_Linear_System;

    // Linear system.
    RCP_Linear_System b_system;

    // Problem database.
    RCP_ParameterList b_db;

  public:

    Solver_Base_Tmpl(RCP_ParameterList db)
    {
        b_db = db;
    }

    virtual ~Solver_Base_Tmpl(){};

    const Linear_System_t &get_linear_system() const {return *b_system;}

    // >>> INHERITED METHODS

    // Write u-vector into the state.
    void write_u_into_state(Teuchos::RCP<const MV>, State &state) const;

};

//===========================================================================//
/*!
 * \brief Write a u-vector into the state.
 */
template <class T>
void Solver_Base_Tmpl<T>::write_u_into_state(Teuchos::RCP<const MV>  u,
                                             State &state) const
{
    REQUIRE(state.num_cells() * b_system->get_dims()->num_equations()
             * state.num_groups() <=
             VectorTraits<T>::local_length(u));
    REQUIRE(state.num_cells() == state.mesh().num_cells());

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

    Teuchos::ArrayView<const double> data = VectorTraits<T>::get_data(u);

    // loop over groups
    for (int g = 0; g < Ng; ++g)
    {
        // get all the fluxes for this group
        View_Field flux = state.flux(g, g);
        CHECK(flux.size() == Nc);

        // loop over cells on this domain
        for (int cell = 0; cell < Nc; ++cell)
        {
            // loop over moments
            for (int n = 0; n < N; ++n)
            {
                CHECK(b_system->index(g, n, cell) <
                      VectorTraits<T>::local_length(u) );

                // get the moment from the solution vector
                u_m[n] = data[b_system->index(g, n, cell)];
            }

            // assign the scalar flux (SPN Phi_0 moment) to the state ->
            // CURRENTLY WE ONLY SUPPORT 1 UNKNOWN PER CELL
            flux[cell] =  Moment_Coefficients::u_to_phi(
                u_m[0], u_m[1], u_m[2], u_m[3]);
        }
    }
}

} // end namespace profugus

#endif // spn_Solver_Base_hh

//---------------------------------------------------------------------------//
//                 end of Solver_Base.hh
//---------------------------------------------------------------------------//
