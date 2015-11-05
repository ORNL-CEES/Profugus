//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/Fission_Matrix_Solver.hh
 * \author Thomas M. Evans
 * \date   Mon Dec 08 17:18:34 2014
 * \brief  Fission_Matrix_Solver class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef mc_Fission_Matrix_Solver_hh
#define mc_Fission_Matrix_Solver_hh

#include <memory>

#include "Teuchos_RCP.hpp"
#include "AnasaziMultiVecTraits.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ArrayView.hpp"

#include "geometry/Cartesian_Mesh.hh"
#include "solvers/LinAlgTypedefs.hh"
#include "solvers/ShiftedOperator.hh"
#include "solvers/LinearSolver.hh"
#include "spn/Linear_System.hh"
#include "spn/Isotropic_Source.hh"
#include "spn/VectorTraits.hh"
#include "Physics.hh"

// Remove this once templated
#include "geometry/RTK_Geometry.hh"

namespace profugus
{

//===========================================================================//
/*!
 * \class Fission_Matrix_Solver
 * \brief Solve the fission matrix acceleration correction equation.
 *
 * The correction equation is solved in SPN space.
 */
//===========================================================================//

template<class T>
class Fission_Matrix_Solver
{
  public:
    //@{
    //! Typedefs.
    typedef ShiftedOperator<T>                          ShiftedOperator_t;
    typedef Teuchos::RCP<ShiftedOperator_t>             RCP_ShiftedOperator;
    typedef Teuchos::RCP<LinearSolver<T>>               RCP_LinearSolver;
    typedef Physics<Core>                               Physics_t;
    typedef typename Physics_t::Fission_Site            Fission_Site;
    typedef typename Physics_t::Fission_Site_Container  Fission_Site_Container;
    typedef profugus::Linear_System<T>                  Linear_System_t;
    typedef Teuchos::RCP<Linear_System_t>               RCP_Linear_System;
    typedef typename Linear_System_t::MV                Vector_t;
    typedef typename Linear_System_t::RCP_MV            RCP_Vector;
    typedef Teuchos::RCP<const Vector_t>                RCP_Const_Vector;
    typedef Teuchos::ParameterList                      ParameterList_t;
    typedef Teuchos::RCP<ParameterList_t>               RCP_ParameterList;
    typedef typename Linear_System_t::RCP_Mesh          RCP_Mesh;
    typedef typename Linear_System_t::RCP_Mat_DB        RCP_Mat_DB;
    typedef typename Linear_System_t::RCP_Indexer       RCP_Indexer;
    typedef Teuchos::ArrayView<const double>            Const_Array_View;
    typedef std::shared_ptr<Cartesian_Mesh>             SP_Cart_Mesh;
    //@}

  private:
    // >>> DATA

    // SPN mesh.
    RCP_Mesh    d_mesh;
    RCP_Indexer d_indexer;

    // SPN materials.
    RCP_Mat_DB d_mat;

    // Original SPN linear system.
    RCP_Linear_System d_system;

    // Acceleration equation operator.
    RCP_ShiftedOperator d_operator;

    // Linear solver for the acceleration equations.
    RCP_LinearSolver d_solver;

    // Global SPN mesh (needed to map partitioned SPN mesh back to domain
    // replicated fission site locations).
    SP_Cart_Mesh d_global_mesh;

    // Correction vector (dimensioned on the partitioned mesh)
    RCP_Vector d_g;

  public:
    //! Constructor.
    Fission_Matrix_Solver(RCP_ParameterList fm_db, RCP_Mesh mesh,
                          RCP_Indexer indexer, RCP_Mat_DB mat,
                          RCP_Linear_System system, SP_Cart_Mesh global_mesh,
                          double keff);

    // Set the eigenvectors of the SPN system.
    void set_eigenvectors(RCP_Vector forward, RCP_Vector adjoint);

    // Set the fission source at the beginning of the cycle.
    void set_u_begin(const Fission_Site_Container &f, double k_l);

    // Set the correction at the end of the cycle.
    void solve(const Fission_Site_Container &f);

    // >>> ACCESSORS

    //! Get the correction vector.
    RCP_Const_Vector get_g() const { return d_g; }

    //! Get the current fission density.
    Const_Array_View current_f() const
    {
        return Teuchos::arrayViewFromVector(d_q_field);
    }

    //! Number of iterations of last solve.
    int iterations() const { return d_solver->num_iters(); }

  private:
    // >>> IMPLEMENTATION

    typedef VectorTraits<T>                                VTraits;
    typedef typename T::MV                                 MultiVector_t;
    typedef typename T::OP                                 Operator_t;
    typedef Anasazi::MultiVecTraits<double, MultiVector_t> ATraits;

    // Setup the solver options.
    void solver_db(RCP_ParameterList mc_db);

    // Build the Bphi mat-vec product from the fission sites.
    RCP_Vector build_Bphi(const Fission_Site_Container &f);

    // External source used to build the action of B\phi.
    std::shared_ptr<Isotropic_Source> d_q;

    // Source shapes (chi), and field.
    Isotropic_Source::Source_Shapes d_q_shapes;
    Isotropic_Source::Source_Field  d_q_field;

    // Global fission source.
    Isotropic_Source::Source_Field d_global_q_field;

    // Forward and adjoint eigenvectors.
    RCP_Vector d_forward, d_adjoint;

    // MC eigenvalue at l.
    double d_k_l;

    //  Work vectors for rhs.
    RCP_Vector d_work;
};

} // end namespace profugus

#endif // mc_Fission_Matrix_Solver_hh

//---------------------------------------------------------------------------//
//                 end of Fission_Matrix_Solver.hh
//---------------------------------------------------------------------------//
