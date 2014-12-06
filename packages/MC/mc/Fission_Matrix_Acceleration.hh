//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/Fission_Matrix_Acceleration.hh
 * \author Thomas M. Evans
 * \date   Wed Nov 12 14:54:34 2014
 * \brief  Fission_Matrix_Acceleration class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef mc_Fission_Matrix_Acceleration_hh
#define mc_Fission_Matrix_Acceleration_hh

#include <memory>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ArrayView.hpp"
#include "AnasaziMultiVecTraits.hpp"

#include "harness/DBC.hh"
#include "xs/Mat_DB.hh"
#include "mesh/Mesh.hh"
#include "solvers/LinAlgTypedefs.hh"
#include "solvers/ShiftedOperator.hh"
#include "solvers/LinearSolver.hh"
#include "spn/Linear_System.hh"
#include "spn/VectorTraits.hh"
#include "spn/Isotropic_Source.hh"
#include "spn_driver/Problem_Builder.hh"
#include "Physics.hh"

namespace profugus
{

//===========================================================================//
/*!
 * \class Fission_Matrix_Acceleration
 * \brief Do Fission Matrix acceleration of fission source.
 */
/*!
 * \example mc/test/tstFission_Matrix_Acceleration.cc
 *
 * Test of Fission_Matrix_Acceleration.
 */
//===========================================================================//

class Fission_Matrix_Acceleration
{
  public:
    // Typedefs.
    typedef spn::Problem_Builder                 Problem_Builder_t;
    typedef Problem_Builder_t::ParameterList     ParameterList;
    typedef Problem_Builder_t::RCP_ParameterList RCP_ParameterList;
    typedef Problem_Builder_t::RCP_Mesh          RCP_Mesh;
    typedef Problem_Builder_t::RCP_Indexer       RCP_Indexer;
    typedef Problem_Builder_t::RCP_Global_Data   RCP_Global_Data;
    typedef Problem_Builder_t::RCP_Mat_DB        RCP_Mat_DB;
    typedef Physics::Fission_Site                Fission_Site;
    typedef Physics::Fission_Site_Container      Fission_Site_Container;

  protected:
    // >>> DATA

    // Problem database.
    RCP_ParameterList b_db;

    // Mesh objects
    RCP_Mesh        b_mesh;
    RCP_Indexer     b_indexer;
    RCP_Global_Data b_gdata;

    // Material database.
    RCP_Mat_DB b_mat;

  public:
    // Constructor.
    Fission_Matrix_Acceleration() { /* ... */ }

    // Destructor.
    virtual ~Fission_Matrix_Acceleration() { /* ... */ }

    //! Build the SPN problem.
    virtual void build_problem(const Problem_Builder_t &builder) = 0;

    //! Initialize the acceleration.
    virtual void initialize(RCP_ParameterList mc_db) = 0;

    //! Start cycle initialization with fission container at \f$l\f$.
    virtual void start_cycle(double k_l, const Fission_Site_Container &f) = 0;

    //! End cycle acceleration to create fission container at \f$l+1\f$.
    virtual void end_cycle(Fission_Site_Container &f) = 0;

    // >>> ACCESSORS

    //! Get SPN mesh.
    const Mesh& mesh() const { return *b_mesh; }

    //! Get SPN materials.
    const Mat_DB &mat() const { return *b_mat; }
};

//===========================================================================//
/*!
 * \class Fission_Matrix_Acceleration_Impl
 * \brief Do implementation of fission matrix acceleration of fission source.
 */
//===========================================================================//

template<class T>
class Fission_Matrix_Acceleration_Impl : public Fission_Matrix_Acceleration
{
  public:
    // Typedefs.
    typedef profugus::Linear_System<T>               Linear_System_t;
    typedef typename Linear_System_t::RCP_Dimensions RCP_Dimensions;
    typedef Teuchos::RCP<Linear_System_t>            RCP_Linear_System;
    typedef typename Linear_System_t::RCP_Operator   RCP_Operator;
    typedef typename Linear_System_t::RCP_Matrix     RCP_Matrix;
    typedef typename Linear_System_t::Vector_t       Vector_t;
    typedef typename Linear_System_t::RCP_Vector     RCP_Vector;
    typedef Teuchos::ArrayView<const double>         Const_Array_View;
    typedef ShiftedOperator<T>                       ShiftedOperator_t;
    typedef Teuchos::RCP<ShiftedOperator_t>          RCP_ShiftedOperator;
    typedef Teuchos::RCP<LinearSolver<T>>            RCP_LinearSolver;
    typedef Teuchos::RCP<const Vector_t>             RCP_Const_Vector;

  private:
    // >>> DATA

    // Problem dimensions.
    RCP_Dimensions d_dim;

    // Original SPN linear system.
    RCP_Linear_System d_system;

    // djoint solutions to the SPN problem.
    RCP_Vector d_adjoint;

    // Keff from the SPN solve.
    double d_keff;

    // Acceleration equation operator.
    RCP_ShiftedOperator d_operator;

    // Linear solver for the acceleration equations.
    RCP_LinearSolver d_solver;

    // Work vectors for rhs.
    RCP_Vector d_g, d_work;

    // MC eigenvalue at the beginning of cycle (l).
    double d_k_l;

  public:
    // Constructor.
    Fission_Matrix_Acceleration_Impl();

    // Build the SPN problem.
    void build_problem(const Problem_Builder_t &builder);

    // Initialize the acceleration.
    void initialize(RCP_ParameterList mc_db);

    // Start cycle initialization.
    void start_cycle(double k_l, const Fission_Site_Container &f);

    // End cycle acceleration.
    void end_cycle(Fission_Site_Container &f);

    // >>> ACCESSORS

    //! Get the linear system of SPN equations.
    const Linear_System_t& spn_system() const { return *d_system; }

    //! Get the eigenvalue of the SPN system.
    double keff() const { return d_keff; }

    //! Get the adjoint eigenvector of the SPN system.
    Const_Array_View adjoint() const
    {
        return VectorTraits<T>::get_data(d_adjoint);
    }

    //! Get the vector at \f$ l\f$.
    RCP_Const_Vector get_Bpsi_l() const { return d_work; }

    //! Get the correction vector.
    RCP_Const_Vector get_g() const { return d_g; }

  private:
    // >>> IMPLEMENTATION

    typedef VectorTraits<T>                                VTraits;
    typedef typename T::MV                                 MultiVector_t;
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
};

} // end namespace profugus

#endif // mc_Fission_Matrix_Acceleration_hh

//---------------------------------------------------------------------------//
//                 end of Fission_Matrix_Acceleration.hh
//---------------------------------------------------------------------------//
