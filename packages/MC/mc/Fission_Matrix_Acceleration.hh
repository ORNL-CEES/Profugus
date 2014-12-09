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
#include <vector>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ArrayView.hpp"
#include "AnasaziMultiVecTraits.hpp"

#include "xs/Mat_DB.hh"
#include "mesh/Mesh.hh"
#include "spn/Linear_System.hh"
#include "spn/VectorTraits.hh"
#include "spn_driver/Problem_Builder.hh"
#include "Fission_Matrix_Solver.hh"
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
    typedef Teuchos::RCP<const Vector_t>             RCP_Const_Vector;
    typedef Teuchos::RCP<Fission_Matrix_Solver<T>>   RCP_Fission_Matrix_Solver;

  private:
    // >>> DATA

    // Problem dimensions.
    RCP_Dimensions d_dim;

    // Original SPN linear system.
    RCP_Linear_System d_system;

    // Fission matrix acceleration equation solver.
    RCP_Fission_Matrix_Solver d_fm_solver;

    // Forward and adjoint solutions to the SPN problem.
    RCP_Vector d_adjoint;
    RCP_Vector d_forward;

    // Eigenvalue of the SPN system.
    double d_keff;

    // MC eigenvalue at the beginning of cycle (l).
    double d_k_l;

    // Damping.
    double d_beta;

    // Multiplicative corrections.
    std::vector<double> d_nu;

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

    //! Get the fission matrix solver.
    const Fission_Matrix_Solver<T>& solver() const { return *d_fm_solver; }

    //! Get the eigenvalue of the SPN system.
    double keff() const { return d_keff; }

    //@{
    //! Get the eigenvectors of the SPN system.
    Const_Array_View adjoint() const
    {
        return VectorTraits<T>::get_data(d_adjoint);
    }
    Const_Array_View forward() const
    {
        return VectorTraits<T>::get_data(d_forward);
    }
    //@}

  private:
    // >>> IMPLEMENTATION

    typedef VectorTraits<T>                                VTraits;
    typedef typename T::MV                                 MultiVector_t;
    typedef typename T::OP                                 Operator_t;
    typedef Anasazi::MultiVecTraits<double, MultiVector_t> ATraits;

    // Convert g to fissions.
    void convert_g(RCP_Const_Vector g);
};

} // end namespace profugus

#endif // mc_Fission_Matrix_Acceleration_hh

//---------------------------------------------------------------------------//
//                 end of Fission_Matrix_Acceleration.hh
//---------------------------------------------------------------------------//
