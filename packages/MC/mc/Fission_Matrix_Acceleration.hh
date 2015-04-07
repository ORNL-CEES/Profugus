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
#include "Teuchos_ArrayRCP.hpp"
#include "AnasaziMultiVecTraits.hpp"

#include "utils/Serial_HDF5_Writer.hh"
#include "xs/Mat_DB.hh"
#include "mesh/Mesh.hh"
#include "mesh/Global_Mesh_Data.hh"
#include "geometry/Cartesian_Mesh.hh"
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
    typedef std::shared_ptr<Cartesian_Mesh>      SP_Cart_Mesh;

  protected:
    // >>> DATA

    // Problem database.
    RCP_ParameterList b_db;

    // Mesh objects
    RCP_Mesh        b_mesh;
    RCP_Indexer     b_indexer;
    RCP_Global_Data b_gdata;

    // Global SPN mesh (needed to map partitioned SPN mesh back to domain
    // replicated fission site locations).
    SP_Cart_Mesh b_global_mesh;

    // Material database.
    RCP_Mat_DB b_mat;

    // Number of global sites after the latest correction.
    int b_current_Np;

    // Update fission sites.
    bool b_update_Np;

  public:
    // Constructor.
    Fission_Matrix_Acceleration() : b_current_Np(0) { /* ... */ }

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

    //! Return the number of fission sites after acceleration.
    int num_sites() const { return b_current_Np; }

    //! Update the number of source particles each cycle.
    bool update_num_particles() const { return b_update_Np; }

    //! Get SPN mesh.
    const Mesh& mesh() const { return *b_mesh; }

    //! Get global mesh.
    const Cartesian_Mesh& global_mesh() const { return *b_global_mesh; }

    //! Get SPN materials.
    const Mat_DB &mat() const { return *b_mat; }

    // >>> DIAGNOSTICS
#ifdef USE_HDF5
    //! Add diagnostics to HDF5 output.
    virtual void diagnostics(Serial_HDF5_Writer &writer) const = 0;
#endif
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
    typedef typename Linear_System_t::MV             Vector_t;
    typedef typename Linear_System_t::RCP_MV         RCP_Vector;
    typedef Teuchos::ArrayView<const double>         Const_Array_View;
    typedef Teuchos::ArrayRCP<const double>          Const_ArrayRCP;
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
    Const_ArrayRCP adjoint() const
    {
        return VectorTraits<T>::get_data(d_adjoint);
    }
    Const_ArrayRCP forward() const
    {
        return VectorTraits<T>::get_data(d_forward);
    }
    //@}

    //! Get the current fission density.
    Const_Array_View multiplicative_correction() const
    {
        return Teuchos::arrayViewFromVector(d_nu);
    }

    // >>> DIAGNOSTICS
#ifdef USE_HDF5
    //! Add diagnostics to HDF5 output.
    void diagnostics(Serial_HDF5_Writer &writer) const;
#endif

  private:
    // >>> IMPLEMENTATION

    typedef VectorTraits<T>                                VTraits;
    typedef typename T::MV                                 MultiVector_t;
    typedef typename T::OP                                 Operator_t;
    typedef Anasazi::MultiVecTraits<double, MultiVector_t> ATraits;

    // Convert g to fissions.
    void convert_g(RCP_Const_Vector g);

    // L2 norm of correction.
    std::vector<double> d_norms;

    // Maximum correction.
    std::vector<double> d_correction;

    // Iteration record of solves.
    std::vector<int> d_iterations;

    // Number of global fission sites after each correction.
    std::vector<int> d_Np;

    // Cycle counters.
    int  d_cycle_ctr, d_cycle_begin, d_cycle_end;
    bool d_accelerate;
};

} // end namespace profugus

#endif // mc_Fission_Matrix_Acceleration_hh

//---------------------------------------------------------------------------//
//                 end of Fission_Matrix_Acceleration.hh
//---------------------------------------------------------------------------//
