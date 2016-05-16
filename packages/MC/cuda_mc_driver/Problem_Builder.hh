//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc_driver/Problem_Builder.hh
 * \author Thomas M. Evans
 * \date   Wed Mar 12 22:25:22 2014
 * \brief  Problem_Builder class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_mc_driver_Problem_Builder_hh
#define cuda_mc_driver_Problem_Builder_hh

#include <string>
#include <vector>
#include <memory>
#include <unordered_map>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_Array.hpp"

#include "cuda_mc/Physics.hh"
#include "cuda_mc/Box_Shape.hh"
#include "cuda_mc/Variance_Reduction.hh"
#include "cuda_mc/Tallier.hh"

#include "cuda_utils/Shared_Device_Ptr.hh"

namespace cuda_mc
{

//===========================================================================//
/*!
 * \class Problem_Builder
 * \brief Read and initialize an mc problem,
 */
//===========================================================================//

template <class Geometry>
class Problem_Builder
{
  public:
    //@{
    //! Typedefs.
    typedef Geometry                                          Geom_t;
    typedef Teuchos::ParameterList                            ParameterList;
    typedef Teuchos::RCP<ParameterList>                       RCP_ParameterList;
    typedef profugus::Physics<Geom_t>                         Physics_t;
    typedef profugus::Tallier<Geom_t>                         Tallier_t;
    typedef cuda::Shared_Device_Ptr<Physics_t>                SDP_Physics;
    typedef cuda::Shared_Device_Ptr<Geom_t>                   SDP_Geometry;
    typedef cuda::Shared_Devoce_Ptr<cuda_profugus::Box_Shape> SDP_Shape;
    typedef profugus::Variance_Reduction<Geom_t>              VR_t;
    typedef std::shared_ptr<VR_t>                             SP_Var_Reduction;
    typedef std::shared_ptr<Tallier_t>                        SP_Tallier;
    typedef typename Tallier_t::SP_Tally                      SP_Tally;
    //@}

  private:
    // >>> DATA

    // Problem-parameterlist (talks to solver components).
    RCP_ParameterList d_db;

    // Physics and geometry.
    SDP_Physics  d_physics;
    SDP_Geometry d_geometry;

    // Variance reduction.
    SP_Var_Reduction d_var_reduction;

    // External source shape.
    SDP_Shape d_shape;

    // Problem talliers.
    SP_Tallier d_tallier;

  public:
    // Constructor.
    Problem_Builder();

    // Setup the problem.
    void setup(RCP_ParameterList master);

    // >>> ACCESSORS

    //! Get problem database.
    RCP_ParameterList problem_db() const { return d_db; }

    //! Get the geometry.
    SDP_Geometry get_geometry() const { return d_geometry; }

    //! Get the physics.
    SDP_Physics get_physics() const { return d_physics; }

    //! Get the external source shape (could be null).
    SDP_Shape get_source_shape() const { return d_shape; }

    //! Get the variance reduction.
    SP_Var_Reduction get_var_reduction() const { return d_var_reduction; }

    //! Get the tallier.
    SP_Tallier get_tallier() const { return d_tallier; }

  private:
    // >>> IMPLEMENTATION

    // Teuchos typedefs.
    typedef Teuchos::Array<int>         OneDArray_int;
    typedef Teuchos::Array<double>      OneDArray_dbl;
    typedef Teuchos::Array<std::string> OneDArray_str;

    // Build implementation.
    void build_geometry(RCP_ParameterList master);
    void build_physics();
    void build_var_reduction();
    void build_source(const ParameterList &source_db);
    void build_tallies();
    void build_spn_problem();

    RCP_ParameterList d_matdb;
};

} // end namespace cuda_mc

#endif // cuda_mc_driver_Problem_Builder_hh

//---------------------------------------------------------------------------//
//                 end of Problem_Builder.hh
//---------------------------------------------------------------------------//
