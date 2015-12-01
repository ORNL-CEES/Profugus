//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc_driver/Problem_Builder.hh
 * \author Thomas M. Evans
 * \date   Wed Mar 12 22:25:22 2014
 * \brief  Problem_Builder class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef mc_driver_Problem_Builder_hh
#define mc_driver_Problem_Builder_hh

#include <string>
#include <vector>
#include <memory>
#include <unordered_map>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_Array.hpp"

#include "mc/Physics.hh"
#include "mc/Shape.hh"
#include "mc/Variance_Reduction.hh"
#include "mc/Tallier.hh"
#include "mc/Fission_Matrix_Acceleration.hh"

namespace mc
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
    typedef Geometry                                      Geom_t;
    typedef Teuchos::ParameterList                        ParameterList;
    typedef Teuchos::RCP<ParameterList>                   RCP_ParameterList;
    typedef profugus::Physics<Geom_t>                     Physics_t;
    typedef profugus::Tallier<Geom_t>                     Tallier_t;
    typedef profugus::Fission_Matrix_Acceleration<Geom_t> FM_Acceleration_t;
    typedef std::shared_ptr<Physics_t>                    SP_Physics;
    typedef std::shared_ptr<Geom_t>                       SP_Geometry;
    typedef std::shared_ptr<profugus::Shape>              SP_Shape;
    typedef profugus::Variance_Reduction<Geom_t>          VR_t;
    typedef std::shared_ptr<VR_t>                         SP_Var_Reduction;
    typedef std::shared_ptr<Tallier_t>                    SP_Tallier;
    typedef typename Tallier_t::SP_Tally                  SP_Tally;
    typedef std::shared_ptr<FM_Acceleration_t>            SP_FM_Acceleration;
    //@}

  private:
    // >>> DATA

    // Problem-parameterlist (talks to solver components).
    RCP_ParameterList d_db;

    // Physics and geometry.
    SP_Physics  d_physics;
    SP_Geometry d_geometry;

    // Variance reduction.
    SP_Var_Reduction d_var_reduction;

    // Fission matrix acceleration.
    SP_FM_Acceleration d_fm_acceleration;

    // External source shape.
    SP_Shape d_shape;

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
    SP_Geometry get_geometry() const { return d_geometry; }

    //! Get the physics.
    SP_Physics get_physics() const { return d_physics; }

    //! Get the external source shape (could be null).
    SP_Shape get_source_shape() const { return d_shape; }

    //! Get the variance reduction.
    SP_Var_Reduction get_var_reduction() const { return d_var_reduction; }

    //! Get the tallier.
    SP_Tallier get_tallier() const { return d_tallier; }

    //! Get the fission matrix acceleration.
    SP_FM_Acceleration get_acceleration() const { return d_fm_acceleration; }

  private:
    // >>> IMPLEMENTATION

    // Teuchos typedefs.
    typedef Teuchos::Array<int>         OneDArray_int;
    typedef Teuchos::Array<double>      OneDArray_dbl;
    typedef Teuchos::Array<std::string> OneDArray_str;

    // Acceleration typedefs.
    typedef typename FM_Acceleration_t::Problem_Builder_t SPN_Builder;

    // Build implementation.
    void build_geometry(RCP_ParameterList master);
    void build_physics();
    void build_var_reduction();
    void build_source(const ParameterList &source_db);
    void build_tallies();
    void build_spn_problem();

    RCP_ParameterList d_matdb;
};

} // end namespace mc

#endif // mc_driver_Problem_Builder_hh

//---------------------------------------------------------------------------//
//                 end of Problem_Builder.hh
//---------------------------------------------------------------------------//
