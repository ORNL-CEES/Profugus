//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc_driver/Geometry_Builder.hh
 * \author Steven Hamilton
 * \date   Wed Nov 25 12:58:58 2015
 * \brief  Geometry_Builder class definition.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef mc_driver_Geometry_Builder_hh
#define mc_driver_Geometry_Builder_hh

#include <memory>
#include <unordered_map>

#include "geometry/RTK_Geometry.hh"
#include "geometry/Mesh_Geometry.hh"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_TwoDArray.hpp"

namespace mc
{

//===========================================================================//
/*!
 * \class Geometry_Builder
 * \brief Build profugus Geometry
 */
//===========================================================================//
template <class Geometry>
class Geometry_Builder
{
  public:

    typedef std::shared_ptr<Geometry>             SP_Geometry;
    typedef Teuchos::RCP<Teuchos::ParameterList>  RCP_ParameterList;

    Geometry_Builder()
    {
        VALIDATE(false,"Missing a specialization");
    }

    SP_Geometry build(RCP_ParameterList master)
    {
        return SP_Geometry();
    }

};

// Specialization for RTK
template <>
class Geometry_Builder<profugus::Core>
{
  public:

    typedef profugus::Core                          Geom_t;
    typedef std::shared_ptr<profugus::Core>         SP_Geometry;
    typedef Teuchos::RCP<Teuchos::ParameterList>    RCP_ParameterList;

    Geometry_Builder(){};

    SP_Geometry build(RCP_ParameterList master);

  private:

    // Geometry typedefs.
    typedef Geom_t::Array_t      Core_t;
    typedef Geom_t::SP_Array     SP_Core;
    typedef Core_t::SP_Object    SP_Lattice;
    typedef Core_t::Object_t     Lattice_t;
    typedef Lattice_t::SP_Object SP_Pin_Cell;
    typedef Lattice_t::Object_t  Pin_Cell_t;

    typedef Teuchos::Array<int>         OneDArray_int;
    typedef Teuchos::Array<double>      OneDArray_dbl;
    typedef Teuchos::Array<std::string> OneDArray_str;
    typedef Teuchos::TwoDArray<int>     TwoDArray_int;
    typedef Teuchos::TwoDArray<double>  TwoDArray_dbl;

    // General typedefs.
    typedef std::unordered_map<int, SP_Lattice>  Lattice_Hash;
    typedef std::unordered_map<int, SP_Pin_Cell> Pin_Hash;

    RCP_ParameterList d_db;
    RCP_ParameterList d_coredb;
    RCP_ParameterList d_assblydb;
    RCP_ParameterList d_pindb;

    SP_Lattice build_axial_lattice(const TwoDArray_int &map, double height);
};

// Specialization for Mesh_Geometry
template <>
class Geometry_Builder<profugus::Mesh_Geometry>
{
  public:

    typedef std::shared_ptr<profugus::Mesh_Geometry>  SP_Geometry;
    typedef Teuchos::RCP<Teuchos::ParameterList>      RCP_ParameterList;
    typedef Teuchos::Array<int>                       OneDArray_int;
    typedef Teuchos::Array<double>                    OneDArray_dbl;

    SP_Geometry build(RCP_ParameterList master);
};

} // end namespace mc

#endif // mc_driver_Geometry_Builder_hh

//---------------------------------------------------------------------------//
//                 end of Geometry_Builder.hh
//---------------------------------------------------------------------------//
