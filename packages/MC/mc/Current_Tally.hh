//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   MC/mc/Current_Tally.hh
 * \author Steven Hamilton
 * \date   Thu Apr 28 20:19:42 2016
 * \brief  Current_Tally class definition.
 * \note   Copyright (C) 2016 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef MC_mc_Current_Tally_hh
#define MC_mc_Current_Tally_hh

#include <memory>

#include "utils/View_Field.hh"
#include "Tally.hh"
#include "Physics.hh"
#include "Particle.hh"

namespace profugus
{

//===========================================================================//
/*!
 * \class Current_Tally
 * \brief Tally of net current on surfaces of Cartesian mesh.
 *
 * This class tallies the net energy-integrated current across the surfaces
 * of a Cartesian mesh.  These currents can be used to compute balance tables
 * or as part of acceleration methods such as CMFD.
 *
 * \sa Current_Tally.t.hh for detailed descriptions.
 */
/*!
 * \example mc/test/tstCurrent_Tally.cc
 *
 * Test of Current_Tally.
 */
//===========================================================================//

template <class Geometry>
class Current_Tally : public Surface_Tally<Geometry>
{
  public:

    typedef Surface_Tally<Geometry>                 Base;
    typedef Physics<Geometry>                       Physics_t;
    typedef Particle<Geometry>                      Particle_t;
    typedef std::shared_ptr<Physics_t>              SP_Physics;
    typedef Teuchos::RCP<Teuchos::ParameterList>    RCP_ParameterList;

    // Constructor
    Current_Tally(RCP_ParameterList          pl,
                  SP_Physics                 physics,
                  const std::vector<double> &x_edges,
                  const std::vector<double> &y_edges,
                  const std::vector<double> &z_edges);

    void tally_surface(const Particle_t &particle) override;

    void end_history() override;

    void finalize(double num_particles) override;

    void reset() override;

    // Get tally results
    profugus::const_View_Field<double> x_current() const
    {
        return profugus::make_view(d_x_current);
    }
    profugus::const_View_Field<double> y_current() const
    {
        return profugus::make_view(d_y_current);
    }
    profugus::const_View_Field<double> z_current() const
    {
        return profugus::make_view(d_z_current);
    }
    profugus::const_View_Field<double> x_current_std_dev() const
    {
        return profugus::make_view(d_x_current_std_dev);
    }
    profugus::const_View_Field<double> y_current_std_dev() const
    {
        return profugus::make_view(d_y_current_std_dev);
    }
    profugus::const_View_Field<double> z_current_std_dev() const
    {
        return profugus::make_view(d_z_current_std_dev);
    }
    profugus::const_View_Field<double> x_flux() const
    {
        return profugus::make_view(d_x_flux);
    }
    profugus::const_View_Field<double> y_flux() const
    {
        return profugus::make_view(d_y_flux);
    }
    profugus::const_View_Field<double> z_flux() const
    {
        return profugus::make_view(d_z_flux);
    }
    profugus::const_View_Field<double> x_flux_std_dev() const
    {
        return profugus::make_view(d_x_flux_std_dev);
    }
    profugus::const_View_Field<double> y_flux_std_dev() const
    {
        return profugus::make_view(d_y_flux_std_dev);
    }
    profugus::const_View_Field<double> z_flux_std_dev() const
    {
        return profugus::make_view(d_z_flux_std_dev);
    }

  private:

    std::string d_problem_name;

    // Edges defining tally mesh
    std::vector<double> d_x_edges;
    std::vector<double> d_y_edges;
    std::vector<double> d_z_edges;

    // Surface areas
    std::vector<double> d_x_areas;
    std::vector<double> d_y_areas;
    std::vector<double> d_z_areas;

    // Current on x, y, and z faces
    std::vector<double> d_x_current;
    std::vector<double> d_y_current;
    std::vector<double> d_z_current;
    std::vector<double> d_x_flux;
    std::vector<double> d_y_flux;
    std::vector<double> d_z_flux;

    // Std deviations
    std::vector<double> d_x_current_std_dev;
    std::vector<double> d_y_current_std_dev;
    std::vector<double> d_z_current_std_dev;
    std::vector<double> d_x_flux_std_dev;
    std::vector<double> d_y_flux_std_dev;
    std::vector<double> d_z_flux_std_dev;

    // Per-history accumulation
    std::vector<double> d_x_current_hist;
    std::vector<double> d_y_current_hist;
    std::vector<double> d_z_current_hist;
    std::vector<double> d_x_flux_hist;
    std::vector<double> d_y_flux_hist;
    std::vector<double> d_z_flux_hist;
};

} // end namespace profugus

#endif // MC_mc_Current_Tally_hh

//---------------------------------------------------------------------------//
//                 end of Current_Tally.hh
//---------------------------------------------------------------------------//
