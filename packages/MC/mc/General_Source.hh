//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/General_Source.hh
 * \author Steven Hamilton
 * \date   Mon Apr 04 20:38:12 2016
 * \brief  General_Source class definition.
 * \note   Copyright (C) 2016 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef mc_General_Source_hh
#define mc_General_Source_hh

#include <memory>

#include "Source.hh"
#include "Physics.hh"
#include "Particle.hh"
#include "rng/RNG_Control.hh"

namespace profugus
{

//===========================================================================//
/*!
 * \class General_Source
 * \brief Defines a space- and energy-dependent source.
 */
/*!
 * \example mc/test/tstGeneral_Source.cc
 *
 * Test of General_Source.
 */
//===========================================================================//

template <class Geometry>
class General_Source : public Source<Geometry>
{
  public:

      typedef Source<Geometry>                      Base;
      typedef Physics<Geometry>                     Physics_t;
      typedef Particle<Geometry>                    Particle_t;
      typedef std::shared_ptr<Geometry>             SP_Geometry;
      typedef std::shared_ptr<Physics_t>            SP_Physics;
      typedef std::shared_ptr<Particle_t>           SP_Particle;
      typedef std::shared_ptr<RNG_Control>          SP_RNG_Control;
      typedef def::Space_Vector                     Space_Vector;
      typedef Teuchos::RCP<Teuchos::ParameterList>  RCP_ParameterList;

      General_Source(RCP_ParameterList  db,
                     SP_Geometry        geometry,
                     SP_Physics         physics,
                     SP_RNG_Control     rng_control);

      void build_source(Teuchos::TwoDArray<double> &source);

      // >> VIRTUAL PUBLIC INTERFACE

      //! Is source finished with all particles
      bool empty() const override
      {
          return (d_np_left == 0);
      }

      //! Number of particles run so far
      def::size_type num_run() const
      {
          return d_np_run;
      }

      //! Number of particles left to run
      def::size_type num_left() const
      {
          return d_np_left;
      }

      //! Number of particles to transport on this domain.
      def::size_type num_to_transport() const
      {
          return d_np_domain;
      }

      //! Total number of particles to transport on all domains.
      def::size_type total_num_to_transport() const
      {
          return d_np_total;
      }

      //! Get particle from source
      SP_Particle get_particle() override;

  private:

      using Base::b_geometry;
      using Base::b_rng_control;
      using Base::b_physics;
      using Base::b_nodes;

      // Total source particles and number for this domain
      def::size_type d_np_total;
      def::size_type d_np_domain;

      // Particle weight
      double d_wt;

      // Number of particles left and number run so far
      def::size_type d_np_left;
      def::size_type d_np_run;

      // CDF to determine cell locations
      std::vector<double> d_cell_cdf;

      // CDF within each cell to determine energy group
      std::vector<std::vector<double>> d_erg_cdfs;
};

} // end namespace profugus

#endif // mc_General_Source_hh

//---------------------------------------------------------------------------//
//                 end of General_Source.hh
//---------------------------------------------------------------------------//
