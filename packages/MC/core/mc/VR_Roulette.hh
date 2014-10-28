//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   core/mc/VR_Roulette.hh
 * \author Thomas M. Evans
 * \date   Fri May 09 13:09:37 2014
 * \brief  VR_Roulette class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef core_mc_VR_Roulette_hh
#define core_mc_VR_Roulette_hh

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Variance_Reduction.hh"

namespace profugus
{

//===========================================================================//
/*!
 * \class VR_Roulette
 * \brief Provides functions for performing variance reduction during MC
 * transport.
 *
 * The following variance reduction methods are defined in this class.
 *
 * \section vr_weight_roulette Weight Roulette
 *
 * When the particle weight drops below a cutoff (\f$W_c\f$) defined such that
 * \f$W_s>W_c\f$, with probability
 * \f[
   \nu = \frac{w}{W_s}
 * \f]
 * the particle survives with new weight \f$W_s\f$.  Otherwise the particle is
 * killed. Obviously, when \f$w_c = 0\f$ weight roulette will be off.
 *
 * \section variance_reduction_db Standard DB entries for VR_Roulette
 *
 * The following entries control variance reduction for roulette in MC.
 *
 * \arg \c weight_cutoff (double) sets a weight cutoff below which particles
 * are rouletted (default: 0.25); a value of 0.0 turns off weight roulette
 *
 * \arg \c weight_survival (double) sets the survival weight for weight
 * roulette (default: 2 * \c weight_cutoff )
 */
/*!
 * \example mc/test/tstVR_Roulette.cc
 *
 * Test of VR_Roulette.
 */
//===========================================================================//

class VR_Roulette : public Variance_Reduction
{
    typedef Variance_Reduction Base;

  public:
    //@{
    //! Useful typedefs.
    typedef Geometry_t::Space_Vector      Space_Vector;
    typedef Particle_t::RNG_t             RNG_t;
    typedef Teuchos::ParameterList        ParameterList_t;
    typedef Teuchos::RCP<ParameterList_t> RCP_Std_DB;
    //@}

  private:
    // >>> DATA

    // Weight-roulette parameters.
    double d_Wc; // weight-cutoff
    double d_Ws; // survival-weight

  public:
    // Constructor.
    explicit VR_Roulette(RCP_Std_DB db);

    // >>> VARIANCE REDUCTION INTERFACE

    //! Do nothing at surfaces
    void post_surface(Particle_t& particle, Bank_t& bank) const { /* * */ }

    // Do weight roulette at collisions
    void post_collision(Particle_t& particle, Bank_t& bank) const;

    // >>> ACCESSORS

    //@{
    //! Weight roulette parameters.
    double weight_survival() const { return d_Ws; }
    double weight_cutoff() const { return d_Wc; }
    //@}
};

} // end namespace profugus

#endif // core_mc_VR_Roulette_hh

//---------------------------------------------------------------------------//
//                 end of VR_Roulette.hh
//---------------------------------------------------------------------------//
