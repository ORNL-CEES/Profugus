//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/VR_Roulette.cuh
 * \author Thomas M. Evans
 * \date   Fri May 09 13:09:37 2014
 * \brief  VR_Roulette class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_mc_VR_Roulette_cuh
#define cuda_mc_VR_Roulette_cuh

#include "Particle.cuh"
#include "cuda_utils/CudaDBC.hh"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

namespace cuda_mc
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

template <class Geometry>
class VR_Roulette
{
  public:
    //@{
    //! Useful typedefs.
    typedef Geometry                            Geometry_t;
    typedef Particle<Geometry_t>                Particle_t;
    typedef cuda_profugus::Space_Vector         Space_Vector;
    typedef Teuchos::ParameterList              ParameterList_t;
    typedef Teuchos::RCP<ParameterList_t>       RCP_Std_DB;
    //@}

  private:

    // Weight-roulette parameters.
    double d_Wc; // weight-cutoff
    double d_Ws; // survival-weight

  public:
    // Constructor.
    explicit VR_Roulette(RCP_Std_DB db);

    // >>> VARIANCE REDUCTION INTERFACE

    //! Do nothing at surfaces
    __device__ void post_surface(Particle_t& particle) const { /* * */ }

    // Do weight roulette at collisions
    __device__ inline void post_collision(Particle_t& particle) const;

    __host__ __device__ bool uses_splitting() const { return false; }

    // >>> ACCESSORS

    //@{
    //! Weight roulette parameters.
    double weight_survival() const { return d_Ws; }
    double weight_cutoff() const { return d_Wc; }
    //@}
};

} // end namespace cuda_mc

#include "VR_Roulette.i.cuh"

#endif // cuda_mc_VR_Roulette_cuh

//---------------------------------------------------------------------------//
//                 end of VR_Roulette.cuh
//---------------------------------------------------------------------------//
