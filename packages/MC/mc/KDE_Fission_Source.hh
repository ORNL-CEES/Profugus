//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   mc/KDE_Fission_Source.hh
 * \author Gregory Davidson
 * \date   Mon Nov 23 15:47:23 2015
 * \brief  KDE_Fission_Source class declaration.
 * \note   Copyright (c) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef MC_mc_KDE_Fission_Source_hh
#define MC_mc_KDE_Fission_Source_hh

#include <memory>

#include "Fission_Source.hh"
#include "KDE_Kernel.hh"

namespace profugus
{

//===========================================================================//
/*!
 * \class KDE_Fission_Source
 * \brief Samples from the KDE_Kernel when sampling a particle.
 */
/*!
 * \example mc/test/tstKDE_Fission_Source.cc
 *
 * Test of KDE_Fission_Source.
 */
//===========================================================================//

template<class Geometry>
class KDE_Fission_Source : public Fission_Source<Geometry>
{
    typedef Fission_Source<Geometry>  Base;

  public:
    //@{
    //! Useful typedefs
    typedef typename Base::RCP_Std_DB        RCP_Std_DB;
    typedef typename Base::SP_Geometry       SP_Geometry;
    typedef typename Base::SP_Physics        SP_Physics;
    typedef typename Base::SP_RNG_Control    SP_RNG_Control;
    typedef typename Base::SP_Particle       SP_Particle;
    typedef typename Base::SP_Fission_Sites  SP_Fission_Sites;
    typedef typename Base::Fission_Site      Fission_Site;
    typedef typename Base::Particle_t        Particle_t;
    typedef typename Base::Space_Vector      Space_Vector;
    typedef KDE_Kernel<Geometry>             KDE_Kernel_t;
    typedef std::shared_ptr<KDE_Kernel_t>    SP_KDE_Kernel;
    typedef typename KDE_Kernel_t::cell_type cell_type;
    //@}

  private:
    // Expose base class data
    using Base::d_fission_sites;
    using Base::d_wt;
    using Base::d_num_left;
    
  public:
    // Constructor
    KDE_Fission_Source(RCP_Std_DB     db,
                       SP_Geometry    geometry,
                       SP_Physics     physics,
                       SP_RNG_Control rng_control);

    // Build a source from a fission site container
    virtual void build_source(SP_Fission_Sites &fission_sites) override;

    // Sample a particle
    virtual SP_Particle get_particle() override;

    //! Get the bandwidth
    double bandwidth(cell_type cellid) const
    {
        return d_kernel->bandwidth(cellid);
    }

  private:
    // >>> IMPLEMENTATION FUNCTION
    Fission_Site& sample_fission_site(RNG &rng) const;

    // >>> IMPLEMENTATION DATA
    // Stores the KDE kernel
    SP_KDE_Kernel d_kernel;
};

//---------------------------------------------------------------------------//
} // end namespace profugus

#endif // MC_mc_KDE_Fission_Source_hh

//---------------------------------------------------------------------------//
// end of MC/mc/KDE_Fission_Source.hh
//---------------------------------------------------------------------------//
