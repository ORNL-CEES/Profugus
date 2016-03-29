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

class KDE_Fission_Source : public Fission_Source
{
    typedef Fission_Source  Base;

  public:
    //@{
    //! Useful typedefs
    typedef Base::RCP_Std_DB            RCP_Std_DB;
    typedef Base::SP_Fission_Sites      SP_Fission_Sites;
    typedef std::shared_ptr<KDE_Kernel> SP_KDE_Kernel;
    typedef KDE_Kernel::cell_type       cell_type;
    //@}

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
