//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/Step_Selector.hh
 * \author Thomas M. Evans
 * \date   Monday May 12 12:11:34 2014
 * \brief  Step_Selector class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_mc_Step_Selector_hh
#define cuda_mc_Step_Selector_hh

#include "cuda_utils/CudaMacros.hh"
#include "Definitions.hh"

namespace cuda_profugus
{

//===========================================================================//
/*!
 * \class Step_Selector
 * \brief Chooses the minimum step size from multiple choices.
 */
//===========================================================================//

class Step_Selector
{
  private:
    // >>> DATA

    // Current step and tag.
    double d_step;
    events::Event    d_tag;

  public:
    //! Constructor.
    PROFUGUS_HOST_DEVICE_FUNCTION
    Step_Selector()
        : d_step(0.0)
        , d_tag(events::DEAD)
    { /* ... */ }

    //! Initialize the selection process with a tag and step.
    PROFUGUS_HOST_DEVICE_FUNCTION
    void initialize(double step, events::Event tag)
    {
        d_step = step;
        d_tag  = tag;
    }

    //! Submit a new step length and compare to the current step.
    PROFUGUS_HOST_DEVICE_FUNCTION
    void submit(double step, events::Event tag)
    {
        if (step < d_step)
        {
            d_step = step;
            d_tag  = tag;
        }
    }

    // >>> ACCESSORS

    //! Return the current step.
    PROFUGUS_HOST_DEVICE_FUNCTION
    double step() const { return d_step; }

    //! Return the tag associated with the current step.
    PROFUGUS_HOST_DEVICE_FUNCTION
    events::Event tag() const { return d_tag; }
};

} // end namespace cuda_profugus

#endif // cuda_mc_Step_Selector_hh

//---------------------------------------------------------------------------//
//              end of cuda_mc/Step_Selector.hh
//---------------------------------------------------------------------------//
