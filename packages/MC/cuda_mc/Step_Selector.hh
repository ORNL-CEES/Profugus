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

namespace cuda_profugus
{

//===========================================================================//
/*!
 * \class Step_Selector
 * \brief Chooses the minimum step size from multiple choices.
 *
 * This class tracks the minimum step size (along with an associated tag)
 * among many choices.  For example,
 * \code
   Step_Selector step;

   while (transport)
   {
       // get distance to next boundary
       step.initialize(d_bnd, DIST_TO_BND);

       // calculate dist to collision
       step.submit(d_col, DIST_TO_COL);

       // calculate dist to tally surface
       step.submit(d_sur, DIST_TO_SUR);

       // now, get the minimum step
       d_step = step.step();

       // process the results
       if (step.tag() == DIST_TO_COL)
           process_collision();
       else if (step.tag() == DIST_TO_BND)
           process_boundary();
       else
           tally_surface();

       // ...
   }
 * \endcode
 *
 * Down the road we may template on operation type so that we can take the
 * minimum or maximum step size.
 */
/*!
 * \example cuda_mc/test/tstStep_Selector.cc
 *
 * Test of Step_Selector.
 */
//===========================================================================//

class Step_Selector
{
  private:
    // >>> DATA

    // Current step and tag.
    double d_step;
    int    d_tag;

  public:
    //! Constructor.
    Step_Selector()
        : d_step(0.0)
        , d_tag(0)
    {
    }

    //! Initialize the selection process with a tag and step.
    PROFUGUS_HOST_DEVICE_FUNCTION
    void initialize(double step, int tag)
    {
        d_step = step;
        d_tag  = tag;
    }

    //! Submit a new step length and compare to the current step.
    PROFUGUS_HOST_DEVICE_FUNCTION
    void submit(double step, int tag)
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
    int tag() const { return d_tag; }
};

} // end namespace cuda_profugus

#endif // cuda_mc_Step_Selector_hh

//---------------------------------------------------------------------------//
//              end of cuda_mc/Step_Selector.hh
//---------------------------------------------------------------------------//
