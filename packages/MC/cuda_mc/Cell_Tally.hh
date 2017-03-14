//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/Cell_Tally.hh
 * \author Stuart Slattery
 * \brief  Cell class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_mc_Cell_Tally_hh
#define cuda_mc_Cell_Tally_hh

#include <Teuchos_Array.hpp>
#include <Teuchos_ParameterList.hpp>

#include "utils/Definitions.hh"

#include "cuda_utils/Shared_Device_Ptr.hh"
#include "cuda_utils/Stream.hh"

#include "Definitions.hh"
#include "Particle_Vector.hh"
#include "Tally.hh"

namespace cuda_profugus
{
//===========================================================================//
/*!
 * \class Cell_Tally
 * \brief Cell_Tally class for MC transport with batch statistics.
 */
/*!
 * \example mc/test/tstCell_Tally.cc
 *
 * Test of Particle.
 */
//===========================================================================//
template <class Geometry>
class Cell_Tally : public Pathlength_Tally<Geometry>
{
  public:
    //@{
    //! Typedefs.
    typedef Particle_Vector<Geometry> Particle_Vector_t;
    typedef typename Geometry::Geo_State_t Geo_State_t;
    typedef events::Event Event_t;
    typedef Teuchos::RCP<Teuchos::ParameterList>    RCP_Std_DB;
    //@}

  private:
    // >>> DATA

    // Geometry over which we are tallying.
    cuda_utils::Shared_Device_Ptr<Geometry> d_geometry;

    // Number of statistical batches in the tally.
    int d_num_batch;

    // Number of cells in the tally.
    int d_num_cells;

    // Cell tallies. Indexed as d_tally[batch][cell]. On-device.
    double* d_tally;

    // Execution stream.
    cuda_utils::Stream<cuda_utils::arch::Device> d_stream;

    // Database
    RCP_Std_DB d_db;

    // Storage for first and second moments of tally result
    std::vector<double> d_first, d_second;

  public:

    // >>> HOST API

    // Constructor
    Cell_Tally( RCP_Std_DB db,
                const cuda_utils::Shared_Device_Ptr<Geometry>& geometry,
		const int num_batch );

    // Destructor.
    ~Cell_Tally();

    // Get the number of cells in the tally.
    int num_cells() const { return d_num_cells; }

    // Get the number of batches in the tally.
    int num_batch() const { return d_num_batch; }

    // Tally the particles in a vector.
    void accumulate(
	const cuda_utils::Shared_Device_Ptr<Particle_Vector_t>& particles ) override;

    // Query if this tally is on during inactive cycles.
    bool inactive_cycle_tally() const override { return false; }

    // Finalize the tally.
    void finalize( double num_particles ) override;

    // Copy the first and second tally moments from the device to the
    // host. The moments are lazy-evaluated in this function and indexed by
    // cell.
    void copy_moments_to_host( Teuchos::Array<double>& first_moment,
			       Teuchos::Array<double>& second_moment ) const
    {
        first_moment = d_first;
        second_moment = d_second;
    }
};

//---------------------------------------------------------------------------//

} // end namespace cuda_profugus

#endif // cuda_mc_Cell_Tally_hh

//---------------------------------------------------------------------------//
//                 end of Cell_Tally.hh
//---------------------------------------------------------------------------//
