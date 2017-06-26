//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/Fission_Source.hh
 * \author Stuart Slattery
 * \date   Mon May 05 14:22:46 2014
 * \brief  Fission_Source class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_mc_Fission_Source_hh
#define cuda_mc_Fission_Source_hh

#include "Teuchos_ArrayView.hpp"

#include "cuda_utils/Stream.hh"
#include "cuda_geometry/Cartesian_Mesh.hh"
#include "Fission_Rebalance.hh"
#include "Source.hh"
#include "Particle_Vector.hh"

#include "cuda_utils/Shared_Device_Ptr.hh"

namespace cuda_profugus
{

//===========================================================================//
/*!
 * \class Fission_Source
 * \brief Defines both an initial and within-cycle fission source for k-code
 * calculations.
 *
 * Through the database passed into the fission source the client specifies
 * the number of fission particles to run each cycle.  On cycles other than
 * the first, the number of fission source particles is determined by the
 * number of sites sampled during the previous generation.  To keep the
 * effective number of particles constant between cycles, fission particles
 * are born with the following weight,
 * \f[
   w = N_p / M\:,
 * \f]
 * where \f$N_p\f$ is the number of requested particles per cycle and \f$M\f$
 * is the total number of particles sampled from the fission source.  In the
 * first cycle \f$M = N_p\f$; however, in various parallel decompositions the
 * number of particles run in the first cycle may be slightly less than the
 * number of particles requested (but never more).  Thus, even in the first
 * cycle the weight of particles may be slightly greater than one because the
 * weight per particle is normalized to the \b requested number of particles
 * per cycle.  Using this weight preserves the constant total weight per
 * cycle.  The requested number of particles can be queried through the Np()
 * member function.
 *
 * \section fission_source_db Standard DB Entries for Fission_Source
 *
 * The database entries that control the cuda_profugus::Fission_Source are:
 *
 * \arg \c init_fission_src (vector<double>) box defined as (lox, hix, loy,
 * hiy, loz, hiz) in which the initial fission source is sampled (default: the
 * entire geometry)
 *
 * \arg \c Np (int) number of particles to use in each cycle (default:
 * 1000)
 */
/*!
 * \example cuda_mc/test/tstFission_Source.cc
 *
 * Test of Fission_Source.
 */
//===========================================================================//

template <class Geometry>
class Fission_Source : public Source<Geometry>
{
  public:
    //@{
    //! Typedefs.
    typedef Geometry                                    Geometry_t;
    typedef Physics<Geometry_t>                         Physics_t;
    typedef Teuchos::ParameterList			            ParameterList_t;
    typedef Teuchos::RCP<ParameterList_t>		        RCP_Std_DB;
    typedef typename Physics_t::Fission_Site            Fission_Site;
    typedef typename Physics_t::Fission_Site_Container  Fission_Site_Container;
    typedef typename Geometry_t::Space_Vector           Space_Vector;
    typedef cuda_utils::Shared_Device_Ptr<Geometry>     SDP_Geometry;
    typedef cuda_utils::Shared_Device_Ptr<Physics_t>    SDP_Physics;
    typedef std::shared_ptr<Fission_Site_Container>     SP_Fission_Sites;
    typedef Fission_Rebalance<Geometry_t>               Fission_Rebalance_t;
    typedef std::shared_ptr<Fission_Rebalance_t>        SP_Fission_Rebalance;
    typedef cuda_utils::Shared_Device_Ptr<Cartesian_Mesh>     SDP_Cart_Mesh;
    typedef Teuchos::ArrayView<const double>            Const_Array_View;
    typedef def::Vec_Dbl                                Vec_Dbl;
    typedef def::Vec_Int                                Vec_Int;
    typedef def::size_type                              size_type;
    //@}

  public:
    // Constructor.
    Fission_Source(const RCP_Std_DB& db,
                   const SDP_Geometry& geometry,
                   const SDP_Physics& physics,
                   const std::vector<double>& volumes,
                   const def::Space_Vector& low_edge,
                   const def::Space_Vector& high_edge);

    // Destructor.
    ~Fission_Source();

    // Build the initial fission source.
    void build_initial_source();

    // Build the initial source from a mesh distribution.
    void build_initial_source( const SDP_Cart_Mesh& mesh,
			       Const_Array_View& fis_dens);

    // Build a source from a fission site container.
    void build_source(SP_Fission_Sites &fission_sites);

    // Create a fission site container.
    SP_Fission_Sites create_fission_site_container() const;

    //! Whether this is the initial source distribution or is unbuilt
    bool is_initial_source() const { return !d_fission_sites; }

    // >>> DERIVED PUBLIC INTERFACE

    //! Get particles from the source.
    void get_particles(
	cuda_utils::Shared_Device_Ptr<Particle_Vector<Geometry> >& particles ) override;

    //! Boolean operator for source (true when source still has particles).
    bool empty() const override { return d_np_left == 0; }

    //! Number of particles to transport in the source on the current domain.
    size_type num_to_transport() const override { return d_np_domain; }

    //! Total number of particles to transport in the entire problem/cycle.
    size_type total_num_to_transport() const override { return d_np_total; }

    //! Total number of requested particles per cycle.
    size_type Np() const override { return d_np_requested; }

    //! Number transported so far on this domain.
    size_type num_run() const override { return d_np_run; }

    //! Number left to transport on this domain.
    size_type num_left() const override { return d_np_left; }

    // >>> CLASS ACCESSORS

    //! Get the current fission site container.
    const Fission_Site_Container& fission_sites() const
    {
        if (!d_fission_sites)
            return d_dummy_container;
        return *d_fission_sites;
    }

    //! Get fission source lower coords for testing purposes
    const Space_Vector& lower_coords() const { return d_lower; }

    //! Get fission source width for testing purposes
    const Space_Vector& width() const { return d_width; }

    //! Set a new number per cycle.
    void update_Np(size_type np) { d_np_requested = np; }

  private:
    // >>> IMPLEMENTATION

    // Build the domain replicated fission source.
    void build_DR( const SDP_Cart_Mesh& mesh, Const_Array_View& fis_dens);

    // Sample the geometry to generate particles
    void sample_geometry(
	cuda_utils::Shared_Device_Ptr<Particle_Vector<Geometry> >& particles,
	const int num_particle,
	const unsigned int num_blocks,
	const unsigned int threads_per_block );

    // Sample the fission mesh to generate particles
    void sample_mesh(
	cuda_utils::Shared_Device_Ptr<Particle_Vector<Geometry> >& particles,
	const int num_particle,
	const unsigned int num_blocks,
	const unsigned int threads_per_block );

    // Sample the fission_sites to generate particles
    void sample_fission_sites(
	cuda_utils::Shared_Device_Ptr<Particle_Vector<Geometry> >& particles,
	const int num_particle,
	const unsigned int num_blocks,
	const unsigned int threads_per_block );

  private:
    // >>> DATA

    // Execution stream.
    cuda_utils::Stream<cuda_utils::arch::Device> d_stream;

    // Geometry
    SDP_Geometry d_geometry;

    // Physics
    SDP_Physics d_physics;

    // Fission site container.
    SP_Fission_Sites d_fission_sites;

    // Fission rebalance (across sets).
    SP_Fission_Rebalance d_fission_rebalance;

    // Initial fission source lower coords and width.
    Space_Vector d_lower;
    Space_Vector d_width;

    // Requested particles per cycle.
    size_type d_np_requested;

    // Number of particles: total, domain
    size_type d_np_total, d_np_domain;

    // Particle weight.
    double d_wt;

    // Number of source particles left in the current domain.
    size_type d_np_left;

    // Number of particles run on the current domain.
    size_type d_np_run;

    // Dummy fission site container.
    Fission_Site_Container d_dummy_container;

    // Mesh-based starting distribution.
    int			d_current_cell;
    Vec_Int		d_fis_dist;
    SDP_Cart_Mesh	d_fis_mesh;

    // Mesh cell volumes
    std::vector<double> d_volumes;

    // Size of the particle vector.
    int d_vector_size;

    // Device fission sites.
    Fission_Site* d_fission_sites_device;

    // Device fission cells.
    int* d_fission_cells_device;

    // number of inactive cycles.
    int d_num_inactive_cycle;

    // active cycle counter.
    int d_active_cycle_counter;
};

} // end namespace cuda_profugus

#endif // cuda_mc_Fission_Source_hh

//---------------------------------------------------------------------------//
//                 end of Fission_Source.hh
//---------------------------------------------------------------------------//
