//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/Fission_Matrix_Tally.hh
 * \author Thomas M. Evans
 * \date   Tue Jul 22 15:09:01 2014
 * \brief  Fission_Matrix_Tally class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef mc_Fission_Matrix_Tally_hh
#define mc_Fission_Matrix_Tally_hh

#include <string>
#include <memory>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "utils/Serial_HDF5_Writer.hh"
#include "geometry/Mesh_Geometry.hh"
#include "Tally.hh"
#include "Fission_Matrix_Processor.hh"

namespace profugus
{

//===========================================================================//
/*!
 * \class Fission_Matrix_Tally
 * \brief Tally the fission matrix during transport.
 */
/*!
 * \example mc/test/tstFission_Matrix_Tally.cc
 *
 * Test of Fission_Matrix_Tally.
 */
//===========================================================================//

class Fission_Matrix_Tally : public Source_Tally,
                             public Pathlength_Tally
{
    typedef Tally Base;

  public:
    //@{
    //! Typedefs.
    typedef Physics_t::SP_Geometry                  SP_Geometry;
    typedef std::shared_ptr<Mesh_Geometry>          SP_Mesh_Geometry;
    typedef Teuchos::ParameterList                  ParameterList_t;
    typedef Teuchos::RCP<ParameterList_t>           RCP_Std_DB;
    typedef Fission_Matrix_Processor::Idx           Idx;
    typedef Fission_Matrix_Processor::Sparse_Matrix Sparse_Matrix;
    typedef Fission_Matrix_Processor::Denominator   Denominator;
    //@}

  private:
    // Geometry.
    SP_Geometry d_geometry;

    // Fission matrix mesh.
    SP_Mesh_Geometry d_fm_mesh;

    // Fission matrix tallies.
    Sparse_Matrix d_numerator;
    Denominator   d_denominator;

    // Fission matrix processor.
    Fission_Matrix_Processor d_processor;

  public:
    // Constructor.
    Fission_Matrix_Tally(RCP_Std_DB db, SP_Physics physics,
                         SP_Mesh_Geometry fm_mesh);

    // Manually build the global fission matrix at the current state.
    void build_matrix();

    // Get the fission matrix processor (results).
    Fission_Matrix_Processor& processor() { return d_processor; }

    //! Query if fission matrix tallying has started.
    bool tally_started() const { return d_tally_started; }

    // >>> INHERITED INTERFACE

    //! Tally events at particle birth.
    void birth(const Particle_t &p);

    //! Track particle, using pre-calculated physics information (multipliers)
    void accumulate(double step, const Particle_t &p);

    //! Accumulate first and second moments
    void end_history() { /* * */ }

    //! Do post-processing on first and second moments
    void finalize(double num_particles) { /* * */ }

    //! Begin active cycles in a kcode calculation (no-op)
    void begin_active_cycles() { /* * */ }

    //! Begin a new cycle in a kcode calculation (no-op)
    void begin_cycle() { /* * */ }

    //! End a cycle in a kcode calculation.
    void end_cycle(double num_particles);

    //! Clear/re-initialize all tally values between solves
    void reset();

  private:
    // >>> IMPLEMENTATION

    // Geometric state on the fission tally mesh geometry.
    Mesh_Geometry::Geo_State_t d_fm_state;

    // Fission matrix birth cell metadata index.
    const unsigned int d_birth_idx;

    // Fission matrix generation options.
    int d_cycle_start;
    int d_cycle_out;
    int d_cycle_ctr;

    // Output file name and hdf5 writer.
    std::string        d_filename;
    Serial_HDF5_Writer d_writer;

    // Tally started flag.
    bool d_tally_started;
};

} // end namespace profugus

#endif // mc_Fission_Matrix_Tally_hh

//---------------------------------------------------------------------------//
//                 end of Fission_Matrix_Tally.hh
//---------------------------------------------------------------------------//
