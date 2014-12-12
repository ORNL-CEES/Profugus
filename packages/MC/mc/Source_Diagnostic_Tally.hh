//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/Source_Diagnostic_Tally.hh
 * \author Thomas M. Evans
 * \date   Tue Dec 09 16:32:04 2014
 * \brief  Source_Diagnostic_Tally class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef mc_Source_Diagnostic_Tally_hh
#define mc_Source_Diagnostic_Tally_hh

#include <memory>
#include <vector>
#include <string>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "utils/Serial_HDF5_Writer.hh"
#include "geometry/Mesh_Geometry.hh"
#include "Tally.hh"

namespace profugus
{

//===========================================================================//
/*!
 * \class Source_Diagnostic_Tally
 * \brief Source diagnostic tally for KCode problems.
 *
 * Tally MC source densities on an input mesh.
 *
 */
/*!
 * \example mc/test/tstSource_Diagnostic_Tally.cc
 *
 * Test of Source_Diagnostic_Tally.
 */
//===========================================================================//

class Source_Diagnostic_Tally : public Source_Tally
{
    typedef Source_Tally Base;

  public:
    //@{
    //! Typedefs.
    typedef Physics_t::SP_Geometry         SP_Geometry;
    typedef std::shared_ptr<Mesh_Geometry> SP_Mesh_Geometry;
    typedef Teuchos::ParameterList         ParameterList_t;
    typedef Teuchos::RCP<ParameterList_t>  RCP_Std_DB;
    //@}

  private:
    // >>> DATA

    // Tally mesh.
    SP_Mesh_Geometry d_mesh;

    // Problem geometry.
    SP_Geometry d_geometry;

    // Geometric state on the fission tally mesh geometry.
    Mesh_Geometry::Geo_State_t d_state;

    // In a KCODE calculation the source density is the particle density,
     //    Np * wt = Nr
    // then
    //    sum (wt) / Nr = (n Nr/Np) / Nr = n / Np
    // in a Kcode calculation we finalize with Nr, so the sum of weights over
    // Nr is equal to the sum of particles (if each had weight 1) divided by
    // the actual number transported
    std::vector<double> d_source_density;

  public:
    // Constructor.
    Source_Diagnostic_Tally(RCP_Std_DB db, SP_Physics physics,
                            SP_Mesh_Geometry mesh, bool inactive);

    // >>> PUBLIC INTERFACE

    // Tally events at particle birth.
    void birth(const Particle_t &p);

    // End a cycle in a kcode calculation.
    void end_cycle(double num_particles);

    // Clear/re-initialize all tally values.
    void reset();

  private:
    // >>> IMPLEMENTATION

    // Output file name and hdf5 writer.
    std::string        d_filename;
    Serial_HDF5_Writer d_writer;

    // Cycle counter.
    int d_cycle_ctr;
};

} // end namespace profugus

#endif // mc_Source_Diagnostic_Tally_hh

//---------------------------------------------------------------------------//
//                 end of Source_Diagnostic_Tally.hh
//---------------------------------------------------------------------------//
