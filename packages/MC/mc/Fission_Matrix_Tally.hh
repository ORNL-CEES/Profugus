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

#include <memory>
#include <unordered_map>
#include <vector>
#include <utility>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "harness/DBC.hh"
#include "geometry/Mesh_Geometry.hh"
#include "Tally.hh"

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

class Fission_Matrix_Tally : public Tally
{
    typedef Tally Base;

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

    //! Hash table for pair of ints.
    struct Idx_Hash
    {
      public:
        std::hash<int> d_hash;
        int            d_N;

        //! Constructor.
        Idx_Hash(int N = 0) : d_N(N) {/*...*/}

        size_t operator()(const std::pair<int, int> &x) const
        {
            REQUIRE(d_N > 0);
            return d_hash(x.first + d_N * x.second);
        }
    };
    
  public:
    // Sparse matrix storage for FM Tally.
    typedef std::pair<int, int>                       Idx;
    typedef std::unordered_map<Idx, double, Idx_Hash> Sparse_Matrix;

  private:
    // Geometry.
    SP_Geometry d_geometry;

    // Fission matrix mesh.
    SP_Mesh_Geometry d_fm_mesh;

    // Fission matrix tallies.
    Sparse_Matrix       d_numerator;
    std::vector<double> d_denominator;

  public:
    // Constructor.
    Fission_Matrix_Tally(RCP_Std_DB db, SP_Physics physics,
                         SP_Mesh_Geometry fm_mesh);

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
    int d_cycle_out;
    int d_cycle_ctr;
};

} // end namespace profugus

#endif // mc_Fission_Matrix_Tally_hh

//---------------------------------------------------------------------------//
//                 end of Fission_Matrix_Tally.hh
//---------------------------------------------------------------------------//
