//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   spn/FV_Gather.hh
 * \author Thomas M. Evans
 * \date   Tue Oct 30 11:04:52 2012
 * \brief  FV_Gather class definition.
 * \note   Copyright (C) 2012 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef spn_FV_Gather_hh
#define spn_FV_Gather_hh

#include "Teuchos_RCP.hpp"

#include "harness/DBC.hh"
#include "comm/global.hh"
#include "utils/Definitions.hh"
#include "mesh/LG_Indexer.hh"
#include "Moment_Coefficients.hh"
#include "SDM_Face_Field.hh"

namespace profugus
{

//===========================================================================//
/*!
 * \class FV_Gather
 * \brief Gather off-processor diffusion matrices.
 */
/*!
 * \example spn/test/tstFV_Gather.cc
 *
 * Test of FV_Gather.
 */
//===========================================================================//

class FV_Gather
{
  public:
    //@{
    //! Typedefs.
    typedef Teuchos::RCP<Moment_Coefficients>        RCP_Moment_Coefficients;
    typedef SDM_Face_Field                           Face_Field_t;
    typedef Teuchos::RCP<Face_Field_t>               RCP_Face_Field;
    typedef profugus::Vector_Lite<RCP_Face_Field, 2> Face_Fields;
    typedef Teuchos::RCP<SDM_Face_Field::Mesh_t>     RCP_Mesh;
    typedef LG_Indexer                               Indexer_t;
    //@}

    //! Low/High face enumeration.
    enum LOHI {LO = 0, HI = 1};

  private:
    // >>> DATA

    // Mesh.
    RCP_Mesh d_mesh;

    // Moment coefficient generator.
    RCP_Moment_Coefficients d_coefficients;

    // Incoming face fields.
    Face_Fields d_incoming_I;
    Face_Fields d_incoming_J;

    // Outgoing face fields.
    Face_Fields d_outgoing_I;
    Face_Fields d_outgoing_J;

  public:
    // Constructor.
    FV_Gather(RCP_Mesh mesh, RCP_Moment_Coefficients coefficients,
              const Indexer_t &indexer);

    // Gather data for a given equation order.
    void gather(int eqn);

    // >>> ACCESSORS

    // Get diffusion matrices on a given face.
    RCP_Face_Field low_side_D(int face) const;
    RCP_Face_Field high_side_D(int face) const;

  private:
    // >>> IMPLEMENTATION

    typedef profugus::Vector_Lite<int, 2>               Tuple;
    typedef profugus::Vector_Lite<profugus::Request, 2> Handles;
    typedef Moment_Coefficients::Serial_Matrix          Serial_Matrix;

    // Neighbor blocks (domains).
    Tuple d_neighbor_I;
    Tuple d_neighbor_J;

    // Request handles.
    Handles d_request_I;
    Handles d_request_J;

    // Number of block meshes in (i,j) directions.
    Tuple d_Nb;

    // The (i,j) index of this block.
    Tuple d_ij;

    // Number of domains.
    int d_domain, d_domains;

    // Convert (i,j) block indices to domain index.
    int convert(int i, int j)
    {
        Require (i < d_Nb[def::I]);
        Require (j < d_Nb[def::J]);
        Ensure (i + j * d_Nb[def::I] < d_domains);
        return i + j * d_Nb[def::I];
    }

    // Post receives.
    void post_receives();

    // Fill i,j face fields.
    void fill_I_face(int eqn, RCP_Face_Field field, int i);
    void fill_J_face(int eqn, RCP_Face_Field field, int j);
};

} // end namespace profugus

#endif // spn_FV_Gather_hh

//---------------------------------------------------------------------------//
//              end of spn/FV_Gather.hh
//---------------------------------------------------------------------------//
