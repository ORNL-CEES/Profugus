//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   xs/XS_Builder.hh
 * \author Thomas M. Evans
 * \date   Wed Feb 05 19:44:36 2014
 * \brief  XS_Builder class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef xs_XS_Builder_hh
#define xs_XS_Builder_hh

#include <vector>
#include <string>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_DefaultComm.hpp"

#include "utils/Static_Map.hh"
#include "XS.hh"

namespace profugus
{

//===========================================================================//
/*!
 * \class XS_Builder
 * \brief Build an XS from input.
 */
/*!
 * \example xs/test/tstXS_Builder.cc
 *
 * Test of XS_Builder.
 */
//===========================================================================//

class XS_Builder
{
  public:
    //@{
    //! Typedefs.
    typedef std::string                 std_string;
    typedef XS                          XS_t;
    typedef Teuchos::RCP<XS_t>          RCP_XS;
    typedef Static_Map<int, std_string> Matid_Map;
    typedef std::vector<std_string>     Vec_Str;
    //@}

  private:
    // >>> DATA

    // Built cross sections.
    RCP_XS d_xs;

  public:
    // Constructor.
    XS_Builder();

    // Open and broadcast an XML cross section file.
    void open_and_broadcast(const std_string &xml_file);

    // Build the cross sections.
    void build(const Matid_Map &map);

    // Build the cross sections over certain groups with a specific Pn order.
    void build(const Matid_Map &map, int pn, int g_first, int g_last);

    //! Get the cross sections.
    RCP_XS get_xs() const { return d_xs; }

    // >>> FILE ACCESSORS

    //! Get the Pn order of scattering data in the file.
    int pn_order() const { return d_pn_order; }

    //! Get the number of groups in the file.
    int num_groups() const { return d_num_groups; }

    //! Get the list of materials defined in the file.
    const Vec_Str& materials() const { return d_matids; }

  private:
    // >>> IMPLEMENTATION

    // Typedefs.
    typedef Teuchos::ParameterList      ParameterList;
    typedef Teuchos::RCP<ParameterList> RCP_ParameterList;
    typedef Teuchos::Comm<int>          Comm;
    typedef Teuchos::RCP<const Comm>    RCP_Comm;
    typedef XS_t::OneDArray             OneDArray;
    typedef XS_t::TwoDArray             TwoDArray;

    // Teuchos communicator.
    RCP_Comm d_comm;

    // Parameterlist of cross sections.
    RCP_ParameterList d_plxs;

    // Pn order and number of groups in the cross section file.
    int d_pn_order, d_num_groups;

    // Group velocites (may not be present; not needed for static problems).
    OneDArray d_velocity;

    // Group bounds (may not be present).
    OneDArray d_bounds;

    // Materials in the file.
    Vec_Str d_matids;
};

} // end namespace profugus

#endif // xs_XS_Builder_hh

//---------------------------------------------------------------------------//
//                 end of XS_Builder.hh
//---------------------------------------------------------------------------//
