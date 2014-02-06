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

#include <string>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_DefaultComm.hpp"

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
    typedef std::string        std_string;
    typedef XS                 XS_t;
    typedef Teuchos::RCP<XS_t> RCP_XS;
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

  private:
    // >>> IMPLEMENTATION

    // Typedefs.
    typedef Teuchos::ParameterList      ParameterList;
    typedef Teuchos::RCP<ParameterList> RCP_ParameterList;
    typedef Teuchos::Comm<int>          Comm;
    typedef Teuchos::RCP<const Comm>    RCP_Comm;

    // Teuchos communicator.
    RCP_Comm d_comm;

    // Parameterlist of cross sections.
    RCP_ParameterList d_plxs;
};

} // end namespace profugus

#endif // xs_XS_Builder_hh

//---------------------------------------------------------------------------//
//                 end of XS_Builder.hh
//---------------------------------------------------------------------------//
