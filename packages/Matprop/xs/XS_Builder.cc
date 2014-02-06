//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   xs/XS_Builder.cc
 * \author Thomas M. Evans
 * \date   Wed Feb 05 19:44:36 2014
 * \brief  XS_Builder member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "harness/DBC.hh"
#include "XS_Builder.hh"

namespace profugus
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
XS_Builder::XS_Builder()
    : d_comm(Teuchos::DefaultComm<int>::getComm())
{
}

//---------------------------------------------------------------------------//
// PUBLIC FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Open an xml file of cross sections and broadcast the data.
 */
void XS_Builder::open_and_broadcast(const std_string &xml_file)
{
    // make the new parameterlist
    d_plxs = Teuchos::rcp(new ParameterList("cross sections"));

    // read the data on every domain
    Teuchos::updateParametersFromXmlFileAndBroadcast(
        xml_file.c_str(), d_plxs.ptr(), *d_comm);

    Ensure (!d_plxs.is_null());
}

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of XS_Builder.cc
//---------------------------------------------------------------------------//
