//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   driver/Problem_Builder.cc
 * \author Thomas M. Evans
 * \date   Wed Mar 12 22:25:22 2014
 * \brief  Problem_Builder member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <vector>

#include "Teuchos_XMLParameterListHelpers.hpp"

#include "harness/DBC.hh"
#include "xs/XS_Builder.hh"
#include "Problem_Builder.hh"

namespace profugus
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
Problem_Builder::Problem_Builder()
    : d_comm(Teuchos::DefaultComm<int>::getComm())
{
    // build validator
    d_validator = Teuchos::rcp(new ParameterList("validator"));

    // build sublists
    ParameterList core, assbly, mat, mesh, problem;

    d_validator->set("CORE", core);
    d_validator->set("ASSEMBLIES", core);
    d_validator->set("MATERIAL", core);
    d_validator->set("MESH", core);
    d_validator->set("PROBLEM", core);
}

//---------------------------------------------------------------------------//
// PUBLIC INTERFACE
//---------------------------------------------------------------------------//
/*!
 * \brief Setup the problem.
 */
void Problem_Builder::setup(const std::string &xml_file)
{
    // make the master parameterlist
    auto master = Teuchos::rcp(new ParameterList);

    // read the data on every domain
    Teuchos::updateParametersFromXmlFileAndBroadcast(
        xml_file.c_str(), master.ptr(), *d_comm);

    // validate the parameter list
    Insist (master->isSublist("CORE"),
            "CORE block not defined in input.");
    Insist (master->isSublist("ASSEMBLIES"),
            "ASSEMBLIES block not defined in input.");
    Insist (master->isSublist("MATERIAL"),
            "MATERIAL block not defined in input.");
    Insist (master->isSublist("MESH"),
            "MESH block not defined in input.");
    Insist (master->isSublist("PROBLEM"),
            "PROBLEM block not defined in input.");

    // store the individual parameter lists
    d_coredb   = Teuchos::sublist(master, "CORE");
    d_assblydb = Teuchos::sublist(master, "ASSEMBLIES");
    d_matdb    = Teuchos::sublist(master, "MATERIAL");
    d_meshdb   = Teuchos::sublist(master, "MESH");
    d_db       = Teuchos::sublist(master, "PROBLEM");

    Check (!d_db.is_null());

    // build mesh
    build_mesh();

    // build material database
    build_matdb();
}

//---------------------------------------------------------------------------//
// PRIVATE IMPLEMENTATION
//---------------------------------------------------------------------------//
/*!
 * \brief Build/partition the mesh.
 */
void Problem_Builder::build_mesh()
{
    Require (d_coredb->isParameter("axial list"));

    // get the axial core map
    const auto &axial_list = d_coredb->get<OneDArray_str>("axial list");
    Check (!axial_list.empty());
}

//---------------------------------------------------------------------------//
/*!
 * \brief Build the material data.
 *
 * For now, we build all the cross sections in the problem on every domain.
 * For the mini-app, this is not expected to be an overburdening cost.
 */
void Problem_Builder::build_matdb()
{
    typedef XS_Builder::Matid_Map Matid_Map;
    typedef XS_Builder::RCP_XS    RCP_XS;

    Require (d_matdb->isParameter("mat list"));
    Validate (d_matdb->isParameter("xs library"),
              "Inline cross sections not implemented yet.");

    // get the material list off of the database
    const auto &mat_list = d_matdb->get<OneDArray_str>("mat list");

    // convert the matlist to a mat-id map
    Matid_Map matids;
    for (int id = 0, N = mat_list.size(); id < N; ++id)
    {
        matids.insert(Matid_Map::value_type(id, mat_list[id]));
    }
    matids.complete();
    Check (matids.size() == mat_list.size());

    // make a cross section builder
    XS_Builder builder;

    // broadcast the raw cross section data
    builder.open_and_broadcast(d_matdb->get<std::string>("xs library"));

    // get the number of groups and moments in the cross section data
    int Ng_data = builder.num_groups();
    int N_data  = builder.pn_order();

    // determine the moment order of the problem
    int pn_order = d_db->get("Pn_order", N_data);
    Validate (pn_order <= N_data, "Requested Pn scattering order of "
              << pn_order << " is greater than available data Pn order of "
              << N_data);

    // get the number of groups required
    int g_first = d_db->get("g_first", 0);
    int g_last  = d_db->get("g_last", Ng_data - 1);
    Validate (1 + (g_last - g_first) <= Ng_data, "Energy group range exceeds "
              << "number of groups in data, 1 + g_last - g_first = "
              << 1 + (g_last - g_first) << " > " << Ng_data);

    // build the cross sections
    builder.build(matids, pn_order, g_first, g_last);
    RCP_XS xs = builder.get_xs();
    Check (xs->num_mat() == matids.size());
    Check (xs->num_groups() == 1 + (g_last - g_first));

    // build the material database
    d_mat = Teuchos::rcp(new Mat_DB);

    // set the cross sections
    d_mat->set(xs, d_mesh->num_cells());
}

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of Problem_Builder.cc
//---------------------------------------------------------------------------//
