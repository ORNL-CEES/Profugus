//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   xs/XS_Builder.cc
 * \author Thomas M. Evans
 * \date   Wed Feb 05 19:44:36 2014
 * \brief  XS_Builder member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Teuchos_XMLParameterListHelpers.hpp"

#include "harness/DBC.hh"
#include "XS_Builder.hh"

namespace profugus
{

//---------------------------------------------------------------------------//
// UNNAMED NAMESPACE

namespace
{

// Names of different total (1D cross sections) - these are ordered the same
// as the XS_Types enumeration in XS.
const char *totals_labels[XS::END_XS_TYPES] = {
    "sigma_t",
    "sigma_f",
    "nu_sigma_f",
    "chi"};

// Names of different Pn scattering cross sections (2D scattering).
const int MAX_PN_ORDER = 8;
const char *scat_labels[MAX_PN_ORDER] = {
    "sigma_s0",
    "sigma_s1",
    "sigma_s2",
    "sigma_s3",
    "sigma_s4",
    "sigma_s5",
    "sigma_s6",
    "sigma_s7"};

}

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

    Check (!d_plxs.is_null());
    Check (d_plxs->isParameter("num groups"));
    Check (d_plxs->isParameter("pn order"));

    // assign number of groups and pn order
    d_pn_order   = d_plxs->get<int>("pn order");
    d_num_groups = d_plxs->get<int>("num groups");

    // get the materials in the file
    d_matids.clear();
    for (ParameterList::ConstIterator itr = d_plxs->begin();
         itr != d_plxs->end(); ++itr)
    {
        if (d_plxs->isSublist(itr->first))
        {
            d_matids.push_back(itr->first);
        }
    }

    Ensure (d_matids.size() > 0);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Build the cross sections.
 *
 * \param map map of material ids to named materials in the xml file
 */
void XS_Builder::build(const Matid_Map &map)
{
    // load all of the cross sections
    build(map, d_pn_order, 0, d_num_groups - 1);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Build the cross sections for a speficied pn-order and over a range
 * of groups.
 *
 * This function is used to load a subset of cross sections on the file into
 * the XS database.  For example, a file may contain P5 cross sections and 23
 * groups.  The user may request a XS containing only P0 and P1 cross sections
 * over groups 3-8.  The number of groups in the generated XS is
 * \f[
   N_g = 1 + g_{\mbox{last}} - g_{\mbox{first}}
 * \f]
 *
 * The group structure of the generated XS will be over the range \c [0,G]
 * such
 * \f[
   g_{\mbox{new}} = g\_{\mbox{old}} - g_{\mbox{first}}
 * \f]
 * Note that the group truncation does \b no renormalizing or corrections.
 * Groups outside the specified range are simply ignored.  Generally, group
 * truncation is used to solve for the downscatter-only part of a problem
 * where the lower energy groups are completely decoupled from the high energy
 * groups.
 *
 * \param pn_order pn order of cross sections to load into the XS database
 *
 * \param g_first first group to include in the database
 * \param g_last last group to include in the database
 *
 * \param map map of material ids to named materials in the xml file
 */
void XS_Builder::build(const Matid_Map &map,
                       int              pn_order,
                       int              g_first,
                       int              g_last)
{
    Require (map.completed());
    Require (1 + g_last - g_first <= d_num_groups);
    Require (pn_order <= d_pn_order);

    // number of groups in the produced XS
    int num_groups = 1 + g_last - g_first;
    Check (num_groups >= 0);

    // create g-end (1 past the last valid entry)
    int g_end = g_last + 1;
    Check (g_end - g_first == num_groups);

    // make a new xs database
    d_xs = Teuchos::rcp(new XS);
    d_xs->set(pn_order, num_groups);

    // iterate through the map and assign the cross sections
    for (Matid_Map::const_iterator itr = map.begin();
         itr != map.end(); ++itr)
    {
        // get the material name and integer id that it will be assigned to in
        // the cross section database
        int matid                 = itr->first;
        const std_string &matname = itr->second;
        Check (d_plxs->isSublist(matname));

        // get the sublist
        const ParameterList &mpl = d_plxs->sublist(matname);

        // loop through existing totals
        for (int t = 0; t < XS::END_XS_TYPES; ++t)
        {
            // check to see if the parameter exists in the xml file
            if (mpl.isParameter(totals_labels[t]))
            {
                // read the 1D cross sections from the parameterlist
                const OneDArray &sig_file = mpl.get<OneDArray>(totals_labels[t]);

                // truncate the range
                OneDArray sigma(sig_file.begin() + g_first,
                                sig_file.begin() + g_end);
                Check (sigma.size() == num_groups);

                // add the cross sections to the xs database
                d_xs->add(matid, t, sigma);
            }
        }

        // loop through Pn moments and add scattering
        for (int n = 0; n <= pn_order; ++n)
        {
            // add all the moments that exist in the xml file
            if (mpl.isParameter(scat_labels[n]))
            {
                // read the 2D scattering cross sections from the
                // parameterList
                const TwoDArray &sig_file = mpl.get<TwoDArray>(scat_labels[n]);

                // truncate the range
                TwoDArray sigma(num_groups, num_groups, 0.0);
                for (int g = 0; g < num_groups; ++g)
                {
                    for (int gp = 0; gp < num_groups; ++gp)
                    {
                        sigma(g, gp) = sig_file(g+g_first, gp+g_first);
                    }
                }

                // add the cross sections to the xs database
                d_xs->add(matid, n, sigma);
            }
        }
    }

    // complete the cross sections
    d_xs->complete();
}

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of XS_Builder.cc
//---------------------------------------------------------------------------//
