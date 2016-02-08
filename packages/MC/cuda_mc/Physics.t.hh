//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/Physics.t.hh
 * \author Thomas M. Evans
 * \date   Thursday May 1 11:14:55 2014
 * \brief  Physics template member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_mc_Physics_t_hh
#define cuda_mc_Physics_t_hh

#include <sstream>
#include <algorithm>

#include "harness/Soft_Equivalence.hh"
#include "harness/Warnings.hh"
#include "utils/Constants.hh"
#include "utils/Vector_Functions.hh"
#include "utils/View_Field.hh"
#include "Physics.hh"

namespace cuda_mc
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor that implicitly creates Group_Bounds
 */
template <class Geometry>
Physics<Geometry>::Physics(RCP_Std_DB db,
                           SP_XS      mat)
    : d_Ng(mat->num_groups())
    , d_Nm(mat->num_mat())
    , d_scatter_vec(mat->num_groups()*mat->num_mat())
    , d_fissionable_vec(d_Nm)
{
    REQUIRE(!db.is_null());

    //auto mat_dev = std::make_shared<XS_Dev_t>(*mat);
    d_mat_host = cuda::Shared_Device_Ptr<XS_Dev_t>(*mat);
    d_mat = d_mat_host.get_device_ptr();

    REQUIRE(d_mat);
    REQUIRE(d_mat_host.get_host_ptr()->num_groups() > 0);
    REQUIRE(d_mat_host.get_host_ptr()->num_mat() > 0);

    d_scatter = d_scatter_vec.data();
    d_fissionable = d_fissionable_vec.data();

    // implicit capture flag
    d_implicit_capture = db->get("implicit_capture", true);

    // check for balanced scattering tables
    d_check_balance = db->get("check_balance", false);

    // turn check balance on if we are not doing implicit capture
    if (!d_implicit_capture) d_check_balance = true;

    // get the material ids in the database
    def::Vec_Int matids;
    mat->get_matids(matids);
    CHECK(matids.size() == d_Nm);

    // calculate total scattering over all groups for each material and
    // determine if fission is available for a given material
    std::vector<double> host_scatter(d_Nm*d_Ng);
    std::vector<int> host_fissionable(d_Nm);
    for (auto m : matids)
    {
        // get the P0 scattering matrix for this material
        const auto &sig_s = mat->matrix(m, 0);
        CHECK(sig_s.numRows() == d_Ng);
        CHECK(sig_s.numCols() == d_Ng);

        // loop over all groups and calculate the in-scatter from other
        // groups and add them to the group OUT-SCATTER; remember, we
        // store data as inscatter for the deterministic code
        for (int g = 0; g < d_Ng; g++)
        {
            // get the g column (g->g' scatter stored as g'g in the matrix)
            const auto *column = sig_s[g];

            // add up the scattering
            for (int gp = 0; gp < d_Ng; ++gp)
            {
                host_scatter[group_mat_index(g,m)] += sig_s(gp,g);
            }
        }

        // check scattering correctness if needed
        if (d_check_balance)
        {
            for (int g = 0; g < d_Ng; g++)
            {
                int ind = group_mat_index(g,m);

                if (host_scatter[ind] >
                    mat->vector(m, XS_t::TOTAL)[g])
                {
                    std::ostringstream mm;
                    mm << "Scattering greater than total "
                       << "for material" << m << " in group " << g
                       << ". Total xs is "
                       << mat->vector(m, XS_t::TOTAL)[g]
                       << " and scatter is " << host_scatter[ind];

                    // terminate if we are running analog
                    if (!d_implicit_capture)
                        VALIDATE(false, mm.str());
                    // else add to warnings
                    else
                        ADD_WARNING(mm.str());
                }
            }
        }

        // see if this material is fissionable by checking Chi
        host_fissionable[m] = mat->vector(m, XS_t::CHI).normOne() > 0.0 ?
                              true : false;
    }

    // Assign host data to device vectors
    d_scatter_vec.assign(profugus::make_view(host_scatter));
    d_fissionable_vec.assign(profugus::make_view(host_fissionable));

    ENSURE(d_Nm > 0);
    ENSURE(d_Ng > 0);
}

} // end namespace cuda_mc

#endif // cuda_mc_Physics_t_hh

//---------------------------------------------------------------------------//
//                 end of Physics.t.hh
//---------------------------------------------------------------------------//
