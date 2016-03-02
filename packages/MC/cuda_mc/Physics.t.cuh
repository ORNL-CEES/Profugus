//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/Physics.t.cuh
 * \author Thomas M. Evans
 * \date   Thursday May 1 11:14:55 2014
 * \brief  Physics template member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_mc_Physics_t_cuh
#define cuda_mc_Physics_t_cuh

#include <sstream>
#include <algorithm>

#include "harness/Soft_Equivalence.hh"
#include "harness/Warnings.hh"
#include "utils/Constants.hh"
#include "utils/Vector_Functions.hh"
#include "utils/View_Field.hh"
#include "Physics.cuh"

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
                           SP_XS      mat_host,
                           SDP_XS_Dev mat)
{
    REQUIRE(!db.is_null());

    d_Ng = mat_host->num_groups();
    d_Nm = mat_host->num_mat();
    REQUIRE( d_Ng > 0 );
    REQUIRE( d_Nm > 0 );

    auto err = cudaMalloc( (void**)&d_scatter, d_Ng*d_Nm*sizeof(double) );
    REQUIRE( err == cudaSuccess );
    err = cudaMalloc( (void**)&d_fissionable, d_Nm*sizeof(int) );
    REQUIRE( err == cudaSuccess );

    d_mat = mat.get_device_ptr();
    CHECK( d_mat );

    // implicit capture flag
    d_implicit_capture = db->get("implicit_capture", true);

    // check for balanced scattering tables
    bool check_balance = db->get("check_balance", false);

    // turn check balance on if we are not doing implicit capture
    if (!d_implicit_capture) check_balance = true;

    // calculate total scattering over all groups for each material and
    // determine if fission is available for a given material
    std::vector<double> host_scatter(d_Nm*d_Ng);
    std::vector<int> host_fissionable(d_Nm);

    // XS_Device matids always run 0->Nm
    for (int m = 0; m < d_Nm; ++m)
    {
        // get the P0 scattering matrix for this material
        const auto &sig_s = mat_host->matrix(m, 0);
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
        if (check_balance)
        {
            for (int g = 0; g < d_Ng; g++)
            {
                int ind = group_mat_index(g,m);

                if (host_scatter[ind] >
                    mat_host->vector(m, XS_t::TOTAL)[g])
                {
                    // terminate if we are running analog
                    if (!d_implicit_capture)
                    {
                        INSIST(false, "Scattering greater than total.");
                    }
                    // else add to warnings
                    else
                    {
                        ADD_WARNING("Scattering greater than total.");
                    }
                }
            }
        }

        // see if this material is fissionable by checking Chi
        host_fissionable[m] = mat_host->vector(m, XS_t::CHI).normOne() > 0.0 ?
                              true : false;
    }

    // Assign host data to device vectors
    err = cudaMemcpy(d_scatter, &host_scatter[0],
                     host_scatter.size()*sizeof(double),
                     cudaMemcpyHostToDevice);
    REQUIRE( err == cudaSuccess );
    err = cudaMemcpy(d_fissionable, &host_fissionable[0],
                     host_fissionable.size()*sizeof(int),
                     cudaMemcpyHostToDevice);
    REQUIRE( err == cudaSuccess );

    cudaDeviceSynchronize();
}

// Destructor
template <class Geometry>
Physics<Geometry>::~Physics()
{
    auto err = cudaFree(d_scatter);
    err = cudaFree(d_fissionable);
}

} // end namespace cuda_mc

#endif // cuda_mc_Physics_t_cuh

//---------------------------------------------------------------------------//
//                 end of Physics.t.cuh
//---------------------------------------------------------------------------//
