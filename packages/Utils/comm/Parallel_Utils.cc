//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Utils/comm/Parallel_Utils.cc
 * \author Thomas M. Evans
 * \date   Fri Jul 13 15:19:06 2007
 * \brief  Parallel utility function definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "harness/DBC.hh"
#include "harness/Soft_Equivalence.hh"
#include "global.hh"
#include "Parallel_Utils.hh"

#include <iostream>
using std::endl;
using std::cout;

namespace profugus
{
//---------------------------------------------------------------------------//
// SPECIAL GLOBAL EQUIVALENCE FUNCTIONALITY
//---------------------------------------------------------------------------//
/*!
 * \brief Function to check the equivalence of an int across all processors.
 *
 * This function is (hopefully) a temporary parallel check function that more
 * properly belongs in C4.  It is used to check the equivalence of a given
 * integer across all processors.  This is used for Design By Contract.
 *
 * \param local_value integer value to check against
 * \return true if equivalent across all processors; false if not
 */
bool check_global_equiv(int local_value)
{
    int node  = profugus::node();
    int nodes = profugus::nodes();

    // passing condition
    bool pass = false;

    // return true if serial, if not then do check on all processors
    if (nodes == 1)
        pass = true;
    else
    {
        // value from processor above local processor
        int neighbors_value;

        if (node > 0 && node < nodes - 1)
        {
            profugus::send(&local_value, 1, node - 1, 600);
            profugus::receive(&neighbors_value, 1, node + 1, 600);
            if (local_value == neighbors_value) pass = true;
        }
        else if (node == nodes - 1)
        {
            profugus::send(&local_value, 1, node - 1, 600);
            pass = true;
        }
        else if (node == 0)
        {
            profugus::receive(&neighbors_value, 1, node + 1, 600);
            if (local_value == neighbors_value) pass = true;
        }
        else
        {
            INSIST(0, "Something is wrong with nodes!");
        }
    }

    // sync everything so we don't leave before all processors are finished
    profugus::global_barrier();

    // return result
    return pass;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Function to check the equivalence of a double across all
 * processors.
 *
 * This function is the same as check_global_equiv(int) except that doubles
 * are compared to precision eps.
 *
 * \param local_value integer value to check against
 * \param eps precision of double, default 1e-8
 * \return true if equivalent across all processors; false if not
 */
bool check_global_equiv(double local_value, double eps)
{
    using profugus::soft_equiv;

    int node  = profugus::node();
    int nodes = profugus::nodes();

    // passing condition
    bool pass = false;

    // return true if serial, if not then do check on all processors
    if (nodes == 1)
        pass = true;
    else
    {
        // value from processor above local processor
        double neighbors_value;

        if (node > 0 && node < nodes - 1)
        {
            profugus::send(&local_value, 1, node - 1, 600);
            profugus::receive(&neighbors_value, 1, node + 1, 600);
            pass = soft_equiv(neighbors_value, local_value, eps);
        }
        else if (node == nodes - 1)
        {
            profugus::send(&local_value, 1, node - 1, 600);
            pass = true;
        }
        else if (node == 0)
        {
            profugus::receive(&neighbors_value, 1, node + 1, 600);
            pass = soft_equiv(neighbors_value, local_value, eps);
        }
        else
        {
            INSIST(0, "Something is wrong with nodes!");
        }
    }

    // sync everything so we don't leave before all processors are finished
    profugus::global_barrier();

    // return result
    return pass;
}

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of Parallel_Utils.cc
//---------------------------------------------------------------------------//
