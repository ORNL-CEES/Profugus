//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   spn/Isotropic_Source.cc
 * \author Thomas M. Evans
 * \date   Fri Feb 14 18:14:05 2014
 * \brief  Isotropic_Source member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <algorithm>
#include "utils/Constants.hh"
#include "Isotropic_Source.hh"

namespace profugus
{

//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 *
 * \param num_cells number of cells in local mesh
 */
Isotropic_Source::Isotropic_Source(int num_cells)
    : d_ids(num_cells, 0)
    , d_source(num_cells, 0.0)
{
}

//---------------------------------------------------------------------------//
/*!
 * \brief Truncate the source.
 *
 * The source gets truncated into a range specified by \c [g_first,g_last].
 * The resulting group will still be defined over the range \c [0,G] but now
 * \c G=1+g_last-g_first.
 *
 * The source shapes are \b not renormalized over the truncated range.  The
 * source strength per group stays the same.
 *
 * \param g_first first group in truncated range \c [0,G]
 * \param g_last last group in truncated range \c [g_first,G]
 */
void Isotropic_Source::truncate(int g_first,
                                int g_last)
{
    REQUIRE(g_first >= 0 && g_first <= g_last);
    REQUIRE(g_last < this->num_groups());

    // number of groups in new source
    int num_groups = 1 + g_last - g_first;
    CHECK(num_groups > 0);

    // return if the number of groups are equivalent
    if (num_groups == this->num_groups()) return;

    // make a new shapes array
    Shape ns(num_groups, 0.0);

    // loop through spectra and make new sources
    for (int s = 0, S = d_shapes.size(); s < S; ++s)
    {
        // loop through spectra and make reassignment
        for (int g = 0, fullg = g_first; g < num_groups; ++g, fullg++)
        {
            ns[g] = d_shapes[s][fullg];
        }

        // assign the new shape
        d_shapes[s] = ns;
    }

    ENSURE(this->num_groups() == num_groups);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Set the source.
 *
 * \param ids cell-field of ids giving the source shape id in a cell
 *
 * \param shapes energy shape of each source; dimensioned
 * [num_shapes][num_groups]
 *
 * \param source source strengths in \f$\mbox{particles}/\mbox{cm}^3\f$ in
 * each cell
 */
void Isotropic_Source::set(const ID_Field      &ids,
                           const Source_Shapes &shapes,
                           const Source_Field  &source)
{
    using profugus::constants::inv_four_pi;

    REQUIRE(!shapes.empty());
    REQUIRE(ids.size() == d_ids.size());
    REQUIRE(source.size() == d_source.size());

    // assign the shapes
    d_shapes = shapes;

    // copy the ids
    std::copy(ids.begin(), ids.end(), d_ids.begin());

    // copy the source strengths
    for (int cell = 0, N = source.size(); cell < N; cell++)
    {
        // convert isotropic source to per steradian
        d_source[cell] = inv_four_pi * source[cell];

        CHECK(d_source[cell] >= 0.0);
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Verify that the sources have been assigned.
 */
bool Isotropic_Source::verify() const
{
    // make sure shapes have been assigned
    if (d_shapes.empty()) return false;

    // make sure the number of groups is consistent
    int num_groups = d_shapes[0].size();
    for (int i = 0, N = d_shapes.size(); i < N; i++)
        if (d_shapes[i].size() != num_groups) return false;

    // check the min and max elements in the ids to make sure they are in
    // range
    if (d_ids.empty())
        return false;

    ID_Field::const_iterator min = std::min_element(d_ids.begin(), d_ids.end());
    ID_Field::const_iterator max = std::max_element(d_ids.begin(), d_ids.end());

    if (*min < 0)
        return false;

    if (*max >= d_shapes.size())
        return false;

    if (d_ids.size() != d_source.size()) return false;

    return true;
}

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of Isotropic_Source.cc
//---------------------------------------------------------------------------//
