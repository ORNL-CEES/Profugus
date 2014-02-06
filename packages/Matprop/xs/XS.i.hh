//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   xs/XS.i.hh
 * \author 9te
 * \date   Wed Jan 29 15:27:36 2014
 * \brief  Member definitions of class XS.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef xs_XS_i_hh
#define xs_XS_i_hh

namespace profugus
{

//---------------------------------------------------------------------------//
/*!
 * \brief Check to see if a given matid exists.
 */
bool XS::has(int matid) const
{
    Require (!d_totals.empty());
    return d_totals[TOTAL].exists(matid);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Return the 1-D data vector for a given matid and type.
 */
const XS::Vector& XS::vector(int matid,
                             int type) const
{
    Require (type < d_totals.size());
    Require (d_totals[type].exists(matid));
    return *d_totals[type][matid];
}

//---------------------------------------------------------------------------//
/*!
 * \brief Return the 2-D data matrix for a given matid and Pn order.
 */
const XS::Matrix& XS::matrix(int matid,
                             int pn) const
{
    Require (pn < d_scatter.size());
    Require (d_scatter[pn].exists(matid));
    return *d_scatter[pn][matid];
}

} // end namespace profugus

#endif // xs_XS_i_hh

//---------------------------------------------------------------------------//
//                 end of XS.i.hh
//---------------------------------------------------------------------------//
