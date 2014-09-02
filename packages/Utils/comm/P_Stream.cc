//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   comm/P_Stream.cc
 * \author Thomas M. Evans
 * \date   Mon Oct 08 10:31:04 2012
 * \brief  P_Stream member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <iomanip>

#include "harness/DBC.hh"
#include "P_Stream.hh"

namespace profugus
{

//---------------------------------------------------------------------------//
// P_OUT CLASS MEMEBER DEFINITIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
P_Out::P_Out(int master)
    : d_master(master)
{
}

//---------------------------------------------------------------------------//
/*!
 * \brief Overloaded stream operator for manipulators that do not require an
 * argument.
 *
 * This overloaded operator takes manipulators that do not require an
 * argument, ie. endl, scientific, fixed, etc.
 */
P_Out& P_Out::operator<<(FP f)
{
    REQUIRE(f);
    f(*this);
    return *this;
}

//---------------------------------------------------------------------------//
// STREAM MANIPULATORS
//---------------------------------------------------------------------------//
/*!
 * \brief endl for P_Out streams.
 */
P_Out& endl(P_Out &p)
{
    if (node() == p.master())
    {
        std::cout << std::endl;
    }
    return p;
}

//---------------------------------------------------------------------------//
/*!
 * \brief scientific for P_Out streams.
 */
P_Out& scientific(P_Out &p)
{
    if (node() == p.master())
    {
        std::cout << std::scientific;
    }
    return p;
}

//---------------------------------------------------------------------------//
/*!
 * \brief fixed for P_Out streams.
 */
P_Out& fixed(P_Out &p)
{
    if (node() == p.master())
    {
        std::cout << std::fixed;
    }
    return p;
}

//---------------------------------------------------------------------------//
/*!
 * \brief left for P_Out streams.
 */
P_Out& left(P_Out &p)
{
    if (node() == p.master())
    {
        std::cout << std::left;
    }
    return p;
}

//---------------------------------------------------------------------------//
/*!
 * \brief right for P_Out streams.
 */
P_Out& right(P_Out &p)
{
    if (node() == p.master())
    {
        std::cout << std::right;
    }
    return p;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Implementation of setw for parallel streams.
 *
 * The protoype for this function must be of type P_Manip::Manipulator;
 */
void setw_implementation(int n)
{
    // apply the standard setw to cout
    std::cout << std::setw(n);
}

//---------------------------------------------------------------------------//
/*!
 * \brief setw for P_Out streams.
 */
P_Manip<int> setw(int n)
{
    // make and return a parallel implementation class that uses the setw
    // implementation
    return P_Manip<int>(setw_implementation, n);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Implementation of setfill for parallel streams.
 *
 * The protoype for this function must be of type P_Manip::Manipulator;
 */
void setfill_implementation(char c)
{
    // apply the standard setfill to cout
    std::cout << std::setfill(c);
}

//---------------------------------------------------------------------------//
/*!
 * \brief setfill for P_Out streams.
 */
P_Manip<char> setfill(char c)
{
    // make and return a parallel implementation class that uses the setfill
    // implementation
    return P_Manip<char>(setfill_implementation, c);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Implementation of setprecision for parallel streams.
 *
 * The protoype for this function must be of type P_Manip::Manipulator;
 */
void setprecision_implementation(int n)
{
    // apply the standard setw to cout
    std::cout << std::setprecision(n);
}

//---------------------------------------------------------------------------//
/*!
 * \brief setprecision for P_Out streams.
 */
P_Manip<int> setprecision(int n)
{
    // make and return a parallel implementation class that uses the
    // setprecision implementation
    return P_Manip<int>(setprecision_implementation, n);
}

//---------------------------------------------------------------------------//
// GLOBAL EXTERNAL DEFINITIONS
//---------------------------------------------------------------------------//

#ifdef UTILS_POUT

//! Use the default master as processor 0.
P_Out pout;

//! Make a local-node pout.
/// It gets the node assigned later in the call to initialize().
P_Out pnout;

#else

//@{
//! Effectively turns off parallel output by setting master to -1.
P_Out pout(-1);
P_Out pnout(-1);
//@}

#endif

//! Always define a pcout object, this cannot be toggled.
P_Out pcout;

//! Always make a pncout local-node outputter.
/// It gets the node assigned later in the call to initialize().
P_Out pncout;

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of P_Stream.cc
//---------------------------------------------------------------------------//
