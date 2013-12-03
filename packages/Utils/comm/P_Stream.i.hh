//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   comm/P_Stream.i.hh
 * \author Thomas M. Evans
 * \date   Mon Oct 08 10:31:04 2012
 * \brief  Member definitions of class P_Stream.
 * \note   Copyright (C) 2012 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef comm_P_Stream_i_hh
#define comm_P_Stream_i_hh

namespace nemesis
{

//---------------------------------------------------------------------------//
// P_OUT MEMBER DEFINITIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Stream out of overloaded types.
 *
 * This will correctly do stream output for any type that can be sent to
 * standard std::ostream objects.
 */
template<class T>
P_Out& P_Out::operator<<(const T &t)
{
    if (node() == d_master)
    {
        std::cout << t;
    }
    return *this;
}

//---------------------------------------------------------------------------//
/*
 * \brief Overloaded stream operator for manipulators that require an
 * argument.
 *
 * This overloaded operator takes manipulators that require an argument,
 * ie. setw(n), setfill('c'), etc.
 */
template<class T>
P_Out& P_Out::operator<<(const P_Manip<T> &p)
{
    if (d_master == node())
    {
        Require (p.action);
        p.action(p.argument);
    }
    return *this;
}

} // end namespace nemesis

#endif // comm_P_Stream_i_hh

//---------------------------------------------------------------------------//
//              end of comm/P_Stream.i.hh
//---------------------------------------------------------------------------//
