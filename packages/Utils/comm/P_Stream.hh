//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Utils/comm/P_Stream.hh
 * \author Thomas M. Evans
 * \date   Mon Oct 08 10:31:04 2012
 * \brief  P_Stream class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef Utils_comm_P_Stream_hh
#define Utils_comm_P_Stream_hh

#include <iostream>
#include "global.hh"

namespace profugus
{

//===========================================================================//
/*!
 * \class P_Manip
 * \brief Manipulator class for parallel stream objects.
 *
 * For the following constructions,
 * \code
 * My_Ostream out;
 * out << setw(12) << x << endl;
 * \endcode
 * The stream manipulator \c setw requires an argument whereas \c endl does
 * not.  By default, the compiler will add an implicit argument of type
 * \c My_Ostream to \c endl. It is unable to do the same for \c setw because
 * it takes an argument.  Hence, this class provides a bridge between the
 * stream output class and the manipulator.  In the example above, we add the
 * following overloaded operator to the \c My_Ostream class:
 * \code
 * template<class T>
 * My_Ostream& My_OStream::operator<<(const P_Manip<T> &p)
 * {
 *     p.action(p.argument);
 * }
 * \endcode
 * Now we need to make the manipulator that returns a P_Manip:
 * \code
 * P_Manip<int> setw(int n)
 * {
 *     return setw(setw_implementation, n);
 * }
 * \endcode
 * Finally, we make the implementation (which the user never sees) of the \c
 * setw operation:
 * \code
 * void setw_implementation(int n)
 * {
 *     std::cout << std::setw(n);
 * }
 * \endcode
 * The implementation functions must correspond to the Manipulator signature,
 * ie.
 * \code
 * void implementation(T value);
 * \endcode
 */
//===========================================================================//

template<class T>
class P_Manip
{
  public:
    //! Function pointer to a manipulator that takes an argument.
    typedef void (*Manipulator)(T);

  public:
    //! Constructor.
    P_Manip(Manipulator m, T v)
        : action(m), argument(v) {}

    //! Action of the manipulator.
    Manipulator action;

    //! Argument of the manipulator.
    T argument;
};

//===========================================================================//
/*!
 * \class P_Out
 * \brief Class for doing parallel stream output.
 *
 * The P_Out class effectively replaces the following constructs:
 * \code
 * if (node() == 0)
 * {
 *     cout << "Here is the value on node " << node() << ": "
 *          << setw(12) << fixed << x << endl;
 * }
 * \endcode
 * with
 * \code
 * P_Out my_pout;
 * // ...
 *
 * my_pout << "Here is the value on node " << node() << ": "
 *         << setw(12) << fixed << x << endl;
 * \endcode
 *
 * The optional argument to the constructor allows output to be directed to
 * any processor.
 *
 * There are four P_Out objects defined in global scope:
 *
 * \arg profugus::pout can be toggled on/off via UTILS_POUT at configure
 * time; only does output on node 0
 *
 * \arg profugus::pcout is always on; only does output on node 0
 *
 * \arg profugus::pnout can be toggled on/off via UTILS_POUT at configure
 * time; does output from the local node
 *
 * \arg profugus::pncout is always on; does output from the local node
 *
 * In effect profugus::pnout and profugus::pout will give identical behavior as
 * std::cout wit the exception that profugus::pnout can be toggled off.
 *
 * \sa P_Manip, UTILS_POUT, pout, pcout, pnout, pncout
 */
/*!
 * \example comm/test/tstP_Stream.cc
 *
 * Test of parallel stream output classes.
 */
//===========================================================================//

class P_Out
{
  private:
    //>>> DATA

    // Function pointers for stream manipulators (ie. endl, setw).
    typedef P_Out& (*FP)(P_Out &);

    // Master node to write output to.
    int d_master;

  public:
    // Constructor.
    P_Out(int master = 0);

    // Overloaded stream operators.
    template<class T> P_Out& operator<<(const T &t);

    // Stream manipulators that do not require an argument.
    P_Out& operator<<(FP f);

    // Stream manipulators that require an argument.
    template<class T> P_Out& operator<<(const P_Manip<T> &p);

    //! Return the master node.
    int master() const { return d_master; }

    //! Set the master node.
    void set_master(int m) { d_master = m; }
};

//---------------------------------------------------------------------------//
// FREE FUNCTION STREAM MANIPULATORS FOR P_OUT
//---------------------------------------------------------------------------//

P_Out& endl(P_Out &p);
P_Out& scientific(P_Out &p);
P_Out& fixed(P_Out &p);
P_Out& left(P_Out &p);
P_Out& right(P_Out &p);

//---------------------------------------------------------------------------//
// ARGUMENT BASED STREAM MANIPULATORS
//---------------------------------------------------------------------------//

P_Manip<int> setw(int n);
P_Manip<char> setfill(char c);
P_Manip<int> setprecision(int n);

//---------------------------------------------------------------------------//
// EXTERNAL, GLOBAL PARALLEL STREAM OUTPUT OBJECTS
//---------------------------------------------------------------------------//
/*!
 * \page UTILS_POUT Using the Utils Parallel Output Facility
 *
 * If defined, then pout will output to processor 0.  Otherwise, output does
 * not occur.
 */
//---------------------------------------------------------------------------//

// Standard parallel stream output.
extern P_Out pout;   // off when UTILS_POUT undefined
extern P_Out pcout;  // always on
extern P_Out pnout;  // off when UTILS_POUT undefined
extern P_Out pncout; // always on

} // end namespace profugus

//---------------------------------------------------------------------------//
// INLINE AND TEMPLATE DEFINITIONS
//---------------------------------------------------------------------------//

#include "P_Stream.i.hh"

#endif // Utils_comm_P_Stream_hh

//---------------------------------------------------------------------------//
//              end of comm/P_Stream.hh
//---------------------------------------------------------------------------//
