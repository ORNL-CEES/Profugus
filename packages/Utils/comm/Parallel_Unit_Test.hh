//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   comm/Parallel_Unit_Test.hh
 * \author Thomas M. Evans
 * \date   Wed Jan  2 15:44:20 2008
 * \brief  Parallel_Unit_Test member definitions.
 * \note   Copyright (C) 2008 Oak Ridge National Laboratory, UT-Battelle, LLC.
 *
 * This file provides a definition for Parallel_Unit_Test.  The purpose of
 * this class is to encapsulate the keywords and behavior of nemesis parallel
 * unit tests.
 */
//---------------------------------------------------------------------------//

#ifndef comm_Parallel_Unit_Test_hh
#define comm_Parallel_Unit_Test_hh

#include <iostream>
#include "harness/Unit_Test.hh"

namespace nemesis
{

//===========================================================================//
/*!
 * \class Parallel_Unit_Test
 * \brief This class encapsulates services for parallel unit tests.
 *
 * This class inherits from Unit_Test.  Much of the documentation for the
 * services of this class is provided in Unit_Test.hh
 *
 * \sa nemesis::Unit_Test for detailed description of all the Unit_Test
 * classes.
 *
 * \par Code Sample:
 * \code
int main(int argc, char *argv[])
{
    try
    {
        nemesis::Parallel_Unit_Test ut( argc, argv, release );
        tstOne(ut);
        ut.status();
    }
    catch( nemesis::assertion &err )
    {
        std::string msg = err.what();
        if( msg != std::string( "Success" ) )
        { cout << "ERROR: While testing " << argv[0] << ", "
               << err.what() << endl;
            return 1;
        }
        return 0;
    }
    catch (exception &err)
    {
        cout << "ERROR: While testing " << argv[0] << ", "
             << err.what() << endl;
        return 1;
    }
    catch( ... )
    {
        cout << "ERROR: While testing " << argv[0] << ", "
             << "An unknown exception was thrown" << endl;
        return 1;
    }
    return 0;
}
 * \endcode
 *
 * \warning The output from this class is closely tied to the nemesis python
 * script \c tools/regression_filter.py that is used during \c gmake \c check.
 * Changing the format or keyword in the output streams from this class should
 * be coordinated with the regular expression matches found in \c
 * tools/regression_filter.py.
 *
 * \warning The differences between Scalar_Unit_Test and Parallel_Unit_Test
 * are correlated to the nemesis m4 macros \c AC_RUNTESTS and \c
 * AC_TEST_APPLICATION.  Changes to these classes should be coordinated with
 * changes to these nemesis m4 macro command
 */
/*!
 * \example comm/test/tstParallel_Unit_Test.cc
 * This unit test demonstrates typical usage for Parallel_Unit_Test.
 */
//===========================================================================//

class Parallel_Unit_Test : public Unit_Test
{
  public:
    //! Default constructor.
    Parallel_Unit_Test(int & argc, char **&argv, string_fp_void release_,
                       std::ostream & out_ = std::cout);

    //!  The copy constructor is disabled.
    Parallel_Unit_Test(const Parallel_Unit_Test &rhs);

    //! Destructor.
    ~Parallel_Unit_Test();

    // >>> MANIPULATORS

    //! The assignment operator is disabled.
    Parallel_Unit_Test& operator=(const Parallel_Unit_Test &rhs);

    // >>> ACCESSORS

    //! Provide a report of the number of unit test passes and fails.
    void status(void);
};

} // end namespace nemesis

#endif // comm_Parallel_Unit_Test_hh

//---------------------------------------------------------------------------//
//              end of comm/Parallel_Unit_Test.hh
//---------------------------------------------------------------------------//
