//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   harness/Scalar_Unit_Test.hh
 * \author Thomas M. Evans
 * \date   Wed Jan  2 11:54:57 2008
 * \brief  Scalar_Unit_Test class definition.
 * \note   Copyright (C) 2008 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef harness_Scalar_Unit_Test_hh
#define harness_Scalar_Unit_Test_hh

#include <iostream>
#include "Unit_Test.hh"

namespace nemesis
{

//===========================================================================//
/*!
 * \class Scalar_Unit_Test
 * \brief This class provides services for scalar unit tests.
 *
 * This class inherits from Unit_Test.  Much of the documentation for the
 * services of this class is provided in Unit_Test.hh
 *
 * \sa nemesis::Unit_Test for detailed description of all the Unit_Test
 * classes.
 *
 * \par Code Sample:
 *
 * Scalar Unit_Tests should have the following syntax.
 * \code
int main(int argc, char *argv[])
{
    try
    {
        nemesis::Scalar_Unit_Test ut( argc, argv, release );
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
 * \test All of the member functions of this class are tested by
 * harness/test/tstScalar_Unit_Test.cc, including the early exit caused by
 * \c --version on the command line.
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
 * changes to these nemesis m4 macro commands.
 */
/*!
 * \example harness/test/tstScalar_Unit_Test.cc
 * The unit test for and example usage of the Scalar_Unit_Test class.
 */
//===========================================================================//

class Scalar_Unit_Test : public Unit_Test
{
  public:

    // CREATORS

    //! Default constructors.
    Scalar_Unit_Test( int & argc, char **&argv, string_fp_void release_,
              std::ostream & out_ = std::cout );

    //! The copy constructor is disabled.
    Scalar_Unit_Test(const Scalar_Unit_Test &rhs);

    //! Destructor.
    ~Scalar_Unit_Test(void){ out << resultMessage() << std::endl; return; }

    // MANIPULATORS

    //! The assignment operator for Scalar_Unit_Test is disabled.
    Scalar_Unit_Test& operator=(const Scalar_Unit_Test &rhs);
};

} // end namespace nemesis

#endif // harness_Scalar_Unit_Test_hh

//---------------------------------------------------------------------------//
//              end of harness/Scalar_Unit_Test.hh
//---------------------------------------------------------------------------//
