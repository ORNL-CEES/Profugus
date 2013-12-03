//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   harness/Unit_Test.hh
 * \author Thomas M. Evans
 * \date   Wed Jan  2 11:54:17 2008
 * \brief  Unit_Test class definition.
 * \note   Copyright (C) 2008 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef harness_Unit_Test_hh
#define harness_Unit_Test_hh

#include <string>
#include <list>
#include <iostream>
#include <map>
#include "DBC.hh"

namespace nemesis
{

//===========================================================================//
/*!
 * \class Unit_Test
 * \brief Object to encapsulate unit testing of Nemesis classes and functions.
 *
 * This is a virtual class.  You should use one of the following Unit_Test
 * classes in your test appliation:
 *
 * \li Scalar_Unit_Test - Used for testing code that does not use parallel
 * communication (comm).
 * \li Parallel_Unit_Test - Used for testing code that does use parallel
 * communications (comm).
 *
 * This unit test classification is tied into the Nemesis Build System.  Unit
 * tests are declared in each package's configure.ac file using the following
 * syntax that corresponds to the unit test's classification:
 *
 * \li \c AC_RUNTESTS( tstName, scalar ) - Run the test tstName as a scalar
 * process (expectes a ScalarUnit_Test object).
 * \li \c AC_RUNTESTS( tstName, 2 5 ) - Run the test tstName under MPI twice.
 * Once with 2 and again with 5 processors (expects a Parallel_Unit_Test
 * object).
 *
 * \sa Unit_Test.cc for additional details.
 *
 * \par Code Sample:
 *
 * Scalar_Unit_Tests should have the following syntax.
 * \code
int main(int argc, char *argv[])
{
    try
    {
        // Test ctor for Scalar_Unit_Test (also tests Unit_Test ctor and member
        // function setTestName).
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
 */
//===========================================================================//

class Unit_Test
{
  public:

    // NESTED CLASSES AND TYPEDEFS

    //! Typedef for function pointer to this package's release function.
    typedef std::string (*string_fp_void)(void);

    // CREATORS

    //! Default constructors.
    Unit_Test( int & argc, char **&argv, string_fp_void release_,
              std::ostream & out_ = std::cout );

    //! The copy constructor is disabled.
    Unit_Test(const Unit_Test &rhs);

    //! Destructor is virtual because this class will be inherited from.
    virtual ~Unit_Test(void){/*empty*/};

    // MANIPULATORS

    //! The assignment operator is disabled.
    Unit_Test& operator=(const Unit_Test &rhs);

    //! Change the target for output
    // void setostream( std::ostream out_ ) { out = out_; return; }

    // ACCESSORS
    bool failure(int line);
    bool failure(int line, const char *file);
    bool failure(const std::string &failmsg);
    bool passes(const std::string &failmsg);

    //! This pure virtual function must be provided by the inherited class.
    /// It should provide output concerning the status of Unit_Test.
    void status(void) const { out << resultMessage() << std::endl; return; }

    //! Reset the pass and fail counts to zero.
    void reset() { numPasses=0;numFails=0; return; }

    // DATA

    //! The number of passes found for this test.
    unsigned numPasses;

    //! The number of failures found for this test.
    unsigned numFails;

    // Features
    static std::map< std::string, unsigned >
    get_word_count( std::ostringstream const & data, bool verbose=false );

  protected:

    // IMPLEMENTATION
    std::string resultMessage(void) const;
    std::string setTestName( std::string const fqName );
    std::string setTestPath( std::string const fqName );

    // DATA

    //! The name of this unit test.
    std::string const testName;
    //! Relative path to the unit test.
    std::string const testPath;

    //! Function pointer this package's release(void) function
    string_fp_void release;

    //! Where should output be sent (default is std::cout)
    std::ostream & out;
};

} // end namespace nemesis

#endif // harness_Unit_Test_hh

//---------------------------------------------------------------------------//
//              end of harness/Unit_Test.hh
//---------------------------------------------------------------------------//
