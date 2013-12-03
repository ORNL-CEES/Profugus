//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   harness/Demangle_Type.hh
 * \author Seth R Johnson
 * \date   Wed Oct 16 08:23:54 2013
 * \brief  Demangle_Type class definition.
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef harness_Demangle_Type_hh
#define harness_Demangle_Type_hh

#include <string>
#include <typeinfo>

namespace nemesis
{

//===========================================================================//
/*!
 * \class Demangle_Type
 * \brief Utility function for demangling C++ types (specifically with GCC).
 *
 * See: http://stackoverflow.com/questions/281818/unmangling-the-result-of-stdtype-infoname
 */
/*!
 * \example harness/test/tstDemangle_Type.cc
 *
 * Test of Demangle_Type.
 */
//===========================================================================//

class Demangle_Type
{
  public:
    // Get the pretty typename of a variable (static)
    template<typename T>
    static std::string demangled_type()
    {
        return Demangle_Type::demangle(typeid(T).name());
    }

    // Get the pretty typename of a variable (dynamic)
    template<typename T>
    static std::string demangled_type(const T& t)
    {
        return Demangle_Type::demangle(typeid(t).name());
    }

    // Demangle the name that comes from typeid
    static std::string demangle(const char* typeid_name);

  private:
    // Prevent construction
    Demangle_Type();
};

} // end namespace nemesis

#endif // harness_Demangle_Type_hh

//---------------------------------------------------------------------------//
//                 end of Demangle_Type.hh
//---------------------------------------------------------------------------//
