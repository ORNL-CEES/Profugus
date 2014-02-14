//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   utils/Definitions.hh
 * \author Thomas M. Evans
 * \date   Tue Sep  4 16:23:33 2007
 * \brief  Definitions that are used by classes in denovo.
 * \note   Copyright (C) 2007 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef utils_Definitions_hh
#define utils_Definitions_hh

#include <cstddef>

#include <vector>
#include <string>

#include <utils/config.h>

// Declare Vector_Lite
namespace profugus
{
template <class T, size_t N> class Vector_Lite;
}

namespace def
{

//! XYZ directions.
enum XYZ {X       = 0,
          Y       = 1,
          Z       = 2,
          END_XYZ = 3};

//! IJK directions.
enum IJK {I       = 0,
          J       = 1,
          K       = 2,
          END_IJ  = 2,
          END_IJK = 3};

//! Direction cosines.
enum Direction_Cosines {MU                    = 0,
                        ETA                   = 1,
                        XI                    = 2,
                        END_DIRECTION_COSINES = 3};

//! Block/Cell categories.
enum Entity_Categories { PROBLEM_BOUNDARY   = -1,
                         INTERNAL           = 0,
                         PROCESSOR_BOUNDARY = 1 };

//! Three-dimensional space vector.
typedef profugus::Vector_Lite<double, 3> Space_Vector;

//@{
//! Common vector typedefs.
typedef std::vector<int>         Vec_Int;
typedef std::vector<double>      Vec_Dbl;
typedef std::vector<float>       Vec_Flt;
typedef std::vector<std::string> Vec_String;
typedef std::vector<char>        Vec_Char;
typedef std::vector<bool>        Vec_Bool;
//@}

//@{
//! PODs.
typedef UTILS_UNSIGNED_INT2 tiny_int;
typedef std::size_t         size_type;
//@}

} // end namespace profugus::def

#endif // utils_Definitions_hh

//---------------------------------------------------------------------------//
//              end of utils/Definitions.hh
//---------------------------------------------------------------------------//
