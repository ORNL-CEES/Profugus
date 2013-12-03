##---------------------------------------------------------------------------##
## ProfugusAddCFlags.cmake
## Seth R Johnson
## Thu Sep 13 14:02:49 EDT 2012
##---------------------------------------------------------------------------##
## Copyright (C) 2012 Oak Ridge National Laboratory, UT-Battelle, LLC.
## --------------------------------------------------------------------------##

INCLUDE(CheckCXXCompilerFlag)
INCLUDE(CheckCCompilerFlag)

# See if the compiler allows the given flags; add them to C and CXX flags if
# so.
#
# This is useful for excluding warnings in GCC without breaking older versions
# of GCC.
MACRO(PROFUGUS_ADD_CXX_FLAGS)
  FOREACH(THEFLAG ${ARGN})
    string(REGEX REPLACE "[^0-9a-zA-Z]" "_" FLAGNAME ${THEFLAG})
    check_cxx_compiler_flag("${THEFLAG}" PROFUGUS_USE_CXX_FLAG_${FLAGNAME})
    if(PROFUGUS_USE_CXX_FLAG_${FLAGNAME})
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${THEFLAG}")
    endif()
  ENDFOREACH()
ENDMACRO()

# See if the compiler allows the given flags; add them to C and CXX flags if
# so.
#
# This is useful for excluding warnings in GCC without breaking older versions
# of GCC.
MACRO(PROFUGUS_ADD_C_FLAGS)
  FOREACH(THEFLAG ${ARGN})
    string(REGEX REPLACE "[^0-9a-zA-Z]" "_" FLAGNAME ${THEFLAG})
    check_c_compiler_flag("${THEFLAG}" PROFUGUS_USE_C_FLAG_${FLAGNAME})
    if(PROFUGUS_USE_C_FLAG_${FLAGNAME})
      set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${THEFLAG}")
    endif()
  ENDFOREACH()
ENDMACRO()

##---------------------------------------------------------------------------##
## end of ProfugusAddCFlags.cmake
##---------------------------------------------------------------------------##
