##---------------------------------------------------------------------------##
## Profugus/cmake/NativeRepositoriesList.cmake
## Thomas M. Evans
## Tuesday July 17 10:40:13 2012
##---------------------------------------------------------------------------##

# Native repositories for Profugus project builds
SET(NATIVE_REPOS .)

# Search to see if MCLS exists in the Profugus source tree
FIND_PATH(MCLS_EXISTS
  NAMES MCLS/PackagesList.cmake
  PATHS ${CMAKE_CURRENT_SOURCE_DIR})
IF (MCLS_EXISTS)
  SET(NATIVE_REPOS MCLS ${NATIVE_REPOS})
ELSE()
  MESSAGE(STATUS "MCLS repository is not available")
ENDIF()

# Search to see if ParaSails exists in the Profugus source tree
FIND_PATH(ParaSails_EXISTS
  NAMES ParaSails/PackagesList.cmake
  PATHS ${CMAKE_CURRENT_SOURCE_DIR})
IF (ParaSails_EXISTS)
  SET(NATIVE_REPOS ParaSails ${NATIVE_REPOS})
ELSE()
  MESSAGE(STATUS "ParaSails repository is not available")
ENDIF()

# Assume the user has already symlinked Trilinos into the current dir
SET(NATIVE_REPOS Trilinos ${NATIVE_REPOS})

# Set the native repos
SET(${PROJECT_NAME}_NATIVE_REPOSITORIES ${NATIVE_REPOS})
UNSET(NATIVE_REPOS)

##---------------------------------------------------------------------------##
## end of Profugus/cmake/NativeRepositoriesList.cmake
##---------------------------------------------------------------------------##
