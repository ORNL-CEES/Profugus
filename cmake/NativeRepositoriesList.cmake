##---------------------------------------------------------------------------##
## Profugus/cmake/NativeRepositoriesList.cmake
## Thomas M. Evans
## Tuesday July 17 10:40:13 2012
##---------------------------------------------------------------------------##

# Native repositories for Profugus project builds
SET(NATIVE_REPOS .)

# Assume the user has already symlinked Trilinos into the current dir
SET(NATIVE_REPOS Trilinos ${NATIVE_REPOS})

SET(${PROJECT_NAME}_NATIVE_REPOSITORIES ${NATIVE_REPOS})
UNSET(NATIVE_REPOS)

##---------------------------------------------------------------------------##
## end of Profugus/cmake/NativeRepositoriesList.cmake
##---------------------------------------------------------------------------##
