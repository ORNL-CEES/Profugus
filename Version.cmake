##---------------------------------------------------------------------------##
## Profugus/Version.cmake
## Thomas M. Evans
## Monday December 2 21:37:2 2013
##---------------------------------------------------------------------------##
# Single file that needs to be changed on a release branch
# or on the development branch in order to configure Trilinos
# for release mode and set the version.

SET(Profugus_VERSION 0.0.0)
SET(Profugus_MAJOR_VERSION 00)
SET(Profugus_MAJOR_MINOR_VERSION 040000)
SET(Profugus_VERSION_STRING "0.0.0 (Dev)")
SET(Profugus_ENABLE_DEVELOPMENT_MODE_DEFAULT ON) # Change to 'OFF' for a release

# Used by testing scripts and should not be used elsewhere
SET(Profugus_REPOSITORY_BRANCH "master" CACHE INTERNAL "")
SET(Profugus_TESTING_TRACK "" CACHE INTERNAL "")

##---------------------------------------------------------------------------##
## end of Profugus/Version.cmake
##---------------------------------------------------------------------------##
