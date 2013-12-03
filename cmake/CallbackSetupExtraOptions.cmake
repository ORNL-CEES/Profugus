##---------------------------------------------------------------------------##
## Profugus/cmake/CallbackSetupExtraOptions.cmake
## Thomas M. Evans
## Monday December 2 21:36:44 2013
##---------------------------------------------------------------------------##

IF (Profugus_SOURCE_DIR)
  # We need to inject the Profugus/cmake directory to find several
  # Profugus-specific macros
  SET(CMAKE_MODULE_PATH  ${CMAKE_MODULE_PATH} "${Profugus_SOURCE_DIR}/cmake")
ENDIF()

##---------------------------------------------------------------------------##

# Enable documentation for this project
MACRO(TRIBITS_REPOSITORY_SETUP_EXTRA_OPTIONS)

  #enable Profugus documentation
  INCLUDE(ProfugusDoc)

ENDMACRO()

##---------------------------------------------------------------------------##
## end of Profugus/cmake/CallbackSetupExtraOptions.cmake
##---------------------------------------------------------------------------##
