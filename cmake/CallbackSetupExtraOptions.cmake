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

  # Add install RPATH when building shared
  IF(BUILD_SHARED_LIBS)
    SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
    SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
    print_var(CMAKE_INSTALL_RPATH)
    IF (${CMAKE_SYSTEM_NAME} STREQUAL "Darwin")
      IF(POLICY CMP0042)
        SET(CMAKE_MACOSX_RPATH ON)
        SET(CMAKE_SKIP_BUILD_RPATH FALSE)
        SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
      ELSE()
        SET(CMAKE_INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/lib"
          CACHE STRING "On MacOSX the location of installed shared libs.")
        print_var(CMAKE_INSTALL_NAME_DIR)
      ENDIF()
    ENDIF()
  ENDIF()

  # Disable EpetraExt HDF5 if we are loading HDF5
  SET(EpetraExt_ENABLE_HDF5 OFF CACHE BOOL "Turn off HDF5 in Trilinos.")

ENDMACRO()

##---------------------------------------------------------------------------##
## end of Profugus/cmake/CallbackSetupExtraOptions.cmake
##---------------------------------------------------------------------------##
