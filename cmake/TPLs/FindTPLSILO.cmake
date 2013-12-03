INCLUDE(TribitsTplDeclareLibraries)

IF(NOT TPL_ENABLE_HDF5)
  MESSAGE(FATAL_ERROR "\nSILO TPL requires that HDF5 support is enabled.\n\n")
ELSE()
  TRIBITS_TPL_DECLARE_LIBRARIES( SILO
    REQUIRED_HEADERS silo.h
    REQUIRED_LIBS_NAMES siloh5
    )
ENDIF()
