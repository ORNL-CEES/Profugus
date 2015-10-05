
IF(NOT TPL_ENABLE_Boost)
  MESSAGE(FATAL_ERROR "\nHPX TPL requires that Boost support is enabled.\n\n")
ELSE()
  TRIBITS_TPL_FIND_INCLUDE_DIRS_AND_LIBRARIES( HPX
    REQUIRED_HEADERS hpx/hpx_main.hpp
    REQUIRED_LIBS_NAMES hpx
    )
ENDIF()
