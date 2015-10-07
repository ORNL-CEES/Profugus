IF(NOT TPL_ENABLE_Boost)
  MESSAGE(FATAL_ERROR "\nHPX TPL requires that Boost support is enabled.\n\n")
ENDIF()

SET(REQUIRED_HEADERS hpx/hpx_main.hpp)
SET(REQUIRED_LIBS_NAMES hpx)

#
# Second, search for HPX components (if allowed) using the standard
# FIND_PACKAGE(HPX ...).
#
TRIBITS_TPL_ALLOW_PRE_FIND_PACKAGE(HPX  HPX_ALLOW_PREFIND)
IF (HPX_ALLOW_PREFIND)

  MESSAGE("-- Using FIND_PACKAGE(HPX ...) ...") 

  FIND_PACKAGE(HPX)

  IF (HPX_FOUND)
    # Tell TriBITS that we found HPX and there no need to look any further!
    SET(TPL_HPX_INCLUDE_DIRS ${HPX_INCLUDE_DIRS} CACHE PATH
      "HPX include dirs")
    SET(TPL_HPX_LIBRARIES ${HPX_LIBRARIES} CACHE FILEPATH
      "HPX libraries")
    SET(TPL_HPX_LIBRARY_DIRS ${HPX_LIBRARY_DIRS} CACHE PATH
      "HPX library dirs")
  ENDIF()

ENDIF()

#
# Third, call TRIBITS_TPL_FIND_INCLUDE_DIRS_AND_LIBRARIES()
#
TRIBITS_TPL_FIND_INCLUDE_DIRS_AND_LIBRARIES( HPX
  REQUIRED_HEADERS ${REQUIRED_HEADERS}
  REQUIRED_LIBS_NAMES ${REQUIRED_LIBS_NAMES}
  )
# NOTE: If FIND_PACKAGE(HPX ...) was called and successfully found HPX, then
# TRIBITS_TPL_FIND_INCLUDE_DIRS_AND_LIBRARIES() will use the already-set
# variables TPL_HPX_INCLUDE_DIRS and TPL_HPX_LIBRARIES and then print them
# out (and set some other standard variables as well).  This is the final
# "hook" into the TriBITS TPL system.
