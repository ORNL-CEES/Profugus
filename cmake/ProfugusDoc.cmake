##---------------------------------------------------------------------------##
## ProfugusDoc.cmake
## Tribits-based Denovo documentation utility macros
## Seth R Johnson
## Tuesday August 7 16:18:45 2012
##---------------------------------------------------------------------------##

find_package(Doxygen)
if(DOXYGEN_FOUND)
  SET(ENABLE_DOCUMENTATION OFF CACHE BOOL
    "Enable documentation for the ${PROJECT_NAME} project." )
  SET(ENABLE_DOCUMENTATION_PDF OFF CACHE BOOL
    "Enable PDF Doxygen processing.")
  SET(ENABLE_DOCUMENTATION_XML OFF CACHE BOOL
    "Enable XML Doxygen processing.")

  # if pdf documentation is on, make sure documentation is on
  if(ENABLE_DOCUMENTATION_PDF AND NOT ENABLE_DOCUMENTATION)
      SET(ENABLE_DOCUMENTATION ON)
  endif()

  # if xml documentation is on, make sure documentation is on
  if(ENABLE_DOCUMENTATION_XML AND NOT ENABLE_DOCUMENTATION)
      SET(ENABLE_DOCUMENTATION ON)
  endif()

  if(ENABLE_DOCUMENTATION)
    add_custom_target(doc
      COMMENT "Creating documentation"
      )
  endif()
endif()

include(CMakeParseArguments)

# This global property takes advantage of the Tribits call sequence to track
# dependencies between documented packages alongside the actual package
# dependencies
set_property(GLOBAL PROPERTY DOCUMENTED_SUBPACKAGE_FULLNAMES "")

################################################################################
# Set up documentation for a particular subpackage
# This assumes SUBPACKAGE_NAME is set by Tribits, along with all the
# accompanying properties.
#
# DENOVO_ADD_DOC(
#   [IS_PACKAGE_LEVEL]
#   [LIST_EXTERNALS]
#   [INPUTS dir|file [...]]
#   [CONF_INPUTS file [...]]
#   [IMAGE_PATH dir]
#   [DEPENDS name1 [name2...] ]
#   )
#
# IS_PACKAGE_LEVEL  Indicate that this is a *package* we're documenting
#                   rather than a subpackage
#
# LIST_EXTERNALS    Show classes from dependent projects (via tag files) in the
#                   current documentation directory's class list
#
# INPUTS            Directories and files to be passed as inputs to the Doxygen
#                   configuration, as paths relative to the current source dir
#
# CONF_INPUTS       Directories and files to be passed as inputs to the Doxygen
#                   configuration, as paths relative to the current *bin* dir
#
# IMAGE_PATH        Directory name where included images can be found
#
# DEPENDS           Extra paths whose modification should force the
#                   documentation to regenerate
#
# If neither INPUTS nor CONF_INPUTS is set, it defaults to searching the current
# directory and the 'doc' subdirectory for valid input.

# Use the true macro depending on whether Doxygen is found.
if(DOXYGEN_FOUND AND ENABLE_DOCUMENTATION)
  macro(DENOVO_ADD_DOC)
    _DENOVO_ADD_DOC(${ARGN})
  endmacro()
else()
  macro(DENOVO_ADD_DOC)
    # no-op
    #message(STATUS "skipping documentation for ${SUBPACKAGE_FULLNAME}")
  endmacro()
endif()


macro(_DENOVO_ADD_DOC)
  # Prefixes of "DDOC" are parsed arguments
  # Prefixes of "EXNIL_DOC" are used in the Doxyfile.in configuration
  cmake_parse_arguments(DDOC
    "IS_PACKAGE_LEVEL;LIST_EXTERNALS"
    "IMAGE_PATH"
    "INPUTS;CONF_INPUTS;DEPENDS"
    ${ARGN})

  if(DDOC_IS_PACKAGE_LEVEL)
    # we're not a subpackage, and therefore ${SUBPACKAGE_FULLNAME} is invalid
    set(SUBPACKAGE_NAME "")
    set(SUBPACKAGE_FULLNAME ${PACKAGE_NAME})
  endif()

  #message(STATUS "Configuring documentation in ${SUBPACKAGE_FULLNAME}")

  # >>> SET UP DOCUMENTED SUBPACKAGES
  # - Get the global list of documented subpackages; add the current subpackage
  #   to that list
  # - Get the Tribits dependency list for the current subpackage
  # - For each dependency, see if it is in the list of documented subpackages
  #
  # result: _SUBPACKAGE_DEPENDS is set to documented subpackages that are
  #         explicit dependencies of the current subpackage
  get_property(_ALL_DOCUMENTED_SUBPACKAGES GLOBAL
    PROPERTY DOCUMENTED_SUBPACKAGE_FULLNAMES)
  set_property(GLOBAL APPEND
    PROPERTY DOCUMENTED_SUBPACKAGE_FULLNAMES ${SUBPACKAGE_FULLNAME})

  set(_SUBPACKAGE_DEPENDS)
  foreach(_SPFN ${${SUBPACKAGE_FULLNAME}_LIB_REQUIRED_DEP_PACKAGES})
    # See if each dependency has documentation
    list(FIND _ALL_DOCUMENTED_SUBPACKAGES ${_SPFN} _FOUND_SUBPACKAGE)
    if(NOT _FOUND_SUBPACKAGE EQUAL -1)
      # The dependency has documentation; add a line in the tag file
      set(_TAGFILE "${PROJECT_BINARY_DIR}/doc/${_SPFN}.tag")
      set(_BINDIR ${${_SPFN}_BINARY_DIR})
      set(EXNIL_DOC_IN_TAGFILES "${EXNIL_DOC_IN_TAGFILES} \\
        ${_TAGFILE}=${_BINDIR}/html")
      list(APPEND _SUBPACKAGE_DEPENDS ${_SPFN})
    endif()
  endforeach()

  # >>> SET UP DOXYFILE INPUTS
  set(EXNIL_DOC_PACKAGE_NAME "${PACKAGE_NAME}: ${SUBPACKAGE_NAME}")
  set(EXNIL_DOC_PACKAGE_VERSION "${${PROJECT_NAME}_VERSION_STRING}")
  set(EXNIL_DOC_WARN_IF_UNDOCUMENTED "YES")
  set(EXNIL_DOC_GENERATE_TAGFILE "${PROJECT_BINARY_DIR}/doc/${SUBPACKAGE_FULLNAME}.tag")
  set(EXNIL_DOC_IMAGE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/doc ${DDOC_IMAGE_PATH} ${PARENT_PACKAGE_SOURCE_DIR}")
  set(EXNIL_DOC_LIST_EXTERNALS ${DDOC_LIST_EXTERNALS})

  # >>> Turn on PDF if requested
  set(EXNIL_DOXYGEN_PDF "NO")
  if(ENABLE_DOCUMENTATION_PDF)
    set(EXNIL_DOXYGEN_PDF "YES")
  endif()

  # >>> Turn on XML if requested
  set(EXNIL_DOXYGEN_XML "NO")
  if(ENABLE_DOCUMENTATION_XML)
    set(EXNIL_DOXYGEN_XML "YES")
  endif()

  # > Set up input list
  # If no inputs are supplied, assume we're looking for stuff in the current
  # source directory and the "doc" subdirectory
  set(EXNIL_DOC_INPUT_FILES)
  if(NOT DDOC_INPUTS AND NOT DDOC_CONF_INPUTS)
    set(DDOC_INPUTS . doc)
  endif()
  # Paths relative to subpackage source directory
  foreach(input_value ${DDOC_INPUTS})
    set(EXNIL_DOC_INPUT_FILES "${EXNIL_DOC_INPUT_FILES} \\
    ${${SUBPACKAGE_FULLNAME}_SOURCE_DIR}/${input_value}"
      )
  endforeach(input_value)
  # Paths relative to current binary directory
  foreach(input_value ${DDOC_CONF_INPUTS})
    set(EXNIL_DOC_INPUT_FILES "${EXNIL_DOC_INPUT_FILES} \\
    ${CMAKE_CURRENT_BINARY_DIR}/${input_value}"
      )
  endforeach(input_value)

  # > Set up tag file list
  set(EXNIL_DOC_IN_TAGFILES) # Series of of Doxygen tag files that we use
  foreach(_SPFN ${_SUBPACKAGE_DEPENDS})
    set(_TAGFILE "${PROJECT_BINARY_DIR}/doc/${_SPFN}.tag")
    set(_BINDIR ${${_SPFN}_BINARY_DIR})
    set(EXNIL_DOC_IN_TAGFILES "${EXNIL_DOC_IN_TAGFILES} \\
      ${_TAGFILE}=${_BINDIR}/html")
  endforeach()

  # >>> CONFIGURE DOXYFILE
  set(_CONFIGURED_DOXYFILE "${CMAKE_CURRENT_BINARY_DIR}/Doxyfile")
  configure_file(
    ${PROJECT_SOURCE_DIR}/doc/Doxyfile.in
    ${_CONFIGURED_DOXYFILE}
    @ONLY #operate on @VAR@, not ${VAR}
    )

  # >>> CREATE TARGET AND DEPENDENCIES
  set(_TARGET "doc_${SUBPACKAGE_FULLNAME}")
  add_custom_command(
    OUTPUT ${EXNIL_DOC_GENERATE_TAGFILE}
    COMMAND ${DOXYGEN_EXECUTABLE} ${_CONFIGURED_DOXYFILE}
    DEPENDS ${_CONFIGURED_DOXYFILE} ${DDOC_DEPENDS}
    COMMENT "Generating ${PACKAGE_NAME} ${SUBPACKAGE_NAME} docs with Doxygen"
    )
  add_custom_target(${_TARGET}
    COMMENT "Created documentation for ${PACKAGE_NAME} ${SUBPACKAGE_NAME}"
    DEPENDS ${EXNIL_DOC_GENERATE_TAGFILE}
    )

  foreach(_SPFN ${_SUBPACKAGE_DEPENDS})
    #message(STATUS "Adding dependency on doc_${_SPFN} to ${_TARGET}")
    add_dependencies(${_TARGET} "doc_${_SPFN}")
  endforeach()

  add_dependencies(doc ${_TARGET})
endmacro()

##---------------------------------------------------------------------------##
## end of ProfugusDoc.cmake
##---------------------------------------------------------------------------##
