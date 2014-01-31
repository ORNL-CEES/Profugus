##---------------------------------------------------------------------------##
## xs/test/CMakeLists.txt
## Thomas M. Evans
## Tuesday December 3 20:32:38 2013
##---------------------------------------------------------------------------##

INCLUDE(TribitsCopyFilesToBinaryDir)
INCLUDE(UtilsTest)

##---------------------------------------------------------------------------##
## TESTS
##---------------------------------------------------------------------------##

ADD_UTILS_TEST(tstTeuchos.cc NP 1)
ADD_UTILS_TEST(tstXS.cc      NP 1)

##---------------------------------------------------------------------------##
## COPY FILES
##---------------------------------------------------------------------------##

TRIBITS_COPY_FILES_TO_BINARY_DIR(XS_CROSS_SECTION_XML
  SOURCE_FILES three_group_xs.xml
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
  EXEDEPS tstTeuchos)

##---------------------------------------------------------------------------##
## end of xs/test/CMakeLists.txt
##---------------------------------------------------------------------------##
