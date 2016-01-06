##---------------------------------------------------------------------------##
## xs/test/CMakeLists.txt
## Thomas M. Evans
## Tuesday December 3 20:32:38 2013
##---------------------------------------------------------------------------##

INCLUDE(UtilsTest)

##---------------------------------------------------------------------------##
## TESTS
##---------------------------------------------------------------------------##

ADD_UTILS_TEST(tstTeuchos.cc         NP 1    )
ADD_UTILS_TEST(tstXS.cc              NP 1    )
ADD_UTILS_TEST(tstMat_DB.cc          NP 1    )
ADD_UTILS_TEST(tstEnergy_Collapse.cc NP 1    )
ADD_UTILS_TEST(tstXS_Builder.cc      NP 1 2 4)

##---------------------------------------------------------------------------##
## COPY FILES
##---------------------------------------------------------------------------##

TRIBITS_COPY_FILES_TO_BINARY_DIR(XS_CROSS_SECTION_XML
  SOURCE_FILES xs3GP0.xml xs5GP1.xml
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
  EXEDEPS tstTeuchos tstXS_Builder tstMat_DB)

##---------------------------------------------------------------------------##
## end of xs/test/CMakeLists.txt
##---------------------------------------------------------------------------##
