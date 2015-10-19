##---------------------------------------------------------------------------##
## mc/test/CMakeLists.txt
## Thomas M. Evans
## Wednesday July 11 14:35:42 2012
##---------------------------------------------------------------------------##
## Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
##---------------------------------------------------------------------------##
## CMAKE for spn
##---------------------------------------------------------------------------##

INCLUDE(UtilsTest)

##---------------------------------------------------------------------------##
## LIBRARIES
##---------------------------------------------------------------------------##

TRIBITS_ADD_LIBRARY(mc_test_lib
  SOURCES
  SourceTestBase.cc
  TransporterTestBase.cc
  TESTONLY )

##---------------------------------------------------------------------------##
## TESTING
##---------------------------------------------------------------------------##

ADD_UTILS_TEST(tstBank.cc                  NP 1              )
ADD_UTILS_TEST(tstSampler.cc               NP 1              )
ADD_UTILS_TEST(tstParticle.cc              NP 1              )
ADD_UTILS_TEST(tstGroup_Bounds.cc          NP 1              )
ADD_UTILS_TEST(tstPhysics.cc               NP 1              )
ADD_UTILS_TEST(tstVR_Roulette.cc           NP 1              )
ADD_UTILS_TEST(tstCell_Tally.cc            NP 1 4            )

ADD_UTILS_TEST(tstTallier.cc                    )

ADD_UTILS_TEST(tstSource_Transporter.cc       DEPLIBS mc_test_lib)
ADD_UTILS_TEST(tstDomain_Transporter.cc  NP 1 DEPLIBS mc_test_lib)
ADD_UTILS_TEST(tstUniform_Source.cc           DEPLIBS mc_test_lib)
ADD_UTILS_TEST(tstFixed_Source_Solver.cc      DEPLIBS mc_test_lib)

##---------------------------------------------------------------------------##
##                    end of mc/test/CMakeLists.txt
##---------------------------------------------------------------------------##
