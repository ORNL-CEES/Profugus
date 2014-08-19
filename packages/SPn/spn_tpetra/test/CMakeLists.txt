##---------------------------------------------------------------------------##
## spn_tpetra/test/CMakeLists.txt
## Thomas M. Evans
## Wednesday July 11 14:35:42 2012
##---------------------------------------------------------------------------##
## Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
##---------------------------------------------------------------------------##
## CMAKE for spn
##---------------------------------------------------------------------------##

INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(UtilsTest)

##---------------------------------------------------------------------------##
## FILES
##---------------------------------------------------------------------------##

# Make testing library
TRIBITS_ADD_LIBRARY(spn_tpetra_test_lib
  SOURCES Test_XS.cc
  TESTONLY
  DEPLIBS SPn)

##---------------------------------------------------------------------------##
## TESTING
##---------------------------------------------------------------------------##

ADD_UTILS_TEST(tstEnergy_Restriction.cc                             )
ADD_UTILS_TEST(tstEnergy_Prolongation.cc                            )
ADD_UTILS_TEST(tstEnergy_Multigrid.cc           DEPLIBS spn_tpetra_test_lib)
ADD_UTILS_TEST(tstLinear_System_FV.cc           DEPLIBS spn_tpetra_test_lib)
ADD_UTILS_TEST(tstFixed_Source_Solver.cc        DEPLIBS spn_tpetra_test_lib)
ADD_UTILS_TEST(tstEigenvalue_Solver.cc          DEPLIBS spn_tpetra_test_lib)

##---------------------------------------------------------------------------##
##                    end of spn_tpetra/test/CMakeLists.txt
##---------------------------------------------------------------------------##
