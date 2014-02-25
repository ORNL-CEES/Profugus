##---------------------------------------------------------------------------##
## mesh/test/CMakeLists.txt
## Thomas M. Evans
## Wednesday July 11 14:35:42 2012
##---------------------------------------------------------------------------##
## Copyright (C) 2012 Oak Ridge National Laboratory, UT-Battelle, LLC.
##---------------------------------------------------------------------------##
## CMAKE for spn
##---------------------------------------------------------------------------##

INCLUDE(TribitsCopyFilesToBinaryDir)
INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(UtilsTest)

##---------------------------------------------------------------------------##
## TESTING
##---------------------------------------------------------------------------##

ADD_UTILS_TEST(tstDecomposition.cc              )
ADD_UTILS_TEST(tstShiftedOperator.cc            )
ADD_UTILS_TEST(tstShiftedInverseOperator.cc     )
ADD_UTILS_TEST(tstStratimikosSolver.cc      NP 1   4)
ADD_UTILS_TEST(tstArnoldi.cc                NP 1 2  )
ADD_UTILS_TEST(tstInverseOperator.cc        NP 1   4)

#ADD_NEMESIS_TEST(tstAztec.cc                        )
#ADD_NEMESIS_TEST(tstBiCGStab.cc                     )
#ADD_NEMESIS_TEST(tstEigenvalueSolverBuilder.cc      )
#ADD_NEMESIS_TEST(tstLinearSolverBuilder.cc          )
#ADD_NEMESIS_TEST(tstPowerIteration.cc               )
#ADD_NEMESIS_TEST(tstPreconditionerBuilder.cc        )
#ADD_NEMESIS_TEST(tstRayleighQuotient.cc             )
#ADD_NEMESIS_TEST(tstRichardson.cc                   )

##---------------------------------------------------------------------------##
## FILES
##---------------------------------------------------------------------------##

TRIBITS_COPY_FILES_TO_BINARY_DIR(XMLSOLVERINP
  SOURCE_FILES
  aztecoo.xml
  belos.xml
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
  EXEDEPS tstStratimikosSolver tstInverseOperator
  )

##---------------------------------------------------------------------------##
##                    end of mesh/test/CMakeLists.txt
##---------------------------------------------------------------------------##
