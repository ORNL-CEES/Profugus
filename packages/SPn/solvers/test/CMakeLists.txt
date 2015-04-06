##---------------------------------------------------------------------------##
## mesh/test/CMakeLists.txt
## Thomas M. Evans
## Wednesday July 11 14:35:42 2012
##---------------------------------------------------------------------------##
## Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
##---------------------------------------------------------------------------##
## CMAKE for spn
##---------------------------------------------------------------------------##

INCLUDE(UtilsTest)

##---------------------------------------------------------------------------##
## TESTING
##---------------------------------------------------------------------------##

ADD_UTILS_TEST(tstDecomposition.cc           )
ADD_UTILS_TEST(tstShiftedOperator.cc         )
ADD_UTILS_TEST(tstShiftedInverseOperator.cc  )
ADD_UTILS_TEST(tstRichardson.cc              )
ADD_UTILS_TEST(tstPowerIteration.cc          )
ADD_UTILS_TEST(tstRayleighQuotient.cc        )
ADD_UTILS_TEST(tstDavidsonEigensolver.cc     )
ADD_UTILS_TEST(tstPreconditionerBuilder.cc   )
ADD_UTILS_TEST(tstLinearSolverBuilder.cc     )
ADD_UTILS_TEST(tstEigenvalueSolverBuilder.cc )
ADD_UTILS_TEST(tstStratimikosSolver.cc       )
ADD_UTILS_TEST(tstBelosSolver.cc             )
ADD_UTILS_TEST(tstArnoldi.cc                 )
ADD_UTILS_TEST(tstInverseOperator.cc         )
ADD_UTILS_TEST(tstModelEvaluatorWrapper.cc   )
ADD_UTILS_TEST(tstAndersonSolver.cc          )

##---------------------------------------------------------------------------##
## FILES
##---------------------------------------------------------------------------##

TRIBITS_COPY_FILES_TO_BINARY_DIR(XMLSOLVERINP
  SOURCE_FILES
  aztecoo.xml
  belos.xml
  mcls.xml
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
  EXEDEPS tstStratimikosSolver tstInverseOperator
  )

##---------------------------------------------------------------------------##
##                    end of mesh/test/CMakeLists.txt
##---------------------------------------------------------------------------##
