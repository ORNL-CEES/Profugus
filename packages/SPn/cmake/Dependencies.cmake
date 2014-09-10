##---------------------------------------------------------------------------##
## SPn/cmake/Dependencies.cmake
## Thomas M. Evans
## Saturday July 14 16:54:41 2012
##---------------------------------------------------------------------------##

##---------------------------------------------------------------------------##
## Define sub-packages
##---------------------------------------------------------------------------##

SET(LIB_REQUIRED_DEP_PACKAGES
  Matprop Utils Teuchos Epetra Thyra Stratimikos AztecOO Belos
  Ifpack Anasazi Tpetra Ifpack2)

SET(LIB_OPTIONAL_DEP_PACKAGES
    ML MCLS MueLu)

SET(TEST_REQUIRED_DEP_PACKAGES)
SET(TEST_OPTIONAL_DEP_PACKAGES)
SET(LIB_REQUIRED_DEP_TPLS)
SET(LIB_OPTIONAL_DEP_TPLS)
SET(TEST_REQUIRED_DEP_TPLS)
SET(TEST_OPTIONAL_DEP_TPLS)

##---------------------------------------------------------------------------##
##                  end of SPn/cmake/Dependencies.cmake
##---------------------------------------------------------------------------##
