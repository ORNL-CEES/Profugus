##---------------------------------------------------------------------------##
## Alea/cmake/Dependencies.cmake
## Thomas M. Evans
## Saturday July 14 16:54:41 2012
##---------------------------------------------------------------------------##

##---------------------------------------------------------------------------##
## Define sub-packages
##---------------------------------------------------------------------------##

SET(LIB_REQUIRED_DEP_PACKAGES
    Utils SPn Teuchos Tpetra Ifpack2 KokkosCore KokkosAlgorithms)

SET(LIB_OPTIONAL_DEP_PACKAGES Belos Temere MCLS)
SET(TEST_REQUIRED_DEP_PACKAGES)
SET(TEST_OPTIONAL_DEP_PACKAGES)
SET(LIB_REQUIRED_DEP_TPLS)
SET(LIB_OPTIONAL_DEP_TPLS)
SET(TEST_REQUIRED_DEP_TPLS)
SET(TEST_OPTIONAL_DEP_TPLS)

##---------------------------------------------------------------------------##
##                  end of Alea/cmake/Dependencies.cmake
##---------------------------------------------------------------------------##
