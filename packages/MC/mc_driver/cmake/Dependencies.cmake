##---------------------------------------------------------------------------##
## MC/cmake/Dependencies.cmake
## Thomas M. Evans
## Saturday July 14 16:54:41 2012
##---------------------------------------------------------------------------##

##---------------------------------------------------------------------------##
## Define sub-packages
##---------------------------------------------------------------------------##

SET(LIB_REQUIRED_DEP_PACKAGES
  Teuchos Epetra Matprop Utils SPn MCmc_geometry MCmc_physics)

SET(LIB_OPTIONAL_DEP_PACKAGES)

SET(TEST_REQUIRED_DEP_PACKAGES)
SET(TEST_OPTIONAL_DEP_PACKAGES)
SET(LIB_REQUIRED_DEP_TPLS Boost HPX)
SET(LIB_OPTIONAL_DEP_TPLS)
SET(TEST_REQUIRED_DEP_TPLS)
SET(TEST_OPTIONAL_DEP_TPLS)

##---------------------------------------------------------------------------##
##                  end of MC/cmake/Dependencies.cmake
##---------------------------------------------------------------------------##
