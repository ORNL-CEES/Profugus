##---------------------------------------------------------------------------##
## Utils/cmake/Dependencies.cmake
## Thomas M. Evans
## Saturday July 14 16:54:41 2012
##---------------------------------------------------------------------------##

##---------------------------------------------------------------------------##
## Define sub-packages
##---------------------------------------------------------------------------##

SET(SUBPACKAGES_DIRS_CLASSIFICATIONS_OPTREQS
  Tools   tools     SS  REQUIRED
  Harness harness   SS  REQUIRED
  Release release   SS  REQUIRED
  Comm    comm      SS  REQUIRED
  Gtest   gtest     SS  REQUIRED
  )

SET(LIB_REQUIRED_DEP_PACKAGES)
SET(LIB_OPTIONAL_DEP_PACKAGES)
SET(TEST_REQUIRED_DEP_PACKAGES)
SET(TEST_OPTIONAL_DEP_PACKAGES)
SET(LIB_REQUIRED_DEP_TPLS)
SET(LIB_OPTIONAL_DEP_TPLS)
SET(TEST_REQUIRED_DEP_TPLS)
SET(TEST_OPTIONAL_DEP_TPLS)

##---------------------------------------------------------------------------##
##                  end of Utils/cmake/Dependencies.cmake
##---------------------------------------------------------------------------##

