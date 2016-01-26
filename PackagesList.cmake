##---------------------------------------------------------------------------##
## Profugus/PackagesList.cmake
## Thomas M. Evans
## Monday December 2 21:24:6 2013
##---------------------------------------------------------------------------##


##---------------------------------------------------------------------------##
## PACKAGES PROVIDED
##---------------------------------------------------------------------------##

TRIBITS_REPOSITORY_DEFINE_PACKAGES(
  Utils     packages/Utils     SS
  CudaUtils packages/CudaUtils SS
  Matprop   packages/Matprop   SS
  SPn       packages/SPn       SS
  MC        packages/MC        SS
  Alea      packages/Alea      SS
  )

##---------------------------------------------------------------------------##
## PLATFORM SUPPORT
##---------------------------------------------------------------------------##

TRIBITS_DISABLE_PACKAGE_ON_PLATFORMS(Profugus Windows)

##---------------------------------------------------------------------------##
## end of PackagesList.cmake
##---------------------------------------------------------------------------##
