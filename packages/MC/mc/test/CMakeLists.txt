##---------------------------------------------------------------------------##
## mc/test/CMakeLists.txt
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

##---------------------------------------------------------------------------##
## TESTING
##---------------------------------------------------------------------------##

ADD_UTILS_TEST(tstBank.cc              NP 1              )
ADD_UTILS_TEST(tstSampler.cc           NP 1              )
ADD_UTILS_TEST(tstParticle.cc          NP 1              )
ADD_UTILS_TEST(tstGroup_Bounds.cc      NP 1              )
ADD_UTILS_TEST(tstPhysics.cc           NP 1              )
#ADD_UTILS_TEST(tstFission_Comm.cc      NP 4        LEGACY)
#ADD_UTILS_TEST(tstFission_Rebalance.cc NP 1 2 4 12       )

##---------------------------------------------------------------------------##
##                    end of mc/test/CMakeLists.txt
##---------------------------------------------------------------------------##
