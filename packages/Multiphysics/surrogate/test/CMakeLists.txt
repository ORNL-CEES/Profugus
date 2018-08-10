##---------------------------------------------------------------------------##
## surrogate/test/CMakeLists.txt
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

##---------------------------------------------------------------------------##
## TESTING
##---------------------------------------------------------------------------##

ADD_UTILS_TEST(tstMulti_Pin_Conduction.cc     NP 1              )
ADD_UTILS_TEST(tstMulti_Pin_Subchannel.cc     NP 1              )
ADD_UTILS_TEST(tstMultiphysics_Driver.cc      NP 1              )
ADD_UTILS_TEST(tstSingle_Pin_Conduction.cc    NP 1              )
ADD_UTILS_TEST(tstSingle_Pin_Subchannel.cc    NP 1              )
ADD_UTILS_TEST(tstTwo_Group_Cross_Sections.cc NP 1              )
ADD_UTILS_TEST(tstTwo_Group_Diffusion.cc      NP 1              )
ADD_UTILS_TEST(tstWater_Properties.cc         NP 1              )

##---------------------------------------------------------------------------##
## FILES
##---------------------------------------------------------------------------##

##---------------------------------------------------------------------------##
##                    end of surrogate/test/CMakeLists.txt
##---------------------------------------------------------------------------##
