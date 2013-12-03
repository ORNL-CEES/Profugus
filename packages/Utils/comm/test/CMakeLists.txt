#---------------------------------------------------------------------------##
## comm/test/CMakeLists.txt
## Thomas M. Evans
## Wednesday July 11 14:35:42 2012
##---------------------------------------------------------------------------##
## Copyright (C) 2011 Oak Ridge National Laboratory, UT-Battelle, LLC.
##---------------------------------------------------------------------------##
## CMAKE for comm
##---------------------------------------------------------------------------##

INCLUDE(NemesisTest)

##---------------------------------------------------------------------------##
## TESTS
##---------------------------------------------------------------------------##


ADD_NEMESIS_TEST(tstP_Stream.cc           LEGACY)
ADD_NEMESIS_TEST(tstAbort.cc              LEGACY)
ADD_NEMESIS_TEST(tstBroadcast.cc          LEGACY)
ADD_NEMESIS_TEST(tstComm_Dup.cc           LEGACY)
ADD_NEMESIS_TEST(tstDefault.cc            LEGACY)
ADD_NEMESIS_TEST(tstParallel_Unit_Test.cc LEGACY)
ADD_NEMESIS_TEST(tstParallel_Utils.cc     LEGACY)
ADD_NEMESIS_TEST(tstPingPong.cc           LEGACY)
ADD_NEMESIS_TEST(tstReduction.cc          LEGACY)
ADD_NEMESIS_TEST(tstRequest.cc            LEGACY)
ADD_NEMESIS_TEST(tstSerial.cc             LEGACY)
ADD_NEMESIS_TEST(tstSpinLock.cc           LEGACY)
ADD_NEMESIS_TEST(tstSplit.cc              LEGACY)
ADD_NEMESIS_TEST(tstGather.cc             LEGACY)

IF(HAVE_SYS_TIMES_H)
  ADD_NEMESIS_TEST(tstTime.cc          LEGACY)
  ADD_NEMESIS_TEST(tstTiming.cc   NP 1 LEGACY)
ENDIF()

# For now we are leaving the FORTRAN bindings off because we don't use them
#ADD_NEMESIS_TEST(tstComm_FC.cc LEGACY)

##---------------------------------------------------------------------------##
##                   end of comm/test/CMakeLists.txt
##---------------------------------------------------------------------------##
