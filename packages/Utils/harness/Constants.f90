!-----------------------------------------------------------------------------!
! \file   harness/Constants.f90
! \author Thomas M. Evans
! \date   Thu Sep  3 16:37:02 2009
! \brief  HARNESS_CONSTANTS module definition.
! \note   Copyright (C) 2009 Oak Ridge National Laboratory, UT-Battelle, LLC.
!-----------------------------------------------------------------------------!
! $Id: template.f90,v 1.1 2009/09/03 01:59:52 9te Exp $
!-----------------------------------------------------------------------------!

!=============================================================================!
! \module harness_constants
! \brief  Defines useful FORTRAN constants.
!=============================================================================!

module harness_constants

  use harness_data_types, only : REAL8
  implicit none

  ! Fractions.
  real(REAL8), parameter :: zero        = 0.0_REAL8
  real(REAL8), parameter :: one         = 1.0_REAL8
  real(REAL8), parameter :: two         = 2.0_REAL8
  real(REAL8), parameter :: three       = 3.0_REAL8
  real(REAL8), parameter :: four        = 4.0_REAL8
  real(REAL8), parameter :: one_half    = one / 2.0_REAL8
  real(REAL8), parameter :: one_third   = one / 3.0_REAL8
  real(REAL8), parameter :: one_fourth  = one / 4.0_REAL8
  real(REAL8), parameter :: one_fifth   = one / 5.0_REAL8
  real(REAL8), parameter :: one_sixth   = one / 6.0_REAL8
  real(REAL8), parameter :: one_seventh = one / 7.0_REAL8
  real(REAL8), parameter :: one_eighth  = one / 8.0_REAL8

  ! Constants.
  real(REAL8), parameter :: pi          = 3.1415926535897931159979635_REAL8
  real(REAL8), parameter :: two_pi      = two * pi
  real(REAL8), parameter :: four_pi     = four * pi
  real(REAL8), parameter :: inv_two_pi  = one / two_pi
  real(REAL8), parameter :: inv_four_pi = one / four_pi

end module harness_constants

!-----------------------------------------------------------------------------!
!                               end of Constants.f90
!-----------------------------------------------------------------------------!

