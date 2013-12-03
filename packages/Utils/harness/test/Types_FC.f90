!-----------------------------------------------------------------------------!
! \file   harness/Types_FC.f90
! \author Thomas M. Evans
! \date   Fri Sep  4 10:35:29 2009
! \brief  Tests of FORTRAN types.
! \note   Copyright (C) 2009 Oak Ridge National Laboratory, UT-Battelle, LLC.
!-----------------------------------------------------------------------------!
! $Id: template.f90,v 1.1 2009/09/03 01:59:52 9te Exp $
!-----------------------------------------------------------------------------!

subroutine test_arrays(int4s, real4s, real8s, fails)

  use harness_data_types, only : INT1, INT2, INT4, INT8, REAL4, REAL8
  implicit none

  integer(INT4), dimension(4), intent(inout) :: int4s

  real(REAL4), dimension(3), intent(inout) :: real4s
  real(REAL8), dimension(3), intent(inout) :: real8s
  
  integer, intent(inout) :: fails

  integer :: i

  ! >>> BODY

  if (int4s(1) /= 2) fails = fails + 1
  if (int4s(2) /= 4) fails = fails + 1
  if (int4s(3) /= 6) fails = fails + 1
  if (int4s(4) /= 8) fails = fails + 1
  
  if (bit_size(1_INT1) /= 8)  fails = fails + 1
  if (bit_size(1_INT2) /= 16) fails = fails + 1
  if (bit_size(1_INT4) /= 32) fails = fails + 1
  if (bit_size(1_INT8) /= 64) fails = fails + 1

  if (real4s(1) /= 2.1_REAL4) fails = fails + 1
  if (real4s(2) /= 4.2_REAL4) fails = fails + 1
  if (real4s(3) /= 6.3_REAL4) fails = fails + 1

  if (real8s(1) /= 5.2_REAL8) fails = fails + 1
  if (real8s(2) /= 6.4_REAL8) fails = fails + 1
  if (real8s(3) /= 7.6_REAL8) fails = fails + 1

end subroutine test_arrays

!-----------------------------------------------------------------------------!
!                               end of Types_FC.f90
!-----------------------------------------------------------------------------!

