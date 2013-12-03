!-----------------------------------------------------------------------------!
! \file   harness/Data_Types.f90
! \author Thomas M. Evans
! \date   Thu Sep  3 16:34:02 2009
! \brief  HARNESS_DATA_TYPES module definition.
! \note   Copyright (C) 2009 Oak Ridge National Laboratory, UT-Battelle, LLC.
!-----------------------------------------------------------------------------!
! $Id: template.f90,v 1.1 2009/09/03 01:59:52 9te Exp $
!-----------------------------------------------------------------------------!

!=============================================================================!
! \module harness_data_types
! \brief  Defines FORTRAN data types.
!=============================================================================!

module harness_data_types

  implicit none

  ! FLOATING-POINT precision control: 4, 8, 16 Byte types.
  integer, parameter :: REAL4  = selected_real_kind(6, 37)
  integer, parameter :: REAL8  = selected_real_kind(13, 300)
  integer, parameter :: REAL16 = selected_real_kind(27, 2400)

  ! INTEGER precision control: 1, 2, 4, 8 Byte types.
  integer, parameter :: INT1 = selected_int_kind(2)
  integer, parameter :: INT2 = selected_int_kind(4)
  integer, parameter :: INT4 = selected_int_kind(9)
  integer, parameter :: INT8 = selected_int_kind(18)

  ! LOGICAL default
  integer, parameter :: LOGIC = kind(.true.)
  
  ! CHARACTER default
  integer, parameter :: CHAR = kind("len=1")

end module harness_data_types

!-----------------------------------------------------------------------------!
!                               end of Data_Types.f90
!-----------------------------------------------------------------------------!

