module CSR_Type

    use :: ISO_Precisions

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!132
!!
!!  Purpose:  Contains complex sparse matrix types for use in Nullarbor.
!!
!! Includes:
!!            CSR - Compressed Sparse Rows
!!
!!    Author:  Edric Matwiejew
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!132

  implicit none

  private

  real :: default_precision

  public :: CSR

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!132

    type :: CSR(working_precision)

!!  Purpose:    Provides a compressed sparse rows (CSR) complex matrix derived type.
!!
!!  Definitions:    working_precision - Precision of non-zero entries stored in values. Defaults to single precision.
!!                  rows - Matrix row dimension.
!!                  columns - Matrix column dimension.
!!                  row_start - start and end indexes for rows in col_index and values.
!!                  col_index - col index of non-zero entries.
!!                  values - Non-zero entries, stored row-wise.

      integer, kind :: working_precision = precision(default_precision)
      integer :: rows, columns
      integer, dimension(:), allocatable :: row_start
      integer, dimension(:), allocatable :: col_index
      complex(working_precision), dimension(:), allocatable :: values
      logical :: symmetric

    end type CSR

end module CSR_Type
