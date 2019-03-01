!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!132
!!
!!    Module:     ISO_Precisions   
!!
!!    Purpose:    Defines single precision (sp), double precision (dp) and quadrupule presision (qp) kind types using the 
!!                iso_fortran_env module types real32, real64 and real128 respectively.
!!
!!    Author:     Edric Matiwejew
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!132

module ISO_Precisions

    use, intrinsic :: iso_fortran_env, only: sp => real32, dp => real64, qp => real128

end module ISO_Precisions


