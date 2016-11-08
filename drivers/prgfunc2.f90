!**************************************************************
!   file: 	        prgfunc2.f90
!   purpose: 	        driver to evaluate
!                       scalar valued function twice
!   creation date:	06 02
!
!   Copyright (C) 2000, 2001
!   FastOpt, Drs. Ralf Giering und Thomas Kaminski GbR
!          , Hamburg, Germany
!
!   Email: info@FastOpt.de
!
!   All rights reserved.
!**************************************************************
program prgfunc2
!*************************************************************
  
  implicit none
  
!=========================================
! declaration
!=========================================
  integer n, m

!==============================================================
! get the number of control variables
!==============================================================
  call setfunc( n, m )
  
!-----------------------------------------
! call the subroutine
!-----------------------------------------
  call doit( n, m )

end program prgfunc2
      
subroutine doit( n, m )

  implicit none
      
  integer n, m
  REAL   x(n),y(m)

!==============================================================
! initialisize the model
! and set the first guess of the control variables
!==============================================================
  call initfunc( n, x )

!===============================================
! fist call of model
!===============================================
  call func( n, x, m, y )
  
  write(*,*)
  write(*,*) ' the function value is : ', y

!===============================================
! second call of model
!===============================================
  call func( n, x, m, y )

  write(*,*)
  write(*,*) ' after the second call'
  write(*,*) ' the function value is : ', y
      
!===============================================
! postprocessing
!===============================================
  call postfunc( n, x, m, y )

  end subroutine doit

