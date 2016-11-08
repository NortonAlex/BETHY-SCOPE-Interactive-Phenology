!**************************************************************
!   file: 	        prgcost2.f90
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
program prgcost2
!*************************************************************
  
  implicit none
  
!=========================================
! declaration
!=========================================
  integer n

!==============================================================
! get the number of control variables
!==============================================================
  call numbmod( n )
  
!-----------------------------------------
! call the subroutine
!-----------------------------------------
  call doit( n )

end program prgcost2
      
subroutine doit( n )

  implicit none
      
  integer n
  REAL    fc, fc2
  REAL   x(n)
  REAL   adx(n)

  real eps, pert

!==============================================================
! initialisize the model
! and set the first guess of the control variables
!==============================================================
  call initmod( n, x )

!===============================================
! fist call of model
!===============================================
  call model( n, x, fc )
  
  write(*,*)
  write(*,*) ' first call of cost function :  fc  = ', fc

!===============================================
! second call of model
!===============================================
  call model( n, x, fc2 )

  write(*,*)
  write(*,*) ' after the second call'
  write(*,*) ' second call of cost function : fc2 = ', fc2
      
!===============================================
! compare
!===============================================
  eps = 10*epsilon(eps)
  write(*,*) 'Comparision of costfunction calculation from first (fc) and second call (fc2) is'
  if (abs(fc-fc2).lt.abs(eps*(1+fc))) then
     write(*,*) 'OK       : |fc-fc2| = ',abs(fc-fc2),' <  |eps*(1+fc)|'
  else
     write(*,*) 'DIFFERENT: |fc-fc2| = ',abs(fc-fc2),' >= |eps*(1+fc)|'     
  endif

!===============================================
! and now a third  call of model, but with x perturbed
!===============================================
  pert = 1.e8*epsilon(eps)
  x = x + pert
  call model( n, x, fc2 )
  write(*,*) 'Comparision of costfunction calculation from first (fc) and third call (fc2) is'
  if (abs(fc-fc2).gt.abs(eps*(1+fc))) then
     write(*,*) 'OK       : |fc-fc2| = ',abs(fc-fc2),' >= |eps*(1+fc)|'
  else
     write(*,*) 'DIFFERENT: |fc-fc2| = ',abs(fc-fc2),' <  |eps*(1+fc)|'     
  endif

!===============================================
! postprocessing
!===============================================
  call postmod( n, x, fc )

  end subroutine doit
