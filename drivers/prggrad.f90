!**************************************************************
!   file:               prggrad.f90
!   purpose:            driver to compute the gradient of the costfunction
!   creation date:      07 02
!   revised:            04 04  parameter info included
!                       04 04  handing of NaNs included
!
!   Copyright (C) 2000, 2001
!   FastOpt, Drs. Ralf Giering und Thomas Kaminski GbR
!          , Hamburg, Germany
!
!   Email: info@FastOpt.de
!
!   All rights reserved.
!**************************************************************
program fgrad
!*************************************************************
  implicit none

!=========================================
! declaration
!=========================================
  integer :: n
  logical, parameter :: info = .true.

!==============================================================
! get the number of control variables
!==============================================================
  call numbmod( n )

!-----------------------------------------
! call the subroutine
!-----------------------------------------
call doit( n, info )

end program fgrad

subroutine doit( n, info )
  implicit none
  integer ::  n
  logical :: info

  REAL    :: x(n)
  REAL    :: adx(n)

  integer :: nn, m, jmin, jmax, ifail
  REAL    :: fc, objf, adobjf
  integer :: i
  logical :: cold
  integer :: ifunc
! in case of +/- NaN we assign +/- huge/factor
  integer, parameter :: factor = 1e6  

!==============================================================
!==============================================================
!txk all relevant initialisation in numbmod  call initmod( n, x )
!HEW-CHG-041006    call instore( nn, ifunc, fc, m, jmin, jmax, cold, ifail )
!HEW-CHG-041201
  call instoref( nn, ifunc, fc, m, jmin, jmax, cold, ifail )
  if (ifail .ne. 0) then
     print *, 'ERROR reading startup file'
     stop 'error'
  end if
  if (nn .ne. n) then
     print *, 'ERROR: Number of params n (=',n,') differs on startup file (=',nn,')'
     stop 'error'
  end if
!HEW-CHG-041201      call dostore( n, x, .false., 1 )
!
     call dostoref( n, x, .false., 1 , 'OPWARMD' )

  adx(:) = 0.0
  adobjf = 1.0

  if (info) print *,' running adjoint model'
  call model_ad( n, x, adx, objf, adobjf )
  if (objf.lt.huge(objf)) then
  else
     if (info) print*, 'Value of cf = ',objf
     if (info) print*, 'Replaced by = ',huge(objf)/factor
     objf=huge(objf)/factor
  endif

  do i =1, n
     if (adx(i).lt.huge(adx(i))) then
!        print*, 'Value of gradient component ',i, ' = ',adx(i)
!        print*, 'less than ',huge(adx(i))
!        print*, 'unchanged '
     else
        if (info) print*, 'Value of gradient component ',i, ' = ',adx(i)
        if (info) print*, 'Replaced by = ',huge(adx(i))/factor
        adx(i)=huge(adx(i))/factor
     endif
     if (adx(i).gt.-huge(adx(i))) then
!        print*, 'Value of gradient component ',i, ' = ',adx(i)
!        print*, 'greater than ',-huge(adx(i))
!        print*, 'unchanged '
     else
        if (info) print*, 'Value of gradient component ',i, ' = ',adx(i)
        if (info) print*, 'Replaced by = ',-1*huge(adx(i))/factor
        adx(i)=-1*huge(adx(i))/factor
     endif
!    print*, 'Final value = ', adx(i)
  enddo

if (info) print *,' gradient ready, start nextx'

!HEW-CHG-041201call outstor( n, ifunc, objf, m, jmin, jmax )
!
    call outstorf( n, ifunc, objf, m, jmin, jmax )
!HEW-CHG-041201        call dostore( n, adx, .true., 2 )
!
     call dostoref( n, adx, .true., 2 , 'OPWARMD' )

end subroutine doit
