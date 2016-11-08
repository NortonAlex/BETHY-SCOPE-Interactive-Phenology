!**************************************************************
!	purpose: 	compute jacobian in reverse mode
!	creation date:	10/02
!
!   Copyright (C) 2000, 2001
!   FastOpt, Drs. Ralf Giering und Thomas Kaminski GbR
!          , Hamburg, Germany
!
!   Email: info@FastOpt.de
!
!   All rights reserved.
!**************************************************************
program prgjacrevpad
  implicit none
  integer :: n,m

!==============================================================
! get the number of the independent and dependent variables
!==============================================================
  call setfunc( n, m )

!==============================================================
! do the actual computation, see subroutine below
!==============================================================
  call doit( n, m )

end program prgjacrevpad

subroutine doit( n, m )
  implicit none
  integer ::  n, m
  REAL   x(n), y(m)
  REAL   adx(m,n)
  REAL   ady(m,m)
  integer j

!==============================================================
! initialise the model and set control variables
!==============================================================
  call initfunc( n, x )

!==============================================================
! initialise the adjoint variables and compute reverse sensitivity
!==============================================================
  adx = 0.0
  ady = 0.0
  do j =1, m
     ady(j,j) = 1.0
  enddo

  call func_pad( m, n, x, adx, m, m, y, ady, m )

!==============================================================
! postprocessing
!==============================================================
  call POSTJACREV(  n, x, m, y, adx)
  call POSTFUNC(  n, x, m, y )

end subroutine doit
