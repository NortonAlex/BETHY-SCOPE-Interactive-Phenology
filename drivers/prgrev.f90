!**************************************************************
!	purpose: 	compute reverse sensitivity
!	creation date:	06 02
!
!   Copyright (C) 2000, 2001
!   FastOpt, Drs. Ralf Giering und Thomas Kaminski GbR
!          , Hamburg, Germany
!
!   Email: info@FastOpt.de
!
!   All rights reserved.
!**************************************************************
program prgrev
  implicit none
  integer :: n

!==============================================================
! get the number of control variables
!==============================================================
  call numbmod( n )

!==============================================================
! do the actual computation, see subroutine below
!==============================================================
  call doit( n )

end program prgrev

subroutine doit( n )
  implicit none
  integer ::  n
  REAL    :: x(n)
  REAL    :: adx(n)
  REAL    :: objf, adobjf
  integer :: ifail

!==============================================================
! initialise the model and set control variables
!==============================================================
  call initmod( n, x )

  OPEN (39,file='x.b',status='old',form='unformatted',iostat=ifail)
  if (ifail.gt.0) then
  else
     read(39) x
     print*,'post: parameter values overwritten by those from x.b'
  endif
  CLOSE(39)


!==============================================================
! initialise the adjoint variables and compute reverse sensitivity
!==============================================================
  adx(:) = 0.0
  adobjf = 1.0

  call model_ad( n, x, adx, objf, adobjf )

!==============================================================
! postprocessing
!==============================================================
  call postadm( n, x, objf, adx )

end subroutine doit
