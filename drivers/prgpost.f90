program post
  !*************************************************************
  implicit none

!=========================================
! declaration
!=========================================
  integer :: n

!==============================================================
! get the number of control variables
!==============================================================
  call numbmod( n )

!-----------------------------------------
! call the subroutine
!-----------------------------------------
  call postrun( n )

end program post


subroutine postrun( n )
!*************************************************************
  implicit none

  integer :: n
  real    :: fc = 0.

  real    ::  x(n), adx(n)

!==============================================================
! declare local variables
!==============================================================
  integer, parameter  :: ioptns = 7
  integer, parameter  :: itmp   = 8

  integer :: ifail, nn
  logical :: cold
  integer :: jmin, jmax, nupdate, ifunc

  x = 0.
  adx = 0.

!==============================================================
! initialisize the model
! and set the first guess of the control variables
!==============================================================
  call initmod( n, x )

!==============================================================
! start the optimisation
!==============================================================
!HEW-CHG-040930  call instore( nn, ifunc, fc, nupdate, jmin, jmax, cold, ifail )
  call instoref( nn, ifunc, fc, nupdate, jmin, jmax, cold, ifail )

  IF (cold) THEN
     print *, 'This is a *POST* run, no optimised vector available'
      PRINT *, 'control variables are set by initfunc'
      PRINT *, 'gradient set to zero'
      adx=0.
  ELSE
     IF (n .NE. nn) THEN
        STOP 'inconsitence number of control variables'
     ENDIF
!HEW-CHG-040830
!  call dostore( n, x, .false., 1 )
     call dostoref( n, x, .false., 1 , 'OPWARMD' )
!HEW-CHG-040830
!  call dostore( n, adx, .false., 2 )
     call dostoref( n, adx, .false., 2 , 'OPWARMD' )
     print *, 'control vector, function, and gradient read'
  ENDIF  
  OPEN (39,file='x.b',status='old',form='unformatted',iostat=ifail)
  if (ifail.gt.0) then
  else
     read(39) x
     print*,'post: parameter values overwritten by those from x.b'
  endif
  CLOSE(39)

  call postadm( n, x, fc, adx )

end subroutine postrun
