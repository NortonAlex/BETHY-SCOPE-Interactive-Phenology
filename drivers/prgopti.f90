!**************************************************************
!   file:               prgopti.f90
!   purpose:            driver to determine next control vector x
!                        in current optimisation iteration step
!   creation date:      07 02
!   revised:            04 04  parameter info included
!
!   Copyright (C) 2000, 2001
!   FastOpt, Drs. Ralf Giering und Thomas Kaminski GbR
!          , Hamburg, Germany
!
!   Email: info@FastOpt.de
!
!   All rights reserved.
!**************************************************************
program main
!*************************************************************
  implicit none

!=========================================
! declaration
!=========================================
  integer :: n
  logical, parameter :: info = .false.

!==============================================================
! get the number of control variables
!==============================================================
  call numbmod( n )
  if (info) print*,' nextx: numbmod completed '

!-----------------------------------------
! call the subroutine
!-----------------------------------------
  call doit( n, info )

end program main

subroutine doit( n, info )
!*************************************************************
  implicit none

  integer :: n
  real    :: objf
  logical :: info

  real    ::  x(n), adx(n)

!==============================================================
! initialisize the model
! and set the first guess of the control variables
!==============================================================
  call initmod( n, x )
  if (info) print*,' nextx: initmod completed '

!==============================================================
! optimize the costfunction using control variables
!==============================================================
  call optimum( n, x, objf, adx, info )

end subroutine doit


subroutine optimum( n, x, objf, adx, info )
!********************************************************************
!
! input  : n      = amount of control variables
!          x      = vector of first guess of control variables
!
! output : x      = optimal control variables
!          objf   = minmal costfunction
!          adx    = gradient vector at optimum
!
!********************************************************************
  implicit none

!==============================================================
! declare call parameter
!==============================================================
  integer :: n
  REAL    :: x(n)
  REAL    :: objf
  REAL    :: adx(n)
  logical :: info

!==============================================================
! declare local variables
!==============================================================
  integer, parameter  :: ioptns = 7
  integer, parameter  :: itmp   = 8

  REAL   :: xx(n)
  REAL   :: rz1(n)
  REAL   :: rz2(n)
  REAL   :: rz3(n)

  REAL    :: wrz(15)
  integer :: iwrz(4)

  integer :: niter, nfunc, iprint
  REAL    :: epsf, epsx, fmin, epsg
  integer :: nupdate
  REAL    :: eps

  namelist /optim/ niter, nfunc, fmin, iprint     &
  &               , epsf, epsx, epsg              &
  &               , nupdate, eps

  REAL    :: adobjf
  integer :: nn, m, jmin, jmax, itask, ifail
  integer :: i
  logical :: cold
  integer :: ifunc
  integer :: istat


!======================================================
! preset the optimisation switches
!======================================================
  niter   = 20
  nfunc   = 40
  fmin    = 0.0
  iprint  = 20

  epsf    = 1.E-6
  epsx    = 1.E-6
  epsg    = 1.E-6
  eps     = -1.E-6

  nupdate = 5

!======================================================
! read new optimisation switches from file
!======================================================
  open( ioptns, file   = 'lsopt.par'       &
  &            , access = 'SEQUENTIAL'     &
  &            , form   = 'FORMATTED'  )

  read( ioptns, nml=optim, IOSTAT=istat )

  if (istat .le. 0) then
     if (eps .gt. 0.0) then
        epsf = eps
        epsx = eps
        epsg = eps
     end if
     if (info) print *, ' nextx: options have been read'
  else
     if (info) print *, ' nextx: cannot read options'
  end if

  close( ioptns )

!==============================================================
! start the optimisation
!==============================================================
!HEW-CHG-050109    call instore( n, ifunc, objf, nupdate, jmin, jmax, cold, ifail )
  call instoref( n, ifunc, objf, nupdate, jmin, jmax, cold, ifail )
  print*,'cold=',cold

  if (.not. cold) then
!HEW-CHG-050109  call dostore( n, x, .false., 1 )
     call dostoref( n, x, .false., 1 , 'OPWARMD' )
!HEW-CHG-050109     call dostore( n, adx, .false., 2 )
     call dostoref( n, adx, .false., 2 , 'OPWARMD' )
     if (info) print *, ' nextx: control vector, function, and gradient read'
     if ((jmin .eq. 2) .and. (jmax .eq. 0)) then
        jmin  = 1
        jmax  = 0
!HEW-CHG-040830       call outstor( n, ifunc, objf, nupdate, jmin, jmax )
        call outstorf( n, ifunc, objf, nupdate, jmin, jmax )
        itask = 0
        if (info) print *, ' nextx: start optimization'
     else
        itask = 1
        if (info) print *, ' nextx: restart optimization'
        open (unit=itmp,file='lsopt.tmp',form='unformatted')
        read (unit=itmp) rz1, rz2, rz3, wrz, iwrz
        close(unit=itmp)
     end if

     call lsopt( itask, n, x, objf, adx                  &
     &         , epsx, fmin, epsg, iprint                &
     &         , niter, nfunc, nupdate, rz1, rz2, rz3    &
     &         , wrz, iwrz                               &
     &         , ifail )

     open (unit=itmp,file='lsopt.tmp',form='unformatted')
     write(unit=itmp) rz1, rz2, rz3, wrz, iwrz
     close(unit=itmp)

     if ((ifail .eq. 0) .and. (itask .eq.1)) then
!HEW-CHG-040830  call dostore( n, x, .true., 1 )
        call dostoref( n, x, .true., 1 , 'OPWARMD' )
        if (info) print *, ' nextx: control vector stored'
        if (info) print *, ' nextx: start fgrad to compute function and gradient'
     else
        call postadm( n, x, objf, adx )
        open (unit=itmp,file='result')
!!MS$        write(unit=itmp,*) objf
!!MS$        write(unit=itmp,*) x
        close(unit=itmp)
     end if

  else
     jmin  = 2
     jmax  = 0
!HEW-CHG-040830       call outstor( n, ifunc, objf, nupdate, jmin, jmax )
     call outstorf( n, ifunc, objf, nupdate, jmin, jmax )
!HEW-CHG-040830  call dostore( n, x, .true., 1 )
     call dostoref( n, x, .true., 1 , 'OPWARMD' )
     if (info) print *, ' nextx: first guess of control vector stored'
     if (info) print *, ' nextx: start fgrad to compute function and gradient'
  end if

end subroutine optimum
