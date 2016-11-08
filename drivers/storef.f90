!**************************************************************
!   file:               storef.f90
!   purpose:            subroutines to store and restore parameters
!   creation date:      04 10
!   Copyright (C) 2000, 2001
!   FastOpt, Drs. Ralf Giering und Thomas Kaminski GbR
!          , Hamburg, Germany
!   modified by Heinrih Widmann, 2004
!
!   Email: info@FastOpt.de
!
!   All rights reserved.
!**************************************************************
subroutine dostoref( n, x, store, j, fname )

  !INPUT
  ! n      number of records
  ! x      vector to store
  ! store  flag for write (true) or read (flase)
  ! fname  file name

  implicit none
  
  integer n, j
  real   x(n)
  logical store
  integer istat
  character*(*) fname
  
  integer    itape     , ntape
  parameter( itape = 91, ntape = 92 )
  
  open( unit   = ntape &
       , file   = TRIM(fname) &
       , status = 'unknown' &
       , form   = 'unformatted' &
       , access = 'direct' &
       , recl   = n * 8 &
       , iostat = istat &
       )
  
  if (istat .gt. 0) then
     print *, 'dostoref: error during open, iostat=', istat
     stop 'dostoref'
  else
     if (store) then
        write( ntape, rec=j, iostat=istat ) x
     else
        read(  ntape, rec=j, iostat=istat ) x
     endif
     if (istat .gt. 0) then
        print *, 'dostoref: error during I/O, iostat=', istat
        stop 'dostoref'
     endif
  endif
  
  close( unit = ntape )
  
end subroutine dostoref

!HEW #include "define.h"

      subroutine instoref( n, ifunc, fc, m, jmin, jmax, cold, ifail )

      implicit none

!
!     arguments
!
      integer n, ifunc, m, jmin, jmax, ifail
!HEW      DREAL   fc
      REAL   fc
      logical cold

!HEW-CHG-040930
!HEW-DEL-040930#include "param.h"
!HEW-ADD-040930
      integer    itape     , ntape
      parameter( itape = 91, ntape = 92 )

      ifail = 0
      cold  = .true.
      ifunc = 0

! one record = two arrays of real*(isize)

!!$      open( unit   = itape
!!$     $    , file   = 'OPWARMI'
!!$     $    , status = 'unknown'
!!$     $    , form   = 'formatted'
!!$     $    , access = 'sequential'
!!$     $    )
      open(unit=itape,file='OPWARMI',status = 'unknown',form='formatted',access = 'sequential')

      read( itape, *, end=800, err=900 ) n, ifunc, fc, m, jmin, jmax

      cold  = .false.

 800  continue
      close(itape)
      return

 900  continue
      close(itape)
      ifail = 1

      end

      subroutine outstorf( n, ifunc, fc, m, jmin, jmax )
      implicit none
!
!     arguments
!
      integer n, ifunc, m, jmin, jmax
!HEW      DREAL   fc
      REAL   fc
!HEW-CHG-040930
!HEW-DEL-040930#include "param.h"
!HEW-ADD-040930
      integer    itape     , ntape
      parameter( itape = 91, ntape = 92 )

      open(unit=itape,file='OPWARMI',access='sequential',form='formatted')
      rewind (itape)

      write( itape, * ) n, ifunc, fc, m, jmin, jmax

      close( itape )

      end
