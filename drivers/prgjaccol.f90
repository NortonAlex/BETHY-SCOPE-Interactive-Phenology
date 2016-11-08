!**************************************************************
!   file: 	        prgjaccol.f90
!   purpose: 	        driver to compute columns of Jacobian in forward mode
!   creation date:	07 02
!
!   Copyright (C) 2000, 2001
!   FastOpt, Drs. Ralf Giering und Thomas Kaminski GbR
!          , Hamburg, Germany
!
!   Email: info@FastOpt.de
!
!   All rights reserved.
!**************************************************************
program prgjaccol
  implicit none
  integer :: n, m
  integer nbeg, nend, nstep
  namelist /jaccol/ nbeg, nend, nstep
  character*(*) cfile
  parameter ( cfile = 'jaccol.par')
  integer ncol, i

!==============================================================
! get the number of the independent and dependent variables
!==============================================================
  call setfunc( n, m )

!==============================================================
! parameters to define columns
!==============================================================
      nbeg  = 1
      nend  = n
      nstep = 1

      open (UNIT=87, FILE=cfile )
      read (UNIT=87, nml=jaccol, END=100, ERR=100 )
      print '(A)',' parameters for columns of Jacobian read'
 100  close(UNIT=87)

      nbeg = max(1,nbeg)
      nbeg = min(n,nbeg)
      nend = max(1,nend)
      nend = min(n,nend)
      ncol = 0
      do i = nbeg,nend,nstep
         ncol = ncol+1
      enddo

!==============================================================
! wrapper for the actual computation, see subroutines below
! (recursive calls are needed for dynamic allocation)
!==============================================================
      call doit( n, m, ncol, nbeg, nend, nstep )

end program prgjaccol

subroutine doit( n, m, ncol, nbeg, nend, nstep )
      implicit none
      integer n, m, ncol, nbeg, nend, nstep
      integer col(ncol)

      REAL   x(n)
      REAL   y(m)
      REAL   g_x(ncol,n)
      REAL   g_y(ncol,m)

      integer nn
      integer ifail
      logical cold
      integer jmin, jmax, nupdate, ifunc
      real    fc

      integer i,j,icol

!==============================================================
! initialise the model and set control variables
!==============================================================
      call initfunc( n, x )

!==============================================================
! option to read control vector after an optimisation
!==============================================================

  ifail=0
  OPEN (39,file='xf.b',status='old',form='unformatted',iostat=ifail)
  if (ifail.gt.0) then
     print*,'parameter values set by initmod'
  else
     read(39) x
     print*,'parameter values overwritten by those from xf.b'
  endif
  CLOSE(39)

!==============================================================
! save index of computed column
!==============================================================
      icol=0
      do i = nbeg, nend, nstep
         icol=icol+1
         col(icol)=i
      enddo

!==============================================================
! define so-called seed matrix 
! for computing columns of the Jacobian
! first index counts the columns
! second index gives the position of the column in the Jacobian
!==============================================================
      g_x = 0.0
      do icol =1, ncol
         g_x(icol,col(icol)) = 1.0
      enddo

!==============================================================
! run derivative code
!==============================================================
      call func_hes( ncol, n, x, g_x, ncol, m, y, g_y, ncol )

!==============================================================
! postprocessing
!==============================================================
      call POSTJACCOL(  n, x, m, y, g_y, col, ncol )
      call POSTFUNC(  n, x, m, y )

end subroutine doit
