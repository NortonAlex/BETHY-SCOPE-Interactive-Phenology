!**************************************************************
!   file: 	        prgtstjaccol.f90
!   purpose: 	        driver to compute columns of Jacobian in forward mode,
!                       check against finite differences 
!                       and print timing information
!   creation date:	07 02
!                       bug - mixing up column indices - removed in 08 02
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
  namelist /pjaccol/ nbeg, nend, nstep
  character*(*) cfile
  parameter ( cfile = 'pjaccol.par')
  integer i,ncol

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

      open (UNIT=87, FILE=cfile)
 
      read (UNIT=87, NML=pjaccol, END=100, ERR=100 )
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
      call doit( n, m, ncol, nbeg, nend, nstep)

end program prgjaccol


subroutine doit( n, m, ncol, nbeg, nend, nstep )

      implicit none
      integer n, m, ncol, nbeg, nend, nstep

      integer col(ncol)

      REAL   x(n)
      REAL   y(m)
      REAL   g_x(ncol,n)
      REAL   g_y(ncol,m)
      REAL   gfd(ncol,m)

      integer i,j,icol,ifail

!==============================================================
! initialise the model and set control variables
!==============================================================
      call initfunc( n, x )

      ifail=0
      OPEN (39,file='x.b',status='old',form='unformatted',iostat=ifail)
      if (ifail.gt.0) then
         print*,'parameter values set by initfunc'
      else
         read(39) x
         print*,'parameter values overwritten by those from x.b'
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
! second index counts the columns
! first index defines the position of the ith column
!==============================================================
      g_x = 0.0
      do icol =1, ncol
         g_x(icol,col(icol)) = 1.0
      enddo
      g_y = 0.0

!==============================================================
! compute Jacobian forward sensitivity
!==============================================================
      call func_ptl( ncol, n, x, g_x, ncol, m, y, g_y, ncol )
        
!==============================================================
! postprocessing
!==============================================================
      call POSTJACCOL(  n, x, m, y, g_y, col, ncol )
      call POSTFUNC(  n, x, m, y )

end subroutine doit
