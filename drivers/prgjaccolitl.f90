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
!HEW-DEL instore -> instoref    
      call instoref( nn, ifunc, fc, nupdate, jmin, jmax, cold, ifail )
      if (cold) then
         print *, 'control variables are set by initfunc'
      else
!         if (n .ne. nn) then
!            stop 'inconsitent number of control variables'
!         endif
!HEW-CHG-040830
!         call dostore( n, x, .false., 1 )
         call dostoref( nn, x, .false., 1 , 'OPWARMD' )
         print *, 'control variables read from OPWARMD'
      endif

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
      call func_itl( ncol, n, x, g_x, ncol, m, y, g_y, ncol )

!==============================================================
! postprocessing
!==============================================================
      call POSTJACCOL(  n, x, m, y, g_y, col, ncol )
      call POSTFUNC(  n, x, m, y )

end subroutine doit
