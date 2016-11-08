!======================================================================
!   file:               prgdfpmin.f90
!   purpose:            driver for optimisation using dpfmin
!   creation date:      09 05
!
!   prepared by:
!   FastOpt, Drs. Ralf Giering und Thomas Kaminski GbR
!          , Hamburg, Germany
!
!   Email: info@FastOpt.de
!======================================================================

      module mo_n
      implicit none
      integer n
      end module mo_n

      module mo_h
      implicit none
      integer, parameter :: mh = 201
      integer, parameter :: mx = 100
      integer ih
      real hx(mx,mh)
      real hf(mh)
      real hg(mx,mh)
      end module mo_h

      program main
!*************************************************************
      use mo_n
      implicit none

!==============================================================
! get the number of control variables
!==============================================================
      call numbmod( n )

!-----------------------------------------
! call the subroutine
!-----------------------------------------
      call opti

      end

      subroutine opti
      use mo_n
      implicit none

      REAL    objf

      REAL   x(n)
      REAL   adx(n)
      integer ifail

!==============================================================
! initialisize the model
! and set the first guess of the control variables
!==============================================================
      call initmod( n, x )
      ifail=0
      OPEN (39,file='x.b',status='old',form='unformatted',iostat=ifail)
      if (ifail.gt.0) then
         print*,'DPFMIN: parameter values set by initmod'
      else
         read(39) x
         print*,'DPFMIN: parameter values read from x.b'
      endif
      CLOSE(39)

!==============================================================
! optimize the costfunction using control variables
!==============================================================
      call optimum( n, x, objf, adx )

!==============================================================
! start of postprocessor
!==============================================================
      call postmod( n, x, objf )

!      call postadm( n, x, objf, adx )

! save optimim
      OPEN (39,file='xf.b',status='unknown',form='unformatted')
      WRITE(39) x
      CLOSE(39)

      end

      subroutine optimum( n, x, fc, adx )
      use mo_h
      implicit none

!==============================================================
! declare call parameters
!==============================================================
      integer n
      REAL    x(n)
      REAL    fc
      REAL    adx(n)

!==============================================================
! declare local variables
!==============================================================
      integer    ioptns, iter, i, unith
      parameter( ioptns = 7 )
      parameter( unith = 8 )
      REAL epsg
      external  func,dfunc

      namelist/ OPTIM/ epsg

      epsg  = 1.E-2

!======================================================
! read new optimisation switches from file
!======================================================
      open( ioptns, FILE   = 'dfpmin.par'&
                 , ACCESS = 'SEQUENTIAL'&
                 , FORM   = 'FORMATTED'  )
      read( ioptns, OPTIM, end=9000, ERR=8000 )
      print *, ' OPTIMI : options have been read'
 8000 continue
 9000 continue
      close( ioptns )

!==============================================================
! perform optimisation
!==============================================================
      ih = 0 ! our own count of iterations
      hg = 0. 
      hf = 0.
      hx = 0.
      call dfpmin(x,n,epsg,iter,fc,func,dfunc)
      print*,'DPFMIN: exit after ',iter,' iterations.'

!==============================================================
! evaluate final gradient
!==============================================================
      call dfunc( x, adx )
      open (unith,file='h.dat',form='formatted')
      do i=1,ih
         write(unith,*)i,hf(i),sqrt(sum(hg(:,i)**2))
      enddo
      close (unith)
      end

!===============================================
! this subroutine evaluates the function
!===============================================
      REAL function func( x )
      use mo_n
      implicit none
      REAL    x(n)
      integer i
      REAL    fc

!===============================================
! call adjoint model
!===============================================
      fc = 0.
      call model( n, x, fc )
      ! ATTENTION
      ! if model returns fc=NaN (e.g. if model is called with
      ! irregular parameters x) we want to return a value that
      ! is a bit less than "huge" to make numerical recipes routine
      ! lnsrch-routine continue
      !
      ! condition ".not.fc.lt.huge(fc)" becomes true for fc=NaN
      !
      if(.not.fc.lt.huge(fc)/1.e50) then
         print* ,'DFP: fc, huge(fc)/1.e50 ',fc, huge(fc)/1.e50
         fc = huge(fc)/1.e50
         print*, 'PRGDFPMIN: avoiding function overflow'
         call flush(6)
      endif
      if(.not.fc.gt.-huge(fc)/1.e50) then
         print* ,'DFP: fc, huge(fc)/1.e50 ',fc, huge(fc)/1.e50
         fc = huge(fc)/1.e50
         print*, 'PRGDFPMIN: avoiding function overflow'
         call flush(6)
      endif
      func = fc
      end

!===============================================
! this subroutine evaluates the gradient only
! use taf -pure -reverse
! to generate subroutine admodel
!===============================================
      subroutine dfunc( x, dfdx )
      use mo_n
      use mo_h
      implicit none
      REAL    x(n), x_tl(n), dfdx(n)
      REAL    fc, fc_tl
      integer i

!===============================================
! call tangent
!===============================================
      ih = ih + 1
      do i = 1,n
         x_tl = 0.
         x_tl(i) = 1.
         print*,'DPFMIN: evaluating gradient component i = ', i
         call model_tl ( n, x, x_tl, fc, fc_tl )
         hx(i,ih)=x(i)
         hg(i,ih)=fc_tl
         dfdx(i)=fc_tl
      enddo
      hf(ih)=fc
      print*,'i = ', ih,' fc = ',fc
      print*,'i = ', ih,' x = ',hx(1:n,ih)
      print*,'i = ', ih,' g = ',hg(1:n,ih)
! save current x
      OPEN (39,file='x.b',status='unknown',form='unformatted')
      WRITE(39) x
      CLOSE(39)
!      print*,'DPFMIN: current parameter values written to x.b'
      end
