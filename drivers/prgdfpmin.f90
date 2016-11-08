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
      integer, parameter :: mh = 1003
      integer, parameter :: mx = 100
      integer ih, ifunc
      real hx(mx,mh)
      real hf(mh)
      real hg(mx,mh)
      real lx(mx)
      real px(mx)
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

      real         :: pert
      character(8) :: pertfile = 'pert.dat'
      logical      :: ex = .false., irregular=.false.
      integer      :: i
      integer, parameter :: nmaxpert=71

!==============================================================
! initialise the model
! and set the first guess of the control variables
!==============================================================
      call initmod( n, x )
      ifail=0
      OPEN (39,file='x.b',status='old',form='unformatted',iostat=ifail)
      if (ifail.gt.0) then
         print*,'DPFMIN: parameter values set by initmod'
         ! allow to change the first guess, e.g. for identical twin experiments
         inquire (file=pertfile,exist=ex)
         if (ex) then
            open (unit=1,file=pertfile,form='formatted')
            read (1,*) pert
            close (1)
            x(1:nmaxpert) = x(1:nmaxpert) + 1 * pert
            print *, 'adding perturbation of ',pert
         else
            print *, 'no perturbation added'
         endif
      else
         read(39) x
         print*,'DPFMIN: parameter values read from x.b'
      endif
      CLOSE(39)

!==============================================================
! optimize the costfunction using control variables
!==============================================================
      call optimum( n, x, objf, adx )

      do i = 1, n
         if(x(i).lt.sqrt(huge(x(i))).and.x(i).gt.-sqrt(huge(x(i)))) then
         else
            irregular=.true.
         endif
      enddo
!FastOpt      if (sqrt(sum(adx(:)**2)).gt.1.e-5) irregular=.true.
!FastOpt      if (sqrt(sum(adx(:)**2)).gt.1.e-5) print*, 'Opti: irregular'
      if (irregular) return

! save optimim
      print*, 'XXX Save Optimum with value ',objf
      OPEN (39,file='xf.b',status='unknown',form='unformatted')
      WRITE(39) x
      CLOSE(39)
! also in lsopt format
      call outstorf( n, 999, objf, 1, 1, 1)
      call dostoref( n, x, .true., 1 , 'OPWARMD' )
      call dostoref( n, adx, .true., 2 , 'OPWARMD' )

!==============================================================
! start of postprocessor
!==============================================================
      call postadm( n, x, objf, adx )
      call postmod( n, x, objf )

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
      ifunc = 0 ! count of function calls 
      hg = 0. 
      hf = 0.
      hx = 0.
      lx = 0.
      px = 0.
      call dfpmin(x,n,epsg,iter,fc,func,dfunc)
      print*,'DPFMIN: exit after ',iter,' iterations.'
      print*,'XXX function ',fc

!==============================================================
! evaluate final gradient
!==============================================================
      call dfunc( x, adx )
!FastOpt      open (unith,file='h.dat',form='formatted')
      open (unith,file='h.dat',form='formatted',position="append")
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
      use mo_h
      use mo_netcdf
      implicit none
      REAL    x(n)
      integer i, status
      REAL    fc
      logical :: opened 
      character*23 outfile,itername
      integer icount
      save icount
      icount=icount+1

! screen output for control
      print'(a5,2a18)', 'i', 'x(i)'
      do i=1,n
         print'(i5,2e18.6)', i,x(i)
      enddo

!===============================================
! call adjoint model
!===============================================
      fc = 0.
      call model( n, x, fc )
      ifunc = ifunc + 1
      print*, 'XXX Function call ',ifunc,' with value ',fc
      ! ATTENTION
      ! if model returns fc=NaN (e.g. if model is called with
      ! irregular parameters x) we want to return a value that
      ! is a bit less than "huge" to make numerical recipes routine
      ! lnsrch-routine continue
      !
      ! condition ".not.fc.lt.huge(fc)" becomes true for fc=NaN
      !
      if(.not.fc.lt.huge(fc)/1.e50) then
         print*, 'FC too large reduced from ',fc,' to ',huge(fc)/1.e50
         fc = huge(fc)/1.e50
      endif
      if(.not.fc.gt.-huge(fc)/1.e50) then
         fc = huge(fc)/1.e50
         print*, 'PRGDFPMIN: avoiding function overflow'
         call flush(6)
      endif
! save optimim
      print*, 'XXX Save previous point from line search '
      OPEN (39,file='lx.b',status='unknown',form='unformatted')
      WRITE(39) lx(1:n)
      CLOSE(39)
      lx(1:n) = x
      func = fc
! save current control vectors
      write(itername,'(i18)') icount
      outfile='output/eval-'//trim(adjustl(itername))//'.b'
      print*, 'outfile = ',outfile
      call flush(6)
      OPEN (39,file=outfile,status='unknown',form='unformatted')
      WRITE(39) x
      CLOSE(39)
      outfile='output/eval-'//trim(adjustl(itername))//'.dat'
      OPEN (39,file=outfile,status='unknown',form='formatted')
      WRITE(39,'(10(1x,e12.6))') x
      CLOSE(39)
! a hack to close units which are still open
      print*, 'Scanning units'
      do i = 1, 1000
         STATUS = NF_CLOSE(i)
         IF (STATUS .NE. NF_NOERR) then
!            print*, 'Scanning units: Unit ',i,' was open'
         else
            print*, 'Scanning units: Unit ',i,' was closed'
         endif
      enddo
    end function func

!===============================================
! this subroutine evaluates the gradient only
! use taf -pure -reverse
! to generate subroutine admodel
!===============================================
    subroutine dfunc( x, adx )
      use mo_n
      use mo_h
      implicit none
      REAL    x(n), adx(n)
      REAL    fc, adfc
      integer i
      logical :: overflow, underflow
      character*23 outfile,itername

!===============================================
! call adjoint model
!===============================================
      adx  = 0.0
      adfc = 1.0
      call model_ad( n, x, adx, fc, adfc )
      do i = 1,n
         overflow = .false.
         underflow = .false.
         overflow = (.not.adx(i).lt.huge(adx(i)))
         underflow = (.not.adx(i).gt.-huge(adx(i)))
         if (underflow) then
            print*, 'PRGDFPMIN: avoiding derivative underflow for component ',i
            print *, 'underflow ', adx(i)
            adx(i) = -huge(adx(i))/1.e50
            print *, 'underflow corrected', adx(i)
            call flush(6)
         endif
         if (overflow) then
            print*, 'PRGDFPMIN: avoiding derivative overflow for component ',i
            print *, 'overflow ', adx(i)
            adx(i) = huge(adx(i))/1.e50
            print *, 'overflow corrected', adx(i)
            call flush(6)
         endif
         call flush(6)
      enddo
! save state of iteration
      ih = ih + 1
      do i = 1,n
         hx(i,ih)=x(i)
         hg(i,ih)=adx(i)
      enddo
      hf(ih)=fc
      write(*,*) 'Opti: ',ih,hf(ih),sqrt(sum(hg(:,ih)**2))
! screen output for control and gradient
      print'(a5,2a18)', 'i', 'x(i)','grad x(i)'
      do i=1,n
         print'(i5,2e18.6)', i,x(i),adx(i)
      enddo
! save current x
      OPEN (39,file='x.b',status='unknown',form='unformatted')
      WRITE(39) x
      CLOSE(39)
! save previous control vector
      print*, 'XXX Save previous control vector'
      print*, ' px = ', px(1:3)
      OPEN (39,file='px.b',status='unknown',form='unformatted')
      WRITE(39) px(1:n)
      CLOSE(39)
      px(1:n) = x
! save current control vector and gradient
      write(itername,'(i18)') ih 
      outfile='output/opti-'//trim(adjustl(itername))//'.b'
      OPEN (39,file=outfile,status='unknown',form='unformatted')
      WRITE(39) x
      CLOSE(39)
      outfile='output/opti-'//trim(adjustl(itername))//'.dat'
      OPEN (39,file=outfile,status='unknown',form='formatted')
      WRITE(39,'(10(1x,e12.6))') x
      CLOSE(39)
      outfile='output/grad-'//trim(adjustl(itername))//'.b'
      OPEN (39,file=outfile,status='unknown',form='unformatted')
      WRITE(39) adx
      CLOSE(39)
      outfile='output/grad-'//trim(adjustl(itername))//'.dat'
      OPEN (39,file=outfile,status='unknown',form='formatted')
      WRITE(39,'(10(1x,e12.6))') adx
      CLOSE(39)
    end subroutine dfunc
