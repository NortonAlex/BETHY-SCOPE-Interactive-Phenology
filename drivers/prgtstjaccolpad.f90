!**************************************************************
!   file: 	        prgtstjaccolpad.f90
!   purpose: 	        driver to compute columns of Jacobian in reverse mode,
!                       check against finite differences 
!                       and print timing information
!                       derivative interface: func_pad
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
program prgtstjaccolpad
  implicit none
  integer :: n, m
  REAL        epsmach
  parameter ( epsmach = 1.E-20 )
  REAL        epsbeg
  parameter ( epsbeg = 1.E-4 )
  REAL    eps, acc
  integer nbeg, nend, mbeg, mend, nstep, neps, epsstep
  namelist /check/ eps, acc, nbeg, nend, mbeg, mend, nstep, neps, epsstep
  character*(*) cfile
  parameter ( cfile = 'tst_prog.par')
  integer i,ncol

!==============================================================
! get the number of the independent and dependent variables
!==============================================================
  call setfunc( n, m )
  print*, 'n,m = ',n,m
!==============================================================
! parameters for check
!==============================================================
      eps   = epsbeg
      nbeg  = 1
      nend  = n
      nstep = 1
      acc   = 1.e-3

      open (UNIT=87, FILE=cfile )
      read (UNIT=87, nml=check, END=100, ERR=100 )
      print '(A)',' parameters for gradient check read'
 100  close(UNIT=87)

      nbeg = max(1,nbeg)
      nbeg = min(n,nbeg)
      nend = max(1,nend)
      nend = min(n,nend)
      mbeg = max(1,mbeg)
      mbeg = min(m,mbeg)
      mend = max(1,mend)
      mend = min(m,mend)
      ncol = 0
      do i = nbeg,nend,nstep
         ncol = ncol+1
      enddo

!==============================================================
! wrapper for the actual computation, see subroutines below
! (recursive calls are needed for dynamic allocation)
!==============================================================
      call doit( n, m, ncol, nbeg, nend, mbeg, mend, nstep, acc, eps, epsmach  )

end program prgtstjaccolpad

module mo_timing
  integer isec0, isec, crate, isecmax
  integer dt
  integer iclock0(8), iclock(8)

contains

  subroutine inittiming
    implicit none
    call system_clock( count=isec0 )
    call system_clock( count=isec )
    dt = isec - isec0
    call system_clock( count_rate=crate, count_max=isecmax )
  end subroutine inittiming

  subroutine btiming
    implicit none
    call date_and_time( values=iclock0 )
    call system_clock( count=isec0 )
  end subroutine btiming
  
  subroutine etiming(deltat)
    implicit none
    real deltat
    call system_clock( count=isec )
    call date_and_time( values=iclock )
    if (iclock(3).eq.iclock0(3)) then
       deltat=float((isec-isec0) - dt)/crate
    elseif (iclock(3).gt.iclock0(3)) then
       isec=isec+(iclock(3)-iclock0(3))*isecmax
       deltat=float((isec-isec0) - dt)/crate   
    else
       deltat=0.0
       print*,'mo_timing: timing problem'
       print*,' maybe a new month has begun overnight?'
       print*,' to avoid crash set deltat = ',deltat
       print*,' isec0 = ', isec0
       print*,' isec = ', isec
       print*,' iclock0(3) = ', iclock0(3)
       print*,' iclock(3) = ', iclock(3)
    endif
  end subroutine etiming

end module mo_timing

subroutine doit( n, m, ncol, nbeg, nend, mbeg, mend, nstep, acc, eps, epsmach )
      use mo_timing
      implicit none
      integer n, m, ncol, nbeg, nend, mbeg, mend, nstep
      real eps, acc, epsmach
      integer col(ncol)

      REAL   x(n)
      REAL   adx(m,n)
      REAL   ady(m,m)
      REAL   gfd(ncol,m)

      integer i,j,icol
      REAL    y(m), yh(m), xmemo, rel(ncol,m), deltay(m), absmax

      REAL    tfunc, tgrad
      REAL    deltat
      integer istat
      integer niter, ncfunc, ncgrad

!======================================================
! initialisations for timing
!======================================================
      tfunc = 0.0
      tgrad = 0.0
      niter = 1
      ncfunc = 0
      ncgrad = 0
      call inittiming

!==============================================================
! initialise the model and set control variables
!==============================================================
      call initfunc( n, x )

!==============================================================
! save index of computed column
!==============================================================
      icol=0
      do i = nbeg, nend, nstep
         icol=icol+1
         col(icol)=i
      enddo

!==============================================================
! compare components of gradient against finite differences
! and do timing for admodel and each call of the model
!==============================================================
      print *
      print *, '**************************************************'
      print '(A,E9.3)', ' Check of columns of Jacobian USING eps = ', eps
      print *, '**************************************************'
      print '(2A)', ' col i row j         x(i) delty(j)/eps  dy(j)/dx(i) ',&
           '    relative'

! run model for initial values of control variables
      call btiming
      call func( n, x, m, y )
      call etiming(deltat)
      ncfunc=ncfunc+1
      tfunc = tfunc + deltat

! initialise adjoint variables
      adx = 0.0
      ady = 0.0
      do j =1, m
         ady(j,j) = 1.0
      enddo

! run Jacobian for initial values of control variables
      call btiming
      call func_pad( m, n, x, adx, m, m, y, ady, m )
      call etiming(deltat)

      ncgrad=ncgrad+1
      tgrad = tgrad + deltat

! loop over components of the function 
      do icol =1, ncol
         
! perturb one control variable
         xmemo = x(col(icol))
         x(col(icol)) = x(col(icol)) + eps

! run model for perturbed values of control variables
         call btiming
         call func( n, x, m, yh )
         call etiming(deltat)
         ncfunc=ncfunc+1
         tfunc = tfunc + deltat
         x(col(icol)) = xmemo
        
! loop over dependent variables
         do j = mbeg, mend
! compute finite difference approximation
            deltay(j) = yh(j) - y(j)
            if (abs(deltay(j)) .lt. epsmach) then
               gfd(icol,j) = 0.0
            else
               gfd(icol,j) = deltay(j) / eps
            end if
! compute normalised difference
            rel(icol,j) = abs(gfd(icol,j)-adx(j,col(icol)))
            absmax = max( abs(gfd(icol,j)), abs(adx(j,col(icol))) )
            if (absmax .lt. epsmach) then
               rel(icol,j) = 0.0
            else
               rel(icol,j) = rel(icol,j) / absmax
            end if

! print result of test
            if (rel(icol,j) .lt. acc) then
               print '(2i6,4(1x,e12.6)  )',col(icol),j,x(col(icol)),gfd(icol,j),adx(j,col(icol)),rel(icol,j)
            else
               print '(2i6,4(1x,e12.6),a)',col(icol),j,x(col(icol)),gfd(icol,j),adx(j,col(icol)),rel(icol,j), &
                    ' DIFFERENT'
            endif
         end do
      enddo

      tfunc=tfunc/ncfunc
      tgrad=tgrad/ncgrad
      print *

!===============================================
! print timing 
!===============================================
      print *
      print *, '**************************************************'
      print *, '  TIMING OF columns of Jacobian '
      print *, '**************************************************'
      print '(a,i8)'   , ' based on function calls      : ', ncfunc
      print '(a,i8)'   , ' based on jac columns calls   : ', ncgrad
      print '(a,i8)'   , ' no. of columns               : ', ncol
      print '(a,f12.5)', ' run time function            : ', tfunc
      print '(a,f12.5)', ' run time jac columns         : ', tgrad

      if (tfunc .gt. 0.0) then
         print '(2a,f12.5)', ' rel. run time'&
           , ' (JACCOLS)/FUNC          = ', tgrad/tfunc
      endif
      print *

!===============================================
! repeat check output 
!===============================================
      print *
      print *, '**************************************************'
      print '(A,E9.3)', ' Check of columns of Jacobian USING eps = ', eps
      print *, '**************************************************'
      print '(2A)', ' col i row j         x(i) delty(j)/eps  dy(j)/dx(i) ',&
           '    relative'
      do icol = 1, ncol
         do j = mbeg, mend
            if (rel(icol,j) .lt. acc) then
               print '(2i6,4(1x,e12.6)  )',col(icol),j,x(col(icol)),gfd(icol,j),adx(j,col(icol)),rel(icol,j)
            else
               print '(2i6,4(1x,e12.6),a)',col(icol),j,x(col(icol)),gfd(icol,j),adx(j,col(icol)),rel(icol,j), &
                    ' DIFFERENT'
            endif
         end do
      end do
      print *

!==============================================================
! postprocessing
!==============================================================
      call POSTJACREV(  n, x, m, y, adx)
      call POSTFUNC(  n, x, m, y )

end subroutine doit

