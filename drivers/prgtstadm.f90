!**************************************************************
!   file: 	        prgtstadm.f90
!   purpose: 	        driver to check gradient from scalar 
!                       valued adjoint against finite 
!                       differences and print timing information
!   creation date:	06 02
!
!   Copyright (C) 2000, 2001
!   FastOpt, Drs. Ralf Giering und Thomas Kaminski GbR
!          , Hamburg, Germany
!
!   Email: info@FastOpt.de
!
!   All rights reserved.
!**************************************************************
program prgtstadm
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

end program prgtstadm

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

      subroutine doit( n )
      use mo_timing
      implicit none
      integer n

      REAL   adfc
      REAL   x(n)
      REAL   adx(n)
      REAL   gfd(n)

      REAL        epsmach
      parameter ( epsmach = 1.E-20 )
      REAL        epsbeg
      parameter ( epsbeg = 1.E-4 )
      integer i
      REAL    fc, fch, xmemo, rel(n), deltafc, absmax
      REAL    eps, acc
      integer nbeg, nend, nstep
      namelist /check/ eps, acc, nbeg, nend, nstep
      character*(*) cfile
      parameter ( cfile = 'tst.par')

      REAL    tfunc, tgrad, tforw, deltat
      integer istat
      integer niter, ncfunc

      INTEGER nn
      INTEGER ifail
      LOGICAL cold
      INTEGER jmin, jmax, nupdate, ifunc


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

!======================================================
! initialisations for timing
!======================================================
      tfunc = 0.0
      tgrad = 0.0
      niter = 1
      call inittiming

!==============================================================
! initialise the model and set control variables
!==============================================================
      call initmod( n, x )

!==============================================================
! option to read control vector after an optimisation
!==============================================================
!HEW-CHG-040930: instore (in lsopt lib) -> instoref (local)
      CALL instoref( nn, ifunc, fc, nupdate, jmin, jmax, cold, ifail )
      
      IF (cold) THEN
         PRINT *, 'control variables are set by initfunc'
      ELSE
         IF (n .NE. nn) THEN
            STOP 'inconsitence number of control variables'
         ENDIF
!HEW-CHG-040830
!  call dostore( n, x, .false., 1 )
         call dostoref( n, x, .false., 1 , 'OPWARMD' )
         PRINT *, 'control variables read from OPWARMD'
      ENDIF
      
  ifail=0
  OPEN (39,file='x.b',status='old',form='unformatted',iostat=ifail)
  if (ifail.gt.0) then
  else
     read(39) x
     print*,'cost: parameter values overwritten by those from x.b'
  endif
  CLOSE(39)

!==============================================================
! compare components of gradient against finite differences
! and do timing for admodel and each call of the model
!==============================================================
      print *
      print *, '**************************************************'
      print '(A,E9.3)', ' CHECK OF ADJOINT USING eps = ', eps
      print *, '**************************************************'
      print '(2A)', '     I         x(i)  delta f/eps', &
           '     grad f   RELATIVE ERR'

! run model for initial values of control variables
      call btiming
      call model( n, x, fc )
      call etiming(deltat)
      ncfunc = 1
      tfunc = tfunc + deltat

! initialise adjoint
      do i = 1, n
        adx(i) = 0.0
      end do
      adfc = 1.0

! run adjoint for initial values of control variables
      call btiming
      call model_ad( n, x, adx, fch, adfc )
      call etiming(deltat)
      tgrad = tgrad + deltat

      print *,'checking function value of adjoint'
      print '(4(1x,a19))','y(row) fu','y(row) ad','deltay'
      print '(4(1x,e19.12)  )',fc,fch,fch-fc

! loop over components of the gradient 
      do i = nbeg, nend, nstep

! perturb one control variable
         xmemo = x(i)
         x(i) = x(i) + eps

! run model for perturbed values of control variables
         call btiming
         call model( n, x, fch )
         call etiming(deltat)
         ncfunc=ncfunc+1
         tfunc = tfunc + deltat
         x(i) = xmemo
         
! compute finite difference approximation
         deltafc = fch - fc
         if (abs(deltafc) .lt. epsmach) then
           gfd(i) = 0.0
         else
           gfd(i) = deltafc / eps
         end if

! compute normalised difference
         rel(i) = abs(gfd(i)-adx(i))
         absmax = max( abs(gfd(i)), abs(adx(i)) )
         if (absmax .lt. epsmach) then
            rel(i) = 0.
         else
            rel(i) = rel(i) / absmax
         end if

! print result of test
         if (rel(i) .lt. acc) then
            print '(i6,4(1x,e12.6)  )',i,x(i),gfd(i),adx(i),rel(i)
         else
            print '(i6,4(1x,e12.6),a)',i,x(i),gfd(i),adx(i),rel(i), &
                ' DIFFERENT'
         endif
      end do

      tfunc=tfunc/ncfunc
      print *, '**************************************************'
      print *

!===============================================
! print timing 
!===============================================
      print *
      print *, '**************************************************'
      print *, '  TIMING OF ADJOINT '
      print *, '**************************************************'
      print '(a,i8)'   , ' based on function calls      : ', ncfunc
      print '(a,f12.3)', ' run time function            : ', tfunc
      print '(a,f12.3)', ' run time function + adjoint  : ', tgrad

      if (tfunc .gt. 0.) then
         print '(2a,f12.3)', ' rel. run time reverse mode'&
           , ' (FUNC+GRAD)/FUNC : ', tgrad/tfunc
      endif
      print *, '**************************************************'
      print *

!===============================================
! repeat check output 
!===============================================
      print *
      print *, '**************************************************'
      print '(A,E9.3)', ' CHECK OF ADJOINT USING eps = ', eps
      print *, '**************************************************'
      print '(2A)', '     I         x(i)  delta f/eps', &
           '     grad f   RELATIVE ERR'
      do i = nbeg, nend, nstep
         if (rel(i) .lt. acc) then
            print '(i6,4(1x,e12.6)  )',i,x(i),gfd(i),adx(i),rel(i)
         else
            print '(i6,4(1x,e12.6),a)',i,x(i),gfd(i),adx(i),rel(i), &
                ' DIFFERENT'
         endif
      end do
      print *, '**************************************************'

! save gradient
     open (39,file='g.b',status='unknown',form='unformatted')
     write(39) adx
     close(39)

! post processing 
     call postmod( n, x, fc )
end subroutine doit
