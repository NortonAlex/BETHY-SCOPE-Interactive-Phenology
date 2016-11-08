!**************************************************************
!   file: 	        prgtsttlm.f90
!   purpose: 	        driver to check gradient of scalar valued function
!                       from multiple runs of product tangent linear * vector 
!                       against finite differences 
!                       and print timing information
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
program prgtsttlm
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

end program prgtsttlm

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
      REAL   g_x(n)
      REAL   g_fc(n)
      REAL   gfd(n)

      REAL        epsmach
      parameter ( epsmach = 1.E-20 )
      REAL        epsbeg
      parameter ( epsbeg = 1.E-4 )
      integer i
      REAL    fc, fc2, fch, xmemo, rel(n), deltafc, absmax
      REAL    eps, acc
      integer nbeg, nend, nstep
      namelist /check/ eps, acc, nbeg, nend, nstep
      character*(*) cfile
      parameter ( cfile = 'tst.par')

      REAL    tfunc, tgrad, tforw, deltat
      integer istat
      integer niter, ncfunc, ncgrad

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
      ncgrad = 0
      call inittiming

!==============================================================
! initialise the model and set control variables
!==============================================================
      call initmod( n, x )

!==============================================================
! option to read control vector after an optimisation
!==============================================================
!HEW-CHG-040930: instore (from library) -> instoref (local)
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
      
!==============================================================
! compare components of gradient against finite differences
! and do timing for admodel and each call of the model
!==============================================================
      print *
      print *, '**************************************************'
      print '(A,E9.3)', ' CHECK OF TLM USING eps = ', eps
      print *, '**************************************************'
      print '(2A)', '     I         x(i)  delta f/eps', &
           '     grad f   RELATIVE ERR'

! run model for initial values of control variables
      call btiming
      call model( n, x, fc )
      call etiming(deltat)
      ncfunc = 1
      tfunc = tfunc + deltat

! loop over components of the gradient 
      do i = nbeg, nend, nstep

! initialise tangent linear: perturb in i-th direction
         g_x = 0.0
         g_x(i) = 1.0

! run TLM for initial values of control variables
         call btiming
         call model_tl( n, x, g_x, fc2, g_fc(i) )
         call etiming(deltat)
         ncgrad=ncgrad+1
         tgrad = tgrad + deltat

         if(i.eq.nbeg) then
            print *,'checking function value of tangent'
            print '(4(1x,a19))','y(row) fu','y(row) ad','deltay'
            print '(4(1x,e19.12)  )',fc,fc2,fc2-fc
         endif

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
         rel(i) = abs(gfd(i)-g_fc(i))
         absmax = max( abs(gfd(i)), abs(g_fc(i)) )
         if (absmax .lt. epsmach) then
            rel(i) = 0.0
         else
            rel(i) = rel(i) / absmax
         end if

! print result of test
         if (rel(i) .lt. acc) then
            print '(i6,4(1x,e12.6)  )',i,x(i),gfd(i),g_fc(i),rel(i)
         else
            print '(i6,4(1x,e12.6),a)',i,x(i),gfd(i),g_fc(i),rel(i), &
                ' DIFFERENT'
         endif
      end do

      

      tfunc=tfunc/ncfunc
      tgrad=tgrad/ncgrad
      print *, '**************************************************'
      print *

!===============================================
! print timing 
!===============================================
      print *
      print *, '**************************************************'
      print *, '  TIMING OF TLM '
      print *, '**************************************************'
      print '(a,i8)'   , ' based on function calls      : ', ncfunc
      print '(a,i8)'   , ' based on tangent linear calls: ', ncgrad
      print '(a,f12.3)', ' run time function            : ', tfunc
      print '(a,f12.3)', ' run time function + tangent  : ', tgrad

      if (tfunc .gt. 0.0) then
         print '(2a,f12.3)', ' rel. run time forw mode'&
           , ' (FUNC+GRAD)/FUNC          : ', tgrad/tfunc
      endif
      print *, '**************************************************'
      print *

!===============================================
! repeat check output 
!===============================================
      print *
      print *, '**************************************************'
      print '(A,E9.3)', ' CHECK OF TLM USING eps = ', eps
      print *, '**************************************************'
      print '(2A)', '     I         x(i)  delta f/eps', &
           '     grad f   RELATIVE ERR'
      do i = nbeg, nend, nstep
         if (rel(i) .lt. acc) then
            print '(i6,4(1x,e12.6)  )',i,x(i),gfd(i),g_fc(i),rel(i)
         else
            print '(i6,4(1x,e12.6),a)',i,x(i),gfd(i),g_fc(i),rel(i), &
                ' DIFFERENT'
         endif
      end do
      print *, '**************************************************'
end subroutine doit
