!**************************************************************
!   file: 	        prgtstvtlm.f90
!   purpose: 	        driver to eval. Jacobian columns in forward mode,
!                       check against FD and do timing 
!   creation date:	10/07
!
!   Copyright (C) 2000-2007
!   FastOpt, GmbH, Hamburg, Germany
!   http://FastOpt.com
!   All rights reserved.
!
!**************************************************************
program prgtstvtlm
  implicit none
  integer :: n, m                 ! dimensions of function
  integer :: nbeg = 1             ! first component to test
  integer :: nend = 1             ! last component to test
  integer :: nstep = 1            ! step size through components
  integer :: neps = 1             ! # of FD intervals
  real    :: eps = 1.e-6          ! first FD interval
  integer :: epsstep = 10         ! factor between two FD intervals
  real    :: epsmach = 1.e-20     ! machine precision
  real    :: acc = 1.e-3          ! required relative accuracy
  integer :: niter = 1            ! # of derivative calls (for timing)
  real    :: tlspeed = 5.         ! minimum relative speed of TLM
  character*(7) :: cfile = 'tst.par'
  integer :: i, ncol
  namelist /check/ nbeg, nend, nstep, neps, eps, epsstep, epsmach, acc, niter, tlspeed

  ! get the dimensions of the function
  call numbmod( n)
  m = 1
  nend  = n
  ! read namelist and overwrite defaults
  open (UNIT=87, FILE=cfile )
  read (UNIT=87, nml=check, END=100, ERR=100 )
  print '(A)',' parameters for read from file ', cfile
100 close(UNIT=87)
  ! set Jacobian dimensions
  nbeg = max(1,nbeg)
  nbeg = min(n,nbeg)
  nend = max(1,nend)
  nend = min(n,nend)
  ncol = 0
  do i = nbeg,nend,nstep
     ncol = ncol+1
  enddo
  ! trick for dynamic allocaton
  call doit( n, m, ncol, nbeg, nend, nstep, neps, eps, epsstep, epsmach, acc, niter, tlspeed )
end program prgtstvtlm

!**************************************************************
module mo_timing
  ! timing module
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
    real(kind=8) deltat
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

!**************************************************************
subroutine doit( n, m, ncol, nbeg, nend, nstep, neps, epsbeg, epsstep, epsmach, acc, niter, tlspeed )
  use mo_timing
  implicit none
  integer :: n, m, ncol, nbeg, nend, nstep, neps, epsstep, niter
  real    :: eps, epsbeg, epsmach, acc, tlspeed

  integer :: col(ncol)
  real    :: x(n)
  real    :: d_x(ncol,n)
  real    :: d_y(ncol,m), d_fc(ncol)
  real    :: gfd(ncol,m)
  real    :: y1(m), y2(m), ytlm(m), xmemo, rel(ncol,m), deltay(m), absmax, fc, fctlm

  logical, parameter :: plot = .true.

  real ::    tfunc, tgrad, deltat
  integer :: ncfunc, ncgrad, i, j, icol, iter, punit(ncol)
  character*30 :: pfile(ncol)

  ! initialisations for timing
  tfunc = 0.0
  tgrad = 0.0
  ncfunc = 0
  ncgrad = 0
  call inittiming

 ! initialise the model and set control variables
  call initmod( n, x )

  ! save index of computed column
  icol=0
  do i = nbeg, nend, nstep
     icol=icol+1
     col(icol)=i
  enddo

  ! run model for initial values of control variables
  do iter=1 , niter 
     fc=0.
     call btiming
     call model( n, x, fc)
     y1(1) = fc
     call etiming(deltat)
     ncfunc=ncfunc+1
     tfunc = tfunc + deltat
  enddo

  ! initialise Jacobian sub matrix of interest to identity
  ! second index counts the columns
  ! first index defines the position of the ith column
  d_x = 0.0
  do icol =1, ncol
     d_x(icol,col(icol)) = 1.0
  enddo

  ! run Jacobian for initial values of control variables
  do iter=1 , niter 
     fctlm=0.
     d_fc=0.
     call btiming
     call model_vtl( ncol, n, x, d_x, ncol, fctlm, d_fc, ncol )
     d_y(:,1) = d_fc(:)
     ytlm(1) = fctlm 
     call etiming(deltat)
     ncgrad=ncgrad+1
     tgrad = tgrad + deltat
  enddo
  print *
  print *, '**************************************************'
  print *, ' Comparing function values function and TLM codes'
  print *, '**************************************************'
  print '(a6,4(1x,a17))','row','y(row) func','y(row) tlm','deltay'
  do j=1,m
     print '(i6,4(1x,e17.10)  )',j,y1(j),ytlm(j),ytlm(j)-y1(j)
  enddo

  if (plot) then
  ! open plot files
     do icol = 1, ncol
        punit(icol) = 900 + icol
        write (pfile(icol),*) col(icol)
        pfile(icol) = 'fd-tlm-component-'//trim(adjustl(pfile(icol)))//'.set'
        open (unit = punit(icol), file = pfile(icol))
        write(punit(icol),'(a1,a65,i10)') '#','|AD - finite differences| per function component for parameter ',icol 
        write(punit(icol),'(a1,a17,10(1x,i17))') '#','eps',(j,j=1,m)
     enddo
  endif

  print *
  print *, '**************************************************'
  print *, ' Finite Differences vs Jacobian columns'
  print *, '**************************************************'
  print '(2a6,7(1x,a17))', 'col i', 'row j','eps','x(i)','y2(j)','y1(j)', 'delty(j)/eps','dy(j)/dx(i) ','rel diff'

  ! loop over independents
  do icol =1, ncol
     xmemo = x(col(icol))
     do iter =1, neps
        eps = epsbeg * epsstep**(iter-1)
        ! perturb one control variable
        x(col(icol)) = xmemo + eps
        ! run model for perturbed values of control variables
        fc=0.
        call btiming
        call model( n, x, fc)
        y2(1) = fc
        call etiming(deltat)
        ncfunc=ncfunc+1
        tfunc = tfunc + deltat
        ! loop over dependent variables
        do j = 1, m
           ! compute finite difference approximation
           deltay(j) = y2(j) - y1(j)
           if (abs(deltay(j)) .lt. epsmach) then
              gfd(icol,j) = 0.0
           else
              gfd(icol,j) = deltay(j) / eps
           end if
           ! compute normalised difference
           rel(icol,j) = abs(gfd(icol,j)-d_y(icol,j))
           absmax = max( abs(d_y(icol,j)), abs(gfd(icol,j))) 
           if (absmax .lt. epsmach) then
              rel(icol,j) = 0.0
           else
              rel(icol,j) = rel(icol,j) / absmax
           end if
           ! print result of test
           if (rel(icol,j) .lt. acc) then
              print '(2i6,7(1x,e17.10)  )',col(icol),j,eps,xmemo,y2(j),y1(j),gfd(icol,j),d_y(icol,j),rel(icol,j)
           else
              print '(2i6,7(1x,e17.10),a)',col(icol),j,eps,xmemo,y2(j),y1(j),gfd(icol,j),d_y(icol,j),rel(icol,j), &
                   ' DIFFERENT'
           endif
        enddo
        if (plot) then
           write(punit(icol),'(6(1x,e17.10))') eps,(max(rel(icol,j),epsmach),j=1,m)
        endif
     enddo
     x(col(icol)) = xmemo
  enddo

  ! timing output
  tfunc=tfunc/ncfunc
  tgrad=tgrad/ncgrad
  print *
  print *, '**************************************************'
  print *, '  Run time '
  print *, '**************************************************'
  print '(a,i8)'   , '   based on function calls      : ', ncfunc
  print '(a,i8)'   , '   based on jac columns calls   : ', ncgrad
  print '(a,i8)'   , '   no. of columns               : ', ncol
  print '(a,f14.5)', '   run time function            : ', tfunc
  print '(a,f14.5)', '   run time jac columns         : ', tgrad

  if (tfunc .gt. 0.0) then
     if (tgrad/tfunc.lt.TLSPEED) then
        print '(a,f14.5)',   '   rel. run time jac/func       : ', tgrad/tfunc
     else
        print '(a,f14.5,a)', '   rel. run time jac/func       : ', tgrad/tfunc, ' TOO SLOW'
     endif
  else
     print*, '  runtime function too short, increase niter'
  endif
  print *

  if (plot) then
  ! open plot files
     do icol = 1, ncol
        close (punit(icol))
     enddo
  endif

end subroutine doit
