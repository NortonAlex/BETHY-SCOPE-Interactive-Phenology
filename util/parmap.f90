!**************************************************************
!   file: 	        testmapping.f90
!   purpose: 	        test the mapping of the prameters
!   creation date:	03/08/07
!
!   Copyright (C) 2003
!   MPI-met , CCDAS consortium, Heinrich Widmann
!          , Hamburg, Germany
!
!   All rights reserved.
!**************************************************************
PROGRAM parmap
!*************************************************************
  
  USE mo_namelist


  IMPLICIT NONE

!=========================================
! declaration
!=========================================
  INTEGER i,nvar

  WRITE(*,*) '-----------------------------------------------------------'
  WRITE(*,*) 'Check parameter mapping between'
  WRITE(*,*) '  physical parameters p and'
  WRITE(*,*) '  natural parameters x (normal Gauss distrib. with sigma 1.'
  WRITE(*,*) 'by disturbations of x0 by x = x0 +/- sigma'
  WRITE(*,*) '-----------------------------------------------------------'

!==============================================================
! .. read namelist variables
!==============================================================
  CALL get_namelist


!==============================================================
! get the number of control variables
!==============================================================
  nvar = count_params()

  CALL doit(nvar)

END PROGRAM parmap

SUBROUTINE doit( n )

  USE mo_namelist
  USE mo_mapping 
  USE costf
  USE mo_vegetation 
  USE mo_beta
  USE mo_ctrl
  USE mo_climate
  USE mo_io

  IMPLICIT NONE

  INTEGER n,i
  REAL    fc,map_error
  REAL, PARAMETER :: eps = 1e-10
  REAL   x(n),xcopy(n),xcheck(n),totsig(n)
!  CHARACTER , DIMENSION(58,10), POINTER :: parmameptr
  CHARACTER (LEN= *), PARAMETER :: chtitle = '(A5, 2X, A16 , 2X, A15, 2X, A17)'
  CHARACTER (LEN= *), PARAMETER :: results = '(A12, 2X, F15.6, 2X, F15.6, 2X, F15.6)'
!  CHARACTER (*,LEN= *), PARAMETER :: parnames = (/'vm'/)
!  CHARACTER parnames(21)*20
!  REAL par(21)

!  parnames(0)='vm' ; par(0)=vm(1)

!  print*,'parnames(0)=',parnames(0)

! .. load options to common block, allocate memory
  CALL initopt

!==============================================================
! read and initialize params
!==============================================================
  WRITE(*,*) '  - Read first guess p0 from parmeter file ',param_file
  CALL init_param ( n, param_file)

  ! redo 
!  DO i=1,n
!      IF (xsf(i)<0 .and. pflag(i)==0) THEN
!        p0(i)= EXP(px0(i))
!        p0sig(i)= EXP(p0s(i))
!     ENDIF
!  ENDDO

! .. initialize model
  CALL initialize (mpot, grid_file)

! mapping from control parameters 
! to bethy parameters according PFT specification
  CALL init_mapping(mapping_file)           

!==============================================================
! initialisize the model
! and set the first guess of the control variables
!==============================================================
  CALL initmod( n, x )

  tspd=1440/dtime

!------------------------------------------------------------------
! allocate arrays
!------------------------------------------------------------------
  CALL dummy

  CALL model_allocate


  xcopy = x
  xcheck = x
  totsig = 0.

! x -> p
  CALL mapping(x,xsf,px0,p0sl,p0su,pflag,lb,ub,a,b,n)

  WRITE(*,*)
  WRITE(*,*) '0. Results after p0 -> x -> p, without shift:'

  CALL printparams(totsig,n)

! px0 -> xcheck
  CALL remapping(xcheck,xsf,px0,p0sl,p0su,pflag,lb,ub,a,b,n)

  map_error = maxval(abs(x-xcheck))
  IF (map_error .gt. eps) THEN
     WRITE(*,*) 'Maximal difference abs(x - x_re)=', map_error,' for parameter ', maxloc(abs(x-xcheck))
     WRITE(*,*) '   should be just numerical noise during mapping !?? '
  ENDIF

  x = xcopy + 1
  x(size(x)) = xcopy(size(x))  ! map for global offset is identity
  CALL mapping(x,xsf,px0,p0sl,p0su,pflag,lb,ub,a,b,n)
  WRITE(*,*)
  WRITE(*,*) '1. Results after p0 -> x -> x+sigma -> p=p+, right shift by sigma=1. '
  CALL printparams(totsig,n)

  x = xcopy - 1
  x(size(x)) = xcopy(size(x))  ! don't map global offset
  CALL mapping(x,xsf,px0,p0sl,p0su,pflag,lb,ub,a,b,n)
  WRITE(*,*)
  WRITE(*,*) '2. Results after p0 -> x -> x-sigma -> p=p-, left shift by sigma=1.'
  CALL printparams(totsig,n)

  WRITE(*,*)
  WRITE(*,*) '3. Compare deviations with uncertainties:'
  WRITE(*,*) '  - a priori errorbar err0 := p0sl+p0su'
  WRITE(*,*) '  - mapped and shift errorbar err_map := abs(p+-p0) + abs(p- -p0)'
  WRITE(*,*) '  - relative difference := err_map/err0'
  WRITE(*,*)

  WRITE(*,*) " Parameter                Errorbars                     Relative difference"
  WRITE(*,*) "                    a priori / mapped and shifted               "

! compute normalised difference
!!$         rel(i) = abs(gfd(i)-adx(i))
!!$         absmax = max( abs(gfd(i)), abs(adx(i)) )
!!$         if (absmax .lt. epsmach) then
!!$            rel(i) = 0.
!!$         else
!!$            rel(i) = rel(i) / absmax
!!$         end if

  WRITE(*, results) "vm(1)", p0sl(12)+p0su(12),totsig(12), totsig(12)/(p0sl(12)+p0su(12))
  WRITE(*, results) "jmf(1)", p0sl(25)+p0su(25),totsig(25), totsig(25)/(p0sl(25)+p0su(25))
  i = 27
  WRITE(*, results) "fautleaf(1)",p0sl(i)+p0su(i) ,totsig(i),  totsig(i)/(p0sl(i)+p0su(i))
  i=i+1
  WRITE(*, results) "ccost(1)", p0sl(i)+p0su(i) ,totsig(i),  totsig(i)/(p0sl(i)+p0su(i))
  i=i+1
  WRITE(*, results) "q10f(1)", p0sl(i)+p0su(i) ,totsig(i),  totsig(i)/(p0sl(i)+p0su(i))
   i=i+1
  WRITE(*, results) "q10s(1)",p0sl(i)+p0su(i) ,totsig(i),  totsig(i)/(p0sl(i)+p0su(i))
   i=i+1
  WRITE(*, results) "tauf(1)", p0sl(i)+p0su(i) ,totsig(i),  totsig(i)/(p0sl(i)+p0su(i))
   i=i+1
  WRITE(*, results) "aw(1)", p0sl(i)+p0su(i) ,totsig(i),  totsig(i)/(p0sl(i)+p0su(i))
   i=i+1
  WRITE(*, results) "fracs(1)   ", p0sl(i)+p0su(i) ,totsig(i),  totsig(i)/(p0sl(i)+p0su(i))
   i=i+1
  WRITE(*, results) "er(1)", p0sl(i)+p0su(i) ,totsig(i),  totsig(i)/(p0sl(i)+p0su(i))
   i=i+1
  WRITE(*, results) "ev(1)", p0sl(i)+p0su(i) ,totsig(i),  totsig(i)/(p0sl(i)+p0su(i))
   i=i+1
  WRITE(*, results) "eo(1)  ", p0sl(i)+p0su(i) ,totsig(i),  totsig(i)/(p0sl(i)+p0su(i))
   i=i+1
  WRITE(*, results) "ec(1)", p0sl(i)+p0su(i) ,totsig(i),  totsig(i)/(p0sl(i)+p0su(i))
   i=i+1
  WRITE(*, results) "ek(1)", p0sl(i)+p0su(i) ,totsig(i),  totsig(i)/(p0sl(i)+p0su(i))
   i=i+1
  WRITE(*, results) "alpha(1)", p0sl(i)+p0su(i) ,totsig(i),  totsig(i)/(p0sl(i)+p0su(i))
   i=i+1
  WRITE(*, results) "alc4(1)", p0sl(i)+p0su(i) ,totsig(i),  totsig(i)/(p0sl(i)+p0su(i))
   i=i+1
  WRITE(*, results) "kc0(1) ", p0sl(i)+p0su(i) ,totsig(i),  totsig(i)/(p0sl(i)+p0su(i))
   i=i+1
  WRITE(*, results) "ko0(1)", p0sl(i)+p0su(i) ,totsig(i),  totsig(i)/(p0sl(i)+p0su(i))
   i=i+1
  WRITE(*, results) "tgam(1)", p0sl(i)+p0su(i) ,totsig(i),  totsig(i)/(p0sl(i)+p0su(i))
   i=i+1
  WRITE(*, results) "beta(1)", p0sl(i)+p0su(i) ,totsig(i),  totsig(i)/(p0sl(i)+p0su(i))

END SUBROUTINE doit

SUBROUTINE printparams(totsig,n)
  USE mo_mapping 
  USE mo_vegetation 
  USE mo_beta
  USE costf
  USE tm

  IMPLICIT NONE

  CHARACTER (LEN= *), PARAMETER :: ptitle = '(A12, 2X, A15 , 2X, A15, 2X, A15, 2X, A9, 2X, A9)'
  CHARACTER (LEN= *), PARAMETER :: params = '(A12, 2X, F15.6, 2X, F15.6, 2X, F15.6, 2X, F15.6, 2X, F15.6)'

  INTEGER i,n
  REAL totsig(n),p(n)



  WRITE(*, ptitle) "Parameter", "p0" , "p", "p-p0","p0sl","p0su"

  totsig(12)=totsig(12)+abs(vm(1)-px0(12))
!+totsig(12)
  WRITE(*, params) "vm(1)", px0(12),vm(1), vm(1)-px0(12),p0sl(12),p0su(12)
  totsig(25)=totsig(25)+abs(jmf(1)-px0(25))
  WRITE(*, params) "jmf(1)", px0(25),jmf(1), jmf(1)-px0(25),p0sl(25),p0su(25)
  i = 27
  totsig(i)=totsig(i)+abs(px0(i)-fautleaf(1))
  WRITE(*, params) "fautleaf(1)", px0(i),fautleaf(1), fautleaf(1)-px0(i),p0sl(i),p0su(i)
  i=i+1
  totsig(i)=totsig(i)+abs(px0(i)-ccost(1))
  WRITE(*, params) "ccost(1)", px0(i),ccost(1), ccost(1)-px0(i),p0sl(i),p0su(i)
  i=i+1
  totsig(i)=totsig(i)+abs(px0(i)-q10f(1))
  WRITE(*, params) "q10f(1)", px0(i),q10f(1), q10f(1)-px0(i),p0sl(i),p0su(i)
   i=i+1
  totsig(i)=totsig(i)+abs(px0(i)-q10s(1))
  WRITE(*, params) "q10s(1)", px0(i),q10s(1), q10s(1)-px0(i),p0sl(i),p0su(i)
   i=i+1
  totsig(i)=totsig(i)+abs(px0(i)-tauf(1))
  WRITE(*, params) "tauf(1)", px0(i),tauf(1), tauf(1)-px0(i),p0sl(i),p0su(i)
   i=i+1
  totsig(i)=totsig(i)+abs(px0(i)-aw(1))
  WRITE(*, params) "aw(1)", px0(i),aw(1), aw(1)-px0(i),p0sl(i),p0su(i)
   i=i+1
  totsig(i)=totsig(i)+abs(px0(i)-fracs(1))
  WRITE(*, params) "fracs(1)", px0(i),fracs(1), fracs(1)-px0(i),p0sl(i),p0su(i)
   i=i+1
  totsig(i)=totsig(i)+abs(px0(i)-er(1))
  WRITE(*, params) "er(1)", px0(i),er(1), er(1)-px0(i),p0sl(i),p0su(i)
   i=i+1
  totsig(i)=totsig(i)+abs(px0(i)-ev(1))
  WRITE(*, params) "ev(1)", px0(i),ev(1), ev(1)-px0(i),p0sl(i),p0su(i)
   i=i+1
  totsig(i)=totsig(i)+abs(px0(i)-eo(1))
  WRITE(*, params) "eo(1)  ", px0(i),eo(1), eo(1)-px0(i),p0sl(i),p0su(i)
   i=i+1
  totsig(i)=totsig(i)+abs(px0(i)-ec(1))
  WRITE(*, params) "ec(1)", px0(i),ec(1), ec(1)-px0(i),p0sl(i),p0su(i)
   i=i+1
  totsig(i)=totsig(i)+abs(px0(i)-ek(1))
  WRITE(*, params) "ek(1)", px0(i),ek(1), ek(1)-px0(i),p0sl(i),p0su(i)
   i=i+1
  totsig(i)=totsig(i)+abs(px0(i)-alpha(1))
  WRITE(*, params) "alpha(1)", px0(i),alpha(1), alpha(1)-px0(i),p0sl(i),p0su(i)
   i=i+1
  totsig(i)=totsig(i)+abs(px0(i)-alc4(1))
  WRITE(*, params) "alc4(1)", px0(i),alc4(1), alc4(1)-px0(i),p0sl(i),p0su(i)
   i=i+1
  totsig(i)=totsig(i)+abs(px0(i)-kc0(1))
  WRITE(*, params) "kc0(1) ", px0(i),kc0(1), kc0(1)-px0(i),p0sl(i),p0su(i)
   i=i+1
  totsig(i)=totsig(i)+abs(px0(i)-ko0(1))
  WRITE(*, params) "ko0(1)", px0(i),ko0(1), ko0(1)-px0(i),p0sl(i),p0su(i)
   i=i+1
  totsig(i)=totsig(i)+abs(px0(i)-tgam(1))
  WRITE(*, params) "tgam(1)", px0(i),tgam(1), tgam(1)-px0(i),p0sl(i),p0su(i)
   i=i+1
  totsig(i)=totsig(i)+abs(px0(i)-beta(1))
  WRITE(*, params) "beta(1)", px0(i),beta(1), beta(1)-px0(i),p0sl(i),p0su(i)
   i=i+13
  totsig(i)=totsig(i)+abs(px0(i)-glob_offset)
  WRITE(*, params) "glob_offset", px0(i),glob_offset, glob_offset-px0(i),p0sl(i),p0su(i)
  
END SUBROUTINE printparams

