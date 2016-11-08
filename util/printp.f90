!**************************************************************
!   file:               printp.f90
!   purpose:            prints parameter of (choosen) OPWARMD
!   creation date:      07 02 
!                       10 03 HEW: modification w.r.t. parameter trafo
!                       
!   Copyright (C) 2000, 2001, 2003
!
!   CCDAS consortium
!
!   All rights reserved.
!**************************************************************
PROGRAM printp

  USE mo_namelist
  USE mo_mapping

  IMPLICIT NONE 

  REAL, ALLOCATABLE, DIMENSION(:) :: x, p0, p, p0sl, p0su, xsf, pflag, lb, ub, a, b
!  REAL, ALLOCATABLE, DIMENSION(:) :: x, p
  REAL :: fc
  INTEGER :: n_var,i,j,k
  CHARACTER (LEN= 20) :: filename, progName
  CHARACTER (LEN= *), PARAMETER :: ptitle = '(A12, 2X, A15 , 2X, A15, 2X, A15, 2X, A9, 2X, A9)'
  INTEGER, EXTERNAL :: IARGC

! Determine name of OPWARMD file
   j=IARGC()
   IF(j.GT.1) THEN
      CALL getarg(0,progName)       ! determine absolute program name, i.e. including path)
      PRINT*,"Usage: ",TRIM(progName)," [ <opwarmd file> ] "
      STOP
   END IF
   IF(j.EQ.0) THEN
      filename="OPWARMD"       ! Default name of opwarmd file
   ELSE
      CALL getarg(1,filename)  ! Name of opwarmd file
   END IF

! .. read namelist variables
  CALL get_namelist

! .. read parameter file
  OPEN(4,file='./input/pdf_params',status='old')
  REWIND 4
  OPEN(3,file=param_file,status='old')
  READ(3,*)
  READ(3,*) n_var
  write(*,*) 'n_var ',n_var
  ALLOCATE(x(n_var), p0(n_var), p(n_var), xsf(n_var), pflag(n_var), lb(n_var), ub(n_var))
  ALLOCATE(a(n_var) , b(n_var))
!  ALLOCATE(x(n_var), p(n_var))
  ALLOCATE(p0sl(n_var), p0su(n_var))
  DO i=1,n_var
     READ(3,*) p0(i),p0sl(i),p0su(i),xsf(i),pflag(i),lb(i),ub(i)
     READ(4,*) a(i),b(i)
!,p0(i),p0sl(i),p0su(i),xsf(i),pflag(i),lb(i),ub(i),a(i),b(i)
!,pflag(i),lb(i),ub(i),a(i),b(i)
  END DO
  CLOSE(3)
  CLOSE(4)

  OPEN(2,file=filename,form='unformatted',access='direct',recl=8*n_var,err=100,status='old')
  READ(2,rec=1) x
  CLOSE(2)

  WRITE(*, ptitle) "No", "initial p0" ,"predicted p", "prior unc. - / +","natural x"

!  DO i=1,n_var
  p=trafo(x,xsf,p0,p0sl,p0su,pflag,lb,ub,a,b,n_var)
  DO i=1,n_var
     if (p(i).le.lb(i)) then
        WRITE(*,*) 'ERROR in printp(): p less then lb after trafo for par ',i
        WRITE(*,*) 'pflag =',pflag(i)
     endif
     WRITE(6,*) i,p0(i),p(i),p0sl(i),p0su(i),x(i)
  END DO

200 CONTINUE
  PRINT*, "End of file ",filename, "reached"
  stop

100 CONTINUE
  PRINT*, "ERROR : Cannot open file ",filename
  stop

END PROGRAM printp
