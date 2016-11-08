!**************************************************************
!   file:               printparams.f90
!   purpose:            prints parameter of (choosen) OPWARMD
!   creation date:      07 02 
!                       10 03 HEW: modification for choose OPWARMD
!                       
!   Copyright (C) 2000, 2001, 2003
!
!   CCDAS consortium
!
!   All rights reserved.
!**************************************************************
PROGRAM printparams

  USE mo_namelist

  IMPLICIT NONE 

  REAL, ALLOCATABLE, DIMENSION(:) :: x0, x0sl,x0su,xsf,x,new_xsf,d,pflag
  REAL :: fc
  INTEGER :: n_var,i,j,k
  CHARACTER (LEN= 20) :: filename, progName


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
  OPEN(1,file=param_file,status='old')
  READ(1,*)
  READ(1,*) n_var
  ALLOCATE(d(n_var), x0(n_var), x(n_var), x0sl(n_var),x0su(n_var))
  ALLOCATE( xsf(n_var), new_xsf(n_var),pflag(n_var))
  DO i=1,n_var
     READ(1,*) x0(i),x0sl(i),x0su(i),xsf(i),pflag(i)
  END DO
  CLOSE(1)


  WRITE(*,*) filename
  OPEN(1,file=filename,form='unformatted',access='direct',recl=8*n_var)
  READ(1,rec=1) x
  CLOSE(1)
  
  DO i=1,n_var
     IF (pflag(i).GT.0) STOP 'Parameter is bounded, use printp!'
     IF (xsf(i)>0.) THEN
        WRITE(6,*)i,x0(i),x0sl(i),x(i)*x0sl(i)/ABS(xsf(i))
     ELSE
        WRITE(6,*)i,x0(i),x0sl(i),EXP((x(i)/ABS(xsf(i)))*LOG(x0sl(i)))
     ENDIF
  END DO

  STOP
END PROGRAM printparams
