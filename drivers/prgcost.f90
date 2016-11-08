!**************************************************************
!   file: 	        prgcost.f90
!   purpose: 	        driver to evaluate
!                       scalar valued function twice
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
PROGRAM prgcost
!*************************************************************
  
  IMPLICIT NONE
  
!=========================================
! declaration
!=========================================
  INTEGER n

!==============================================================
! get the number of control variables
!==============================================================
  CALL numbmod( n )
  
!-----------------------------------------
! call the subroutine
!-----------------------------------------
  CALL doit( n )

END PROGRAM prgcost
      
SUBROUTINE doit( n )


  IMPLICIT NONE
      
  INTEGER n
  REAL    fc
  REAL   x(n)
  REAL   adx(n)

  INTEGER nn
  INTEGER ifail
  LOGICAL cold
  INTEGER jmin, jmax, nupdate, ifunc


!==============================================================
! initialisize the model
! and set the first guess of the control variables
!==============================================================

  CALL initmod( n, x )

!==============================================================
! option to read control vector after an optimisation
!==============================================================
  CALL instoref( nn, ifunc, fc, nupdate, jmin, jmax, cold, ifail )
  
  IF (cold) THEN
     PRINT *, 'control variables are set by initmod'
  ELSE
     IF (n .NE. nn) THEN
        STOP 'inconsitence number of control variables'
     ENDIF
     CALL dostoref( n,   x, .FALSE., 1 ,'OPWARMD' )
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

!===============================================
! fist call of model
!===============================================
  CALL model( n, x, fc )
  
  WRITE(*,*)
  WRITE(*,*) ' the cost function is : ', fc
      
!===============================================
! postprocessing
!===============================================
  CALL postmod( n, x, fc )

  END SUBROUTINE doit
