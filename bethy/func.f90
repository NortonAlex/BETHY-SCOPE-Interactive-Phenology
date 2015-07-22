SUBROUTINE setfunc( n, m )
  IMPLICIT NONE
  INTEGER n, m
  CALL numbmod( n )
  m = n
END SUBROUTINE setfunc


SUBROUTINE initfunc( n, x )
  IMPLICIT NONE
  INTEGER n
  REAL    x(n)
  CALL initmod( n, x )
END SUBROUTINE initfunc


SUBROUTINE func( n, x, m, y )
!********************************************************************
!
! input  : N      amount of parameters
!          X      array of parameters
!
! output : FC     value of cost function
!          Y      gradient of cost function
!
!********************************************************************
  IMPLICIT NONE

  INTEGER :: n, m
  REAL    :: x(n), y(m)
  REAL    :: fc, adfc
  INTEGER :: i

!==============================================
! initialize adjoint variables
!==============================================
  y(:) = 0.0
  adfc = 1.0

!==============================================
! compute gradient and cost function
!==============================================
  CALL model_ad( n, x, y, fc, adfc )

END SUBROUTINE func


SUBROUTINE postfunc( n, x, m, y )
!********************************************************************
! postprocessing for Hessian run
!********************************************************************
  IMPLICIT NONE

  INTEGER :: n, m
  REAL    :: x(n), y(n)
  INTEGER :: j

      PRINT *,'**************************************'
      PRINT *,'***         Gradient VALUE          **'
      PRINT *,'**************************************'
      DO J = 1,M
         PRINT '(X,I4,20(X,E12.6))', J, Y(J)
      ENDDO

END SUBROUTINE postfunc
SUBROUTINE POSTJACCOL  (nvar, x, m, y, g_y, col, ncol )
  USE mo_namelist, ONLY: outdir
  IMPLICIT NONE
  INTEGER nvar, m, ncol, col(ncol)
  REAL   g_y(ncol,m)
  REAL   x(nvar)
  REAL   y(m)

  ! local variables
  INTEGER i_stat, i_year, i_month , i
  CHARACTER (len=80) :: ofname
  INTEGER icol ,j

  PRINT *
  PRINT *,    '**************************************************'
  PRINT '(A)',' Hessian columns'
  PRINT *,    '**************************************************'
  PRINT '(A)', ' col row       x(col)       y(row)  d2fc(row)/dx2(col)'
  DO icol = 1, ncol
     DO j = 1, m
        PRINT '(2i4,2(1x,e12.6),(8x,e12.6)  )',&
             & col(icol),j,x(col(icol)),y(j),g_y(icol,j)
     END DO
  END DO
  PRINT *
! save hessian
   ofname=TRIM(outdir)//'hesscol.bin'
   OPEN (39,file=ofname,status='unknown',form='unformatted')
   WRITE(39) g_y
   CLOSE(39)
END SUBROUTINE POSTJACCOL
