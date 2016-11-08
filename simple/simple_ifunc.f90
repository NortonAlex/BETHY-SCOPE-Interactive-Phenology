SUBROUTINE setfunc( nbgr, m )
  IMPLICIT NONE

  INTEGER nbgr, m, n

  CALL numbmod( n )
  m = n
  nbgr = 1
END SUBROUTINE setfunc

SUBROUTINE initfunc( nbgr, xbgr )
  IMPLICIT NONE
  INTEGER nbgr
  REAL    xbgr(nbgr)
  xbgr = 0.
END SUBROUTINE initfunc

SUBROUTINE func( nbgr, xbgr, m, y )
  use coef
  IMPLICIT NONE

  INTEGER :: nbgr
  REAL    :: xbgr(nbgr)

  INTEGER :: m
  REAL    :: x(m), y(m)
  REAL    :: fc, adfc
  INTEGER :: i,j,ibgr
  REAL, PARAMETER :: sfac = 1.e12! background flux pert in GtC/gridcell/year

!==============================================
! initialize model and parameters
!==============================================
  CALL initmod( m, x )

!==============================================
! mapping
!==============================================
  a = a + xbgr(1)

!==============================================
! initialize adjoint variables
!==============================================
  y(:) = 0.0
  adfc = 1.0

!==============================================
! compute gradient and cost function
!==============================================
  CALL model_ad( m, x, y, fc, adfc )

!==============================================
! reset f_bgr (to allow multiple fd tests)
!==============================================
  a = a - xbgr(1)

END SUBROUTINE func


SUBROUTINE postfunc( nbgr, x, m, y )
!********************************************************************
! postprocessing for Hessian run
!********************************************************************
  IMPLICIT NONE

  INTEGER :: nbgr, m
  REAL    :: x(nbgr), y(m)
  INTEGER :: j

      PRINT *,'**************************************'
      PRINT *,'***         Gradient VALUE          **'
      PRINT *,'**************************************'
      DO J = 1,M
         PRINT '(X,I4,20(X,E12.6))', J, Y(J)
      ENDDO

END SUBROUTINE postfunc
SUBROUTINE POSTJACCOL  (nvar, x, m, y, g_y, col, ncol )
  IMPLICIT NONE
  INTEGER nvar, m, ncol, col(ncol)
  REAL   g_y(ncol,m)
  REAL   x(nvar)
  REAL   y(m)

  ! local variables
  INTEGER i_stat, i_year, i_month , i
  INTEGER icol ,j

  PRINT *
  PRINT *,    '**************************************************'
  PRINT '(A)',' Derivatives'
  PRINT *,    '**************************************************'
  PRINT '(A)', ' col row       b(col)     dJdx(row)  d2Jdx(row)/db(col)'
  DO icol = 1, ncol
     DO j = 1, m
        PRINT '(2i4,2(1x,e12.6),(8x,e12.6)  )',&
             & col(icol),j,x(col(icol)),y(j),g_y(icol,j)
     END DO
  END DO
END SUBROUTINE POSTJACCOL

SUBROUTINE POSTJACREV(  n, x, m, y, adx)
  IMPLICIT NONE
  INTEGER n, m
  REAL   adx(m,n)
  REAL   x(n)
  REAL   y(m)
  
  INTEGER i,j

  PRINT *
  PRINT *,    '**************************************************'
  PRINT '(A)',' Jacobian'
  PRINT *,    '**************************************************'
  PRINT '(A)', ' col row       x(col)       y(row)  dy(row)/dx(col)'
  DO i = 1, n
     DO j = 1, m
        PRINT '(2i4,2(1x,e12.6),(5x,e12.6)  )',&
             &           i,j,x(i),y(j),adx(j,i)
     END DO
  END DO
END SUBROUTINE POSTJACREV
