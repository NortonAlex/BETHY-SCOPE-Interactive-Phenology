!====================================================================
! The subroutine evaluates the cost function "FC"
! depending on the control vector "X(N)".
!====================================================================
      SUBROUTINE FUNC( N, X, M, Y )
      USE costf
      implicit none
      INTEGER N, M
      REAL    X(N), fc, Y(M)

      CALL model( n, x, fc)
      CALL getlsq( n, x, m, y)

      END

!$taf subroutine func adname = func_gtb
      subroutine vecjac( n, x, res, m, y, v)
      implicit none
      integer n, m
      real x(n)
      real y(m)
      real v(m)
      real res(n)

      res = 0.
      call func_gtb( n, x, res, m, y, v )
      end

SUBROUTINE setfunc( n, m )
  use tm
  IMPLICIT NONE
  INTEGER n, m
  CALL numbmod( n )
  m = size(conc) + n
END SUBROUTINE setfunc
SUBROUTINE initfunc( n, x )
  IMPLICIT NONE
  INTEGER n
  REAL x(n)
  CALL initmod( n, x )
END SUBROUTINE initfunc

SUBROUTINE postfunc( n, x, m, y )
!********************************************************************
! postprocessing for Hessian run
!********************************************************************
  IMPLICIT NONE

  INTEGER :: n, m
  REAL :: x(n), y(m)
  INTEGER :: j

      PRINT *,'**************************************'
      PRINT *,'***         Function VALUE          **'
      PRINT *,'**************************************'
      DO J = 1,M
         PRINT '(X,I4,20(X,E12.6))', J, Y(J)
      ENDDO

END SUBROUTINE postfunc
