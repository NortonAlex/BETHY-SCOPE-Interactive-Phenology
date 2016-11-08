MODULE coef
  IMPLICIT NONE 
  REAL    A, B, C, D, E, F
END MODULE coef

!================================================
! This subroutine sets the number
! of control variables
!================================================
      SUBROUTINE NUMBMOD( N )
      use coef
      INTEGER N
  A =  1.E-3
  B = 1.E-2
  C =  1.E-1
  D =  1.E0
  E = 1.E1
  F =  1.E2

      N = 3
      END

!================================================
! The subroutine "MODEL" is called by the
! optimization procedure.
! It has to calculate the cost function "FC"
! depending on the control vector "X(N)".
!================================================
      SUBROUTINE MODEL( N, X, FC )
      use coef
      implicit none
      INTEGER N
      REAL    X(N), FC

      real fc1, fc2

      call bethy (n, x, 3, 10, fc1)
      call bethy (n, x, 4, 2, fc2)
      fc = fc1 + fc2

      END

      SUBROUTINE BETHY( N, X, nt, NG, FC )
      use coef
      implicit none
      INTEGER N, ng, nt
      REAL    X(N), FC

      integer:: ig, it
      real, allocatable :: z(:)

      allocate(z(ng))
      z = x(3)+3

      fc = 0.
      do it = 1, nt
         do ig = 1, ng
            z(ig) = 0.1 * z(ig) ** 2
         enddo
      enddo
      fc = a * (sum(z**2) + (x(1)-2)**2) + x(2)**2

!FastOpt      fc = abs(x(3))  ! just to produce a V-shaped cost function for testing

      deallocate(z)
      END

!================================================
! The subroutine "INITMOD" is called before the
! optimization. It must set a first guess
! of the parameter vector.
! It may also contain the initialization of
! the model.
!================================================
      SUBROUTINE INITMOD( N, X )
      implicit none
      INTEGER N
      REAL    X(N)
      integer i

      do i=1,n
         X(i) = float(i)
      enddo

      END

!================================================
! The subroutine "POSTMOD" is called after the function is evalutated
! It should contain the output of the function value
!================================================
      SUBROUTINE POSTMOD( N, X, FC)

      INTEGER N
      REAL    X(N), FC
      INTEGER I

      PRINT *,'**************************************'
      PRINT *,'***     FUNCTION VALUE              **'
      PRINT *,'**************************************'
      PRINT 9010, FC
      PRINT 9020
      DO I = 1,N
        PRINT 9030, I, X(I)
      ENDDO
 9010 FORMAT(1X,'The value of the cost function is : ',E15.9)
 9020 FORMAT(1X,'The parameter values are :')
 9030 FORMAT(1X,I5,1X,E15.9,1X,E15.9)

      END

!================================================
! The subroutine "POSTADM" is called after
! the optimization.
! It should contain the output of the results.
!================================================
      SUBROUTINE POSTADM( N, X, FC, ADX )

      INTEGER N
      REAL    X(N), FC, ADX(N)
      INTEGER I

      PRINT *,'**************************************'
      PRINT *,'***     RESULT OF OPTIMISATION      **'
      PRINT *,'**************************************'
      PRINT 9010, FC
      PRINT 9020
      DO I = 1,N
        PRINT 9030, I, X(I), ADX(I)
      ENDDO
 9010 FORMAT(1X,'The value of the cost function is : ',E15.9)
 9020 FORMAT(1X,'The optimal parameter values, gradients are :')
 9030 FORMAT(1X,I5,1X,E15.9,1X,E15.9)

      END
