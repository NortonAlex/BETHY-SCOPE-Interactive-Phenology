!================================================
! The subroutine "POSTJACCOL" is called
! after computation of Hessian by prgjaccol
!================================================
      SUBROUTINE POSTJACCOL(  n, x, m, y, g_y, col, ncol )
      implicit none
      integer n, m, ncol, col(ncol)
      REAL   g_y(ncol,m)
      REAL   x(n)
      REAL   y(m)

      integer icol ,j
      print *
      print *,    '**************************************************'
      print '(A)',' Jacobian columns'
      print *,    '**************************************************'
      print '(A)', ' col row       x(col)       y(row)  dy(row)/dx(col)'
      do icol = 1, ncol
         do j = 1, m
            print '(2i4,2(1x,e12.6),(5x,e12.6)  )',&
     &           col(icol),j,x(col(icol)),y(j),g_y(icol,j)
         end do
      end do
      print *
      END

      subroutine setfunc( n, m )
      implicit none
      integer n, m
      call numbmod( n )
      m = n
      end

      subroutine initfunc( n, x )
      implicit none
      integer n
      REAL    x(n)
      call initmod( n, x )
      end

      subroutine func( n, x, m, y )
      implicit none

      integer  n, m
      REAL     x(n), y(n)
      REAL     fc, adfc
      integer  i

!==============================================
! initialize adjoint variables
!==============================================
      do i=1,n
         y(i) = 0.0
      enddo
      adfc = 1.0

!==============================================
! compute gradient and cost function
!==============================================
      call model_ad( n, x, y, fc, adfc )

      end

!================================================
! The subroutine "POSTFUNC" is called at last
! It should output the function value
!================================================
      SUBROUTINE POSTFUNC( N, X, M, Y)

      INTEGER N, M
      REAL    X(N), Y(M)
      INTEGER I, J

      PRINT *,'**************************************'
      PRINT *,'***         Gradient VALUE          **'
      PRINT *,'**************************************'
      DO J = 1,M
         PRINT '(X,I4,20(X,E12.6))', J, Y(J)
      ENDDO

      END
