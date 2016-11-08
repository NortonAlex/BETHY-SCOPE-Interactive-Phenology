!================================================
! The subroutine "POSTJACREV" is called
! after computation of full Jacobian in reverse mode
!================================================
      SUBROUTINE POSTJACREV(  n, x, m, y, adx)
      implicit none
      integer n, m
      REAL   adx(m,n)
      REAL   x(n)
      REAL   y(m)

      integer i,j
      print *
      print *,    '**************************************************'
      print '(A)',' Jacobian'
      print *,    '**************************************************'
      print '(A)', ' col row       x(col)       y(row)  dy(row)/dx(col)'
      do i = 1, n
         do j = 1, m
            print '(2i4,2(1x,e12.6),(5x,e12.6)  )',&
     &           i,j,x(i),y(j),adx(j,i)
         end do
      end do
      print *
      END

      subroutine setfunc( n, m )
      implicit none
      integer n, m
      call numbmod( n )
      m = 1
      end subroutine setfunc

      subroutine initfunc( n, x )
      implicit none
      integer n
      REAL    x(n)
      call initmod( n, x )
      end subroutine initfunc

      subroutine func( n, x, m, y )
      implicit none

      integer :: n, m
      REAL    :: x(n), y(m)
      REAL    :: fc

!==============================================
! initialize adjoint variables
!==============================================
      y(:) = 0.0

!==============================================
! compute gradient and cost function
!==============================================
      call model( n, x, fc )
      y(1) = fc

      end subroutine

!================================================
! The subroutine "POSTFUNC" is called at last
! It should output the function value
!================================================
      SUBROUTINE POSTFUNC( N, X, M, Y)

      INTEGER N, M
      REAL    X(N), Y(M)
      INTEGER I, J

      PRINT *,'**************************************'
      PRINT *,'***         Function VALUE          **'
      PRINT *,'**************************************'
      DO J = 1,M
         PRINT '(X,I4,20(X,E12.6))', J, Y(J)
      ENDDO

      END

      SUBROUTINE POSTJACCOL  (nvar, x, m, y, g_y, col, ncol )
      IMPLICIT NONE
      integer nvar, m, ncol, col(ncol)
      REAL   g_y(ncol,m)
      REAL   x(nvar)
      REAL   y(m)
      
                                ! local variables
      INTEGER i_stat, i_year, i_month , i
      integer icol ,j
      
      print *
      print *,    '**************************************************'
      print '(A)',' Jacobian'
      print *,    '**************************************************'
      print '(A)', ' col row       x(col)       y(row)  dy(row)/dx(col)'
      do icol = 1, ncol
         do j = 1, m
            print '(2i4,2(1x,e12.6),(8x,e12.6)  )',&
     &           col(icol),j,x(col(icol)),y(j),g_y(icol,j)
         end do
      end do
      print *
      END 
