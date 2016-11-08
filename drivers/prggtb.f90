PROGRAM gtb_fortran
  IMPLICIT NONE
  integer n,m
  CALL setfunc(n,m)
! using subroutine structure for dynamic allocation
  call doit(n,m)
end PROGRAM gtb_fortran

subroutine doit(n,m)
  IMPLICIT NONE

  INTEGER :: i, j, m, n, check, itercount
!  PARAMETER (n = 3, m=109)
  REAL :: alpha, dirderiv, delta, normrquadr, wert, norm
  REAL :: x(n), xnew(n), y(m), ynew(m), ydummy(m), gradphi(n)
  REAL :: s1(n), p(n), As(m), r(m), help1(n+1), help(n), mu(n)
  REAL :: rho(n), y_k(m)

  REAL :: R0(n+1,n), Q1(m,n+1)
  REAL :: eps, max_iter
  PARAMETER (eps = 0.00001, max_iter =2000 )

 integer :: n_p_sec, ia, ie
  real :: t
  call system_clock(count_rate=n_p_sec)!once
  call system_clock(count=ia)

! FastOpt missing init
  r0 = 0.
  q1 = 0.

!  CALL setfunc(n,m)
  DO i = 1, n     !initialize A_0 = I
   R0(i,i) = 1.  
   Q1(i,i) = 1.
  END DO

  PRINT *, 'm : ', m, 'n : ', n
  CALL initfunc(n, x)
!  PRINT * , 'x0    : ', x 

  CALL func(n,x,m,y)
!  PRINT * ,  'y0    : ', y 
  DO i = 1, m 
     ydummy(i) = y(i)
  END DO
  CALL vecjac(n,x,gradphi,m,ydummy,y)
!  CALL vecjac(n,x,m,ydummy,gradphi)

  DO i = 1, m 
     y(i) = ydummy(i)
  END DO
  itercount = 0

  CALL norm2(gradphi,n,norm)
!  PRINT * ,  'norm gradphi    : ', norm

  DO WHILE (norm > eps .AND. itercount <= max_iter)
! Compute step from R^TR*p p=-F'(x)^TF(x)
     DO i = 1, n 
        gradphi(i) = -1.* gradphi(i)
     END DO
     CALL fbsolv(R0,p,gradphi,n)
!     PRINT * ,  'Schritt p : ', p 

!-gradphi^Tp for linesearch
     CALL scalar(gradphi, p , n, wert)
     dirderiv = -1.*wert 

     alpha = 1.

!Linesearch
     CALL quadratic_interpolation(alpha, x,dirderiv, p, m, n, check)
   
     IF(check==-1) THEN
        WRITE (*, *) 'bad quadratic interpolation - baling out \n'
     ENDIF
!     WRITE (*,*) 'Calculated step length is ', alpha

!Making the step
     DO i = 1, n
        s1(i) =  alpha*p(i)
        xnew(i) = x(i) + s1(i)
     END DO

! Update of A=QR (transposed Broyden)
     CALL func(n,xnew,m,ynew)
     CALL vecjac(n,xnew,gradphi,m,ydummy,ynew)

     DO i = 1, m
        ydummy(i) = ynew(i)
     END DO
!     CALL vecjac(n,xnew,m,ydummy,gradphi)
 
     DO i = 1, m
        y_k(i) = (ynew(i) - y(i))/ alpha
     END DO

     CALL ttrimatvec(R0,n,s1,n+1,help1)
     CALL matvec(Q1,n+1,help1,m,As)

     DO i = 1, m
        r(i)=y_k(i) - As(i)
        ydummy(i) = r(i)
     END DO

     CALL vecjac(n,xnew,mu,m,ynew,ydummy)
!     CALL vecjac(n,xnew,m,ydummy,mu)
  
     CALL  mattransvec(Q1,n+1,r,m,help1)
     CALL ttrimattransvec(R0,n,help1,n+1,help)

     CALL scalar(r,r,m, normrquadr)
     DO i = 1, n
        rho(i) = (mu(i) - help(i))/ normrquadr
     END DO
 
     CALL qrrankone(Q1, R0, r, rho, m, n)
 !    PRINT *, 'Q1 ', Q1
!     PRINT *, 'R0 ', R0

! Update 
     DO i = 1, n
        x(i) = xnew(i)
     END DO
     DO i = 1, m
        y(i) = ynew(i)
     END DO

     itercount = itercount + 1
     CALL norm2(gradphi,n,norm)
  END DO
  call system_clock(count=ie)
  WRITE (*,*) 'A solution is found afer ', itercount, 'iterations.'
  PRINT *, 'Norm gradphi = ', norm
  PRINT *, 'It is x = '
  PRINT *, x
  CALL eval_penalty(x, wert, m, n)
  PRINT *, 'And 0.5*|F(x)|^2 = ', wert
  t = (ie-ia)/real(n_p_sec)
  write(unit=*,fmt=*)  "CPU time in sec:",t
END subroutine doit
