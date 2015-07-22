MODULE fluo_func 

IMPLICIT NONE 


CONTAINS 

! f corresponds to Psofunction in the matlab tool of Van der Tol 
FUNCTION func(K,km,LAI,q,dso,xl)
IMPLICIT NONE

! Input variables
REAL, INTENT(IN)                 :: xl
REAL, INTENT(IN)                 :: K,km,LAI,q,dso

! Ouput variables
double precision                 ::  func

! Local variables
REAL                             :: alf

if (dso/=0) then
    alf          =   (dso/q) *2/(km+K)
   func  =   exp((K+km)*LAI*xl + sqrt(K*km)*LAI/(alf  )*(1-exp(xl*(alf )))) !% [nl+1]  factor for correlation of Ps and Po
else
   func  =   exp((K+km)*LAI*xl - sqrt(K*km)*LAI*xl)  !% [nl+1]  !factor for correlation of Ps and Po
endif

return
END FUNCTION func


SUBROUTINE  simpson(k1,k2,lai,q,dso,func,a,b,integral,n)
!==========================================================
! Integration of f(x) on [a,b]
! Method: Simpson rule for n intervals  
! written by: Alex Godunov (October 2009)
!----------------------------------------------------------
! IN:
! func   - Function to integrate (supplied by a user)
! a       - Lower limit of integration
! b       - Upper limit of integration
! n   - number of intervals
! OUT:
! integral - Result of integration
!==========================================================
implicit none
double precision func, a, b, integral,s
double precision h, x
REAL  :: k1,k2,lai,q,dso
integer nint
integer n, i,j

do j=1,16
   n = n*2

! if n is odd we add +1 to make it even
if((n/2)*2.ne.n) n=n+1

! loop over n (number of intervals)
s = 0.0
h = (b-a)/dfloat(n)

do i=2, n-2, 2
   x   = a+dfloat(i)*h
   s = s + 2.0*func(k1,k2,lai,q,dso,x) + 4.0*func(k1,k2,lai,q,dso,x+h)
end do
integral = (s + func(k1,k2,lai,q,dso,a) + func(k1,k2,lai,q,dso,b) + &
         &  4.0*func(k1,k2,lai,q,dso,a+h))*h/3.0
end do 

return
END SUBROUTINE simpson


SUBROUTINE  expint(ind,nk,x,intexp) 
!  
!   ROUTINE:   expint
!   PURPOSE:   compute the exponential integral of x.
!              the exponential integral of index n is given by
!              integral( exp(-tx)/x^n dx), x ranges from 1 to infinity
!  
!   USEAGE:    result=expint(ind,x,nz)
!  
!   INPUT:
!     ind      type of exponential integral, for example use ind=1
!              to get first exponential integral, E1
!     x        argument to exponential integral ( 0 < x < infinity)
!
!   OUTPUT:
!     result   exponential integral
!  
     IMPLICIT NONE 
    
   ! Input variables  
     INTEGER,INTENT(IN)               :: ind,nk 
     REAL, INTENT(IN), DIMENSION(nk)  :: x

  ! output variables 
    REAL,INTENT(OUT), DIMENSION(nk)   :: intexp 

  ! Local variables  
     REAL, DIMENSION(4)               :: a
     REAL, DIMENSION(4)               :: b
     REAL, DIMENSION(6)               :: c
     REAL,   DIMENSION(nk)            :: expx,a1,a2 
     !REAL, ALLOCATABLE, DIMENSION(nk) :: expx,a1,a2 
     INTEGER                          :: n 



       a(1) = 0.2677737343
       a(2) = 8.6347608925
       a(3) = 18.0590169730
       a(4) = 8.5733287401

       b(1) = 3.9584969228
       b(2) = 21.0996530827
       b(3) = 25.6329561486
       b(4) = 9.5733223454
 
       c(1) = -0.57721566
       c(2) = 0.99999193
       c(3) = -0.24991055
       c(4) = 0.05519968
       c(5) = -0.00976004
       c(6) =  0.00107857
!
      if(ind .le. 0) then 
        write(*,*) 'illegal value of IND in EXPINT'
        stop
      endif
! 
!   print*, ' nk ', nk
!   print*, ' size x ', size(x)
!   print*, ' x ', minval(x), maxval(x), sum(x)
!   print*, ' ind ', ind 

   expx=exp(-x)

  !  print*, ' expx ', minval(expx), maxval(expx), sum(expx)

      WHERE (x .le. 1.) 
        intexp = ((((c(6)*x+c(5))*x+c(4))*x+c(3))*x+c(2))*x+c(1)-log(x)
      ELSEWHERE
        a1 = (((x+a(4))*x+a(3))*x+a(2))*x+a(1)
        a2 = (((x+b(4))*x+b(3))*x+b(2))*x+b(1)
        intexp = expx*a1/(a2*x)
      ENDWHERE 

  !print*, ' bio_func intexp ', minval(intexp), maxval(intexp), sum(intexp) 

! use recursion to get higher order exponential integrals
      if(ind.eq.1) return
      do 10 n = 1,ind-1 
        intexp = (expx-x*intexp)/n 
 10   continue
      return
  END SUBROUTINE expint 


SUBROUTINE mask_ind(n,mask,mask_indices)
! returns a packed list  of subscripts for which mask is true
  implicit none
  integer                           :: n
  logical, intent(in), dimension(n) :: mask
  !integer, dimension(count(mask)) :: mask_ind
  integer, dimension(count(mask))  :: mask_indices
  integer :: i_mask, i_ind
  integer, dimension(:), allocatable, save  :: indices
  logical, save :: first = .true.
! we keep a list of indices to avoid calculating  them each time, was causing
! memory  leak


  if( first) then ! allocate and assign indices
     first = .false.
     allocate( indices( size( mask)))
     indices = (/ ( i_ind, i_ind= 1, size( mask) ) /)
  endif
  ! keep biggest index list we need
  if( size(mask) > size( indices)) then ! we need to reallocate and expand
     deallocate( indices)
     allocate( indices( size( mask)))
     indices = (/ ( i_ind,i_ind=1, size( mask) ) /)
  endif
 mask_indices = pack( indices, mask)

RETURN
END SUBROUTINE mask_ind


! The routine interpolates the data y at xi with output yi  
SUBROUTINE interpol(n,ni,x,y,xi,yi)

IMPLICIT NONE 

! Input variables 
INTEGER, INTENT(IN)                             :: n, ni
REAL, INTENT(IN), DIMENSION(n)                  :: x,y 
REAL, INTENT(IN), DIMENSION(ni)                 :: xi 

! Ouput variables 
REAL, INTENT(OUT), DIMENSION(ni)                 :: yi

! Local variables 
INTEGER                                         :: jdeb 
INTEGER                                         :: i,j 
REAL                                            :: a, b 


jdeb=1
DO i=1, ni 


   ! Search for the two points that bound the point xi 
   DO j=jdeb, n-1
        if ((x(j).le.xi(i)).and.(x(j+1).ge.xi(i)) ) then 
 !      print*, ' xj xi xj+1 ', i,j,x(j), xi(i), x(j+1)
        ! We linearly interpolate the data y  at xi  ==> yi = a*xi +b 
        ! The slope  a of the line  
         a = (y(j+1)-y(j))/(x(j+1)-x(j))

        ! The ordinate of the origin        
         b = y(j) - a*x(j)
        
        ! yi is given by :
          yi(i) = a*xi(i) + b          
          
        ! We start the next j loop at jdeb 
          jdeb =  j 
         exit 
        endif
   END DO 

 
END DO 

!DO i=1,ni 
!print*, yi(i)
!END DO 

END SUBROUTINE interpol 

!*****************************************************
!* Sorts an array ARR of length N in ascending order *
!* by straight insertion.                            *
!* ------------------------------------------------- *
!* INPUTS:                                           *
!*          N     size of table ARR                  *
!*          ARR   table to be sorted                 *
!* OUTPUT:                                           *
!*          ARR   table sorted in ascending order    *
!*                                                   *
!* NOTE: Straight insertion is a N² routine and      *
!*       should only be used for relatively small    *
!*       arrays (N<100).                             *
!*****************************************************         
SUBROUTINE unique(N,ARR,res,ks)
  INTEGER, INTENT(IN)                   :: N 
  REAL, INTENT(INOUT),DIMENSION(N)      :: ARR
  REAL, INTENT(OUT), DIMENSION(N)       :: res
  INTEGER, INTENT(OUT)                  :: ks 
  REAL                                  :: a 
  INTEGER                               :: i,j,k 

  do j=2, N
    a=ARR(j)
    do i=j-1,1,-1
      if (ARR(i)<=a) goto 10
      ARR(i+1)=ARR(i)
    end do
        i=0
10  ARR(i+1)=a

  end do

! Revove duplicates values  
  k = 1
  res(1) = ARR(1)
  outer: do i=2,N
     do j=1,k
        if (res(j) == ARR(i)) then
           ! Found a match so start looking again
           cycle outer
        end if
     end do
     ! No match found so add it to the output
     k = k + 1
     res(k) = ARR(i)
  end do outer

  ks =  k 

  return

END SUBROUTINE unique 

SUBROUTINE Sint(nx,y,x,yi)

!    % Simpson integration
!    % x and y must be any vectors (rows, columns), but of the same length
!    % x must be a monotonically increasing series

!    % WV Jan. 2013, for SCOPE 1.40
! Translate in Fortran 2013, April, E Koffi

IMPLICIT NONE

! Input variables
INTEGER, INTENT(IN)                :: nx
REAL, DIMENSION(nx), INTENT(IN)    :: x,y

! Ouput variables
REAL, INTENT(OUT)                  :: yi

! Local variables
REAL,DIMENSION(nx-1)               :: step, mean

!    nx   = length(x);
!    if (size(x,1) == 1) THEN
!        x = transpose(x)
!    endif

!    if (size(y,1).ne.1) then
!        y = transpose(y)
!    endif

    step = x(2:nx) - x(1:nx-1)
    mean = .5 * (y(1:nx-1) + y(2:nx))
    yi  = dot_product(mean, step)

END SUBROUTINE Sint


SUBROUTINE satvap(npts,T,es,s)
!function [es,s] = satvap(T)
!%% function [es,s]= satvap(T)
!% Author: Dr. ir. Christiaan van der Tol
!% Date: 2003
!%
!% calculates the saturated vapour pressure at
!% temperature T (degrees C)
!% and the derivative of es to temperature s (kPa/C)
!% the output is in mbar or hPa. The approximation formula that is used is:
!% es(T) = es(0)*10^(aT/(b+T));
!% where es(0) = 6.107 mb, a = 7.5 and b = 237.3 degrees C
!% and s(T) = es(T)*ln(10)*a*b/(b+T)^2

! Inputs
INTEGER, INTENT(IN)               :: npts
REAL, DIMENSION(npts), INTENT(IN) :: T

! Outputs
REAL, INTENT(OUT),DIMENSION(npts) :: es,s

! Local
REAL                              :: a,b

!%% constants
a  = 7.5
b  = 237.3              !; %degrees C

!%% calculations
es  = 6.107*10.**(7.5*T/(b+T))
s   = es*log(10.)*a*b/(b+T)**2

RETURN
END SUBROUTINE satvap


END MODULE fluo_func 
