MODULE mo_helper
  
  IMPLICIT NONE
  
  ! PARAMETERS FOR AUXILLIARY FUNCTIONS
  REAL, PARAMETER :: zmin = 1e-18
!  REAL, PARAMETER :: eta = 0.99999            ! curvature parameter for mins/maxs
  
CONTAINS
    
  !*********************************************************
  !* FUNCTION errf
  !* the (cumulative) error function
  !* numerical recipes in Fortran 77, Chapter 6.2
  !*********************************************************
  
  REAL FUNCTION errf (x)
    REAL :: x, z, t
    z=ABS(x) 
    t=1./(1.+0.5*z) 
    errf=1.-0.5*t*EXP(-z*z-1.26551223+t*(1.00002368+t*(.37409196+ &
         t*(.09678418+t*(-.18628806+t*(.27886807+t*(-1.13520398+ &
         t*(1.48851587+t*(-.82215223+t*.17087277)))))))))
    IF (x<0) errf=1.-errf
  END FUNCTION errf

  !*********************************************************
  !*  FUNCTION mins
  !*  smoothed minimum function
  !*********************************************************

!  REAL FUNCTION mins (x, y)
!    REAL :: x, y
!    REAL :: z
!    z = (x+y)**2 - 4.*eta*x*y
!    IF (z.GE.zmin) THEN
!       mins = (x + y - SQRT(z)) / (2.*eta)
!    ELSE
!       mins = 0.
!    ENDIF
!  END FUNCTION mins
  REAL FUNCTION mins (x, y, eta)
    REAL :: x, y, eta
    REAL :: z
    z = (x+y)**2 - 4.*eta*x*y
    z = max (z, zmin)
    mins = (x + y - SQRT(z)) / (2.*eta)
  END FUNCTION mins

  !*********************************************************
  !*  FUNCTION maxs
  !*  smoothed maximum function
  !*********************************************************

!  REAL FUNCTION maxs (x, y)
!    REAL :: x, y
!    REAL :: z
!    z = (x+y)**2 - 4.*eta*x*y
!    IF (z.GE.zmin) THEN
!       maxs = (x + y + SQRT(z)) / (2.*eta)
!    ELSE
!       maxs = 0.
!    ENDIF
!  END FUNCTION maxs
  REAL FUNCTION maxs (x, y, eta)
    REAL :: x, y, eta
    REAL :: z
    z = (x+y)**2 - 4.*eta*x*y
    z = max (z, zmin)
    maxs = (x + y + SQRT(z)) / (2.*eta)
  END FUNCTION maxs

  !*********************************************************
  !*  FUNCTION minx
  !*  minimum function with exponential transition
  !*********************************************************

  REAL FUNCTION minx (x, y, x0)
    REAL :: x, y, x0
    IF (x.LE.y+x0) THEN
       minx = x - x0*EXP((x-y)/x0-1.)
    ELSE
       minx = y
    ENDIF
  END FUNCTION minx

  !*********************************************************
  !*  FUNCTION maxx
  !*  maximum function with exponential transition
  !*********************************************************

  REAL FUNCTION maxx (x, y, x0)
    REAL :: x, y, x0
    IF (x.GE.y-x0) THEN
       maxx = x + x0*EXP(-(x-y)/x0-1.)
    ELSE
       maxx = y
    ENDIF
  END FUNCTION maxx

END MODULE mo_helper
