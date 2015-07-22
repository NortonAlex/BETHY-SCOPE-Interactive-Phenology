!------------------------------------------------------------------
! initialize first guess vector of control variables
!------------------------------------------------------------------
SUBROUTINE initmod( nvar, x )
  USE costf
  USE mo_trafo

  IMPLICIT NONE

  INTEGER, INTENT(in)  :: nvar
  REAL   , INTENT(out) :: x(nvar)
  INTEGER i  ! HEW-ADD

  x = p2x(px0, xpf, xpfa, xpfb, px0, p0su, a, nvar)

! MAS-RM-100407 removed xsf as scaling flag
!  i = size(xsf)
!  x(1:i) = x(1:i) * abs(xsf)

END SUBROUTINE initmod
