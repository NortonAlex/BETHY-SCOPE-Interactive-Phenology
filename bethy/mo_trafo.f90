MODULE mo_trafo
!------------------------------------------------------------------
! transform unitless parameters x to physical parameters p and vice versa
!------------------------------------------------------------------
    
CONTAINS
  
!------------------------------------------------------------------
! ************** x2p 
!    maps a natural control variable x onto its corr. physical model parameter p
!    depending on its bounds
!  05/05: Txk introducing cap on x to avoid overflow for lognormal trafo
!------------------------------------------------------------------
FUNCTION x2p(x, x0, xpf, xpfa, xpfb, p0, ps0u, a, nx)

  IMPLICIT NONE
  
  INTEGER nx,nvar,k
  REAL, DIMENSION(:) ::  x,x0,p0,ps0u,a,xpfa,xpfb
  INTEGER, DIMENSION(:) ::  xpf
  REAL, DIMENSION(nx) :: x2p
  
  
  ! INPUT:
  ! x    actual value of control variable
  ! p0     first guess of physical parameter
  ! ps0u    a priori uncertainties
  
  nvar = SIZE(p0)

  DO k = 1,nvar    ! loop over all parameters, PFT dependents implicit
     
     IF (abs(xpf(k)) == 3) THEN
        x2p(k) = EXP(a(k)*x(k))
     ELSE IF (abs(xpf(k)) == 2) THEN
        x2p(k) = a(k) * x(k)**2 + xpfa(k)
     ELSE IF (abs(xpf(k)) == 5) THEN
        x2p(k) = ps0u(k) * x0(k) + xpfa(k)        
     ELSE IF (abs(xpf(k)) == 4) THEN
        x2p(k) = (xpfb(k)-xpfa(k))/(1+EXP(-x(k)*a(k))) + xpfa(k)
     ELSE
        x2p(k) = ps0u(k) * x(k)                  
     ENDIF
  ENDDO

  IF (nx > nvar) THEN
     DO k=nvar+1,nx
        x2p(k) = x(k)
     ENDDO
  ENDIF

  RETURN

END FUNCTION x2p

!------------------------------------------------------------------
! ************** p2x 
!    maps a physical model parameter p to its natural control variable x
!    depending on its bounds
!------------------------------------------------------------------
FUNCTION p2x(p,xpf,xpfa,xpfb,p0,ps0u,a,nx)

  IMPLICIT NONE

  INTEGER nx,nvar,k
  REAL, DIMENSION(:) ::  p,p0,ps0u,a,xpfa,xpfb
  INTEGER, DIMENSION(:) :: XPF
  REAL, DIMENSION(nx) :: p2x

  nvar = SIZE(p)

  DO k = 1,nvar    ! loop over all parameters, PFT dependents implicit     

     IF (abs(xpf(k)) == 3) THEN
        p2x(k) =  LOG(p(k)-xpfa(k))/a(k)
     ELSE IF (abs(xpf(k)) == 2) THEN
        p2x(k) =  SQRT((p(k)-xpfa(k))/a(k))
     ELSE IF (abs(xpf(k)) == 5) THEN
        p2x(k) =  p0(k) / ps0u(k)
     ELSE IF (abs(xpf(k)) == 4) THEN
        p2x(k) =  -LOG((xpfb(k)-xpfa(k))/(p(k)-xpfa(k)) - 1)/a(k)
     ELSE
        p2x(k) =  p(k) / ps0u(k) 
     ENDIF

  ENDDO

  IF (nx > nvar) THEN
     DO k=nvar+1,nx
        p2x(k) = 0.
     ENDDO
  ENDIF

END FUNCTION p2x

END MODULE mo_trafo
