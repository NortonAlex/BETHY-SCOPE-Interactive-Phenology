!*********************************************************
!*  SUBROUTINE budget
!*  calculates net flux and grid cell het. respiration
!*********************************************************

SUBROUTINE budget (ng, vp, nrun, outint, scale, cs, lfac, rnpp, resfv, ressv, ress, resf, rnep, rnepp, rnppv)

! .. Use statements
  USE mo_grid, ONLY: gridp, frac

  IMPLICIT NONE

! .. Arguments
  Integer, intent(in) :: ng, vp, nrun, outint, scale
  REAL, DIMENSION (vp), INTENT(inout) :: cs, lfac
  REAL, DIMENSION(0:nrun,outint,vp), INTENT(inout) :: resfv, ressv
  REAL, DIMENSION(0:nrun,outint,ng,3), INTENT(inout) :: rnepp
  REAL, DIMENSION(0:nrun,outint,ng), INTENT(in) :: rnpp
  REAL, DIMENSION(0:nrun,outint,vp), INTENT(in) :: rnppv
  REAL, DIMENSION(0:nrun,outint,ng), INTENT(out) :: rnep, ress, resf

! .. Local variables

  INTEGER :: k, j, np, jold

!$taf store cs, resfv, ressv  = scale_tape, rec = scale
! ... rescaling and sums over sub-grid cells
  resf = 0.
  ress = 0.

!$taf loop = parallel
  np=1
  jold=gridp(1)  
  DO k = 1,vp
     j=gridp(k)
     
     resfv(:,:,k) = resfv(:,:,k) * lfac(k)
     ressv(:,:,k) = ressv(:,:,k) * cs(k)
     
     resf(:,:,j) = resf(:,:,j) + resfv(:,:,k) * frac(k)
     ress(:,:,j) = ress(:,:,j) + ressv(:,:,k) * frac(k)
     
     IF (j/=jold) np=1
     rnepp(:,:,j,np) = rnppv(:,:,k) - ressv(:,:,k) - resfv(:,:,k)
     np=np+1
     jold=j
  
  ENDDO

  rnep = rnpp - ress - resf

end SUBROUTINE budget
