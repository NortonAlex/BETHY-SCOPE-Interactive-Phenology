!*********************************************************
!*  SUBROUTINE io_deallocate
!*********************************************************
SUBROUTINE io_deallocate
  USE mo_constants
!HEW-ADD-050307:
  USE mo_namelist, ONLY : wiflg,vflg
  USE mo_io
  USE mo_pheno, only: mlai
  USE mo_hydro, only: mpasm

  IMPLICIT NONE

  DEALLOCATE (mtran, mpp, mrt, mtmp)

  Deallocate (mlai,mpasm)

  IF (vflg==1)  DEALLOCATE(mvp)
  IF (wiflg==1) DEALLOCATE(mwind)

END SUBROUTINE io_deallocate


