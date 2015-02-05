!*********************************************************
!*  SUBROUTINE io_allocate
!*********************************************************
SUBROUTINE io_allocate
  USE mo_constants
!HEW-ADD-050307:
  USE mo_namelist, ONLY : nrun,wiflg,vflg
  USE mo_io
  USE mo_pheno, only: mlai
  USE mo_hydro, only: mpasm

  IMPLICIT NONE
  INTEGER :: nm

  nm=(nrun+1)*12 ! number of simulated months, spin-up climatology only once

  ALLOCATE (mtran(nm,ng), mpp(nm,ng),mrt(nm,ng), mtmp(nm,ng))

  allocate ( mpasm(nm,ng,nv),mlai(nm,ng,nv))

  IF (vflg==1)  ALLOCATE(mvp(nm,ng))
  IF (wiflg==1) ALLOCATE(mwind(nm,ng))

END SUBROUTINE io_allocate
