!*********************************************************
!*  SUBROUTINE climate_allocate
!*********************************************************
SUBROUTINE climateold_allocate
  USE mo_constants
!HEW-ADD-050307:
  USE mo_namelist, ONLY : intime,vflg,wiflg
  USE mo_climate
 

! .. allocate daily arrays
  IF (intime==1) THEN 
     ALLOCATE (dtmp(31,ng),dtran(31,ng),dpp(31,ng), drt(31,ng))
     IF (wiflg==1) THEN
        ALLOCATE (dwind(31,ng))
     ENDIF
     IF (vflg==1) THEN
        ALLOCATE (dvp(31,ng))
     ENDIF
  ELSE
     ALLOCATE(dtmp(jdpyear,ng),dtran(jdpyear,ng),dpp(jdpyear,ng),drt(jdpyear,ng))
     IF (wiflg==1) THEN
        ALLOCATE (dwind(jdpyear,ng))
     ENDIF
     IF (vflg==1) THEN
        ALLOCATE (dvp(jdpyear,ng))
     ENDIF
  ENDIF

! .. allocate subdaily arrays
  ALLOCATE (htmp(ng,tspd))

! .. allocate daily pft specific arrays
  ALLOCATE (zrhos(vp))

END SUBROUTINE climateold_allocate
