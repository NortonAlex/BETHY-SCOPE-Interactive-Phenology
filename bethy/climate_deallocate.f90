!*********************************************************
!*  SUBROUTINE climate_deallocate
!*********************************************************
SUBROUTINE climateold_deallocate
  USE mo_constants
!HEW-ADD-050307:
  USE mo_namelist, ONLY : intime,vflg,wiflg
  USE mo_climate

  IF (intime==1) THEN
     DEALLOCATE (dtmp,dtran,dpp,drt)
     IF (wiflg==1) THEN
        DEALLOCATE (dwind)
     ENDIF
     IF (vflg==1) THEN
        DEALLOCATE (dvp)
     ENDIF
  ELSE
     DEALLOCATE (dtmp,dtran,dpp,drt)
     IF (wiflg==1) THEN
        DEALLOCATE (dwind)
     ENDIF
     IF (vflg==1) THEN
        DEALLOCATE (dvp)
     ENDIF
  ENDIF

  DEALLOCATE (htmp,zrhos)

END SUBROUTINE climateold_deallocate

