!*********************************************************
!*  SUBROUTINE diagnostics
!*  organizes diagnostic ouput
!*********************************************************

SUBROUTINE diagnostics (ng,vp,zassc,zraut,zgrowth,zmaint,ztrans,zptrans,zpcevp,zpsevp,dapar)

! .. Use Statements ..
  USE mo_constants
  USE mo_diagnostics
!WOK-ADD-070625
!WOK-CHG-070711 dptr -> dptrp
  USE mo_hydro, ONLY: dtrp, dptrp, dpcevp, dpsevp, dievp
  USE mo_climate, ONLY: zrhos
  USE mo_grid, ONLY : gridp

  IMPLICIT NONE

! .. Arguments
  INTEGER, INTENT(in) :: ng, vp
  REAL, DIMENSION(vp), INTENT(in) :: zassc,zraut,zgrowth,zmaint
  REAL, DIMENSION(vp), INTENT(in) :: ztrans,zptrans,zpcevp,zpsevp
  REAL, DIMENSION(vp), INTENT(in) :: dapar

! .. Local Scalars ..
   INTEGER :: j, k, np, jold
   


! units: kg C/m2 per month for rgpp and rnpp
  
!WOK-ADD-070625 sum up daily carbon fluxes
!WOK-ADD-070626 potential soil evaporation from 'zpevp', renamed 'dpet' into 'dptrp'

  np=1
  jold=gridp(1)  


!FastOpt working around a TAF bug
  k = 2*vp
!taf store k = top_tape, rec = 1
  IF (k.GT.0) THEN
     DO k= 1,vp
        j=gridp(k)

        dnpp(k) = dnpp(k) + (zassc(k)-zraut(k)) * zdiagt
! correct transpiration for times of wet canopy
!        dtrp(k) = dtrp(k) + ztrans(k) * (1.-dievp(k) / max (dpcevp(k), 1e-6)) * zdiagt 
        dtrp(k) = dtrp(k) + ztrans(k) * zdiagt 
        dpcevp(k) = dpcevp(k) + zpcevp(k) * zdiagt
        dpsevp(k) = dpsevp(k) + zpsevp(k) * zdiagt
        dptrp(k) = dptrp(k) + zptrans(k) * zdiagt



        dgpp(k) = dgpp(k) + zassc(k) * zdiagt
        drgrw(k) = drgrw(k) + zgrowth(k) * zdiagt
        drmnt(k) = drmnt(k) + zmaint(k) * zdiagt

        ! WOK-CHG-071019 'rrhos' moved here
        drhos(k) = drhos(k) + zrhos(k) * zdiagt / (24. * 3600.)

        nppv(outyear,rmonth,k) = nppv(outyear,rmonth,k) + &
             & (zassc(k)-zraut(k))*zdiagt * rdays(rmonth)
        IF (j/=jold) np=1
        nppp(outyear,rmonth,j,np) = nppp(outyear,rmonth,j,np) + &
             & (zassc(k)-zraut(k))*zdiagt * rdays(rmonth)
        np=np+1
        jold=j


     ENDDO
  ENDIF



END SUBROUTINE diagnostics


!WOK-ADD-070622 new subroutine
!*********************************************************
!*  SUBROUTINE mdiagnostics
!*  records monthly fluxes of carbon and water
!*  from the diurnal cycle calculations representative
!*  of a period of 'dayint' days
!*********************************************************
!WOK-CHG-071031 renamed variables for actual and potential evapotranspiration

SUBROUTINE mdiagnostics (ng,vp,dayint,rday, scale)

! .. Use Statements ..
  USE mo_constants
  USE mo_diagnostics
  USE mo_grid, ONLY : gridp, frac
  USE mo_hydro, ONLY: dtrp, dsevp, dievp, dptrp, dpcevp
  USE mo_calendar
  USE mo_climate

  IMPLICIT NONE
  
  ! .. Arguments ..
  INTEGER, INTENT(in) :: ng, vp, dayint, rday, scale

  ! .. Local Variables ..
  INTEGER :: j, k,days,dnm,dim
  INTEGER :: oy, om
  INTEGER :: daysrmonth, daysom, dayhelp
  INTEGER :: jold, np

! number of days in current month
  daysrmonth = lmonth(rday)-fmonth(rday)+1
  days=dayint

  days=dayint
  IF (scale == 1) THEN

     IF ( (rday+dayint-1) <= lmonth(rday)) THEN

        ! units: kg C/m2 per month for rgpp and rnpp
        DO k = 1, vp
           j = gridp(k)
           rnppv(outyear,rmonth,k) = rnppv(outyear,rmonth,k) + dnpp(k) * days
           rnpp(outyear,rmonth,j) = rnpp(outyear,rmonth,j) + dnpp(k) * frac(k) * days

           rgpp(outyear,rmonth,j) = rgpp(outyear,rmonth,j) + dgpp(k) * frac(k) * days
           rgrw(outyear,rmonth,j) = rgrw(outyear,rmonth,j) + drgrw(k) * frac(k) * days
           rmnt(outyear,rmonth,j) = rmnt(outyear,rmonth,j) + drmnt(k) * frac(k) * days
           raet(outyear,rmonth,j) = raet(outyear,rmonth,j) + &
                & (dtrp(k)+dievp(k)+dsevp(k)) * frac(k) * days
           rpet(outyear,rmonth,j) = rpet(outyear,rmonth,j) + &
                & (dptrp(k) + dpcevp(k)) * frac(k) * days
           rrhos(outyear,rmonth,j) = rrhos(outyear,rmonth,j) + drhos(k) * frac(k) * days / daysrmonth

        ENDDO

     ELSE
        dim = lmonth(rday) - rday +1
        DO k = 1, vp
           j = gridp(k)
           rnppv(outyear,rmonth,k) = rnppv(outyear,rmonth,k) + dnpp(k) * dim
           rnpp(outyear,rmonth,j) = rnpp(outyear,rmonth,j) + dnpp(k) * frac(k) * dim

           rgpp(outyear,rmonth,j) = rgpp(outyear,rmonth,j) + dgpp(k) * frac(k) * dim
           rgrw(outyear,rmonth,j) = rgrw(outyear,rmonth,j) + drgrw(k) * frac(k) * dim
           rmnt(outyear,rmonth,j) = rmnt(outyear,rmonth,j) + drmnt(k) * frac(k) * dim
           raet(outyear,rmonth,j) = raet(outyear,rmonth,j) + &
                & (dtrp(k)+dievp(k)+dsevp(k)) * frac(k) * dim
           rpet(outyear,rmonth,j) = rpet(outyear,rmonth,j) + &
                & (dptrp(k) + dpcevp(k)) * frac(k) * dim
           rrhos(outyear,rmonth,j) = rrhos(outyear,rmonth,j) + drhos(k) * frac(k) * dim / daysrmonth

        ENDDO

        IF (rday+dayint > tdays) THEN
           dayhelp=tdays
        ELSE
           dayhelp=rday+dayint
        ENDIF

        IF (rday.NE.spin(rday)) THEN
           dnm = (dayhelp) - lmonth(rday) -1
           oy = ayear(dayhelp)
           om = amonth(dayhelp)
           daysom = lmonth(dayhelp)-fmonth(dayhelp)+1

           DO k = 1, vp
              j = gridp(k)
              rnppv(oy,om,k) = rnppv(oy,om,k) + dnpp(k) * dnm
              rnpp(oy,om,j) = rnpp(oy,om,j) + dnpp(k) * frac(k) * dnm

              rgpp(oy,om,j) = rgpp(oy,om,j) + dgpp(k) * frac(k) * dnm
              rgrw(oy,om,j) = rgrw(oy,om,j) + drgrw(k) * frac(k) * dnm
              rmnt(oy,om,j) = rmnt(oy,om,j) + drmnt(k) * frac(k) * dnm
              raet(oy,om,j) = raet(oy,om,j) + &
                   & (dtrp(k)+dievp(k)+dsevp(k)) * frac(k) * dnm
              rpet(oy,om,j) = rpet(oy,om,j) + &
                   & (dptrp(k) + dpcevp(k)) * frac(k) * dnm
              rrhos(oy,om,j) = rrhos(oy,om,j) + drhos(k) * frac(k) * dnm / daysom

           ENDDO

        ENDIF
     ENDIF
  ELSEIF (scale == 2) THEN

! Assuming that runs on site scale are always simulating every day of the year
     np=1
     jold=gridp(1)  
!$taf loop = dependent
     DO k = 1, vp
        j = gridp(k)
        rnpp(outyear,doy(rday),j) = rnpp(outyear,doy(rday),j) + dnpp(k) * frac(k)
        rnppv(outyear,doy(rday),k) = dnpp(k)


        IF (j/=jold) np=1
        rnppp(outyear,rmonth,j,np) =  rnppp(outyear,rmonth,j,np) + dnpp(k)*days 
        rgppp(outyear,rmonth,j,np) = rgppp(outyear,rmonth,j,np) + dgpp(k)*days 
        np=np+1
        jold=j
        rgpp(outyear,doy(rday),j) = dgpp(k) * frac(k) 
        rgrw(outyear,rmonth,j) = rgrw(outyear,rmonth,j) + drgrw(k) * frac(k) * days
        rmnt(outyear,rmonth,j) = rmnt(outyear,rmonth,j) + drmnt(k) * frac(k) * days
        raet(outyear,rmonth,j) = raet(outyear,rmonth,j) + &
             & (dtrp(k)+dievp(k)+dsevp(k)) * frac(k) * days
        rpet(outyear,rmonth,j) = rpet(outyear,rmonth,j) + &
             & (dptrp(k) + dpcevp(k)) * frac(k) * days
        rrhos(outyear,rmonth,j) = rrhos(outyear,rmonth,j) + drhos(k) * frac(k) * days / daysrmonth

     ENDDO
     
  ENDIF
  
END SUBROUTINE mdiagnostics


!WOK-ADD-070725 new subroutine
!*********************************************************
!*  SUBROUTINE ddiagnostics
!*  records monthly averages and sums of status variables
!*  from hydrology and phenology
!*********************************************************

SUBROUTINE ddiagnostics (nrun, outint, ng,vp,daysinmonth, aday, scale, fapar)

! .. Use Statements ..
  USE mo_constants
  USE mo_diagnostics
  USE mo_grid, ONLY: gridp, frac
!  USE mo_climate, ONLY: zrhos
  USE mo_hydro, ONLY: pasm, runoff, dsevp, dsnmelt
  USE mo_vegetation, ONLY: lai, zfpar, dapar

  IMPLICIT NONE
  
  ! .. Arguments ..
  INTEGER, INTENT(in) :: nrun, outint, ng, vp, daysinmonth, scale, aday
  REAL, DIMENSION(nrun,outint,ng), INTENT(inout)  :: fapar

  ! .. Local Variables ..
  INTEGER :: j, k, np, jold

  np=1
  jold=gridp(1)  
  DO k = 1, vp
    j = gridp(k)

    IF (outyear > 0) THEN
       IF (j/=jold) np=1
       IF (scale == 1) THEN

          rfparp(outyear,rmonth,j,np) = rfparp(outyear,rmonth,j,np) + zfpar(k) / daysinmonth

          fapar(outyear,rmonth,j) = fapar(outyear,rmonth,j)+ zfpar(k) * frac(k) / daysinmonth
       ELSE

          rfparp(outyear,aday,j,np) = rfparp(outyear,aday,j,np) + zfpar(k)

          fapar(outyear,aday,j) = fapar(outyear,aday,j)+ zfpar(k) * frac(k)
       ENDIF
       np=np+1
       jold=j
    ENDIF

    rlai(outyear,rmonth,j) = rlai(outyear,rmonth,j) + lai(k) * frac(k) / daysinmonth


    rpasm(outyear,rmonth,j) = rpasm(outyear,rmonth,j) + pasm(k) * frac(k) / daysinmonth
    rrunoff(outyear,rmonth,j) = rrunoff(outyear,rmonth,j) + runoff(k) * frac(k)
    rsevp(outyear,rmonth,j) = rsevp(outyear,rmonth,j) + dsevp(k) * frac(k)
    rsnmelt(outyear,rmonth,j) = rsnmelt(outyear,rmonth,j) + dsnmelt(k) * frac(k)
    rfpar(outyear,rmonth,j) = rfpar(outyear,rmonth,j) + zfpar(k) * frac(k) / daysinmonth
    rdapar(outyear,rmonth,j) = rdapar(outyear,rmonth,j) + dapar(k) * frac(k) / daysinmonth

  ENDDO

  
END SUBROUTINE ddiagnostics


!MAS-ADD-070725 new subroutine
!*********************************************************
!*  SUBROUTINE dailysave
!*  records various daily variable
!*********************************************************

SUBROUTINE dailysave (ng,vp,day)

! .. Use Statements ..
  USE mo_constants
  USE mo_diagnostics
  USE mo_grid, ONLY: gridp
  USE mo_pheno, ONLY: leafshed
  USE mo_hydro, ONLY: psoilst
  USE mo_namelist, ONLY: nspin

  IMPLICIT NONE
  
  ! .. Arguments ..
  INTEGER, INTENT(in) :: ng,vp,day

  ! .. Local Variables ..
!  REAL :: div
  INTEGER :: j, k


  DO k = 1, vp
    j = gridp(k)

    dsleafshed(day,k) = leafshed(k)
    dspsoilst(day,k) = psoilst(k)
    
  ENDDO
  
END SUBROUTINE dailysave


!*********************************************************
!*  SUBROUTINE diagdayreset
!*  resets daily totals of carbon fluxes
!*********************************************************
  
SUBROUTINE diagdayreset
  
  USE mo_diagnostics, ONLY: dgpp, dnpp, drgrw, drmnt, drhos
  USE mo_hydro, ONLY: dtrp, dptrp, dpcevp, dpsevp

  IMPLICIT NONE

  dnpp = 0.
  dtrp = 0.
  dpcevp = 0.
  dpsevp = 0.
  dptrp = 0.


  dgpp = 0.
  drgrw = 0.
  drmnt = 0.
  drhos = 0.


END SUBROUTINE diagdayreset



