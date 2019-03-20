! compute net terrestrial carbon flux
!------------------------------------------------------------------ 

SUBROUTINE bethy( nchk, dayint, ng, vp, nrun, outint, scale, flux, fapar)

  ! .. Use Statements ..
  USE mo_constants
  !HEW-ADD-050307:
  USE mo_namelist, ONLY : dtime, nspin, &
       & p1start, p1end, p2start, p2end
  USE mo_carparams
  USE mo_carvar
  USE mo_carbon
  USE mo_climate
  USE mo_diagnostics
  USE mo_grid, ONLY : gridp
  USE mo_io
  USE mo_surface
  USE mo_vegetation
  USE mo_beta
  USE mo_taf
  USE mo_pheno
  USE mo_hydro
  USE mo_calendar
  USE mo_prog
  USE mo_helper, ONLY : minx, maxx
  !ANorton-ADD-150107. Fluorescence
  USE fluo
  
  IMPLICIT NONE

  INTEGER, INTENT(in)  :: nchk, dayint, ng ,vp, nrun, outint, scale
  REAL, DIMENSION(nrun,outint,ng), INTENT(out)  :: fapar, flux


! pjr added extra declaration
  REAL                  :: fluxdiff ! difference in NEP between two periods 

! .. Local Variables ..
  INTEGER             :: oday, iday, rday, dstep, aday
  INTEGER             :: ts, day1p, k, t, yr, mo
  INTEGER :: iflgph, i, j, a, b, c
  REAL :: fx
  INTEGER ::  mday, nyears, its, iday0, iday1, ryear0, adayint
  REAL, DIMENSION(ng) :: area
  REAL, DIMENSION(nrun,12) :: diagh1,diagh2,diagh3,diagh4,diagh5
  REAL,  DIMENSION(vp) :: sfac
  REAL, DIMENSION(vp)  :: lfac

  INTEGER :: jl, jj !FastOpt: new loop counter from diagnostics now in model
!  print *, 'entered the bethy code' 
!------------------------------------------------------------------
! some initializations
! Note, remove inho dependency
!------------------------------------------------------------------

!Output fluo and gppfluo for every day run, not just per rmonth. 
  OPEN(unit=93,file='fluoro.bin',form='unformatted')
  OPEN(unit=94,file='gppfluo.bin',form='unformatted')
  OPEN(unit=95,file='gridp.bin',form='unformatted')
  OPEN(unit=96,file='vg_nv.bin',form='unformatted')
  OPEN(unit=97,file='t_.bin',form='unformatted')
  jdlast = 0
  ragain = .TRUE.

  diagh1 = 0.
  diagh2 = 0.

! init phenology and hydrology 
  CALL inithydro
  CALL initpheno
 
  rbe    = 0.
  fe     = 1. 
  cs     = 0.
  ptv    = 0.
  pgs    = 0.
  pci    = 0.
  ztrans = 0.
  zptrans= 0.
  zpcevp  = 0.
  zpsevp  = 0.
  zrhos  = 0.
  !vm(:) = MAX( vm(:), 1.e-20 )
  rnpp      = 0.
  rnppv      = 0.
  rnep      = 0.
  ress      = 0.
  resf      = 0.
  ressv      = 0.
  resfv      = 0.
  dnpp =0.
  fapar = 0. ! dummy setting
  zgc = 0.
  dsevp = 0.

print *, 'nchk =', nchk
print *, 'tdays =', tdays

!$taf store pasmmax = scale_tape, rec = scale
!------------------------------------------------------------------
! .. outer daily time loop (for checkpointing)
!------------------------------------------------------------------
  DO oday = 1,nchk
     dstep=INT(tdays/nchk)
     IF (MOD(tdays,nchk).GT.0) THEN
        dstep = dstep + 1
     ENDIF
     IF (dstep*oday.GT.tdays) THEN
        dstep = tdays - dstep * (oday-1) 
     ENDIF

!------------------------------------------------------------------
! .. inner daily time loop (for checkpointing)
!------------------------------------------------------------------
     DO iday = 1, dstep
        keyday = iday+ (oday-1) * (dstep) + (scale-1) * (nchk*dstep)
!        print *, 'keyday = ', keyday
!FastOpt !$taf store dnpp,dpcevp,dpsevp,dtrp,esum,lai,laihi,lintw = day_tape, rec = keyday 
!$taf store esum,lai = day_tape, rec = keyday 
!FastOpt !$taf store pasm,psoilst,rbe,rhosn,snow,snowh,zfc,zlai = day_tape, rec = keyday 
!$taf store psoilst,rbe,rhosn,snowh,zfc = day_tape, rec = keyday 
!FastOpt !$taf store daylen, tmpm = day_tape, rec = keyday 
!$taf store daylen = day_tape, rec = keyday 
        rday = iday + dstep*(oday-1)

!------------------------------------------------------------------
! .. set timings
!------------------------------------------------------------------
        ftspd  = REAL(tspd)
        zdiagt = dtime * 60 ! time step in seconds
        tspm   = tspd * ndays
        atspm  = tspd * jdlst
        tspy   = tspd * jdpyear
        nyears = nrun + nspin
        maxts  = nyears * tspd * jdpyear
        ts     = 0

        IF (spin(rday)==rday) THEN
           outyear=0
        ELSE
           outyear=ayear(rday)
        ENDIF
        ryear=ayear(rday)

        IF (spin(rday)==0) THEN
           iday0=rday-dspin
           iday1=rday-dspin+dayint-1
           aday=rday-dspin
        ELSE
           IF (dspin>sdays) THEN
              aday = MOD(rday,sdays)
              iday0 = MOD(rday,sdays)
              iday1=MOD(rday,sdays)+dayint-1
           ELSE
              aday = rday
              iday0 = rday
              iday1 = rday+dayint-1
           ENDIF
        ENDIF

        rmonth=amonth(rday)

        ! Use for monthly mean climate forcing (one diurnal cycle per month) for photosynthesis calcs
	! Calculates the first and last days of the month.
!        iday0=SUM(rdays(1:rmonth))-rdays(rmonth)+1
!        iday1=SUM(rdays(1:rmonth))
!        adayint=iday1-iday0+1

        IF (iday1>sdays) iday1=sdays

        adayint=dayint
        IF (idayint(rday)+adayint>tdays) adayint=tdays-idayint(rday)+1

   print*,"rday,iday0,iday1,sdays,aday"
   PRINT*,rday,iday0,iday1,sdays,aday

        ! Option to use monthly prescribed LAI.
        !  - lai is forced for simulated period, not spin-up.
        !IF (rday > sdays) lai = prescribed_lai(:,rmonth)

        IF (rday == idayint(rday)) THEN
           ryear0 = outyear
           IF (outyear<1) ryear0 = ryear
           print *, '..daycount = ', daycount(iday), iday
           print*,'   rday,aday,iday0,iday1'
           print*,rday,aday,iday0,iday1
           inho=13

           CALL diagdayreset

           DO k=1,vp
              j=gridp(k)
              fe(k) = maxx (psoilst(k), 0., 1e-3)
              fe(k) = minx (fe(k), 1., 1e-3)
              fx = maxx ((psoilst(k)-0.9)/(1-0.9), 0., 1e-3)
              fx = minx (fx, 1., 1e-3)
              zrhos(k) = fx*rhosw(j) + (1-fx)*rhosd(j)
           ENDDO

           CALL climsubday1 (ng, vp, fe, iday0, iday1)
           keydayint = daycount(iday)+ (oday-1) * (dstep) + (scale-1) * (nchk*dstep)
           print *,'    daycount(iday) = ', daycount(iday)
           print *,'        diurnal simulation with mean state from ',iday0,' to ',iday1  
!           print*, 'iday, keydayint = ',iday, keydayint
!$TAF store mu, tmp, pair, cloudf  = dayint_tape, rec = keydayint

!------------------------------------------------------------------
! diurnal timestep loop 
!------------------------------------------------------------------
           DO its = 1, tspd
              keydiurnal = its + (daycount(iday)-1) * (tspd) + (oday-1) * (tspd*dstep) + (scale-1) * (tspd*nchk*dstep)
!              print *, 'keydiurnal = ', keydiurnal, daycount(iday), iday, dayint
!              print *,'[keydiurnal, daycount(iday), iday, dayint] =', keydiurnal, daycount(iday), iday, dayint
              inho = MOD(its+11,24) + 1
              day1p = its-1
              ts=ts+1
              DO k = 1, vp
                 j = gridp (k)
                 IF (snowh(k)>0.) THEN
                    IF (mu(its,j)>0.) THEN
                       zrhos(k) = rhosn(k) + (1.-rhosn(k)) * rhosn(k)**3 * ( &
                            cloudf(j)**2 + (1. - 1.3*cloudf(j)**2) &
                            * EXP (1. - (1.-mu(its,j))**(-2)))
                    ELSE
                       zrhos(k) = rhosn(k)
                    ENDIF
                 ENDIF
              ENDDO
!              print *,'Inside diurnal timestep loop, inho = ', inho,' its = ', its
              CALL climsubday2 (ng, REAL(inho))
!              print *, 'swdown shape:', shape(swdown)
!              print *, 'SWDOWN ARRAY:', swdown
!$TAF store swdown  = diurnal_tape, rec = keydiurnal
              ! .. calculate stomatal conductance
              zlai = lai
!$TAF store zrhos = diurnal_tape, rec = keydiurnal
!FastOpt !$TAF store pardown, zgc, pgs, rbe = diurnal_tape, rec = keydiurnal
              CALL photo1 (ng,vp,swdown,tmp(inho,:),pair, &
                   & zrhos,pardown,fdirpar,zlai,zgc,pgs,zfpar,zfc,zassc,zraut,inho, &
                   & c4flg,ph,class,vm,jmf,zrphc,fautleaf,ccost, &
                   & EC,EO,EV,ER,EK,tgam,alpha,alc4,kc0,ko0) 
!$TAF store pgs,zfpar,eamin  = diurnal_tape, rec = keydiurnal
              ! .. calculate energy balance and evapotranspiration fluxes
              CALL surface (ng,vp,tmp(inho,:),cloudf,wspeed,pair,psoilst,lintw,zfpar, &
                   & zfc,zlai,swdown,zgc,pgs,ptv,eamin,fe,ztrans,zptrans,zpcevp, &
                   & zpsevp,inho,hv,day1p,zrhos,keydiurnal)     
!write (9,'(2I5,9E20.12)') iday, its, ztrans(429), zptrans(429)
!$TAF store ptv  = diurnal_tape, rec = keydiurnal
              ! .. calculate C assimilation
              CALL photo2 (ng,vp,swdown,ptv,pair, &
                   & zrhos,pardown,fdirpar,zlai,zgc,pgs,zfpar,zfc,zassc,zraut,inho, &
                   & c4flg,ph,class,vm,jmf,zrphc,fautleaf,ccost, &
                   & EC,EO,EV,ER,EK,tgam,alpha,alc4,kc0,ko0,zgrowth,zmaint)
              ! .. do diurnal diagnostics 
              ! option to give prescribed lai (from file) to fluorescence calculations
              IF ( inho == 13 ) THEN
              CALL fluorescence (ryear,rmonth,iday,inho,iday0,iday1,swdown,pardown,&
                                & tmp(inho,:),pair,eamin,ca,OX, & 
                                & zlai, &
                                & jmf,vms,EC,EO,EV,ER,EK,kc0s,ko0s,vomf,rdf,&
                                & rfluo,rgppfluo,PAR_scope,PAR_scope_cab,&
                                & rfluo_diurnal,rgppfluo_diurnal,&
                                & rlai_diurnal,rapar_diurnal,raparcab_diurnal,&
                                & rpar_diurnal,rswdown_diurnal)             
              zassc = zgppfluo               ! ANorton. To allow SCOPE-GPP to pass onto subsequent c-balance equations
!              print *,'SCOPE FLUO::', rfluo
!              print *,'SCOPE GPP::', rgppfluo
              ENDIF         ! for selected time of fluo computation
              CALL diagnostics (ng,vp,zassc,zraut,zgrowth,zmaint,ztrans,zptrans,zpcevp,zpsevp)

	   ENDDO ! end diurnal timestep loop 
           !WOK-CHG-100327: save daily transpiration after diurnal cycle
           zdtrp = dtrp
           !MAS-CHG-070719 : moved here to provide rnppv to cbalance
!$taf store dnpp = day_tape, rec = keyday 
           CALL mdiagnostics (ng,vp,adayint,rday,scale)
        ENDIF

        !WOK-CHG-070611 : call to hydrology moved to after daily loop with photosynthesis and energy balance
        ! new subroutine to calculate hydrology, currently only copies global smoisture field on smoist pft*grid vector
!$TAF store dpcevp,dpsevp = day_tape, rec = keyday 
        CALL hydrology (ng,vp,aday,keyday)
        !WOK-CHG-070611 : phenology moved down here and implemented as daily scheme
!FastOpt !$TAF store pasm,zlai,dptrp,lai,tmpm = day_tape, rec = keyday 
!$TAF store pasm,zlai,dptrp,tmpm = day_tape, rec = keyday 
        CALL phenology (aday,keyday)
        !mas-add-070724: save some daily fields (leafshed,psoilst) for cbalance spin-up
        IF (outyear>0)   CALL dailysave(ng,vp,aday)
        !WOK-CHG-070625 need to change 'cbalance' to daily calls, but this is not yet working
        !WOK-CHG-070625 pass daily NPP addition to monthly budget 'nppv' (now used only for spin-up)
        !WOK-CHG-070625 also pass time step as last argument
!FastOpt !$taf store zfpar = day_tape, rec = keyday 
        CALL ddiagnostics (nrun,outint,ng,vp,lmonth(aday)-fmonth(aday)+1,doy(rday),scale,fapar)
!        print *,'..end of inner daily loop'
!        RETURN
     ENDDO ! inner daily loop
!     print *,'.end of outer daily loop' 
  ENDDO ! outer daily loop
!print *,'Outside of outer daily loop'

  ! .. update carbon balance (inc. soil respiration)
!$TAF store dsleafshed,dspsoilst,rnppv = scale_tape, rec = scale
  CALL cbalance (ng, vp, nrun, scale, lfac, cs)



!$TAF store lfac = scale_tape, rec = scale
!  do i=1,nrun
!  do iday=1,outint
!  do k=1,vp
!  write (9,'(4I4,9E20.12)') 0,i,iday,k,resfv(i,iday,k)
!  enddo
!  enddo
!  enddo
  CALL budget(ng, vp, nrun, outint, scale, cs, lfac, rnpp, resfv, ressv, ress, resf, rnep, rnepp, rnppv)
!  do i=1,nrun
!  do iday=1,outint
!  do k=1,vp
!  write (9,'(4I4,9E20.12)') 1,i,iday,k,resfv(i,iday,k)
!  enddo
!  enddo
!  enddo
!  do i=1,nrun
!  do iday=1,outint
!  do k=1,ng
!  write (9,'(3I4,9E20.12)') i,iday,k,rnep(i,iday,k),rnpp(i,iday,k),ress(i,iday,k),resf(i,iday,k)
!  enddo
!  enddo
!  enddo
  !------------------------------------------------------------------
  ! .. write output  
  !------------------------------------------------------------------
  CALL diagout (ng,vp,scale,outint)
  IF (scale == 1) THEN
     flux = rnep(1:,:,:)/1000.
     prog_global(1,:,:,:) = rnep(1:,:,:)/1000.
     prog_global(2,:,:,:) = rnpp(1:,:,:)/1000.
  ELSE
     flux = rnep(1:,:,:)
     prog_sites(:,:,:) = rnpp(1:,:,:)
  ENDIF

CLOSE(93)
CLOSE(94)
CLOSE(95)
CLOSE(96)
CLOSE(97)

END SUBROUTINE bethy
