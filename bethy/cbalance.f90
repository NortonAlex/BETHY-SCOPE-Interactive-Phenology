!*********************************************************
!*  SUBROUTINE cbalance
!*  calculates heterotrophic respiration
! 05/05: Txk revises handling negative pool sizes etc
! 11/06: MAS/Txk avoid division dellai = 0.
!*********************************************************
SUBROUTINE cbalance (ng, vp, nrs, scale, lfac, cs)

! .. Use statements
  USE mo_constants
  USE mo_namelist
  USE mo_grid, ONLY: gridp, frac
  USE mo_vegetation, ONLY: aw, q10f, q10s, fracs, tauf, pft
  USE mo_beta, ONLY: beta
  USE mo_diagnostics, ONLY: ressv, resfv, rnppv, dsleafshed, dspsoilst
  USE mo_calendar
  USE mo_climate, ONLY: dtmin, dtmax

  IMPLICIT NONE

! .. Arguments
  INTEGER, INTENT(in) :: ng, vp, nrs, scale
  REAL, DIMENSION(vp), intent(out) :: lfac, cs

! .. Local variables

  INTEGER ::  k, evflg,j,cspin, d, yr ,mn
  REAL :: taul, fraca, npp, dellai
  REAL :: litterprod, litterloss, cfold, alpha
  REAL :: fw, kf, ks, c, t
  REAL :: minsw
  REAL, DIMENSION(vp) :: cf, nppsum, ressumf, ressums, litterprodsum
  REAL, DIMENSION(ng) :: dsmtmp

 ! .. Parameters
!  REAL, PARAMETER :: eps = 1.e-6 ! a small number
  REAL, PARAMETER :: eps = 1.e-12 ! a small number
  REAL, PARAMETER :: fphox = 1e-3 ! fractional rate of photo-oxidation of litter

!WOK-070723 changed jl->k and jj->j for better legibility and compatibility with mo_pheno & mo_hydro
!           passing 'leafshed' instead of 'lai' (see 1)
!WOK-070723 the following problems with 'cbalance' were identified:
! (1) the leaf shedding rate should be calculated directly by phenology (done)
! (2) litter production is here equated with leaf shedding rate, which is not true,
!     because there is also root litter production and other dead plant material.
!     in the current implementation, respiration from the fast pool
!     can become extremely small in dry areas (alpha<<1). This does not apply to the slow pool,
!     because this pool is scaled so that RES(slow) = NPP - RES(fast) on time average.
!     However: because the fast soil pool is spun up, low respiration (e.g. alpha<<1)
!     will lead to large pool size over a long-enough time. In this way, if
!     litterpod = const. * leafshed, then this error will be compensated by in creasing
!     'cl' by a factor of 'const.'. However, this method may be extremely inefficient
!     at some place where decomposition is extremely slow.
!
!     Possible solution: sum up 'litterprod' over simulation period and require it to be
!     equal to NPP. In this way, leaf shedding is used only to drive the relative
!     temporal changes in the leaf production rate, not the absolute rate; this is probably
!     a reasonable first-order assumption.

  taul = 2.
  minsw=0.0001
  yr = 0
  mn = 0
  litterprodsum = 0.
  cf = 0.
  nppsum  = 0.
  ressums = 0.
  ressumf = 0.





!FastOpt !FastOpt !$taf init cbalance_spin_tape = static,  maxdpy*(aspin+maxyears)*2 
!$taf init cbalance_spin_tape = static,  maxdpy*(aspin+maxyears)
!FastOpt !FastOpt !$taf init cbalance_tape = static,  maxdpy*(maxyears)*2 
!$taf init cbalance_tape = static,  maxdpy*(maxyears)


  DO cspin = aspin/nrs+1, 1, -1 

     DO d = 1, sdays ! sdays is the number of simulation days 

        dsmtmp(:)=(dtmin(:,d)+dtmax(:,d))/2.

!FastOpt !$taf store cf = cbalance_spin_tape, rec = d + (cspin - 1) * sdays + (scale-1) * sdays * (aspin/nrs+1)
!$taf store cf = cbalance_spin_tape, rec = d + (cspin - 1) * sdays 
!$taf loop = parallel
        DO k = 1,vp
           j=gridp(k)
!           alpha = dspsoilst(d,k)+minsw
           alpha = (1-2*minsw)*dspsoilst(d,k)+minsw
           IF (aw(k).GT.0.) THEN
              fw = fphox + (1.-fphox)*alpha**aw(k)
           ELSE
              fw = 1.
           ENDIF
           t = dsmtmp(j)
           kf = fw * q10f(k) ** (t/10.) * tauf(k) / 365
           litterprod = dsleafshed(d,k)
           litterloss = litterprod - (litterprod/kf - cf(k))*(1-exp(-kf))
           ! update state variable 'litter carbon' cf
           cf(k) = litterprod/kf+(cf(k)-litterprod/kf)*exp(-kf)
        ENDDO
     ENDDO
  ENDDO  ! end of spin up

  DO d = 1, sdays ! sdays is the number of simulation days 
     
     dsmtmp(:)=(dtmin(:,d)+dtmax(:,d))/2. 
     
!FastOpt !$taf store cf = cbalance_tape, rec = d  + (scale-1) * sdays 
!$taf store cf = cbalance_tape, rec = d 
!$taf store mn,yr = cbalance_tape, rec = d 
!$taf loop = parallel
     DO k = 1,vp
        j=gridp(k)
!        alpha = dspsoilst(d,k)+minsw
        alpha = (1-2*minsw)*dspsoilst(d,k)+minsw
        IF (aw(k).GT.0.) THEN
           fw = fphox + (1.-fphox)*alpha**aw(k)
        ELSE
           fw = 1.
        ENDIF
        t = dsmtmp(j)
        kf = fw * q10f(k) ** (t/10.) * tauf(k) / 365
        litterprod = dsleafshed(d,k)
        litterloss = litterprod - (litterprod/kf - cf(k))*(1-exp(-kf))
       ! update state variable 'litter carbon' cf
        cf(k) = litterprod/kf+(cf(k)-litterprod/kf)*exp(-kf)

        ks = fw * q10s(k) ** (t/10.)
        if (k==1) then 
           if (fyear(d+dspin)==d+dspin) yr=yr+1
           if (fmonth(d+dspin)==d+dspin) mn=mn+1
           if (mn==13) mn = 1
        endif
        if (scale == 1) then
           resfv(yr,mn,k) = resfv(yr,mn,k) + litterloss * (1.-fracs(k))
           ressv(yr,mn,k) = ressv(yr,mn,k) + ks
        elseif (scale == 2) then
           resfv(yr,doy(d+dspin),k) = litterloss * (1.-fracs(k))
           ressv(yr,doy(d+dspin),k) = ks
        endif
        litterprodsum(k) = litterprodsum(k) + litterprod              
     ENDDO
     
  ENDDO

  nppsum = SUM(SUM(rnppv(1:,:,:),dim=1),dim=1)

!$taf store litterprodsum,nppsum = scale_tape, rec = scale
  DO k=1,vp
     IF ((abs(nppsum(k))) < litterprodsum(k)/eps ) THEN
        lfac(k) = nppsum(k)/ litterprodsum(k)    
     ELSE
        lfac(k) = 1./eps
     ENDIF
  ENDDO

  ressumf = SUM(SUM(resfv(1:,:,:),dim=1),dim=1)
  ressums = SUM(SUM(ressv(1:,:,:),dim=1),dim=1)

! ... scalar for slow pool, which equals slow pool size divided by slow pool turnover time
!     (here conveniently called 'cs')
!$taf store ressums, ressumf = scale_tape, rec = scale
  DO k=1,vp 
     j=gridp(k)
     IF (ABS((nppsum(k)*beta(k) - ressumf(k)*lfac(k))) < ressums(k)/eps ) THEN
        cs(k) = (nppsum(k)*beta(k) - ressumf(k)*lfac(k)) / ressums(k)    
!     IF (ABS((nppsum(k)/beta(k) - ressumf(k)*lfac(k))) < ressums(k)/eps ) THEN
!        cs(k) = (nppsum(k)/beta(k) - ressumf(k)*lfac(k)) / ressums(k)    
     ELSE
!        cs(k) = 1. / eps
        cs(k) = sign (1. / eps, nppsum(k)*beta(k) - ressumf(k)*lfac(k))
     ENDIF

  ENDDO

END SUBROUTINE cbalance



