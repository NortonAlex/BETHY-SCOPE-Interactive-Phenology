MODULE mo_pheno
  
!***********************************************************
!* WOK, 2008-07-30
!* IMPLIFIED PHENOLOGY MODEL
!* simplified phenology model
!***********************************************************

  IMPLICIT NONE
  
!  REAL, ALLOCATABLE, DIMENSION (:,:,:) :: mlai ! monthly LAI fields from external data
  REAL, ALLOCATABLE, DIMENSION (:) :: tmpm  ! air-temperature memory [deg C]
  REAL, ALLOCATABLE, DIMENSION (:) :: laim  ! water limited LAI memory
  REAL, ALLOCATABLE, DIMENSION (:) :: laihi ! highest recorded LAI (with a decay rate, for setting 'zfc')
! WOK-ADD-070723 litter production to be calculated directly in phenology (not in cbalance indirectly)
  REAL, ALLOCATABLE, DIMENSION (:) :: tmpmmult, laimmult ! auxiliary fields
  REAL, ALLOCATABLE, DIMENSION (:) :: leafshed ! output field
  REAL :: laihimult
  REAL, PARAMETER :: taulaihi = 5.0           ! memory time for updating fractional cover
!  REAL, PARAMETER :: laimin = 1e-6            ! minimum LAI for pot. transpiration per LAI estimates
  REAL, PARAMETER :: eta = 0.99999            ! curvature parameter for mins/maxs
! WOK-ADD-070723 the list of controlling parameters
! FREE PARAMETERS
  REAL, ALLOCATABLE, DIMENSION (:) :: PLAIMAX ! maximum LAI
  REAL, ALLOCATABLE, DIMENSION (:) :: PTPHEN  ! leaf onset temperature [deg C]
  REAL, ALLOCATABLE, DIMENSION (:) :: PTPHENR ! range of leaf onset temperature [1/deg C]
  REAL, ALLOCATABLE, DIMENSION (:) :: PDPHEN  ! leaf shedding daylength [hours]
  REAL, ALLOCATABLE, DIMENSION (:) :: PDPHENR ! range of leaf shedding daylength [hours]
!  REAL, ALLOCATABLE, DIMENSION (:) :: PTSHD   ! leaf shedding temperature [deg C]
!  REAL, ALLOCATABLE, DIMENSION (:) :: PTSHDS  ! spread of leaf shedding temperature [1/deg C]
  REAL, ALLOCATABLE, DIMENSION (:) :: PLGR    ! leaf growth factor [1/days]
  REAL, ALLOCATABLE, DIMENSION (:) :: PKL     ! inverse leaf longevity from start of senescense [1/days]
  REAL, ALLOCATABLE, DIMENSION (:) :: PTAUW   ! target survival time at current soil moisture [days]
! PARAMETERS LEFT FIXED
  REAL, ALLOCATABLE, DIMENSION (:) :: PKS     ! inverse memory time for soil moisture-limited LAI [1/days]
  REAL, ALLOCATABLE, DIMENSION (:) :: PKM     ! inverse memory time for air temperature [1/days]

! WOK-090309 'ph' is now only used in nscale
! xph     1: warm-evergreen; 2: cold-evergreen; 3: summergreen; 4: raingreen; 5: grass; 6: annual crop;
  INTEGER, DIMENSION (0:13), PARAMETER :: xph= &
!PFT:0  1  2  3  4  5  6  7  8  9 10 11 12 13
   (/5, 1, 4, 1, 3, 2, 3, 1, 4, 5, 5, 2, 5, 6/)

!  LIST OF PFTs:
!  1:  tropical broadleaf evergreen tree 
!  2:  tropical broadleaf deciduous tree 
!  3:  temperate broadleaf evergreen tree 
!  4:  temperate broadleaf deciduous tree 
!  5:  evergreen coniferous tree 
!  6:  deciduous coniferous tree 
!  7:  evergreen shrub 
!  8:  deciduous shrub 
!  9:  C3 grass 
! 10:  C4 grass 
! 11:  tundra 
! 12:  swamp 
! 13:  arable crop 

CONTAINS
  
  !*********************************************************
  !*  SUBROUTINE initpheno
  !*  reads in LAI from file initializes and allocates memory
  !*********************************************************
  
  SUBROUTINE pheno_allocate(vp)
    
    ! .. Use Statements ..
    USE mo_constants


    IMPLICIT NONE

    INTEGER, INTENT(in) :: vp

!.. auxiliary fields and status variables
    ALLOCATE (tmpm(vp), laim(vp), laihi(vp))
    ALLOCATE (tmpmmult(vp), laimmult(vp))
!.. phenology parameters
    ALLOCATE (plaimax(vp), ptphen(vp), ptphenr(vp), pdphen(vp), pdphenr(vp))
    ALLOCATE (plgr(vp), pkl(vp), ptauw(vp), pks(vp), pkm(vp))
!.. output field
    ALLOCATE (leafshed(vp))

  END SUBROUTINE pheno_allocate


  SUBROUTINE initpheno

    USE mo_grid, ONLY: vp, gridp
    USE mo_vegetation, ONLY : zfc, lai !, ph
    USE mo_climate, ONLY :  dtmin, dtmax
    USE mo_carparams, ONLY : fcmax0
    USE mo_calendar

! .. locals
    INTEGER :: j, k, iday

    pkm = 1. / 30.
    pks = 1. / 30.
!    pks = 1. / 90.
!    write (9,*)'plaimax:  ',plaimax
!    write (9,*)'ptphen:   ',ptphen
!    write (9,*)'ptphenr:  ',ptphenr
!    write (9,*)'pdphen:   ',pdphen
!    write (9,*)'pdphenr:  ',pdphenr
!    write (9,*)'plgr:     ',plgr
!    write (9,*)'pkl:      ',pkl
!    write (9,*)'ptauw:    ',ptauw
!    write (9,*)'pkm:      ',pkm
!    write (9,*)'pks:      ',pks

!   multiplier for advancing temperature memory by one day
    tmpmmult = exp (-pkm)
!   multiplier for advancing soil-water limited LAI memory by one day
    laimmult = exp (-pks)
!   decay multiplier for evergreen LAI
!    laimult = exp (-pkl)
!   the air-temperature memory
    tmpm = 0.
!   the water stress index memory
    laim = 0.
!   decay multiplier for maximum LAI used to set fractional cover
    laihimult = exp (-1./(taulaihi*365.))
!   initialize LAI
    lai = plaimax
!   control for fractional cover 'zfc'    
    laihi = 0.
    zfc = fcmax0
!   output field
    leafshed = 0.
!   spin-up of temperature memory
    DO iday=1,sdays
       DO k = 1, vp
          j = gridp (k)
          tmpm(k) = (dtmin(j,iday)+dtmax(j,iday))/2. * (1. - tmpmmult(k)) + tmpm(k) * tmpmmult(k)
       ENDDO
    ENDDO
    
  END SUBROUTINE initpheno


  SUBROUTINE phenology (iday,keyday)

!------------------------------------------------------------------
! advances LAI and fractional cover by one day
! from its current state to the state at day 'iday'
!------------------------------------------------------------------

! .. use statements
    USE mo_constants
    USE mo_grid, ONLY: gridp, vp
    USE mo_climate, ONLY :  dtmin, dtmax, daylen
    USE mo_hydro, ONLY :  pasm, dptrp
    USE mo_helper, ONLY : errf, mins, maxs, minx, maxx
    USE mo_vegetation, ONLY : zfc, lai, zlai, sla, pft
    USE mo_carparams, ONLY : fcmax0, lailim0, cdrm
  
    IMPLICIT NONE

! .. arguments
    INTEGER, INTENT(in) :: iday, keyday
    
! .. locals
    INTEGER :: j, k
    REAL    :: xdtmp, lait, laiw, lailast, fx, t0, ts, ft, fd, fg
    REAL    :: laimaxw, laimax, r, lailim, wai

!FastOpt !$TAF store tmpm, lai, laim  = day_tape, rec = keyday
!$TAF store lai, laim  = day_tape, rec = keyday
!$taf loop = parallel
    DO k = 1, vp
      j = gridp (k)

      lailast = lai(k)

!      IF (ph(k)==1.or.ph(k)==4) THEN ! warm-evergreen and warm-deciduous phenology
      IF (pft(k)==1.or.pft(k)==2.or.pft(k)==3.or.pft(k)==7) THEN ! warm-evergreen and warm-deciduous phenology
        ! effective maximum LAI, taking into account structural limiations
!        laimax = plaimax(k) * (1. - exp(-laimaxw/plaimax(k)))
        laimax = pasm(k) * zlai(k) / ptauw(k) / maxx (dptrp(k), 1e-3, 2e-2) 
        laimax = mins (laimax, plaimax(k), 0.9)
        ! update water limited LAI memory
        laim(k) = laimax * (1. - laimmult(k)) + laim(k) * laimmult(k)
        ! rate of change of LAI towards limit
        r = plgr(k)
        ! limit LAI
        lailim = laim(k)
        ! update LAI
!        lai(k) = lailim - (lailim - lai(k)) * exp (-r)

!      ELSE IF (ph(k)==2.or.ph(k)==3) THEN ! cold-evergreen and cold-deciduous phenology
      ELSE IF (pft(k)==4.or.pft(k)==5.or.pft(k)==6.or.pft(k)==8.or.pft(k)==11) THEN ! cold-evergreen and cold-deciduous phenology
        ! update memory of daily mean temperature
        tmpm(k) = (dtmin(j,iday)+dtmax(j,iday))/2. * (1. - tmpmmult(k)) + tmpm(k) * tmpmmult(k)
        ! fraction of vegetation above temperature threshold
        ft = errf((tmpm(k)-ptphen(k))/ptphenr(k))
        ! fraction of vegetation above daylength threshold
        fd = errf((daylen(j)-pdphen(k))/pdphenr(k))
        r = ft * fd * plgr(k) + (1. - ft * fd) * pkl(k) + 1e-9
        lailim = maxx (ft * fd * plgr(k) * plaimax(k) / r, 1e-9, 5e-3)
!        lai(k) = plaimax(k) - (plaimax(k) - lai(k)) * exp (-r)
!        lai(k) = lailim - (lailim - lai(k)) * exp (-r)
            
      ELSE ! grass and annual crop phenology
        tmpm(k) = (dtmin(j,iday)+dtmax(j,iday))/2. * (1. - tmpmmult(k)) + tmpm(k) * tmpmmult(k)
        ft = errf((tmpm(k)-ptphen(k))/ptphenr(k))
        laimax = pasm(k) * zlai(k) / ptauw(k) / maxx (dptrp(k), 1e-3, 2e-2) 
        laimax = mins (laimax, plaimax(k), 0.9)
        laim(k) = laimax * (1. - laimmult(k)) + laim(k) * laimmult(k)
        r = ft * plgr(k) + (1. - ft) * pkl(k) + 1e-9
        lailim = maxx (ft * plgr(k) * laim(k) / r, 1e-9, 5e-3)
!        lai(k) = lailim - (lailim - lai(k)) * exp (-r)
        
      ENDIF

!      leafshed(k) = maxx (lailast - lai(k), 0., 1e-3) / sla(k) * 1000. * cdrm
      leafshed(k) = maxx ((lailim-lai(k))*(1.-exp(-r)), 0., 1e-3) / sla(k) * 1000. * cdrm
      lai(k) = lailim - (lailim - lai(k)) * exp (-r)

    ENDDO

!$TAF store lai,laihi  = day_tape, rec = keyday 
!$taf loop = parallel
    DO k = 1, vp
      ! set fractional cover
      laihi(k) = maxs (lai(k), laihi(k), eta) * laihimult
      zfc(k) = maxs (laihi(k) / lailim0, lai(k) / lailim0, eta)
      zfc(k) = mins ( zfc(k), 1., eta) * fcmax0
    ENDDO

  END SUBROUTINE phenology

  SUBROUTINE pheno_deallocate(vp)
  
    IMPLICIT NONE
  
    INTEGER, INTENT(in) :: vp

!$taf next required  = vp
    DEALLOCATE (tmpm, laim, laihi)
!$taf next required  = vp
    DEALLOCATE  (tmpmmult, laimmult)
!$taf next required  = vp
    DEALLOCATE (plaimax, ptphen, ptphenr, pdphen, pdphenr)
!$taf next required  = vp
    DEALLOCATE (plgr, pkl, ptauw, pks, pkm)
!$taf next required  = vp
    DEALLOCATE (leafshed)

  END SUBROUTINE pheno_deallocate


END MODULE mo_pheno
