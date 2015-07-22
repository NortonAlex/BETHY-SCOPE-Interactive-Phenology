MODULE mo_hydro
  
!***********************************************************
!* NECESSARY CHANGES TO SUBROUTINE surface
!* WOK, 2007-06-12
!* (1) existing 'ztrans' as additional output field, containing transpiration only
!* (2) additional input field 'lintw', containing leaf interception water [kg/m^2]
!*     which supercedes 'zlcontent', which is currently a scalar and =0
!* (3) replaced ZLCONTENT by LINTW(JL)
!* WOK, 2007-06-26
!* (4)  added 'zpcevp' & 'zpsevp' and removed 'zctrans' and other unused arguments;
!*      revised entire computetion of canopy and soil net radiation and potential evapotranspiration
!* ASSOCIATED CHANGES ELSEWHERE
!* WOK, 2007-06-12
!* (1) added 'ztrans' to mo_vegetation, veg_allocate and veg_deallocate
!* (2) added line 'ztrans = 0.' to initialization part of SUBROUTINE MODEL
!* WOK, 2007-06-26
!* (3) added 'zpcevp' & 'zpsevp' to mo_vegetation, veg_allocate and veg_deallocate, and as
!*     arguments to 'diagnostics', removed 'zctrans' from the same
!* (4) revised daily water flux variables and diagnostics, 'raet' and 'rpet' are now
!*     actual and potential (total soil and canopy) evapotranspiration
!* WOK, 2007-07-12
!* (5) minimum function for demand/supply equation re-introduced
!* (6) calcuation of transpiration was not conistent with stomatal parameter -> corrected
!*     by recalculation of 'ztrans' after setting final canopy conductance
!* (7) error in calculation of canopy temperature 'ptv' -> need to take into account the
!*     effective radiating canopy area; corrected, allowing for cases of LAI->0
!***********************************************************


!***********************************************************
!* OTHER CHANGES TO FILES OUTSIDE MO_HYDRO & MO_PHENO
!* WOK, 2007-06-25
!* (1) added daily carbon fluxes to mo_diagnostics and diagnostics_allocate, diagnostics_deallocate
!* (2) sum up daily C fluxes in 'diagnostics', reset them in new subroutine 'diagdayreset'
!* (3) renamed daily water fluxes and move handling to 'diagnostics' and new routine 'mdiagnostics'
!* (4) removed explicit diagnostics from 'model', now handled in 'diagnostics' and 'mdiagnostics'
!* (5) soubroutine diagnostics uses 'ndayspm'; loading soil moisture separated in
!*     'hydrology_read' ; strangely, subsequent diurnal cycles within a month almost but
!*     exactly reproduce GPP, NPP, and don't reproduce at all transpiration and int. evaporation
!* WOK, 2007-06-27
!* (6) update surface albedo with snow as function of cloudiness and solar angle in 'model'
!* (7) new calls to 'hydro_deallocate' and 'pheno_deallocate' in model.f90
!* WOK, 2007-06-28
!* (8) redefined phenology classes;
!*     moved definition of phenology classed to mo_pheno
!* WOK, 2007-06-29
!* (9) modified Makefile in subdirectory bethy, moving mo_hydro ahead of mo_pheno
!*     (because mo_pheno contains USE mo_hydro statement)
!*(10) specific leaf-area data in mo_carparams out-of-date, fill in numbers from Diss.
!*     => but need to go back to Schulze et al. (1994)
!* WOK, 2007-07-01
!*(11) moved 'zlai' from mo_pheno to mo_vegetation
!* WOK, 2007-07-04
!* (12) now fully moved the phenology table 'xph' to phenology and deleted it in mo_carparms;
!*      requires USE statement in initialize.f90.
!* (13) removed setting of 'zfc' in mo_carbon
!* WOK, 2007-07-05
!* (14) Correction for time of wet canopy moved to 'diagnostics'
!* WOK, 2007-07-12
!* (15) moved 'lai' from mo_pheno to mo_vegetation, veg_alloacte, veg_deallocate
!* (16) changed 'dptr' to 'dptrp' in m_hydro and diagnostics
!* WOK, 2007-07-14
!* (17) changed meaning of 'zlai', to record LAI at time of calling diurnal cycle, in 'model'
!* WOK, 2007-07-19
!* (18) installed calls to subroutines in mo_climsubday within 'model', this is running in
!*      parallel with the old file handling involving 'getmonth', 'climatereadday', 'climinday'
!*      'climinsubday', 'daytemp', and 'rad', which it should later supplant
!***********************************************************


!***********************************************************
!* NOTES ON WHAT TO CLEAN UP LATER
!* (1) in climate_allocate, daily climate fields are defined
!*     this replicates daily fields that are now read in mo_hydro;
!*     eventually, move away from monthly input data altogether
!*     but maybe allow reading in separate daily input data for spinup
!*     (as currently done in 'loadmfield' called from 'getmonth' in mo_io;
!*      there is also the routine 'loaddfield' in mo_io which is not used)
!* (2) replace the old handling of diurnal data with the new one in mo_climsubday, removing
!*     the obsolete routines (see point 17 above), and moving mo_climsubday into its own file mo_climsubday.f90
!*     clean up a lot of fields (like 'drt', 'dpp', etc.)
!* (3) treatment of 'inkind' -> this should also become compatible with the new mo_climsubday
!* (4) remove wgen.f90, as it is not used any more
!* (5) mysterious use of 'jday' in 'model': even though passed on to 'rad' via 'climinsubday',
!*     this variable, originally meant to contain Julian day, now counts time steps of
!*     the diurnal cycle calls; strangely, in 'rad' (called 'inday' here), it is dysfunctional;
!*     instead, the day within the year is known from the computation of 'spds' and 'cpds'
!*     in 'daytemp', which uses 'midmoday(rmonth)'
!*     => need to clean up computation of diurnal cycle to middle of the period
!*        through which it is called by modifying 'climinday', passing to it
!*        the middle Julian day of the period, i.e. pass (iday1-iday0)/2 to 'climinday' in 'model'
!* (6) the use of 'gday' (documented as 'Julian day') is equally strange: it was set to
!*     iday in 'model', a loop running only over 1, and then used in 'climinday' to read
!*     daily temperature and temperature range into an array
!* (7) check and adjust the accelerated spin-up of carbon pools in 'cbalance', in particular
!*     currently, it uses only monthly (not daily) NPP from the spin-up runs stored in
!*     sub-array 'nppv(0,:,:)'; will need daily NPP values from spin-up period
!* (8) 'bethy1', 'bethy2' & 'surface' do not use any check on executing grid cells that have
!*     no vegetation activity (as has full BETHY), this could be changed for efficiency.
!*     (Indeed, the variable name 'nveglist' used there is misleading, as 'veglist' has been
!*      abandoned.)
!* (9) fully substitute old routines by mo_climsubday - remove a whole range of subroutine, like
!*     climinday, climinsubday etc.
!*(10) clean up the parameter list of bethy1/2 and surface
!*(11) clean up use of inyear0, inyear1 vs. year0, year1 (in control, mo_namelist, mo_hydro)
!*(12) moved mo_climsubday into mo_hydro.f90 for easier compilation -> later move back
!*     into own file mo_climsubday.f90 and change bethy/Makefile accordingly
!***********************************************************
!* MAS/TXK, 2007-07-19
!* slight temporary changes fixes for TAF compliance
!* split off allocation, use separate subroutine
!***********************************************************


  IMPLICIT NONE
  
! this one in 'veg_allocate'
  REAL, ALLOCATABLE, DIMENSION (:) :: psoilst

! newly added hydrological fields, allocated here in mo_hydro:
! another field for soil moisture, this time calculate here locally
! => LATER TO REPLACE THE FIELDS THAT ARE READ FROM EXTERNAL DATA
!    INSTALL A SWITCH SO THAT FIELDS FROM STEP ONE CAN STILL BE USED?
! plant available soil moisture [mm]
  REAL, ALLOCATABLE, DIMENSION (:) :: pasm
! leaf interception water reservoir, max. on-leaf water [mm]
  REAL, ALLOCATABLE, DIMENSION (:) :: lintw
! runoff, snowmelt, throughfall, soil evaporation [mm/day]
  REAL, ALLOCATABLE, DIMENSION (:) :: runoff, dsnmelt, thruf
!WOK-CHG-070625 renamed the daily sums
!WOK-CHG-070711 dptr -> dptrp
! daily sums of transpiration, potential transpiration and intercept evaporation [mm]
  REAL, ALLOCATABLE, DIMENSION (:) :: dtrp, dptrp, dievp
! dtrp stored after call of seasonal cycle:
  REAL, ALLOCATABLE, DIMENSION (:) :: zdtrp
! daily sums of actual and potential soil evaporation, and canopy evaporation [mm]
  REAL, ALLOCATABLE, DIMENSION (:) :: dsevp, dpsevp, dpcevp
!WOK-ADD-070626
! declarations for soil evaporation model:
! Reference: Ritchie, J.T. 1972. Model for predicting evaporation from a
!            row crop with incomplete cover. Water Resour. Res. 8, 1204-1213.
! local variables
  REAL, ALLOCATABLE, DIMENSION (:) :: esum

!WOK-ADD-070626
!WOK-ADD-070627
! declarations for snow model
! daily snow evaporation, snowfall and snowmelt, snow water equiv. [mm], snow height [m], snow albedo
  REAL, ALLOCATABLE, DIMENSION (:) :: dsnevp, snow, snowh, rhosn
!WOK-ADD-070627
!WOK-ADD-070629 removed again, use Christian Reick's 'pseudo soil temperature' instead
! declarations for soil temperature model: soil temp. and 1.5m depth, 
! averge temp. 0-59 days ago, and annual mean temp. as average 0-364 days ago [deg C]
!  REAL, ALLOCATABLE, DIMENSION (:) :: tsoil15, temp30, tempann
! counters over how many consecutive days average has been taken
!  INTEGER :: ntemp30, ntempann

!MAS, 071005: moved to mo_climate & mo_calendar
!!MS$! declarations for daily climate input data
!!MS$! precipitation in mm/day
!!MS$  REAL, ALLOCATABLE, DIMENSION (:,:) :: dprecip, dtmin, dtmax, dswdown
!!MS$! total number of days of model run, first and last year of daily input data
!!MS$  INTEGER :: ddays, yearin0, yearin1
!!MS$! translation table from year (starting with year0 as 1), month, giving first day
!!MS$! and last day (used for the for daily input data)
!!MS$  INTEGER, ALLOCATABLE, DIMENSION (:,:) :: dayp0, dayp1
!!MS$  INTEGER, ALLOCATABLE, DIMENSION (:) :: dayspy  ! holds # of days per year
 
!-----------------------------------------------------------------------------------
! PARAMETERS (to optimize)
!-----------------------------------------------------------------------------------
! PARAMETERS AT EVERY SUB-GRID CELL
! maximum plant avail. soil moisture = "bucket size" [mm]
  REAL, ALLOCATABLE, DIMENSION (:) :: pasmmax

! GRID POINT-WISE PARAMETERS
! declarations for soil evaporation model:
! desorptivity [mm/sqrt(day)], phase 1 soil evaporation [mm]
  REAL, ALLOCATABLE, DIMENSION (:) :: des, evap1
! wet and dry soil albedo
  REAL, ALLOCATABLE, DIMENSION (:) :: rhosw, rhosd

! GLOBAL PARAMETERS
! maximum infiltration rate [mm/day]
  REAL, PARAMETER :: inflmax0 = 100.
! effective leaf water storage per LAI at daily time step [mm]
  REAL, PARAMETER :: lwmax0 = 0.1
!  REAL, PARAMETER :: lwmax0 = 1.0
! width of transition zone where bucket empties or fills only partially
  REAL, PARAMETER :: pasmcrit = 0.1

! AUXILLIARY FIELDS
! initial value of 'pasmmax' read from external data (on different grid!)
  REAL, ALLOCATABLE, DIMENSION (:) :: pasmmax0
! uncertainty of 'pasmmax', currently set to constant 100mm
  REAL, ALLOCATABLE, DIMENSION (:) :: pasmmax_u

CONTAINS
  
  !*********************************************************
  !*  SUBROUTINE inithydro
  !*  reads in soil moisture from file, initializes and allocates memory
  !*********************************************************
  
  SUBROUTINE hydro_allocate (ng, vp, scale)
    
    ! .. Use Statements ..
    USE mo_netcdf
    USE mo_constants
    USE mo_grid, ONLY: gridp, nlon, nlat
    USE mo_namelist, ONLY : nrun, grid_file, site_file

    IMPLICIT NONE
    INTEGER, INTENT(in) :: vp, ng, scale

! .. Local Scalars ..
    TYPE(ncfile) :: infile
    TYPE(ncvar) :: ncwmax, ncrhosw, ncrhosd, ncdes, ncevap1
    TYPE(ncvar) :: ncprecip, nctmax, nctmin, ncswdown
! .. Local Arrays ..
!   this one is new, for bucket size:
    REAL, ALLOCATABLE, DIMENSION(:,:) :: iload1, iload1b, iload1c, iload1d, iload1e
!   and this one for the complete daily meteorological data:
    REAL, ALLOCATABLE, DIMENSION(:,:,:) :: iload2

    INTEGER :: i, j, k, t, ms, iday
    INTEGER :: n, m, yr, rdays0    

!   write (*,*) 'hydro_allocate: ng, vp, scale: ', ng, vp, scale
! ... allocate memory for permanent fields
    ALLOCATE (pasm(vp), pasmmax(vp), lintw(vp))
    ALLOCATE (runoff(vp), thruf(vp), dsnmelt(vp))
    ALLOCATE (dtrp(vp), zdtrp(vp), dptrp(vp), dievp(vp), dsevp(vp), dpsevp(vp), dpcevp(vp))
    ALLOCATE (pasmmax0(ng))
    ALLOCATE (pasmmax_u(vp))
    ALLOCATE (rhosw(ng), rhosd(ng))
    ALLOCATE (des(ng), evap1(ng))
    ALLOCATE (esum(vp))
    ALLOCATE (dsnevp(vp), snow(vp), snowh(vp), rhosn(vp))
    dtrp = 0.
    zdtrp = 0.
    des = 0.
    evap1 = 0.
    dsevp = 0.

!***********************************************************
!* ... load bucket size and soil evaporation parameters from grid file in NetCDF
!* the bucket size is the same for all sub-grid cells within each
!* spatial 2x2 grid cell, but the array 'pasmmax' runs over all
!* sub-grid cells; hence the value is repeated for each sub-grid cell
!* within one spatial grid cell. For the optimisation, these value will
!* be priors. During the optimisation, it will be possible
!* to adjust the value for each sub-grid cell independently. Also, it might
!* later be possible to adjust even the individual priors of the sub-grid cells
!* within each grid cell to allow for differences in rooting depth of the
!* vegetation types present
!***********************************************************

    IF (scale==1) THEN

      ! .. allocate temporary memory
      ALLOCATE (iload1(nlon,nlat),iload1b(nlon,nlat),iload1c(nlon,nlat))
      ALLOCATE (iload1d(nlon,nlat),iload1e(nlon,nlat))
      infile%name=grid_file
      CALL ncopen(infile)
      ! get bucket size
      ncwmax%name='wmax'
      CALL ncread (infile,ncwmax,iload1)
      ! get wet-soil albedo
      ncrhosw%name='rhosw'
      CALL ncread (infile,ncrhosw,iload1b)
      ! get dry-soil albedo
      ncrhosd%name='rhosd'
      CALL ncread (infile,ncrhosd,iload1c)
      ! get desorptivity
      ncdes%name='des'
      CALL ncread (infile,ncdes,iload1d)
      ! get phase 1 evaporation
      ncevap1%name='evap1'
      CALL ncread (infile,ncevap1,iload1e)
      CALL ncclose(infile)
      n=0
      pasmmax0=0
!$taf loop = parallel
      DO i = nlat,1,-1
!$taf loop = parallel
         DO j = 1,nlon  
           IF (iload1(j,i)/=-9999.) THEN
             n=n+1
             pasmmax0(n) = iload1(j,i)
             rhosw(n) = iload1b(j,i)
             rhosd(n) = iload1c(j,i)
             des(n) = iload1d(j,i)
             evap1(n) = iload1e(j,i)
           ENDIF
         ENDDO
      ENDDO
    
      DEALLOCATE (iload1, iload1b, iload1c, iload1d, iload1e)

    ELSE ! (scale==2)

      ! .. read site information from file
      OPEN (120,file=trim(site_file),status='unknown',form='formatted')
      READ (120,*) !header

      DO i = 1, ng
         READ(120,'(60X,F7.1,2F5.2,F6.2,F7.1)') &
        &         pasmmax0(i), rhosw(i), rhosd(i), des(i), evap1(i)
         !write (*,'(F7.1,2F5.2,F6.2,F7.1)') pasmmax0(i), rhosw(i), rhosd(i), des(i), evap1(i)
      ENDDO
  
      CLOSE(120)

    ENDIF

    ! uncertainty for pasmmax, set to 100mm 
    DO k = 1, vp
      n = gridp(k)
      pasmmax_u(k) = 100.
    END DO
  
  END SUBROUTINE hydro_allocate

  SUBROUTINE inithydro
    USE mo_grid, ONLY: gridp
    USE mo_climate, ONLY : zrhos
    
!***********************************************************
!* ... initializations
!***********************************************************

! ... set water reservoirs and various diagnostic fluxes to 0
    pasm = 0.
    psoilst = 0.
    lintw = 0.
    runoff = 0.
    thruf = 0.
    dsnmelt = 0.
    dsevp = 0.
    dpsevp = 0.
    dievp = 0.
    dpcevp = 0.
    dsnevp = 0.
! ... initialize soil albedo (consistent with 0 soil moisture)
    zrhos = rhosd(gridp)
! ... reset state variables of soil evaporation model
    esum = 0.  ! re-use as small, overlapping bucket for soil evaporation
! ... reset state variables of snow model
    snow = 0.
    snowh = 0.
    rhosn = zrhos

  END SUBROUTINE inithydro

 
  SUBROUTINE hydrology (ng,vp,iday,keyday)

    ! .. Use Statements ..
    USE mo_constants
    USE mo_climate
    USE mo_grid, ONLY : gridp
    USE mo_vegetation, ONLY : zfc, lai
    USE mo_helper, ONLY : errf, maxx, minx

    IMPLICIT NONE

    ! .. Arguments
    INTEGER, INTENT(in) :: ng, vp, iday

    ! .. Local variables
    ! infiltration, canopy interception [mm/day]
    REAL :: infl, intcp, fintcp, flow, pasm_save, pasm_next, corr
    REAL :: tevap, evapmax, esum0, esum1, xpp, xpp_help, xe1, xe2
    REAL :: dpsn, xdtmp, fthw, fx, fsnh, densn, densn0
    REAL, DIMENSION (vp) :: mtemp
    INTEGER :: k, j, iday0, keyday


    !------------------------------------------------------
    !     ... snow balance and snow and soil albedo
    !------------------------------------------------------
!$TAF store snowh, snow, rhosn = day_tape, rec = keyday 
    dsnevp = 0.
!$taf loop = parallel
    DO k = 1, vp
       j = gridp (k)
       xdtmp = (dtmin(j,iday)+dtmax(j,iday))/2.
       !     set soil albedo (set it twice, also in 'model')
!       fx=MIN(MAX((psoilst(k)-0.9)/(1-0.9),0.),1.)
       fx = maxx ((psoilst(k)-0.9)/(1-0.9), 0., 1e-3)
       fx = minx (fx, 1., 1e-3)
       zrhos(k) = fx*rhosw(j) + (1-fx)*rhosd(j)
       IF (xdtmp<=5.0.OR.snow(k)>0.0) THEN
          !       snowfall and snowmelt
          dpsn = dprecip(j,iday) * MIN (MAX ((3.3 - xdtmp) / 4.4, 0.), 1.)
!          dpsn = dprecip(j,iday) * minx (maxx ((3.3 - xdtmp) / 4.4, 0., 1e-3), 1., 1e-3)
          !       snow sublimation rate
          dsnevp(k) = MIN (dpsevp(k), snow(k) + dpsn)
!          dsnevp(k) = MIN (dpsevp(k), snow(k) + dpsn, 1e-3)
          !       compute compaction due to gravity
          densn = snow(k) / MAX (snowh(k), 1e-6)
!          densn = snow(k) / maxx (snowh(k), 1e-6, 1e-2)
          densn = densn * (1. + 0.0229 * EXP (-0.021*densn+0.08*xdtmp)*snow(k)/2.)
          !       set lower bound of snow density
          densn = MAX (densn, 30.)
!          densn = MAX (densn, 30., 1e-1)
          !       update snow amount [kg / m^2] and snow height [m] without the new snow
!          dsnmelt(k) = MAX (3.22 * xdtmp, 0.)
          dsnmelt(k) = maxx (3.22 * xdtmp, 0., 1e-2)
          snow(k) = snow(k) - dsnmelt(k) - dsnevp(k)
!          dsnmelt(k) = dsnmelt(k) + MIN (snow(k), 0.)
          dsnmelt(k) = dsnmelt(k) + minx (snow(k), 0., 1e-3)
!          snow(k) = MAX (snow(k), 0.)
          snow(k) = maxx (snow(k), 0., 1e-1)
!          dsnevp(k) = dsnevp(k) + MIN (dsnmelt(k), 0.)
          dsnevp(k) = dsnevp(k) + minx (dsnmelt(k), 0., 1e-3)
!          dsnmelt(k) = MAX (dsnmelt(k), 0.)
          dsnmelt(k) = maxx (dsnmelt(k), 0., 1e-3)
!          dsnevp(k) = MAX (dsnevp(k), 0.)
          dsnevp(k) = maxx (dsnevp(k), 0., 1e-3)
          snowh(k) = snow(k) / densn 
!          fx = SIGN (0.5, xdtmp + 15.) + 0.5
          fx = errf ((xdtmp + 15.) / 1e-2)
          densn0 = (10. + 2.667 * (xdtmp + 30.)) * (1.-fx) + &
!               (50. + 1.7 * MAX (xdtmp + 15.,0.)**1.5) * fx
               (50. + 1.7 * maxx (xdtmp + 15., 0., 1e-2)**1.5) * fx
!          densn0 = MAX (densn0, 30.)
          densn0 = maxx (densn0, 30., 1e-1)
          !       add the new snow
          snow(k) = snow(k) + dpsn
          snowh(k) = snowh(k) + dpsn / densn0
          !       update snow albedo
!          rhosn(k) = MIN (rhosn(k) + dpsn/densn0*10., 0.80)
          rhosn(k) = minx (rhosn(k) + dpsn/densn0*10., 0.80, 1e-2)
!          fthw = SIGN (0.5, xdtmp) + 0.5
          fthw = errf (xdtmp / 1e-2)
!          fsnh = SIGN (0.5, snowh(k)-0.25) + 0.5
          fsnh = errf ((snowh(k)-0.25) / 1e-2)
          rhosn(k) = rhosn(k)  &
!               - MAX (0.107 - 0.214 * rhosn(k), 0.) * fthw * fsnh  &
               - maxx (0.107 - 0.214 * rhosn(k), 0., 1e-3) * fthw * fsnh  &
               - 0.071 * fthw * (1.-fsnh) &
               - 0.006 * (1.-fthw)
!          rhosn(k) = MAX (rhosn(k), zrhos(k))
          rhosn(k) = maxx (rhosn(k), zrhos(k), 1e-3)
!        if (k==321) write (9,'(I4,16F10.3)') iday,xdtmp,dprecip(j,iday),dpsn,dsnmelt(k),snow(k),&
!           dsnevp(k),dpsevp(k),snowh(k),densn,densn0,rhosn(k),fthw,fsnh,zrhos(k)
       ELSE
         dsnmelt(k) = 0.
         dsnevp(k) = 0.
       ENDIF
    ENDDO


!------------------------------------------------------
!     ... canopy water balance
!------------------------------------------------------
!$TAF store lintw = day_tape, rec = keyday 
    DO k = 1, vp
       j = gridp (k)  ! grid cell number as function of sub-grid cell number
       !     canopy interception
       ! WOK-ADD-10-03-29 Scaling potential canopy evaporation by LAI intercept factor
       fintcp = zfc(k)*EXP (-0.5*lai(k)/maxx(zfc(k),1e-6,1e-5))
       intcp = fintcp*dprecip(j,iday)
       !     throughfall (="drip-over" of leaf interception water)
       thruf(k) = maxx (lintw(k) + intcp - lwmax0*lai(k), 0., 1e-3)
       !     re-evaporation from the canopy
       dievp(k) = minx (dpcevp(k)*fintcp, intcp - thruf(k) + lintw(k), 1e-3)
       !     canopy water update
       lintw(k) = lintw(k) + intcp - thruf(k) - dievp(k)
       !     add non-intercepted rain to throughfall (=total throughfall)
       thruf(k) = thruf(k) + dprecip(j,iday) - intcp
!if(k==438) write (9,'(i5,9e20.12)') iday,dprecip(j,iday), thruf(k),fintcp,fintcp*dpcevp(k),dievp(k)
ENDDO

!!$TAF store dsnevp, esum, thruf = day_tape, rec = keyday 
!    DO k = 1, vp
!       j = gridp (k)  ! grid cell number as function of sub-grid cell number
!------------------------------------------------------
!     ... soil evaporation
!------------------------------------------------------
! WOK 2009-01-27 replaced three times SIGN (implemented as step function) by error function with range 0.01
!!       tevap = (MAX(esum(k) - evap1(j), 0.) / des(j)) ** 2
!       tevap = ( maxx (esum(k) - evap1(j), 0., 1e-3) / des(j)) ** 2
!       !txk temp change for differentibility      dsevp(k) = des(j) * (sqrt (tevap+1.) - sqrt(tevap)) + thruf(k)
!       !WOK-070723 small number in second 'sqrt' decreases to make it numerically harmless
!       dsevp(k) = des(j) * (SQRT (tevap+1.) - SQRT(tevap+1e-3)) + thruf(k)
!       !FastOpt; aux array introduced to work around a TAF bug
!       xpp_help = thruf(k)-dsevp(k)
!       !xpp = SIGN (0.5,1.8*xpp_help) + 0.5
!       xpp = errf (1.8*xpp_help/1e-2)
!       dsevp(k) = dsevp(k) * (1.-xpp) + 0.8 * thruf(k) * xpp
!!       evapmax = MAX (dpsevp(k) - dsnevp(k), 0.)
!       evapmax = maxx (dpsevp(k) - dsnevp(k), 0., 1e-3)
!!       esum0 = MAX (esum(k) - thruf(k), 0.)
!       esum0 = maxx (esum(k) - thruf(k), 0., 1e-3)
!       !xe1 = SIGN (0.5, evap1(j)-esum0) + 0.5
!       xe1 = errf ((evap1(j)-esum0)/1e-2)
!       dsevp(k) = dsevp(k) * (1.-xe1) + evapmax * xe1
!!       dsevp(k) = MIN (dsevp(k), evapmax)
!       dsevp(k) = minx (dsevp(k), evapmax, 1e-3)
!       esum1 = esum0 + dsevp(k)
!       !xe2 = (SIGN (0.5, esum1-evap1(j)) + 0.5) * xe1
!       xe2 = errf((esum1-evap1(j))/1e-2) * xe1
!       dsevp(k) = dsevp(k) * (1.-xe2) + &
!            (evapmax - 0.4 * (esum1 - evap1(j))) * xe2
!!       if(k==2.AND.abs(iday-1269).LE.1) write (*,'(I5,9E20.12)') &
!!       iday, pasm(k), pasm(k)/pasmmax(k), dsevp(k), xe2, esum1, evapmax, evap1(j)
!       esum(k) = esum0 + dsevp(k)
!!       dsevp(k) = MAX (dsevp(k), 0.)
!       dsevp(k) = maxx (dsevp(k), 0., 1e-3)
!!       if(k==2.AND.abs(iday-1269).LE.1) write (*,'(I5,9E20.12)') &
!!       iday, pasm(k), pasm(k)/pasmmax(k), dsevp(k), esum(k), dpsevp(k)
!!      really dumb model:
!       dsevp(k) = dpsevp(k)
!!      limit soil evaporation to 10% of total soil water pool
!       dsevp(k) = dsevp(k) * errf (pasm(k)-0.1*pasmmax(k))
!    END DO

!------------------------------------------------------
!     ... soil water balance
!------------------------------------------------------
!FastOpt !$TAF store pasm, dievp, dtrp, dsnmelt = day_tape, rec = keyday 
!$TAF store pasm, dievp, zdtrp, dsnmelt = day_tape, rec = keyday 
!$taf loop = parallel
    DO k = 1, vp
!------------------------------------------------------
!     ... correction of daily transpiration
!------------------------------------------------------
       dtrp(k) = maxx (zdtrp(k) - dievp(k), 1e-6, 1e-3)

!------------------------------------------------------
!    ... soil evaporation
!------------------------------------------------------
       dsevp(k) = esum(k)/min(5.,pasmmax(k)) * dpsevp(k)
       esum(k) = esum(k) + thruf(k) + dsnmelt(k) - dsevp(k) - dtrp(k)
       esum(k) = minx (esum(k), min(5.,pasmmax(k)), 1e-1) ! the size of the soil evaporation bucket
       esum(k) = maxx (esum(k), 0., 1e-1)
!WOK2009-28-01 replaced hard bounds at 0 and pasmmax soil moisture with a statistical approach
!       if(k==438) write(9,'(i5,9E24.16)') iday,esum(k),min(5.,pasmmax(k)),dpsevp(k),dsevp(k)
!       if(k==438) write(9,'(i5,9E24.16)') iday,dptrp(k),dtrp(k),zdtrp(k)

!------------------------------------------------------
!    ... infiltration
!------------------------------------------------------
       !infl = MIN (thruf(k) + dsnmelt(k), pasmmax(k) - pasm(k), inflmax0)
       infl = minx (pasmmax(k) - pasm(k), thruf(k) + dsnmelt(k), pasmcrit)
       infl = minx (infl, inflmax0, pasmcrit)

!------------------------------------------------------
!    ... runoff
!------------------------------------------------------
       runoff(k) = thruf(k) + dsnmelt(k) - infl

!------------------------------------------------------
!    ... daily soil water balance
!------------------------------------------------------
       flow = infl - dsevp(k) - dtrp(k) ! the net in/out flow
       pasm_save = pasm(k)              ! last time step's soil moisture
       pasm_next = pasm(k)+flow         ! tentative next time step's soil moisture
!       if(k==2.AND.abs(iday-1269).LE.2) write (*,'(I5,9E20.12)') &
!       iday, pasm(k), pasm(k)/pasmmax(k), flow, infl, dsevp(k), dtrp(k), thruf(k), dsnmelt(k)
       pasm(k) = maxx (pasm_next, 0., pasmcrit)
       corr = pasm(k) - pasm_save - flow   ! the flow balance correction (always negative)
       dsevp(k) = dsevp(k) + corr          ! add to (i.e. subtract from) soil evaporation
!       if(k==1) write (9,'(I5,9E20.12)') iday, pasm(k), runoff(k), thruf(k), dsnmelt(k), infl
!       if(k==1) write (9,'(I5,9E20.12)') iday, pasm(k), flow, infl, dsevp(k), dtrp(k), thruf(k), dsnmelt(k)
!if(k==2) write (9,'(I5,9E20.12)') iday, pasm(k), pasm(k)/pasmmax(k), flow, infl, dsevp(k), dtrp(k), thruf(k), dsnmelt(k)
!if(k==498) write (9,'(I5,99E20.12)') iday, dprecip(gridp(k),iday), pasm(k), pasm(k)/pasmmax(k), flow, infl, dsevp(k), dtrp(k), thruf(k), dsnmelt(k), runoff(k),lai(k)


!------------------------------------------------------
! end of loop over the grid cells
!------------------------------------------------------
    END DO


! ACTIVATE THIS LINE TO PASS FRACTIONAL SOIL MOISTURE TO REST OF MODEL
!FastOpt !$TAF store pasm = day_tape, rec = keyday 
    psoilst = pasm / pasmmax

  END SUBROUTINE hydrology
  

  SUBROUTINE hydro_deallocate(ng,vp)
  
  IMPLICIT NONE
  INTEGER, INTENT(in) :: ng, vp  

!MAS, 071005: moved to mo_climate  
!!MS$    DEALLOCATE (dayp0,dayp1)
!!MS$    DEALLOCATE (dprecip, dtmin, dtmax, dswdown)
!$taf next required = vp
    DEALLOCATE (pasm, pasmmax, lintw)
!$taf next required = vp
    DEALLOCATE (runoff, thruf, dsnmelt)
!$taf next required = vp
    DEALLOCATE (dtrp, zdtrp, dptrp, dievp, dsevp, dpsevp, dpcevp)
!$taf next required = vp
    DEALLOCATE (pasmmax_u)
!$taf next required = ng
    DEALLOCATE (pasmmax0)
!$taf next required = ng
    DEALLOCATE (rhosw, rhosd)
!$taf next required = ng
    DEALLOCATE (des, evap1)
!$taf next required = ng
    DEALLOCATE (esum)
!$taf next required = ng
    DEALLOCATE (dsnevp, snow, snowh, rhosn)
!    DEALLOCATE (tsoil15, temp30, tempann)

  END SUBROUTINE hydro_deallocate
  
END MODULE mo_hydro
