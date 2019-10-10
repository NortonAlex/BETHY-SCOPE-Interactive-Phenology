! module for  namelist, pjr
MODULE mo_namelist
  
  IMPLICIT NONE

  ! beginning and end of run, will normalize names with Marko
  INTEGER :: year0   = 1979 
  INTEGER :: year1     = 1999

!HEW-ADD-BEG: added options from adbethy.options, now only read in namelist from control
  INTEGER :: irestart  ! restart mode :   0 initial run, 1 restart run
  INTEGER :: day0      ! initial julian day for simulation
  INTEGER :: nrun      ! number of years in this run (change for restart run if app.)
  INTEGER :: nspin     ! number of spin-up years
  INTEGER :: intime    ! time resolution of input data: 
                       !    1 :monthly, 2:daily, 3:6-hourly
  INTEGER :: trans     ! interannual mode ?? :
                       !    0:climatological, 1: transient run
  INTEGER :: inkind    ! radiation flag :  
                       !    1 :PAR, 2:cloudcover, 3:shortwave radiation
  INTEGER :: dtime     ! time step in minutes
  INTEGER :: dayint    ! interval in days the diurnal cycle routines (bethy1/2, 
                       ! surface) are invoked
  INTEGER :: outtime   ! output time step: 
                       !    0 yearly, 1 monthly, 2 daily
  ! grid resolution and region
  REAL :: grdwid
  REAL :: lone, lonw, latn, lats
  INTEGER :: wetflg, wiflg, vflg, rflg
!!HEW-DEL: Never used ???  INTEGER :: rcmode 
!FastOpt use ./input instead of ../input
  ! master flux data file, contains list of flux sites
  ! CHARACTER(len=80) :: fluxfile      = "./control_bethy/eddy_data"    
  ! master station data file, contains station list
  CHARACTER(len=80) :: datfile      = "./control_bethy/gv2_data"    
  ! interpolation file
  CHARACTER(len=80) :: interp_file  = "./input/bethy2tm.mat" 
  ! list of parameters and probably priors
  CHARACTER(len=80) :: param_file   = "./input/test_params"  
  ! file to describe mapping from params to gridcells
  CHARACTER(len=80) :: mapping_file_global = "./input/test_mapping" 
  ! file to describe bethy 2.0*2.0 grid
  CHARACTER(len=80) :: grid_file = "./input/adbethy_grid_79-93.nc" 
  ! file to describe site specific data
  CHARACTER(len=80) :: site_file = "./input/site_specs.dat" 
  ! first part of fapar file names for sites
  CHARACTER(len=80) :: faparfile = "fapar_redres_" 
  ! global fapar file name
  CHARACTER(len=80) :: faparfile_global = "fapar_lores.nc" 
!WOK-ADD-070725 daily climate input files
  ! the following files contain daily climate input data
  INTEGER :: yearin0 = 1979
  INTEGER :: yearin1 = 2005
  CHARACTER(len=80) :: dprecip_file = "./input/climate/forcing.2x2.precip.1979-2005.nc" 
  CHARACTER(len=80) :: dtmax_file = "./input/climate/forcing.2x2.tmax.1979-2005.nc" 
  CHARACTER(len=80) :: dtmin_file = "./input/climate/forcing.2x2.tmin.1979-2005.nc" 
  CHARACTER(len=80) :: dswdown_file = "./input/climate/forcing.2x2.swdown.precip.1979-2005.nc" 
  ! file to read in dummyfluxes for costfunction
  CHARACTER(len=80) :: dummyfluxes = "./input/co2fluxes.dat" 
  ! directory for jacobians
  CHARACTER(len=80) :: jacdir = "./input/jactd/" 
  ! directory for flask station observations
  CHARACTER(len=80) :: datdir = "./input/gv04/" 
  ! directory for eddy flux site observations
  CHARACTER(len=80) :: fluxdir = "./input/eddy_sites/" 
  ! directory for background fluxes
  CHARACTER(len=80) :: bgrdir = "./input/background/" 
  ! patterns for directly added fluxes, basis functions in synthesis inversion  language
  CHARACTER( len=80) :: pattern_file = "flux_patterns.nc"
  CHARACTER(len=80) :: flux_temp_file = "flux_temp.bin" ! storing flux patterns

  ! directory for output files
  CHARACTER(len=80) :: outdir = "../output/" 

  ! fluorescence data
  CHARACTER(len=80) :: fsfile = "input/fsdata/fs_to_bethy-09-11.nc"

  ! data for SCOPE model
  CHARACTER(len=80) :: scopedir = "input/scope/"

  ! prescribed lai
  ! - global simulations 
  CHARACTER(len=80) :: plai_file = "no_file"
  ! - site simulations (default is no_file)
  CHARACTER(len=80) :: site_file_lai = "no_file"

  ! Test in the use of fs data: 0 = do not use, 1 = use
  INTEGER :: ifs = 1

  ! the following variables contain start and end year for site data
  INTEGER :: year0_site = 2000
  INTEGER :: year1_site = 2001
  ! file to describe mapping from params to sites
  CHARACTER(len=80) :: mapping_file_site = "./input/test_mapping" 

  ! flag option to provide forcing directly (e.g. hourly) for site runs
  LOGICAL :: siteforceflag = .FALSE.

  ! flags for pft frac optimisation
  LOGICAL :: optpftg = .FALSE.
  LOGICAL :: optpftl = .FALSE.

  ! flags for bucket size optimisation
  LOGICAL :: optbsmg = .FALSE.
  LOGICAL :: optbsml = .FALSE.

  ! beginning of first flux period
  INTEGER :: p1start = 1
  ! end of first flux period
  INTEGER :: p1end = 1
  ! start of second flux period
  INTEGER :: p2start = 1
  ! end of second flux period
  INTEGER :: p2end = 1

  ! vp block splitting parallelisation
  INTEGER :: nblocks = -1
  INTEGER :: iblock = -1
  CHARACTER(len=80) :: blockvpfile = "no_file"    ! default is no_file i.e. run without block parallelisation. This file specifies the veg-points to run. 

!HEW-CHG : all options read in here on namelist control:
!HEW-CHG  NAMELIST /control/ firstdat, lastdat, datfile, param_file, mapping_file, &
!HEW-CHG       & grid_file, interp_file, spinup_file, climate_file, dummyfluxes, jacdir, outdir, &
!HEW-CHG       & datdir, bgrdir, p1start, p1end, p2start, p2end
!HEW-ADD: added all options from adbethy.options here:
!WOK-ADD-070725 daily climate input files
  NAMELIST /control/ year0, year1, irestart, day0, nspin, &
       & intime, trans, inkind, dtime, dayint, outtime, &
       & grdwid, lonw, lone, latn, lats, &
       & wetflg, rflg, wiflg, vflg, &
       & datfile, param_file, mapping_file_global, &
       & grid_file, interp_file, faparfile,faparfile_global, &
       & yearin0, yearin1, dprecip_file, dtmax_file, dtmin_file, dswdown_file, plai_file, &
       & dummyfluxes, jacdir, outdir, year0_site, year1_site, mapping_file_site, &
       & site_file, site_file_lai, siteforceflag, pattern_file, flux_temp_file,&
       & datdir, fluxdir, bgrdir, optpftg, optpftl, optbsmg, optbsml, &
       & p1start, p1end, p2start, p2end, nblocks, iblock, blockvpfile


CONTAINS

  SUBROUTINE get_namelist

!HEW-ADD: 050304:
    USE mo_grid, ONLY: ng   !,vp,i1,i2
    USE mo_constants

    INTEGER :: nruns
    character(len=200) :: control_filename
    logical :: file_exists

    if( command_argument_count() .eq. 0) then
        control_filename = 'control'
    else if(command_argument_count() .eq. 1) then
        call get_command_argument(1, control_filename)
        inquire(file = trim(control_filename), exist = file_exists)
        if( .not. file_exists ) then
            write(0,*) 'control file "'//trim(control_filename)//'" not found ...'
            stop 
        end if
    else
        stop 'wrong number of arguments ...'
    end if

    WRITE(6,*) '# control file:                ',control_filename

! to read  namelist
!HEW-CHG: 050304: all control files now in dir control_bethy
    OPEN(unit=1, file=trim(control_filename), status='old')
    REWIND 1 
    
    READ(1,NML=control)
    CLOSE(1)
!HEW-DEL    WRITE(6,nml=control)

    WRITE(6,*) '# parameter file:              ',param_file
    WRITE(6,*) '# daily climate input files:   ',dprecip_file
    WRITE(6,*) dtmax_file
    WRITE(6,*) dtmin_file
    WRITE(6,*) dswdown_file
    WRITE(6,*) '# prescribed lai file (global run): ',plai_file
    WRITE(6,*) '# prescribed lai file (site run): ',site_file_lai
    WRITE(6,*) '# grid file:                   ',grid_file
    WRITE(6,*) '# initial year input data:     ',yearin0
    WRITE(6,*) '# end year input data:         ',yearin1
    WRITE(6,*) '# number of spin-up years:     ',nspin


    IF (datfile/='./control_bethy/no_stations') THEN
       WRITE(6,*) 
       WRITE(6,*) 'GLOBAL SIMULATION'
       WRITE(6,*) '# station data file:              ',datfile
       WRITE(6,*) '# mapping file global:            ',mapping_file_global

       print *, 'Grid points ng = ', ng
! Check grid settings
       IF (grdwid==0.5 .AND.  ng /= 62483) THEN
          STOP '0.5 degree reg grid must have 62483 land grid points'
       ELSEIF (grdwid==1. .AND.  ng /= 11069) THEN
          STOP '1 degree equ area grid must have 11069 land grid points'
       ELSEIF (grdwid==2. .AND.  ng /= 3462) THEN
          STOP '2 degree reg area grid must have 3462 land grid points'
       ELSEIF (grdwid==10. .AND.  ng /= 170) THEN
          STOP 'TM2 grid must have 170 land grid points'
       ENDIF

! Split veg-points into blocks 
!       IF ((iblock==-1) .OR. (nblocks==-1)) THEN
!           print*,'  *  vp block par check 1 * :: vp=',vp 
!           i1 = 1
!           i2 = vp
!       else if( (nblocks .ge. 1) .and. (iblock .ge. 1) .and. (iblock .le. nblocks)) then
!           print*,'  *  vp block par check 2 * :: vp=',vp
!           call get_splits(vp, nblocks, iblock, i1, i2)
!       else
!           stop 'iblock and/or nblocks not defined correctly...'
!       end if
!       write(*,*) 'veg-point block parallelisation: (i1,i2) =',i1,i2
!       write(*,*) 'nblocks, iblock =', nblocks, iblock

! Check timings
       IF (year0 < yearin0 .OR. year0 > yearin1) THEN
          PRINT*,'Start of model run, ',year0,', must be within input data time interval ',yearin0,yearin1
          STOP
       ELSEIF (year1 < yearin0 .OR. year1 > yearin1) THEN
          PRINT*,'End of model run, ',year1,', must be within input data time interval ',yearin0,yearin1
          STOP 
       ELSEIF (year0 > year1) THEN
          PRINT*,'End of model run, ',year1,', must be after start of model run, ',year0
          STOP 
       ELSEIF (yearin0 > yearin1) THEN
          PRINT*,'End of input data time interval, ',yearin1,', must be after start of input data time interval, ',yearin1
       ENDIF
       
       nrun = year1-year0+1
       
       IF (nrun<nspin) THEN
          PRINT*,'ATTENTION'
          PRINT*,'Number of spin-up years, ',nspin,', is larger than  number of simulation years, ',nrun
          PRINT*,' The simulation years will be recycled for the spin-up period.'
          PRINT*
       ENDIF

       ! Report options
       WRITE (*,*) '# initial year for simulation:   ', year0
       WRITE (*,*) '# end year for simulation:       ', year1
       WRITE (*,*) '# initial j. day for simulation: ', day0
       IF (day0<1.OR.day0>365) STOP '0 < day0 < 366'
       WRITE (*,*) '# number of simulation years:    ', nrun 

       SELECT CASE (intime)
       CASE (1)
          WRITE (*,*) '# climate input data, t-res.:     monthly'               
       CASE (2)
          WRITE (*,*) '# climate input data, t-res.:     daily'
       CASE (3)
          WRITE (*,*) '# climate input data, t-res.:     6-hourly'
       END SELECT
       SELECT CASE (inkind)
       CASE (1)
          WRITE (*,*) '# climate input data, kind:       PAR' 
       CASE (2)
          WRITE (*,*) '# climate input data, kind:       cloudcover'
       CASE (3)
          WRITE (*,*) '# climate input data, kind:       shortwave rad'
       END SELECT
       WRITE (*,*) '# model time step (min):         ', dtime
       WRITE (*,*) '# diurnal cycles interval:       ', dayint
       WRITE (*,*) '# transient input data:          ', trans
       SELECT CASE (outtime)
       CASE (0)
          WRITE (*,*) '# output time step:               yearly'
          outt=1
       CASE (1)
          WRITE (*,*) '# output time step:               monthly'
          outt=12
       CASE (2)
          WRITE (*,*) '# output time step:               daily'
          outt=jdpyear
       CASE DEFAULT
          outt=jdpyear        
       END SELECT
       WRITE (*,*) '# north/south latitude:          ', latn, lats
       WRITE (*,*) '# west/east longitude:           ', lonw, lone
       WRITE (*,*) '# grid width:                    ', grdwid
       WRITE (*,*) '# temp. scaling with rainfall:    ', wetflg
       WRITE (*,*) '# stochastic rainfall gen.:       ', rflg
       IF (wiflg==1) THEN
          WRITE (*,*) '# loading windfields as input'
       ENDIF
       IF (vflg==1) THEN
          WRITE (*,*) '# loading vapour pressure as input'
       ENDIF
    ELSE
       nrun=0
    ENDIF

    IF (site_file/='./control_bethy/no_sites') THEN
       WRITE(6,*) 
       WRITE(6,*) 'SITE SIMULATION'
       WRITE(6,*) '# site data file:              ',site_file
       WRITE(6,*) '# mapping file site:           ',mapping_file_site
       WRITE(6,*) '# initial year for simulation: ', year0_site
       WRITE(6,*) '# end year for simulation:     ', year1_site
       WRITE(6,*) '# siteforceflag:               ',siteforceflag

       nrun = year1_site-year0_site+1
       IF (nrun<nspin) THEN
          PRINT*,'ATTENTION'
          PRINT*,'Number of spin-up years, ',nspin,', is larger than  number of simulation years, ',nrun
          PRINT*,' The simulation years will be recycled for the spin-up period.'
          PRINT*
       ENDIF
       WRITE(6,*) '# number of simulation years:  ', nrun
       outt=12
    ENDIF


  END SUBROUTINE get_namelist
     
END MODULE mo_namelist
