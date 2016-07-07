!------------------------------------------------------------------
! set the number of control variables
!------------------------------------------------------------------
SUBROUTINE numbmod( nx )
  USE mo_namelist
  USE mo_constants
  USE mo_vegetation
  USE mo_climate
  USE mo_io
  USE mo_beta
  USE costf
  USE tm ! transport module
  USE bgr ! background or specified fluxes
  USE bethy2tm
  USE mo_surface
  USE mo_diagnostics
  USE mo_mapping
  USE mo_pheno
  USE mo_hydro
  USE mo_grid, ONLY: sp, vp, nlon, nlat, nrs
  USE mo_netcdf
  USE mo_prog
  USE mo_config
  USE mo_taf

  !ANorton. Fluorescence
  USE fluo_param, ONLY: fluo_initparam 
  USE fluo

  IMPLICIT NONE

  INTEGER, INTENT(out) :: nx

  ! local variables
  INTEGER :: i
  INTEGER :: nvar


  WRITE(*,*) ''
  WRITE(*,*) '***************************************************'
  WRITE(*,*) '*                     CCDAS V2                    *'
  WRITE(*,*) '* incl full BETHY for inverse modelling studies   *'
  WRITE(*,*) '**************************************************'
  WRITE(*,*) '' 

  nrs = 0.

! .. read namelist variables
  CALL get_namelist

  nvar = count_params(param_file)

  CALL init_param ( nvar, param_file)

! initialize ccn data
  n_stats = init_conc(datfile, year0, year1, datdir)
  IF (n_stats>0) THEN

! .. set output files 
     CALL setio(outdir, 1)

     CALL config_global(grid_file)
     help=0
     WHERE (vtype > 0) help=1
     vp = SUM(help)

     print*,'  ** check: vp is set (numbmod;) ** '

     IF (optpftg) THEN
        IF (optbsmg) THEN
           nx = nvar + vp + vp
        ELSE
           nx = nvar + vp 
        ENDIF
     ELSE
        IF (optbsmg) THEN
           nx = nvar + vp
        ELSE
           nx =  nvar
        ENDIF
     ENDIF

     CALL init_tm( datfile, year0, year1, jacdir, n_stats) ! initialize transport
     CALL init_bgr
     CALL init_bethy2tm( interp_file)
     CALL init_global_mapping(mapping_file_global) ! read in mapping file for global run
     
     ! ANorton.
     ! Here we initialise the parameters for the calculation of the fluoresecence...
     ! They are linked to SCOPE model. 
     CALL fluo_initparam

     ! Here we read the fluorescence data and allocate the array for the BETHY
     ! computed fluorescence 
     CALL  init_fluobethy( year0, year1)

     ! allocate prog for jacobians here to make TAF happy I hope
     ALLOCATE( prog_global(2,nrun,12,ng))
     prog_global = 0.0

  ELSE ! dimensions of conc and f_tm needed by TAF
     ALLOCATE(conc(1))
     ALLOCATE(f_tm(1,1,1))
     nx = nvar
     optpftg = .FALSE.
     optbsmg = .FALSE.
  ENDIF

! initialize eddy data
  n_sites = init_flux(site_file, year0_site, year1_site, fluxdir)
  IF (n_sites>0) THEN

     ! .. set output files 
     CALL setio(outdir, 2)

     nrs = year1_site - year0_site + 1

     CALL config_sites(site_file,n_sites)
     help=0
     WHERE (vtype > 0) help=1
     sp = SUM(help)

     IF (optpftl) THEN
        IF (optbsml) THEN
           nx = nx + sp + sp
        ELSE
           nx = nx + sp 
        ENDIF
     ELSE
        IF (optbsml) THEN
           nx = nx + sp
        ELSE
           nx =  nx
        ENDIF
     ENDIF

     CALL init_site_mapping(mapping_file_site) ! read in mapping file for site runs 

     ! allocate prog for jacobians here to make TAF happy I hope
     ALLOCATE( prog_sites(nrs,366,n_sites))
     prog_sites = 0.0

  ELSE
     optpftl = .FALSE.
     optbsml = .FALSE.
  ENDIF

  tspd   = 1440/dtime

  IF (maxyears < (nrun+nspin)) THEN
     WRITE(6,*) 'numbmod: inconsistent maxyears value for global run'
     WRITE(6,*) 'maxyears = nrun + nspin: ', maxyears,' = ',nrun,' + ', nspin
     STOP
  ELSEIF (maxyears < (nrs+nspin)) THEN
     WRITE(6,*) 'numbmod: inconsistent maxyears value for site runs'
     WRITE(6,*) 'maxyears = nrs + nspin: ', maxyears,' = ',nrs,' + ',nspin
     STOP
  ENDIF

  IF (n_stats.eq.0) THEN ! only sites
     maxkeyday = (nrs+nspin)*maxdpy*2
     maxkeydayint = (nrs+nspin)*maxdpy*2
     maxkeydiurnal = (nrs+nspin)*maxdpy*maxtspd*2
  ELSEIF (n_sites.eq.0) THEN ! only stats
     maxkeyday = (nrun+nspin)*maxdpy
     maxkeydayint = (nrun+nspin)*(1+maxdpy/dayint)
     maxkeydiurnal = (nrun+nspin)*(1+maxdpy/dayint)*maxtspd
  ELSE
     maxkeyday = (max(nrs,nrun)+nspin)*maxdpy*2
     maxkeydayint = (max(nrs,nrun)+nspin)*maxdpy*2
     maxkeydiurnal = (max(nrs,nrun)+nspin)*(1+maxdpy/dayint)*maxtspd*2
  ENDIF
  print*, 'dimensions for storing = ', maxkeyday, maxkeydiurnal

END SUBROUTINE numbmod
