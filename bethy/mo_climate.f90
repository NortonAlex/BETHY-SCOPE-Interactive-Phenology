MODULE mo_climate

!CCCC climate parameter
! txk, 09/01/08 : 
!      removes fdirsw in order to avoid alloc problem in TLM

! .. Use Statements ..
  USE mo_constants, ONLY :  pi

  IMPLICIT NONE

  INTEGER n_sites ! number of sites including sub_sites
  INTEGER r_sites ! number of sites
  INTEGER force_lai_flag   ! flag that switches on (=1) if a site LAI file is given
  INTEGER, ALLOCATABLE, DIMENSION(:) :: sub_sites ! number of sub_sites per site
  INTEGER n_stats ! number of stations
  INTEGER, ALLOCATABLE, DIMENSION(:) :: site_clim
  CHARACTER(len=20), ALLOCATABLE, DIMENSION(:) :: stat_names
  CHARACTER(len=20), ALLOCATABLE, DIMENSION(:) :: site_names

  REAL, ALLOCATABLE, DIMENSION (:,:) :: motmp

  
 ! .. Variables set by subroutine climsubday1, intended for output
  REAL, ALLOCATABLE, DIMENSION (:,:) :: tmp, mu
! original names:                       ---, coszen
  REAL, ALLOCATABLE, DIMENSION (:) :: eamin, pair, wspeed, swratio, cloudf, daylen
! original names:         ea0,   p,    u,      drt,     nc,     dayl

! .. Variables set by subroutine climsubday2, intended for output
  REAL, ALLOCATABLE, DIMENSION (:) :: pardown, swdown, fdirpar!FastOpt, fdirsw
  REAL, ALLOCATABLE, DIMENSION (:) :: lwdown    ! Only used in fluorescence
! original names:         par,     irrin,  pdir,   ,fdir

! Variables for optional sub-daily forcing in site runs
! shape = (tspd,n_sites) = (time-steps per day, number of sites i.e. ng)
! time series of forcing data
  REAL, ALLOCATABLE, DIMENSION (:,:) :: hlwdown,hswdown,hta,hea,hpair  
! the current diurnal time step loop forcing data
  REAL, ALLOCATABLE, DIMENSION (:)   :: ctlwdown,ctswdown,ctta,ctea,ctpair

! .. Variables set by subroutine climsubday1, for module-internal use
  REAL, PRIVATE :: dbodsq

! .. Parameters for model-internal use
  REAL, PRIVATE, PARAMETER :: mutiny = 1e-3
  REAL, PRIVATE, PARAMETER :: nodays = 365.
  REAL, PRIVATE, PARAMETER :: pio180 = pi/180.
  REAL, PRIVATE, PARAMETER :: s0 = 600.0
  REAL, PRIVATE, DIMENSION(5), PARAMETER :: a = &
  (/1.00011, 0.34221E-1, 0.128E-2, 0.719E-3, 0.77E-4/)


  REAL, ALLOCATABLE, DIMENSION (:,:) :: htmp
!  REAL, DIMENSION (ng) :: dayl
  REAL, ALLOCATABLE, DIMENSION (:) :: zrhos
  REAL, ALLOCATABLE, DIMENSION (:,:) :: dtmp, dtran, dpp, drt
  REAL, ALLOCATABLE, DIMENSION (:,:) :: dvp, dwind

! variables for daily input, used in climin and getmonth

  REAL, ALLOCATABLE, DIMENSION (:) :: coszen, spds, cpds

!  variables for zenith ankle, used in climin, rad, daytemp

! declarations for daily climate input data, with dimensions (ng, time) 
! precipitation in mm/day
  REAL, ALLOCATABLE, DIMENSION (:,:) :: dprecip, dtmin, dtmax, dswdown, prescribed_lai
! ! LAI input (sites only), with dimensions (vp, time)
  REAL, ALLOCATABLE, DIMENSION (:,:) :: dlai

CONTAINS


  SUBROUTINE get_global_climate

  ! .. Use Statements ..
    USE mo_netcdf
    USE mo_constants
    USE mo_namelist, ONLY : dprecip_file, dtmax_file, dtmin_file, dswdown_file, plai_file
    USE mo_calendar
    USE mo_grid, ONLY : nlon, nlat, vp

    IMPLICIT NONE

! .. Local Scalars ..
    TYPE(ncfile) :: infile
    TYPE(ncvar) :: ncprecip, nctmax, nctmin, ncswdown

! .. Local Arrays ..
    REAL, ALLOCATABLE, DIMENSION(:,:,:) :: iload
    INTEGER, DIMENSION(3) :: ncstart
    INTEGER :: i, j, n

!***********************************************************
!* ... load daily precipitation in NetCDF
!* later on, tmin, tmax, swdown and lwdown will also have to be
!* loaded and used with the surface routine, replacing the current
!* monthly input data
!***********************************************************

    ! .. allocate temporary memory
    ALLOCATE (iload(nlon,nlat,sdays))
    ncstart=(/1,1,nds0+1/)

    ! get daily precipitation
    infile%name=dprecip_file
    CALL ncopen(infile)
    ncprecip%name='precip'
    CALL ncread (infile,ncprecip,iload,start=ncstart)
    CALL ncclose(infile)
    n=0
!$taf loop = parallel
    DO i = nlat,1,-1
!$taf loop = parallel
       DO j = 1,nlon      
          IF (iload(j,i,1).gt.-1000) THEN 
             n=n+1
             dprecip(n,:) = iload(j,i,:)
          ENDIF
       ENDDO
    ENDDO


    ! get daily maximum temperature
    infile%name=dtmax_file
    CALL ncopen(infile)
    nctmax%name='tmax'
    CALL ncread (infile,nctmax,iload,start=ncstart)
    CALL ncclose(infile)
    n=0
!$taf loop = parallel
    DO i = nlat,1,-1
!$taf loop = parallel
       DO j = 1,nlon      
          IF (iload(j,i,1).gt.-1000) THEN        
             n=n+1
             dtmax(n,:) = iload(j,i,:)
          ENDIF
       ENDDO
    ENDDO

    ! get daily minimum temperature
    infile%name=dtmin_file
    CALL ncopen(infile)
    nctmin%name='tmin'
    CALL ncread (infile,nctmin,iload,start=ncstart)
    CALL ncclose(infile)
    n=0
!$taf loop = parallel
    DO i = nlat,1,-1
!$taf loop = parallel
       DO j = 1,nlon      
          IF (iload(j,i,1).gt.-1000) THEN        
             n=n+1
             dtmin(n,:) = iload(j,i,:)
          ENDIF
       ENDDO
    ENDDO

    ! get daily shortwave downward radiation
    infile%name=dswdown_file
    CALL ncopen(infile)
    ncswdown%name='swdown'
    CALL ncread (infile,ncswdown,iload,start=ncstart)
    CALL ncclose(infile)
    n=0
!$taf loop = parallel
    DO i = nlat,1,-1
!$taf loop = parallel
       DO j = 1,nlon      
          IF (iload(j,i,1).gt.-1000) THEN        
             n=n+1
             dswdown(n,:) = iload(j,i,:)
          ENDIF
       ENDDO
    ENDDO

    DEALLOCATE (iload)

    ! Read global-level LAI data if provided.
    IF (plai_file .ne. 'no_file') THEN
        print*,'# Using forced LAI data from plai_file'
        force_lai_flag = 1
        OPEN(unit=79,file=plai_file,form='formatted',status='old')
        REWIND 79
        DO i = 1,vp
           READ(79,*) (prescribed_lai(i,j),j=1,12)
        END DO
        CLOSE(79)
    ELSE    ! no site LAI file provided in control namelist, so setting it to zero
        force_lai_flag = 0
    END IF

  END SUBROUTINE get_global_climate

  SUBROUTINE climate_allocate(ng,vp,sdays)

    USE mo_constants, ONLY: tspd
    
    IMPLICIT NONE
    
    INTEGER, INTENT(in) :: ng, vp, sdays

!$taf next required = ng, sdays
    ALLOCATE (dprecip(ng,sdays), dtmin(ng,sdays), dtmax(ng,sdays), dswdown(ng,sdays))
    ALLOCATE (dlai(vp,sdays))
    ALLOCATE (motmp(12,ng))
    ALLOCATE (eamin(ng), pair(ng), wspeed(ng), swratio(ng), cloudf(ng), daylen(ng))
    daylen = 0.
    ALLOCATE (pardown(ng), swdown(ng), fdirpar(ng))!FastOpt, fdirsw(ng))
    ALLOCATE (lwdown(ng))
    ALLOCATE (coszen(ng), spds(ng), cpds(ng))
    ALLOCATE (htmp(ng,tspd))
    ALLOCATE (zrhos(vp))
    ALLOCATE (prescribed_lai(vp,12))
    ALLOCATE (hta(sdays*tspd,ng),hswdown(sdays*tspd,ng))
    ALLOCATE (hea(sdays*tspd,ng),hpair(sdays*tspd,ng))
    ALLOCATE (hlwdown(sdays*tspd,ng))
    ALLOCATE (ctswdown(ng),ctta(ng))
    ALLOCATE (ctea(ng),ctpair(ng))
    ALLOCATE (ctlwdown(ng))

    CALL climsubday_allocate(ng)

  END SUBROUTINE climate_allocate



  SUBROUTINE climate_deallocate(ng,vp,sdays)
    
    USE mo_constants, ONLY: tspd
    
    IMPLICIT NONE

    INTEGER, INTENT(in) :: ng, vp, sdays

    DEALLOCATE (dprecip, dtmin, dtmax, dswdown, dlai)
!$taf next required = ng
    DEALLOCATE (motmp)
!$taf next required = ng
    DEALLOCATE (eamin, pair, wspeed, swratio, cloudf, daylen)
!$taf next required = ng
    DEALLOCATE (pardown, swdown, fdirpar)!FastOpt, fdirsw)
    DEALLOCATE (lwdown)
!$taf next required = ng
    DEALLOCATE (coszen, spds, cpds)
!$taf next required = tspd
    DEALLOCATE (htmp) 
!$taf next required = vp
    DEALLOCATE (zrhos) 
    DEALLOCATE (prescribed_lai)
    DEALLOCATE (hta,hswdown,hea,hpair,hlwdown)
    DEALLOCATE (ctswdown,ctta,ctea,ctpair,ctlwdown)

    CALL climsubday_deallocate
    
  END SUBROUTINE climate_deallocate
  
!*********************************************************
!*  SUBROUTINE climsubday1
!*  generates sub-daily temperatures, base vapour pressure
!*  and daily constant variables for sub-daily calculations
!*
!* INPUT: (defined in this module or as arguments)
!*        dtmax, dtmin: daily minimum and maximum temperature [deg C]
!*        dswdown: daily mean shortwave downward radiation [W/m^2]
!*        (direct arguments)
!*        fe: actual / potential evapotranspiration from last time step (at sub-grid)
!*        iday0, iday1: Julian days of first and last day in the series
!*
!* FURTHER OPTIONAL INPUT (to be implemented later):
!*        diurnal (hourly or half-hourly) temperature
!*        diurnal vapour pressure deficit
!*        wind speed (daily mean or diurnal)
!*
!* OUTPUT: (all through variables declared within the module)
!*         representative average diurnal cycle of
!*         tmp: near-surface air temperature [deg C]
!*         mu: solar zenith angle for every hour, for use in climsubday2
!*         and daily values of
!*         eamin: vapour pressure at minimum temperature [Pa]
!*         pair: air pressure [Pa]
!*         wspeed: wind speed above canpopy [m/s]
!*         swratio: ratio of actual to potential shortwave incoming radiation
!*         cloudf: cloud fraction
!*         dayl: day length [hours]
!*********************************************************

  SUBROUTINE climsubday1 (ng, vp, fe, iday0, iday1)

! .. Use Statements ..
    USE mo_constants, ONLY: tspd, relh0, g, amd, ar, u0, pi
!    USE mo_namelist, ONLY : wiflg
    USE mo_grid, ONLY: elev, lat, gridp, frac
    USE mo_carparams, ONLY: laps, p0

    IMPLICIT NONE

! .. Arguments
    INTEGER, INTENT(in) :: ng, vp, iday0, iday1
    REAL, DIMENSION (vp), INTENT(in) :: fe

! .. Local Arrays
    REAL, DIMENSION (ng) :: tfe, atmean, atrange, aswdown, pswdown

! .. Local Scalars ..
    INTEGER :: i, j, ih
    REAL :: r
    REAL :: td, fthw, a0, b0, esta, estb
    REAL :: rdaymid, delta, arg, h0, h1, sd, sd1, hour, dhour, tmin, tmp1
    REAL :: alpha1, alpha2, rtop, tdir, ttot
 
! .. Intrinsic Functions ..
    INTRINSIC ACOS, COS, SIN, MOD, REAL

! .. compute average conditions over period iday0 ... iday1
    r = 1. / REAL (iday1-iday0+1)
    atmean = r * SUM (dtmax(:,iday0:iday1)+dtmin(:,iday0:iday1),2) / 2.
    atrange = r * SUM (dtmax(:,iday0:iday1)-dtmin(:,iday0:iday1),2)
    aswdown = r * SUM (dswdown(:,iday0:iday1),2)
  
! .. average actual / potential evapotranspiration over grid cell
    tfe=0
    DO i=1,vp
       j=gridp(i)
       tfe(j)=tfe(j)+fe(i)*frac(i)
    ENDDO

! .. compute daily course of temperature and daylength
    rdaymid = REAL (iday0+iday1) / 2.
    delta = -23.4*COS(2.*pi*(rdaymid+10.)/365.)
    DO i = 1,ng
      spds(i) = SIN(lat(i)*pio180)*SIN(delta*pio180)
      cpds(i) = COS(lat(i)*pio180)*COS(delta*pio180)
      arg = -spds(i)/cpds(i)
      IF (arg>1.) THEN
        !        polar night:
        daylen(i) = 0.
      ELSE IF (arg<-1) THEN
        !        polar day:
        daylen(i) = 24.
      ELSE
        !        normal day / night:
        daylen(i) = ACOS(arg)/pi*24.
      END IF
      IF (daylen(i)>=4..AND.daylen(i)<=20.) THEN
        !        sunrise
        h0 = 12. - daylen(i)/2.
        !        sundown
        h1 = 12. + daylen(i)/2.
        !        at sundown:
        sd1 = SIN(pi*(2.*h1+(daylen(i)-52.)/2.)/(daylen(i)+4.))
!! unroll zum vektorisieren
        DO ih = 1, tspd
          hour = REAL(ih)
          IF (hour>h0 .AND. hour<h1) THEN
            sd = SIN(pi*(2.*ih+(daylen(i)-52.)/2.)/(daylen(i)+4.))
            tmp(ih,i) = atmean(i) + atrange(i)/2.*sd
          ELSE
            ! temperature at sundown
            tmp1 = atmean(i) + atrange(i)/2.*sd1
            ! hours since sundown
            dhour = MOD(hour-h1+24.,24.)
            tmin = atmean(i) - atrange(i)/2.
            tmp(ih,i) = tmp1 - (tmp1-tmin)*(dhour/(24.-daylen(i)))
          END IF
        END DO
      ELSEIF (daylen(i)>20.) THEN
        DO ih = 1, tspd
          hour = REAL(ih)
          sd = COS(pi*(hour-14.)/(daylen(i)/2.+2.))
          tmp(ih,i) = atmean(i) + atrange(i)/2.*sd
        END DO
      ELSE
        tmp(:,i) = atmean(i)
      END IF
    ENDDO
!    i = 10
!    write (7,'(i5,24f12.5)') iday0, tmp(:,i)

! .. compute air pressure and saturation vapour pressure at diurnal minimum temperature
    DO i = 1,ng
      ! minimum temperature
      td = atmean(i) - atrange(i)/2.
      fthw = SIGN(0.5,atmean(i)) + 0.5
      a0 = 17.269*fthw + 22.33*(1.-fthw)
      b0 = 237.3*fthw + 271.15*(1.-fthw)
      esta = 610.78*EXP(a0*td/(b0+td))
!     base vapour pressure
      eamin(i) = (relh0+(1.-relh0)*tfe(i))*esta
      ! air pressure
      pair(i) = p0*(1./(1.+REAL(elev(i))*laps/(atmean(i)+273.15))) &
                **(g*amd*0.001/ar/laps)        
    ENDDO

! .. set wind speed 
    wspeed=u0 ! constant value (default)
!    IF (wiflg==1) u(:)=dwind(gday,:)

! compute the square of the inverse relative distance between the Sun
! and the Earth (cf. Knapp, C. et al. 1980. Insolation Data Manual,
! Solar Energy Research Institute, p. 252; Paltridge, G. W. and C. M.
! Platt. 1976. Radiative Processes in Meteorology and Climatology,
! Elsevier, p. 57.)
    alpha1 = 2.0*pi*(rdaymid-1.)/nodays
    alpha2 = alpha1*2.0
!   compute the Sun-Earth inverse squared distance:
    dbodsq = a(1)+a(2)*COS(alpha1)+a(3)*SIN(alpha1)+a(4)*COS(alpha2)+ &
          &  a(5)*SIN(alpha2)

! .. compute potential diurnal average incident solar radiation
    pswdown = 0.
    DO ih = 1, tspd
      DO i = 1, ng
        hour = REAL (ih)
        mu(ih,i) = spds(i) - cpds(i)*COS(hour*pi/12.)            
        ! compute terrestrial solar irradiance            
        !     check if solar zenith angle > 89 degrees
        IF (mu(ih,i)>=mutiny) THEN
          ! compute terr. solar irradiance and direct frac for PAR and S
          rtop = s0*dbodsq*mu(ih,i)         
          tdir = EXP(-.185/mu(ih,i)*pair(i)/p0)
          ttot = (.4+.6*tdir)
          pswdown(i) = pswdown(i) + rtop * ttot / 0.43
        ENDIF
      ENDDO
    ENDDO
    pswdown = pswdown / 24.

! .. compute ratio of actual to potential shortwave downward radiation and cloud fraction
    WHERE (aswdown>1e-3.AND.pswdown>=1e-3)
      swratio = aswdown / pswdown
    ELSEWHERE
      swratio = 0.
    END WHERE
    cloudf = MIN((0.9-MIN(swratio,0.9))/0.4,1.)
    
!    i = 10
!    write (8,'(i5,10f12.5)') iday0, daylen(i), eamin(i), pair(i), wspeed(i), cloudf(i), swratio(i), aswdown(i), pswdown(i)

    
    !*         tmp: near-surface air temperature [deg C]
!*         mu: solar zenith angle for every hour, for use in climsubday2
!*         and daily values of
!*         eamin: vapour pressure at minimum temperature [Pa]
!*         pair: air pressure [Pa]
!*         wspeed: wind speed above canpopy [m/s]
!*         swratio: ratio of actual to potential shortwave incoming radiation
!*         cloudf: cloud fraction



!    i = 80
!    write (9,'(a,i5,10f12.5)') 'iday0: ', iday0, atmean(i), atrange(i), daylen(i), cloudf(i), swratio(i), aswdown(i), pswdown(i)
!    write (9,'(a,i5,10f12.5)') 'iday0: ', iday0, atmean(i)-0.5*atrange(i), relh0, tfe(i), eamin(i)
!    write (9,'(a,i5,10f14.5)') 'iday0: ', iday0, atmean(i), REAL(elev(i)), p0, pair(i)

!    i = 10
!    do j = iday0, iday1
!      write (*,*) j, (dtmin(i,j)+dtmax(i,j))/2.
!      write (*,*) j, (dtmax(i,j)-dtmin(i,j))
!      write (*,*) j, dswdown(i,j)
!    enddo
!    write (*,*) atmean(i)
!    write (*,*) atrange(i)
!    write (*,*) aswdown(i)

END SUBROUTINE climsubday1


!*********************************************************
!*  SUBROUTINE climsubday2
!*  generates sub-daily variables and is called
!*  within the sub-daily (usually hourly) loop
!*
!* INPUT: (through variable declared in module)
!*        swratio: ratio of actual to potential shortwave downward radiation
!*
!* FURTHER OPTIONAL INPUT (to be implemented later):
!*        instantaneous photosynthetically active radiation (PAR)
!*        instantaneous shortwave downward radiation
!*
!* OUTPUT: (through variables declared in module)
!*         instantaneous sub-daily (hourly) values of
!*         swdown: incident solar radiation [W/m^2]
!*         pardown: incident photosynthetically active radiation (PAR) [W/m^2]
!*         fdirsw: direct (as opposed to diffuse) fraction of incident solar radiation
!*         fdirpar: direct fraction of incident PAR
!*********************************************************

  SUBROUTINE climsubday2 (ng, hour)

! .. Use Statements ..
    USE mo_carparams, ONLY: p0

    IMPLICIT NONE

! .. Arguments
    REAL, INTENT(in) :: hour ! current time in hours
    INTEGER, INTENT(in) :: ng ! number of grid points/sites

! .. Local Arrays
    REAL, DIMENSION (ng) :: tfe, mtmean, mtrange, aswdown

! .. Local Scalars ..
    INTEGER :: i, ih
    REAL :: rtop, tdir, ttot, ppot, rpx, cf, cfdir, fx
 
! .. compute various radiation components
    ih = REAL (hour + 0.5)
    DO i = 1,ng
      !     check if solar zenith angle > 89 degrees or not
      IF (mu(ih,i)<mutiny) THEN
        pardown(i) = 0.
        swdown(i) = 0.
        fdirpar(i) = 0.
!FastOpt        fdirsw(i) = 0.        
      ELSE               
        ! compute terrestrial solar irradiance            
        rtop = s0*dbodsq*mu(ih,i)
        tdir = EXP(-.185/mu(ih,i)*pair(i)/p0)
       ! compute terr. solar irradiance and direct frac for PAR and S
        ttot = (.4+.6*tdir)
        ppot = rtop*ttot
        pardown(i) = swratio(i)*ppot
        rpx = MIN(swratio(i),0.9)
        cf = .43 + .25 * (1. - rpx) * mu(ih,i)
        cfdir = .34 +  .094 * mu(ih,i)
        swdown(i) = pardown(i)/cf
        fx = (1.-((0.9-rpx)/0.7)**.6667)
        fx = MAX (MIN (fx, 1.), 0.)
        fdirpar(i) = tdir/ttot*fx
!FastOpt        fdirsw(i) = fdirpar(i) * cf / cfdir
      ENDIF
    ENDDO


  END SUBROUTINE climsubday2


  SUBROUTINE climsubday_allocate(ng)
  
    USE mo_constants, ONLY: tspd

    IMPLICIT NONE

    INTEGER, INTENT(in) :: ng ! number of grid points/sites    

    ALLOCATE ( tmp(tspd,ng), mu(tspd,ng) )  

  END SUBROUTINE climsubday_allocate


  SUBROUTINE climsubday_deallocate
  
    IMPLICIT NONE

    DEALLOCATE (tmp, mu)  

  END SUBROUTINE climsubday_deallocate


  SUBROUTINE get_local_climate
    ! .. Use Statements ..
    USE mo_netcdf
!txk    USE costf, ONLY: site_names, n_sites, site_clim
    USE mo_constants
    USE mo_namelist, ONLY : year0_site, year1_site, site_file_lai
    USE mo_calendar
    USE mo_grid, ONLY : vp

    IMPLICIT NONE

! .. Local Arrays ..
    REAL, DIMENSION(48,13) :: in
    REAL, DIMENSION(48,3) :: infn
    REAL, ALLOCATABLE, DIMENSION(:) :: iload

! .. Local Characters ..    
    CHARACTER(len=80) :: site_file, header
    CHARACTER (len=4) :: cy
    CHARACTER(len=80) :: in_fluxdir

! .. Local Scalars ..
    INTEGER :: i, j, ios, ny, y, n, lt, row, col, endyr
    REAL, ALLOCATABLE, DIMENSION(:,:) :: lrec

    TYPE(ncfile) :: infile
    TYPE(ncvar) :: ncprecip, nctmax, nctmin, ncswdown

    in_fluxdir = "input/eddy_sites/"


    ny = year1_site - year0_site +1

    OPEN(unit=10,file=TRIM(in_fluxdir)//'forcing.US-NR1.2015.txt',form='formatted')
    READ(10,*) endyr
    close(10)

    If (year1_site <= endyr) then

       DO n = 1,n_sites

          i = 1
          j = 1

!!MS$          IF (site_clim(n)==1) THEN
             OPEN(unit=10,file=TRIM(in_fluxdir)//'forcing.US-NR1.2015.txt',form='formatted')
             read(10,*) endyr
             READ(10,*) header
             DO WHILE (header .NE. site_names(n))
                READ(10,*)
                READ(10,*)
                READ(10,*)
                READ(10,*)
                READ(10,*) header
             END DO
             READ(10,'(365f8.2)') dtmax(n,:)
             READ(10,'(365f8.2)') dtmin(n,:)
             READ(10,'(365f8.2)') dprecip(n,:)
             READ(10,'(365f8.2)') dswdown(n,:)
             CLOSE(10)

             ! Read site-level LAI data if provided.
             IF (site_file_lai .ne. 'no_file') THEN
                 print*,'# Using forced LAI data from site_file_lai'
                 force_lai_flag = 1
                 OPEN(unit=80,file=TRIM(site_file_lai),form='formatted')
                 READ(80,*) endyr
                 READ(80,*) header
                 DO WHILE (header .NE. site_names(n))
                    READ(80,*)
                    READ(80,*)
                    READ(80,*)
                    READ(80,*)
                    READ(80,*) header
                 END DO
                 DO i = 1,vp
                    READ(80,'(365f8.2)') dlai(i,:)
                 END DO
                 CLOSE(80)
             ELSE    ! no site LAI file provided in control namelist, so setting it to zero
                 force_lai_flag = 0
                 dlai = 0.
             END IF

!!MS$          ELSE
!!MS$             IF (site_names(n)=='hainich') THEN
!!MS$                DO y=1,ny          
!!MS$                   WRITE(cy,'(i4)') (year0_site - 1 + y)
!!MS$                   site_file= TRIM(in_fluxdir)//TRIM(site_names(n))//'_'//TRIM(cy)//'.txt'       
!!MS$                   OPEN(unit=10,file=TRIM(site_file),form='formatted')
!!MS$                   READ(10,*) header
!!MS$
!!MS$                   DO         
!!MS$                      READ(10,*,iostat=ios) in(i,:)
!!MS$                      IF (ios<0) EXIT    
!!MS$
!!MS$                      IF (in(i,5)==-9999..AND.i/=1) THEN
!!MS$                         WRITE(99,*) TRIM(site_names(n)),': missing value in R_g, used value from previous half hour'
!!MS$                         in(i,5) = in(i-1,5)
!!MS$                      ELSEIF (in(i,5)==-9999..AND.i==1) THEN
!!MS$                         in(i,5) = in(48,5)
!!MS$                      ENDIF
!!MS$
!!MS$                      IF (in(i,7)==-9999..AND.i/=1) THEN
!!MS$                         WRITE(99,*) TRIM(site_names(n)),': missing value in T_air, used value from previous half hour'
!!MS$                         in(i,7) = in(i-1,7)
!!MS$                      ELSEIF (in(i,7)==-9999..AND.i==1) THEN
!!MS$                         in(i,7) = in(48,7)
!!MS$                      ENDIF
!!MS$
!!MS$                      IF (in(i,13)==-9999.) THEN
!!MS$                         WRITE(99,*) TRIM(site_names(n)),': missing value in Precip, half hour value set to zero'
!!MS$                         in(i,13) = 0.
!!MS$                      ENDIF
!!MS$
!!MS$                      IF (MOD(i,48)==0) THEN
!!MS$                         dtmax(n,j) = MAXVAL(in(:,7))
!!MS$                         dtmin(n,j) = MINVAL(in(:,7))
!!MS$                         dswdown(n,j) = SUM(in(:,5))/48.
!!MS$                         dprecip(n,j) = SUM(in(:,13))
!!MS$                         i = 0
!!MS$                         j = j + 1
!!MS$                      ENDIF
!!MS$                      i = i + 1
!!MS$                   ENDDO
!!MS$                   CLOSE(10)
!!MS$                ENDDO
!!MS$             ELSE
!!MS$                DO y=1,ny
!!MS$                   WRITE(cy,'(i4)') (year0_site - 1 + y)
!!MS$                   IF ((year0_site - 1 + y) == 2000) THEN
!!MS$                      lt=17568
!!MS$                      ALLOCATE(lrec(lt,255))
!!MS$                   ELSE
!!MS$                      lt=17520
!!MS$                      ALLOCATE(lrec(lt,255))
!!MS$                   ENDIF
!!MS$
!!MS$                   site_file= TRIM(in_fluxdir)//TRIM(site_names(n))//'.'//TRIM(cy)//'.flux.hourly.bin'
!!MS$
!!MS$                   OPEN(unit=10,file=TRIM(site_file),form='unformatted')
!!MS$
!!MS$                   READ(unit=10) lrec
!!MS$                   CLOSE(10)
!!MS$
!!MS$                   DO col=1,lt             
!!MS$                      IF (col>1 .AND. lrec(col,11)==-9999.) THEN
!!MS$                         lrec(col,11)=lrec(col-1,11)
!!MS$                         WRITE(99,*) TRIM(site_names(n)),': missing value in T_air, used value from previous half hour'
!!MS$                      ENDIF
!!MS$                      infn(i,1)=lrec(col,11)
!!MS$                      IF (col>1 .AND. lrec(col,12)==-9999.) THEN
!!MS$                         lrec(col,12)=lrec(col-1,12)
!!MS$                         WRITE(99,*) TRIM(site_names(n)),': missing value in R_g, used value from previous half hour'
!!MS$                      ENDIF
!!MS$                      infn(i,2)=lrec(col,12)
!!MS$                      IF (col>1 .AND. lrec(col,153)==-9999.) THEN
!!MS$                         lrec(col,153) = 0
!!MS$                         WRITE(99,*) TRIM(site_names(n)),': missing value in Precip, half hour value set to zero'
!!MS$                      ENDIF
!!MS$                      infn(i,3)=lrec(col,49)
!!MS$
!!MS$                      IF (MOD(i,48)==0) THEN
!!MS$                         dtmax(n,j) = MAXVAL(infn(:,1))
!!MS$                         dtmin(n,j) = MINVAL(infn(:,1))
!!MS$                         dswdown(n,j) = SUM(infn(:,2))/48.
!!MS$                         dprecip(n,j) = SUM(infn(:,3))
!!MS$                         i = 0
!!MS$                         j = j+1
!!MS$                      ENDIF
!!MS$
!!MS$                      i = i + 1                
!!MS$                   ENDDO
!!MS$                   DEALLOCATE(lrec)
!!MS$                ENDDO
!!MS$             ENDIF
!!MS$          ENDIF
       ENDDO
    else
       
       do n = 1, n_sites

          ! .. allocate temporary memory
          ALLOCATE (iload(ndays_in))
          infile%name= TRIM(in_fluxdir)//'HadCM3_A1_daily_1979_2039_'//TRIM(site_names(n))//'.nc'
          CALL ncopen(infile)
          ! get daily precipitation
          ncprecip%name='precip'
          CALL ncread (infile,ncprecip,iload)
          dprecip(n,:) = iload(1+nds0:sdays+nds0)
          nctmax%name='tmax'
          CALL ncread (infile,nctmax,iload)
          dtmax(n,:) = iload(1+nds0:sdays+nds0)
          nctmin%name='tmin'
          CALL ncread (infile,nctmin,iload)
          dtmin(n,:) = iload(1+nds0:sdays+nds0)
          ncswdown%name='swdown'
          CALL ncread (infile,ncswdown,iload)
          dswdown(n,:) = iload(1+nds0:sdays+nds0)
          CALL ncclose(infile)
          
          DEALLOCATE (iload)
       enddo
    endif
    
  END SUBROUTINE get_local_climate


  SUBROUTINE get_local_climate_subday

  !----------------------------------------------------------------------
  !    READS SUB-DAILY FORCING FOR SITE LEVEL RUNS FROM NETCDF FILES
  !---------------------------------------------------------------------- 

    ! .. Use Statements ..
    USE mo_netcdf
    USE mo_constants
    USE mo_namelist, ONLY : year0_site, year1_site, site_file_lai
    USE mo_calendar

    IMPLICIT NONE

! .. Local Scalars ..
    TYPE(ncfile) :: infile
    TYPE(ncvar) :: ncta, ncsw, ncea, ncps, nclw
    INTEGER :: i, j, n, ny, endyr

! .. Local Arrays ..
    REAL, ALLOCATABLE, DIMENSION(:,:) :: iload
    INTEGER, DIMENSION(2) :: ncstart

! .. Local Characters ..
    CHARACTER(len=80) :: site_file, header
    CHARACTER(len=80) :: in_fluxdir

    in_fluxdir = "input/eddy_sites/"
    ny = year1_site - year0_site +1

    ! READ FROM NETCDF FILE
    ! .. allocate temporary memory
    ALLOCATE (iload(n_sites,sdays*tspd))
    !ncstart=(/1,1/)

    ! Air temperature
    infile%name=TRIM(in_fluxdir)//'forcing-hourly-ta.US-NR1.2015.nc'
    CALL ncopen(infile)
    ncta%name='ta'
    !CALL ncread (infile,ncta,iload,start=ncstart)
    CALL ncread (infile,ncta,iload)
    CALL ncclose(infile)
    n=0

    DO i = 1,n_sites
       hta(:,i) = iload(i,:)
    END DO

    ! Shortwave down radiation
    infile%name=TRIM(in_fluxdir)//'forcing-hourly-sw.US-NR1.2015.nc'
    CALL ncopen(infile)
    ncsw%name='sw'
    CALL ncread (infile,ncsw,iload)
    CALL ncclose(infile)
    n=0

    DO i = 1,n_sites
       hswdown(:,i) = iload(i,:)
    END DO

    ! Longwave down radiation
    infile%name=TRIM(in_fluxdir)//'forcing-hourly-lw.US-NR1.2015.nc'
    CALL ncopen(infile)
    nclw%name='lw'
    CALL ncread (infile,nclw,iload)
    CALL ncclose(infile)
    n=0

    DO i = 1,n_sites
       hlwdown(:,i) = iload(i,:)
    END DO

    ! Actual vapor pressure
    infile%name=TRIM(in_fluxdir)//'forcing-hourly-ea.US-NR1.2015.nc'
    CALL ncopen(infile)
    ncea%name='ea'
    CALL ncread (infile,ncea,iload)
    CALL ncclose(infile)
    n=0

    DO i = 1,n_sites
       hea(:,i) = iload(i,:)
    END DO

    ! Air pressure
    infile%name=TRIM(in_fluxdir)//'forcing-hourly-ps.US-NR1.2015.nc'
    CALL ncopen(infile)
    ncps%name='ps'
    CALL ncread (infile,ncps,iload)
    CALL ncclose(infile)
    n=0

    DO i = 1,n_sites
       hpair(:,i) = iload(i,:)
    END DO


    DEALLOCATE (iload)

    END SUBROUTINE

  SUBROUTINE climsubdayforce (n_sites, iday0, iday1, inho)

  !----------------------------------------------------------------------
  !   GETS SUB-DAILY FORCING FOR SITE LEVEL RUNS FOR CURRENT TIME-STEP
  !----------------------------------------------------------------------
  !     CALCULATE THE AVERAGE FORCING OF THE CURRENT TIME-STEP (inho)
  !     ACROSS THE CURRENT RUN PERIOD.
  !     E.G. IF CURRENT HOUR TO RUN IS 10AM AND DAYINT IS 5 DAYS,
  !     TAKE THE AVERAGE FORCING OF ALL 10AM VALUES OVER THE 5-DAY 
  !     PERIOD. IF DAYINT=1 DAY, THEN JUST GET FORCING FOR THAT HOUR.
  !-----------------------------------------------------------------

! .. Use Statements ..
    USE mo_constants, ONLY: tspd

    IMPLICIT NONE

! .. Arguments
    INTEGER, INTENT(in) :: n_sites, iday0, iday1, inho

! .. Local scalars ..
    ! indexes for forcing arrays for current dayint run period
    INTEGER :: isubday0, isubday1
    INTEGER :: i, j, it
    REAL    :: ntsteps    ! number of days in current dayint period
    REAL    :: xlwdown, xswdown, xta, xea, xpair  ! sum of current time steps climate variables

    isubday0=(iday0-1) * tspd + 1
    isubday1=(iday0-1) * tspd + (iday1-iday0+1)*tspd

    ntsteps = REAL(iday1 - iday0 + 1)

    ! Loop over sites (ng)
    DO i = 1,n_sites
       xlwdown = 0.
       xswdown = 0.
       xta     = 0.
       xea     = 0.
       xpair   = 0.
       ! Loop over days in current dayint period 
       ! and sum them up into x* variables
       DO j = 1,ntsteps
          ! Get index of all curent hours
          it = (24*(iday0-1)+inho) + 24*(j-1)
          xlwdown = xlwdown + hlwdown(it,i)
          xswdown = xswdown + hswdown(it,i)
          xta     = xta + hta(it,i)
          xea     = xea + hea(it,i)
          xpair   = xpair + hpair(it,i)
       END DO
       ! Calculate average of x* variable
       ! i.e. divide sum by number of values
       ctlwdown(i) = xlwdown/ntsteps
       ctswdown(i) = xswdown/ntsteps
       ctta(i)     = xta/ntsteps
       ctea(i)     = xea/ntsteps
       ctpair(i)   = xpair/ntsteps
    END DO

  END SUBROUTINE

END MODULE mo_climate
