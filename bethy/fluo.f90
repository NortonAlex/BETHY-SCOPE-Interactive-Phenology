MODULE fluo

IMPLICIT NONE 

  REAL,DIMENSION(:),ALLOCATABLE               :: dxyp

  REAL,ALLOCATABLE,DIMENSION(:)               :: zgppfluo    !ANorton. SCOPE-GPP over vp 

  ! FS data  
  REAL                                        :: afs    ! Slope between simulated fluo and gpp 
  REAL, ALLOCATABLE, DIMENSION(:), SAVE       :: fs    !fluorescence  simulated 
  REAL, ALLOCATABLE, DIMENSION(:), SAVE       :: fs_obs !fluorescence  data 
  REAL, ALLOCATABLE, DIMENSION(:), SAVE       :: fs_unc !fluorescence  data uncertainty 
  REAL, ALLOCATABLE, DIMENSION(:), SAVE       :: date_fs
  INTEGER                                     :: nobs_fs
  INTEGER, ALLOCATABLE, DIMENSION(:), SAVE    :: indm_fs
  INTEGER, ALLOCATABLE, DIMENSION(:), SAVE    :: indy_fs
  INTEGER, ALLOCATABLE, DIMENSION(:), SAVE    :: indng_fs
  INTEGER, ALLOCATABLE, DIMENSION(:,:),SAVE   :: nbfs
  INTEGER, ALLOCATABLE, DIMENSION(:,:,:)      :: ngfs

 ! Radiation data 
  REAL, ALLOCATABLE, DIMENSION(:,:),SAVE      :: lwd   ! Longwave down radiation 
  REAL, ALLOCATABLE, DIMENSION(:,:),SAVE      :: swd   ! Shortwave down radiation 

CONTAINS

  SUBROUTINE   init_fluobethy(year0, year1)
    USE constants_pjr
    USE mo_constants
    USE mo_namelist, ONLY : nrun
    USE mo_grid, ONLY : ng, vp

    IMPLICIT NONE

    INTEGER, INTENT(in)                     :: year0,year1
    INTEGER                                 :: firstdat, lastdat

    ! local variables   
    INTEGER                                 :: i ! index variables
    
    INTEGER, PARAMETER                      :: inunit=2

    REAL                                    :: dlon,dlat,lat,lata,lato
 
     firstdat=year0 
     lastdat=year1 
    
    !calculate grid cell area
    ALLOCATE (dxyp(n_lats))
    dlon = 360./n_lons
    lata=90.
    DO i=1,n_lats       
       IF (i==1 .OR. i==n_lats) THEN
          dlat = 180./((n_lats-2)*2+2)
       ELSE
          dlat = 2*180./((n_lats-2)*2+2)
       ENDIF       
       lato=lata-dlat
       lat=(lato+lata)/2.
       dxyp(i)=510108933.5e6* dlat / 360. * COS (lat * pi / 180) &
            & * SIN (dlon / 2. * pi / 180)
       lata=lato
    ENDDO

    ! ANorton. SCOPE-GPP over vegetation points
    ALLOCATE ( zgppfluo (vp) )
    zgppfluo = 0

  ! This array allows to have the number of bethy pixels for which observed fs
  ! are available 
    ALLOCATE(nbfs(nrun,outt))
    ALLOCATE(ngfs(nrun,outt,ng))
    nbfs =  0
    ngfs =  0 

 ! We read the long wave radiation and also the shortwave radiation here  
 ! We read here since the long wave radiation is not used by BEHTY .... 
   CALL read_radiation 

! We read the data from fs 
!    call read_obs_fs(lastdat,firstdat)
 !   print*, ' nobs_fs ', nobs_fs 
    
!    CALL fs_allocate


    RETURN
  END SUBROUTINE init_fluobethy
  

  SUBROUTINE fs_allocate
    IMPLICIT NONE 
    ALLOCATE(fs(nobs_fs))
    fs=0.0
  END SUBROUTINE fs_allocate

  SUBROUTINE fs_deallocate
  IMPLICIT NONE 
  DEALLOCATE(fs)
  END SUBROUTINE fs_deallocate

! Routine that reads the radiation data (longwave down and shortwave down)
SUBROUTINE read_radiation
    USE mo_constants
    USE mo_netcdf
    USE mo_namelist
    USE mo_grid

    IMPLICIT NONE

    TYPE (ncfile)                         :: infile_lw,infile_sw
    TYPE (ncvar)                          :: nc_lw,nc_sw
    TYPE (ncvar)                          :: nc_lon,nc_lat 
    INTEGER, DIMENSION( NF_MAX_DIMS)      :: theLens
    REAL, ALLOCATABLE,DIMENSION(:,:,:)    :: field1, field2
    REAL, ALLOCATABLE,DIMENSION(:)        :: lonr,latr 
    
    INTEGER                               :: i, ilon,ilat,ing 

! Test if the lon/lat of BETHY grid cells are here 
!    print*, ' lat lon ', size(lat,1), size(lon,1)
  
   ! Long wave radiation
    infile_lw%name = 'input/climate/lwdown_mean_79-05.nc'
    CALL ncopen(infile_lw)
    nc_lw%name = 'lwdown'
    CALL ncvarinfo(infile_lw, nc_lw, dimlens=thelens)
    ALLOCATE( field1( thelens(1), thelens(2), thelens(3)) )
    CALL ncread( infile_lw, nc_lw, field1)


  ! Short wave radiation (mean value)
    infile_sw%name = 'input/climate/shortdown_mean.nc'
    CALL ncopen( infile_sw)
    nc_sw%name = 'shortwave'
    nc_lon%name = 'lon'
    nc_lat%name = 'lat'
    CALL ncvarinfo(infile_sw, nc_sw, dimlens=thelens)
    ALLOCATE( field2( thelens(1), thelens(2), thelens(3)) )
    ALLOCATE( lonr(thelens(1)))
    ALLOCATE( latr( thelens(2)))
    CALL ncread( infile_sw, nc_sw, field2)
    CALL ncread( infile_sw, nc_lon, lonr)
    CALL ncread( infile_sw, nc_lat, latr)

    print*, ' lon ', minval(lonr), maxval(lonr)
    print*, ' lat ', minval(latr), maxval(latr)
    
! Monthly fields of longwave down and shortwave down
    ALLOCATE( lwd(ng,12))
    ALLOCATE( swd(ng,12))

   ! We condider only the data for BETHY grid cells 
      DO ilat = 1,thelens(2)
        DO ilon = 1,thelens(1) 
          DO ing = 1,ng
             IF (lon(ing)==lonr(ilon).AND.(lat(ing)==latr(ilat))) THEN 
              lwd(ing,:)= field1(ilon,ilat,:)
              swd(ing,:)= field2(ilon,ilat,:)
            !if (ing==31)  print*, ' RADLON ',ing, lon(ing),lonr(ilon)
            !if (ing==31)  print*, ' RADLAT ',ing, lat(ing),latr(ilat)
            !if (ing == 31)  print*, ' VAL ',ing, field1(ilon,ilat,6),field2(ilon,ilat,6)
             !print*, ' SWD ',ing, minval(field2(ilon,ilat,:)),maxval(field2(ilon,ilat,:))
             ENDIF 
          END DO 
        END DO 
      END DO 

  !   print*,' lwd ',  minval(lwd), maxval(lwd)
  !   print*, 'swd ', minval(swd), maxval(swd)

  DEALLOCATE(field1,field2)
  DEALLOCATE(lonr,latr)

END SUBROUTINE read_radiation 


!  ---------------------------------------------------------------------------
! E Koffi May 2012 
! This reads fluorescence data and selected ony data that cover the 
! period of the run 
! ---------------------------------------------------------------------------

 SUBROUTINE  read_obs_fs(lastdat,firstdat)
    USE mo_constants
    USE mo_netcdf
    USE mo_namelist 
    USE mo_grid 
    
    IMPLICIT NONE 

    INTEGER,INTENT(in)                 :: firstdat, lastdat

    ! local variables
     REAL, ALLOCATABLE, DIMENSION(:,:) :: obs_fs !fluorescence  data
     REAL, ALLOCATABLE, DIMENSION(:,:) :: unc_fs !fluorescence  data uncertainty
     REAL, ALLOCATABLE, DIMENSION(:)   :: time_fs 

     REAL                              :: unc_mod 
     REAL                              :: date1
     INTEGER                           :: yrci,yrce,yr,imonth 
     INTEGER                           :: iyr,nyr,nm,m
     INTEGER                           :: ng_bethy, ig 
     INTEGER                           :: iobs_fs
     INTEGER                           :: iyear,im,ing 
     TYPE (ncfile)                     :: infile
     TYPE (ncvar)                      :: nc_fs,nc_fs_sig,nc_fs_time
     INTEGER, DIMENSION( NF_MAX_DIMS)  :: theLens

     INTEGER                           :: nb
     INTEGER                           :: j,jp,k,jl,np,img  
     
!     print*, ' ifs ', ifs 

     if (ifs /=0 ) then 
    ! Uncertainties in BETHY when simulating fs ... To evalutate this 
     unc_mod =  0.2 

!!! FS data 
    print*, 'fsfile: ', fsfile 
    infile %name = fsfile 
    print*, infile%name 
    CALL ncopen( infile)
    nc_fs %name = 'fs'
    CALL ncvarinfo( infile, nc_fs, dimlens=thelens)
    ALLOCATE( obs_fs( thelens(1), thelens(2)) )
    CALL ncread( infile, nc_fs, obs_fs)

    nm=thelens(1)
    ng_bethy=thelens(2)

    nc_fs_sig %name = 'fs_sig'
    CALL ncvarinfo( infile, nc_fs_sig, dimlens=thelens)
    ALLOCATE( unc_fs( thelens(1), thelens(2) ) )
    CALL ncread( infile, nc_fs_sig, unc_fs)

    nc_fs_time%name='date'
    CALL ncvarinfo( infile, nc_fs_time, dimlens=thelens)
    ALLOCATE( time_fs( thelens(1)))
    CALL ncread( infile, nc_fs_time, time_fs)
    CALL ncclose( infile) 


   ! print*, time_fs(1), int(time_fs(1)), firstdat, lastdat
   ! print*, time_fs(thelens(1)), int(time_fs(thelens(1)))

   ! print*,'min max obs_fs ', minval(obs_fs), maxval(obs_fs)
   ! print*,'min max unc_fs ',  minval(unc_fs), maxval(unc_fs)


    ! Search for the number of fs observtions that are in the selected common
    ! period 
      m=0
      iobs_fs=0

    !  print*, ' nm ', nm 
    !  print*, ' ng ', ng 
      ! start loop over all files with observations
       DO imonth= 1,nm

          date1=time_fs(imonth) 
!          print*, ' firsdtat date1  ', firstdat, date1, lastdat+1
!          print*, ' iyear imonth ', imonth 
         
          if ((date1>=firstdat).and.(date1 <lastdat+1)) then 
          DO  ig = 1,ng_bethy 
              ! We consider only the observations when the test > 0 Otherwise we
              ! consider all the BETHY pixels 
             !if (obs_fs(imonth,ig) > 0)  then 
              if (obs_fs(imonth,ig) > -9999)  then 
               iobs_fs = iobs_fs +1
              endif 
          END DO 
         endif 

     END DO 

    nobs_fs=iobs_fs 
    !print*, 'nobs_fs', nobs_fs 

    ! We put the fs data for the selected period 
    ALLOCATE( fs_obs(nobs_fs) )
    ALLOCATE( fs_unc(nobs_fs) )
    ALLOCATE( date_fs(nobs_fs) )
    ALLOCATE( indy_fs(nobs_fs) )
    ALLOCATE( indm_fs(nobs_fs) )
    ALLOCATE( indng_fs(nobs_fs) )

! We fill in the data: observations with uncertainties 
      m=0
      iobs_fs=0

      ! start loop over all files with observations
       DO imonth= 1,nm
          date1 = time_fs(imonth)
          im    = 1+ (date1 - int(date1)) *12 
          iyear = int(date1) - firstdat +1
          
          
          img   = im + (iyear -1)*12   
          
          nb = 0

          if ((date1>=firstdat).and.(date1 < lastdat+1)) then
!          print*, ' firsdtat date1  ', firstdat, date1,  lastdat+1
          !print*, ' iyear imonth ', iyear, im 
          DO  ig = 1,ng_bethy
              ! We consider only the observations when the test > 0 Otherwise we
              ! consider all the BETHY pixels 
              !if (obs_fs(imonth,ig) > 0)  then
              if (obs_fs(imonth,ig) > -9999)  then

                       iobs_fs    = iobs_fs +1
                       nb         = nb +1
 
                fs_obs(iobs_fs)   = obs_fs(imonth,ig)
                fs_unc(iobs_fs)   = unc_fs(imonth,ig) + unc_mod
               indy_fs(iobs_fs)   = iyear  
               indm_fs(iobs_fs)   = im 
              indng_fs(iobs_fs)   = ig 
              date_fs(iobs_fs)    = date1 

           ngfs(iyear,img,nb)     = ig 
              

              endif
          END DO
                   nbfs(iyear,img) = nb

    !      print*,'PIXEL NB ', iyear, imonth, im,img,nbfs(iyear,img)
         endif
 

     END DO

     !  do  iyear = 1, 1
     !   do imonth = 1,12  
     !    np =  0
     ! do j =1, nbfs(iyear,imonth) 
     !    jp =  ngfs(iyear,imonth,j)
     !    do k =1,3 
     !    jl = gridvp(jp,k)
     !    if(jl /=0) np = np +1 
     !    if(jl /=0) print*,' GGG ',jl,frac(jl),vg_nv(jl) 
     !    end do 
     !  end do  
     !     print*, ' VERIF iyear imonth np ', iyear, imonth,np 
     ! end do 
     ! end do 
     
    else 
    nobs_fs=1 

    ALLOCATE( fs_obs(nobs_fs) )
    ALLOCATE( fs_unc(nobs_fs) )
    ALLOCATE( indy_fs(nobs_fs) )
    ALLOCATE( indm_fs(nobs_fs) )
    ALLOCATE( indng_fs(nobs_fs) )
    ALLOCATE( date_fs(nobs_fs) )

    fs_obs   = 0.
    fs_unc   = 999. 
    indy_fs  = 1 
    indm_fs  = 1 
    indng_fs = 1
    date_fs  = -999.
   endif 

   print*, 'nobs_fs :  ', nobs_fs 

    RETURN 

  END SUBROUTINE read_obs_fs 


SUBROUTINE fluorescence(iyear,imonth,iday,ihour,iday0,iday1,irrin,par,&
                  & temp,pres,ea0,Cca,COa,zlai,&
                  & jmf,vm,EC,EO,EV,ER,EK,kc0,ko0,vomf,rdf,&
                  & rfluo,rgppfluo,PAR_scope,PAR_scope_cab,&
                  & rfluo_diurnal,rgppfluo_diurnal,&
                  & rlai_diurnal,rapar_diurnal,raparcab_diurnal,&
                  & rpar_diurnal,rswdown_diurnal,rta_diurnal,&
                  & rfsyieldu_toc_diurnal,rfsyieldh_toc_diurnal,&
                  & rpyieldu_toc_diurnal,rpyieldh_toc_diurnal,&
                  & rnpqyieldu_toc_diurnal,rnpqyieldh_toc_diurnal)

!CALL fluorescence(ryear,rmonth,iday,its,irrin,par,&
!                    & temp,p,ea0,ca,OX,zlai,fracs, &
!                    & jmf,vm,EC,EO,EV,ER,EK,kc0,ko0,shwd,&
!                    & rfluo,rgppfluo,PAR_scope,PAR_scope_cab)
!

! Bethy  
USE mo_constants
USE mo_namelist
USE mo_grid

! Fluo 
USE fluo_param, ONLY : aggreg,fluspect,leafangles,jatmos_file,atmos_file,spectral_nreg,spectral_start,spectral_end, &
                      & spectral_res,tts,tto,Rin,Rli,Ta,pa,ea,LAI,Vcmo,Oa,Cab,option, &
                      & leafbio,Jmo,Cdm,Cw,Csm,N,fqe1,fqe2,avovc,ardvc, &
                      & rho_thermal,tau_thermal,nlazi,nli,nl,psi,Ps,Po,Pso, &
                      & km,Kext,Esun_,Esky_,P,fEsuno,fEskyo,fEsunt,fEskyt,Eplu_,Emin_, &
                      & Lo_,Eout_,Eouto,Eoutt,Rnhs,Rnus,Rnhc,Rnuc,Pnhc,Pnuc,Pnhc_Cab, &
                      & Pnuc_Cab,Fc,tempcor,LoF,Fhem,Fiprof,ifreq_sat,ial,wlf,wle,nwl,nwlP, &
                      & LIDFa,LIDFb,q 
USE fluo_func 
USE mo_rtmo, ONLY : rtmo 
USE chemical, ONLY : biochemical_faq, biochemical
USE mo_rtmf, ONLY : rtmf  
USE mo_vegetation, ONLY : Chl,Cdm_arr,Csm_arr,LIDFa_arr,LIDFb_arr,hc_arr,leafwidth_arr
                      
USE mo_config, ONLY : vps,block_vps

!% Input:
!% Esun_     [W m-2 um]          Vector of incoming shortwave radiation (=<2.5 um)
!% Esky_     [W m-2 um]          Vector of incoming longwave radiation (>2.5 um)
!% tts       [deg]               solar zenith angle
!% Ca        [umol m-3]          Atmospheric CO2 concentration at measurement height
!% Oa        [mmol m-3]          Atmospheric O2 concentration at measurement height
!% Ta        [oC]                Air temperature at measurement height
!% ea        [hPa or mbar]       Atmospheric vapour pressure at measurement height 
!% t         [days]              Time of the year in decimal days
!% pa        [hPa]               Atmospheric vapour pressue
!% u         [m s-1]             Wind speed at measurement height
!% GAM       [J K-1 m-2 s-1/2]   Soil thermal inertia
!% LAI      [m-2 m-2]           One-sided leaf area index
!% Vcmo      [umol m-2 s-1]      maximum carboxylation capacity
!% Jmo       [umol m-2 s-1]      maximum electron transport rate
!% lam       []                  Cowan's stomatal parameter
!%
!% Globals:
!% MH2O      [g mol-1]           molecular mass of water (18)
!% Mair      [g mol-1]           molecular mass of dry air (29.96)
!% rwc       [s m-1]             Within canopy aerodynamic resistance
!% kappa     []                  Von Karman coefficient (0.4)
!% omega     [rad s-1]           frequency of diurnal cycle (2*pi/86400)
!% Tslta     [C]                 seasonal mean temperature of the soil
!% rhoa      [kg m-3]            specific mass of the air (kg m-3)
!% cp        [J kg-1 K-1]        specific heat of dry air
!% g         [m s-2]             gravity acceleration
!% Type      []                  photosynthetic pathway (1: C4, ~1: C3)
!% tto       [deg]               viewing angle
!% psi       [deg]               Observer Azimuth angle
!% nl        []                  number of leaf layers
!%                               given to the new update <0:1]
!%
!% Output:
!% Actot     [umol m-2 s-1]      canopy net CO2 uptake
!% Pntot     [umol m-2 s-1]      total net PAR
!% Lo        [W m-2]             total outgoing radiation in viewing direction
!% Loo       [W m-2]             total outgoing shortwave radiation in viewing direction (<=2.5 um)
!% Lout_     [W m-2 um-1]        outgoing hemispherical radiation per wavelength interval
!% Lo_       [W m-2 um-1]        outgoing radiation per wavelength interval in viewing direction
!% Fc        []                  fraction of sunlit leaves per leaf layer
!
!%% 1. Solar and sky radiation
!% SAIL model for solar and sky radiation, without emission by canopy and
!% soil itself
!% subscript '0' indicates: shaded leaf/soil, '1' indicates: sunlit
!% Rn is net radiation (W m-2), E a spectrally integrated flux in W m-2,
!% and E_ a flux in W m-2 um-1

IMPLICIT NONE 
! Input variables
! Date 
INTEGER            , INTENT(in)           :: iyear,imonth,iday,ihour
INTEGER            , INTENT(in)           :: iday0,iday1

! Radiation short wave and PAR computed from bethy  
REAL, DIMENSION(ng), INTENT(in)           :: irrin,par       

!REAL, DIMENSION(24,ng), INTENT(in)        :: shwd    ! solar radiation as a function of the time
!REAL, DIMENSION(24,ng), INTENT(in)        :: temp_, pres_,vap_pres 


! Meteo and chemistry at the boundary of the leaves  
REAL, DIMENSION(ng), INTENT(in)           :: temp, pres,ea0        ! Climate data
REAL, INTENT(in)                          :: Cca     ! Concentration of CO2 at the leave interface (i.e., air concentration )
REAL, INTENT(in)                          :: COa     ! Concentration of O2 at the leave interface (i.e., air concentration )

! Biochemical 
REAL, DIMENSION(vp), INTENT(in)            :: EC, EO, EV, ER, EK  ! Energy activations of each  PFT data
REAL, DIMENSION(vp), INTENT(in)            :: kc0, ko0
REAL, DIMENSION(vp), INTENT(in)            :: vm                  ! Vcmax for the various PFTs
REAL, DIMENSION(vp), INTENT(in)            :: jmf                 ! jmf --> ajv
REAL, DIMENSION(vp), INTENT(in)            :: vomf                ! Vomax as fraction of Vcmax
REAL, DIMENSION(vp), INTENT(in)            :: rdf                 ! Rd as fraction of Vcmax

! PFT characteristics
REAL, DIMENSION(vp), INTENT(in)            :: zlai       ! LAI data
!REAL, DIMENSION(vp), INTENT(in)            :: fracs      ! fraction of the PFT in the grid cell


! Chlorophyl concentration 
!REAL,INTENT(IN)                             :: Cab,Cdm,Cw,Csm,NN

!Output variable from fluspect
REAL, ALLOCATABLE, DIMENSION(:,:)            :: MfI, MbI, MfII, MbII
REAL, ALLOCATABLE, DIMENSION(:)              :: rho, tau, rs
REAL, ALLOCATABLE, DIMENSION(:)              :: kChlrel

!Output variable from leafangles
REAL, ALLOCATABLE, DIMENSION(:)              :: lidf

! Output variable from fluoresence routine
!REAL, DIMENSION(vp), INTENT(out)            :: zfluo, zgppfluo
REAL, DIMENSION(0:nrun,outt,ng), INTENT(out) :: PAR_scope
REAL, DIMENSION(0:nrun,outt,ng), INTENT(out) :: PAR_scope_cab
REAL, DIMENSION(0:nrun,outt,ng), INTENT(out) :: rfluo,rgppfluo
REAL, DIMENSION(vp)                          :: daygpp, dayfluo
REAL, DIMENSION(0:nrun,365,tspd,vp), INTENT(out) :: rfluo_diurnal,rgppfluo_diurnal
REAL, DIMENSION(0:nrun,365,tspd,vp), INTENT(out) :: rlai_diurnal,rapar_diurnal,raparcab_diurnal
REAL, DIMENSION(0:nrun,365,tspd,vp), INTENT(out) :: rpar_diurnal,rswdown_diurnal
REAL, DIMENSION(0:nrun,365,tspd,vp), INTENT(out) :: rta_diurnal
! Output for top-of-canopy PSII quantum yields
REAL, DIMENSION(0:nrun,365,tspd,vp), INTENT(out) :: rfsyieldu_toc_diurnal,rfsyieldh_toc_diurnal
REAL, DIMENSION(0:nrun,365,tspd,vp), INTENT(out) :: rpyieldu_toc_diurnal,rpyieldh_toc_diurnal
REAL, DIMENSION(0:nrun,365,tspd,vp), INTENT(out) :: rnpqyieldu_toc_diurnal,rnpqyieldh_toc_diurnal

! Local fields 
REAL, DIMENSION(2)                           :: Fs_mat                      ! matrix containing values for probabilities of viewing sunlit/shaded leaves/soil 

REAL, ALLOCATABLE,DIMENSION(:)               :: Agh,Ah,rcwh,Fh,Cih,Cch       ! brown leaves
REAL, ALLOCATABLE,DIMENSION(:)               :: Agu,Au,rcwu,Fu,Ciu,Ccu       ! green leaves

REAL, ALLOCATABLE,DIMENSION(:)               :: phi_p_h,phi_fs_h,phi_npq_h   ! shaded
REAL, ALLOCATABLE,DIMENSION(:)               :: phi_p_u,phi_fs_u,phi_npq_u   ! sunlit

REAL, ALLOCATABLE,DIMENSION(:)               :: A0,Ag0,rcw0,F0a,W0,F0            ! brown leaves
REAL, ALLOCATABLE,DIMENSION(:)               :: A1,Ag1,rcw1,F1a,W1,F1            ! green leaves


REAL, ALLOCATABLE,DIMENSION(:)               :: Tch       ! brown leaves
REAL, ALLOCATABLE,DIMENSION(:)               :: Tcu       ! green leaves
REAL, ALLOCATABLE,DIMENSION(:)               :: Ccs
REAL, ALLOCATABLE,DIMENSION(:)               :: Fout,Fout_profile

REAL                                         :: Actot,Agtot 
REAL                                         :: Actot_faq,Agtot_faq 
REAL                                         :: Pntot, Pntot_Cab   !% Radiative fluxes, spectrally integrated
REAL                                         :: fs_yieldu_toc,fs_yieldh_toc !Fluorescence quantum yield (sunlit and shaded) at top-of-canopy
REAL                                         :: p_yieldu_toc,p_yieldh_toc   !Photochemical quantum yield (sunlit and shaded) at top-of-canopy
REAL                                         :: npq_yieldu_toc,npq_yieldh_toc !NPQ (basal+regulated) quantum yield (sunlit and shaded) at top-of-canopy
REAL                                         :: LoF_jl
REAL                                         :: Cc
REAL                                         :: Ccc,Occ

INTEGER                                      :: jl,jj
INTEGER                                      :: gcmethod

CHARACTER(len=30)                            :: type_integration
INTEGER                                      :: nlh,nlu,np

REAL                                         :: year  ! Year
REAL                                         :: doy   ! Day of the year
REAL                                         :: tm   ! the time of the day when the sun reaches its highest angle
REAL                                         :: t     ! Local time of the day in hours
REAL                                         :: ttsR  ! Sun zenith angle in rad

REAL                                         :: Omega_g  !% slope azimuth angle (deg)
REAL                                         :: Fi_gm   !% slope of the surface (deg)

REAL                                         :: Long, Lati

REAL                                         :: ttsR0, tts0
INTEGER                                      :: n0, n1, ns,kt

INTEGER                                      :: pft 
INTEGER                                      :: j, k, im 
INTEGER                                      :: ok
REAL                                         :: sum_frac_tot,sum_frac_ok
REAL                                         :: sum_fluo
REAL                                         :: vmax
REAL                                         :: frac1 


REAL                                         :: frac_max,v_max
INTEGER                                      :: pft_dominant_cell, jl_max
REAL                                         :: LAI_max   ! Maximum LAI in the grid cell

INTEGER                                      :: faq 
INTEGER                                      :: jd 
REAL                                         :: hm(24)

INTEGER                                      :: it   ! time index for diurnal output arrays

!INTEGER, DIMENSION (12)                      ::  rdays
!DATA rdays /31, 28, 31, 30, 31, 30, 31, 31, 30, 31,30, 31/

! When using Faquahar model to compute GPP and Fluo at leaf level 
gcmethod = 1      ! method gcmethod = 1 (Cowan, 1997?), = 0 (Leuning)
     faq = 0      ! 1: In addition to Collatz model, we use Faquahar model to compute GPP  0:  Only Collatz 



! To refine this 
!               t = 12.                    ! This is the local time at the selected pixel
!              im =  mod(imonth,12)
!if (im == 0 ) im =  12
!             doy =  sum(rdays(1:im)) +15.   ! We consider the middle of the month

              doy = INT((iday0+iday1)/2)        ! Computes mid-point doy of simulated period (i.e. period over which forcing is averaged)

! We have a pb with the hour to fix this 
CALL pb_hour_bethy(hm)

! Read in modtran atmosphere transmittance files to variable array
!CALL read_modtran_files

!$OMP PARALLEL DO DEFAULT(FIRSTPRIVATE) SHARED(vp,gridp,ok,lat,lon,imonth) &
!$OMP& SHARED(spectral_nreg,spectral_start,spectral_end,spectral_res,hm,doy) & 
!$OMP& SHARED(irrin,lwd,temp,pres,ea0,zlai,vg_nv,vm,frac,Cca,COa) &
!$OMP& SHARED(jmf,Cw,N,fqe1,fqe2) &
!$OMP& SHARED(nl,nli,nlazi,faq,EC,EO,EV,ER,EK,kc0,ko0,gcmethod,rfluo,iyear) &
!$OMP& SHARED(rgppfluo,zgppfluo,PAR_scope,PAR_scope_cab,ifreq_sat,psi,tto) & 
!$OMP& SHARED(nwl,wlf,nwlP,Chl,Cdm_arr,Csm_arr,LIDFa_arr,LIDFb_arr,hc_arr,leafwidth_arr) &
!$OMP& SHARED(ihour,iday,rfluo_diurnal,rgppfluo_diurnal,rlai_diurnal,rapar_diurnal,raparcab_diurnal) & 
!$OMP& SHARED(rpar_diurnal,rswdown_diurnal,rta_diurnal) &
!$OMP& SHARED(rfsyieldu_toc_diurnal,rfsyieldh_toc_diurnal) &
!$OMP& SHARED(rpyieldu_toc_diurnal,rpyieldh_toc_diurnal) &
!$OMP& SHARED(rnpqyieldu_toc_diurnal,rnpqyieldh_toc_diurnal) &
!$OMP& SHARED(vps,block_vps,vomf,rdf) &
!$OMP& PRIVATE(MfI,MbI,MfII,MbII,rho,tau,rs,kChlrel,lidf,Agh,Ah,rcwh,Fh,A0,Ag0,rcw0) &
!$OMP& PRIVATE(F0a,F0,W0,Cih,Cch,Tch,Fout,Fout_profile,Agu,Au,rcwu,Fu,A1,Ag1,rcw1,F1a,F1,W1) &
!$OMP& PRIVATE(Ciu,Ccu,Tcu,jl,jj,j,P) &
!$OMP& PRIVATE(phi_p_h,phi_fs_h,phi_npq_h,phi_p_u,phi_fs_u,phi_npq_u)

  DO j = 1,vps
      jl = block_vps(j)   ! jl = vp  
        jj=gridp(jl)      ! jj = gp

 ! do jj = 1, ng

! We verif the frac of the FTs over the selected grid cells. The maximum has to
! be 1
      sum_frac_tot = 0.
      sum_frac_ok  = 0.
      sum_fluo     = 0.

                ok = 0

! Latitude/longitude of the grid cell
             Lati  =  lat(jj)             ! Latitude of the BETHY pixel
             Long  =  lon(jj)             ! Longitude of the BETHY pixel

! Search for the right modtran file 
CALL modtran_ifile (imonth, Lati, jatmos_file) 

! We select the appropriate spectrum  ...  for the mid-latitudes , we have two
! files (winter and summer) and one file. We select the file according to the
! latitude 
CALL aggreg (jatmos_file,spectral_nreg,spectral_start,spectral_end,spectral_res)


! Diurnal variations. We only consider the data at 12 h, which correspond to
! BEHTY its of 24 
  !DO jd = 1, 24 
!  DO jd = 1, 1

!     t = hm(jd)*1. 
!     t = hm(ihour)*1.

! Removed a dependency of SCOPE sub-day timestep on BETHY looping var "its". 
!   BETHY var "inho" (fluor ihour argument) corresponds to: actual_time (in 24hr time) = inho+1
     t = MOD(ihour,24)+1     

! Computation of the sun zenith angle in radians
! We have the time at which the sun reaches its maximum
! We consider the appropriate radiation for the computation
 CALL  calczenithangle(doy,t,tm,Fi_gm,Omega_g,Long,Lati,ttsR)   !%sun zenith angle in rad
tts    = min(85.,ttsR/deg2rad)                                  ! %sun zenith angle in deg


! Radiation..... 
             Rin  = irrin(jj)  ! Shortwave radiation ... from BETHY shortwave at the

             !Rin = swd(jj,imonth)  ! Shortwave radiation. In case we consider the mean
                           !short wave radiation as performed by the routine red_radiation in this module
             Rli  = lwd(jj,imonth)  ! Longwave radiation. We consider the monthly mean
                           ! radiation as performed in the routine read_radiation

! Interface of the leaves  
              Ta  = temp(jj) 
              pa  = pres(jj)/100.  ! In hPa
              ea  = ea0(jj) /100.   ! In hPa 


! We initialise these variables for check
         frac_max = -999
          LAI_max = -999
pft_dominant_cell = -999
           jl_max = -999
            v_max = -999

! We go through up to the 3 PFTs that are in the selected BETHY grid cell
 !     do k =1,3

 !             jl = gridvp(jj,k)

               P = 0.
           Pntot = 0.
       Pntot_Cab = 0.
           Agtot = 0.
       Agtot_faq = 0.
           Actot = 0.
       Actot_faq = 0.
          LoF_jl = 0.
       fs_yieldu_toc = 0.
       fs_yieldh_toc = 0.
        p_yieldu_toc = 0.
        p_yieldh_toc = 0.
        npq_yieldu_toc = 0.
        npq_yieldh_toc = 0.
            
             LAI = -1.
             Cab = -1. 
            Vcmo = -1. 
             pft = -1
           frac1 = 0. 
           LIDFa = 0.
           LIDFb = 0.
           avovc = -1.
           ardvc = -1.

     ALLOCATE(MfI(size(wlf),size(wle)))
     ALLOCATE(MbI(size(wlf),size(wle)))
     ALLOCATE(MfII(size(wlf),size(wle)))
     ALLOCATE(MbII(size(wlf),size(wle)))
     ALLOCATE(rho(nwl))
     ALLOCATE(tau(nwl))
     ALLOCATE(rs(nwl))
     ALLOCATE(kChlrel(nwlP))

     ALLOCATE(lidf(nli))

! We consider only grid cells with PFTs having LAI greater than 0.
  IF (zlai(jl) >0.) THEN

               ok = ok +1

     sum_frac_tot = sum_frac_tot +frac(jl)

              LAI = zlai(jl)    ! LAI for the selected PFT in the grid cell
              pft = vg_nv(jl)
             Vcmo = vm(jl)*1.e6
            frac1 = frac(jl)
              Cab = Chl(jl)

            avovc = vomf(jl)    ! ratio of Vomax to Vcmax
            ardvc = rdf(jl)     ! ratio of Rd to Vcmax

! Concentration of CO2 and O2 in the atmosphere, i.e, at the boundary of the
! leaves

              Cc  = Cca*1e6
              Oa  = COa*1e3 

!   print*, ' VERIF ', Cc, Oa, Vcmo 

! -----------------------------------------------------
 ! Wne using biocchemical base on Faquhar for GPP 
            CCc  = Cca*1e6 *1E3*rhoa/Mair
            Occ  = Coa*1e3 *1E3 *rhoa/Mair

! Jmax : maximum electron transport rate table .. Maybe this needs to be refined
!Jmo             = 2.5*Vm(jl) * 1e6
            Jmo = jmf(jl)*Vcmo


! ----------------------------------------------------------

!  print*, ' VERIF ', Cc, Oa, CCc,Occ,Vcmo 

! We determine the LAI_max for the selected PFT over the selected grid cell for
! check reason
if (zlai(jl).gt.LAI_max) LAI_max = zlai(jl)


! option = 0: C3 plant    1: C4 plant
               option = 0       ! it is put to 0 for C3 plant
if (pft == 10) option = 1  ! C4 crop plant only
!if (pft == 13) option = 1  ! C3 crop plant only


! Concentration of the fluorescence of the leaves We assume this for the moment
!                   Cab = 10.
!      if (pft==1)  Cab = 40.
!      if (pft==2)  Cab = 15.
!      if (pft==3)  Cab = 15.
!      if (pft==13) Cab = 20.
!      if (pft==9)  Cab = 10.
!      if (pft==10) Cab = 5.


! First we consider the same Cab for all the PFTs 
!                 Cab  = 40. 


! Leaf fluorescence concentration  and related parameters... We can only
! consider Cab , but for the moment we leave as follows
           leafbio(1) = Cab
           leafbio(2) = Cdm_arr(jl)
           leafbio(3) = Cw
           leafbio(4) = Csm_arr(jl)
           leafbio(5) = N
           leafbio(6) = fqe1
           leafbio(7) = fqe2
           leafbio(8) = rho_thermal
           leafbio(9) = tau_thermal

!! Option for seasonal variation in parameters.
!! - we give these parameters seasonal variability based on a sine function with
!!   the maximum on the summer solstice and minimum on winter solstice.
!IF (pft.eq.1.or.pft.eq.3.or.pft.eq.5.or.pft.eq.7.or.pft.eq.11.or.pft.eq.12) THEN
!   IF (Lati.gt.0) THEN
!      leafbio(1) = Cab + Cab * 0.1 * SIN(2*3.14159/365*(doy-81))
!      Vcmo = Vcmo + Vcmo * 0.1 * SIN(2*3.14159/365*(doy-81))
!   ELSE
!      leafbio(1) = Cab + Cab * 0.1 * SIN(2*3.14159/365*(doy-264))
!      Vcmo = Vcmo + Vcmo * 0.1 * SIN(2*3.14159/365*(doy-264))
!   ENDIF
!ELSE
!   IF (Lati.gt.0) THEN
!      leafbio(1) = Cab + Cab * 0.5 * SIN(2*3.14159/365*(doy-81))
!      Vcmo = Vcmo + Vcmo * 0.5 * SIN(2*3.14159/365*(doy-81))
!   ELSE
!      leafbio(1) = Cab + Cab * 0.5 * SIN(2*3.14159/365*(doy-264))
!      Vcmo = Vcmo + Vcmo * 0.5 * SIN(2*3.14159/365*(doy-264))
!   ENDIF
!ENDIF

! Computation of the fluorescence matrices 
  CALL fluspect(leafbio,MfI,MbI,MfII,MbII,rho,tau,rs,kChlrel)
! ALLOCATE ARRAYS DEPENDING ON LAI THAT THEY ARE OK FOR THE SELECTED LAYERS
! WHICH DEPEND ON LAI
!  CALL Nlayers(LAI)

!print*,' MfI',minval(MfI),maxval(MfI),size(MfI)
!print*,' MbI',minval(MbI),maxval(MbI),size(MbI)
!print*,' MfII',minval(MfII),maxval(MfII),size(MfII)
!print*,' MbII',minval(MbII),maxval(MbII),size(MbII)
!print*,' rho',minval(rho),maxval(rho),size(rho)
!print*,' tau',minval(tau),maxval(tau),size(tau)
!print*,' rs ',minval(rs),maxval(rs),size(rs)
!print*,' kChlrel ',minval(kChlrel),maxval(kChlrel),size(kChlrel)
! Deallocate field with layers 
!CALL fieldlayer_deallocate

! ALLOCATE THESE FIELDS WITH THE NEW nl  
!CALL fieldlayer_allocate



! We calculate the xlay and dx the thickness of the layers
!  CALL layers

! We calculate leaf hotspot parameter (q) from vegetation height and leafwidth
                q = leafwidth_arr(jl)/hc_arr(jl)

! We call leafangles subroutine to determine leaf-angle distribution function (lidf)
            LIDFa = LIDFa_arr(jl)
            LIDFb = LIDFb_arr(jl)
CALL leafangles(LIDFa,LIDFb,lidf)

! ALLOCATE LOCAL FIELDS 
! For shaded leaves
                 nlh = nl
ALLOCATE(Agh(nlh))
ALLOCATE(Ah(nlh),rcwh(nlh),Fh(nlh))
ALLOCATE(Cih(nlh),Cch (nlh),Tch(nlh))
ALLOCATE(A0(nlh),Ag0(nlh),rcw0(nlh),F0a(nlh),W0(nlh),F0(nlh))
ALLOCATE(phi_p_h(nlh),phi_fs_h(nlh),phi_npq_h(nlh))

! For sunlight leaves
                nlu = nli*nlazi*nl
ALLOCATE(Agu(nlu))
ALLOCATE(Au(nlu),rcwu(nlu),Fu(nlu))
ALLOCATE(Ciu(nlu),Ccu(nlu),Tcu(nlu))
ALLOCATE(A1(nlu),Ag1(nlu),rcw1(nlu),F1a(nlu),W1(nlu),F1(nlu))
ALLOCATE(phi_p_u(nlu),phi_fs_u(nlu),phi_npq_u(nlu))


!%% 1. Solar and sky radiation
!% SAIL model for solar and sky radiation, without emission by canopy and
!% soil itself
!% subscript '0' indicates: shaded leaf/soil, '1' indicates: sunlit
!% Rn is net radiation (W m-2), E a spectrally integrated flux in W m-2,
!% and E_ a flux in W m-2 um-1


! We put the values of the temperture
Tch             = Ta+.1 ! %   Leaf temperature (shaded leaves)
Tcu             = Ta+.3 ! %   Leaf temperature (sunlit leaves)

! Initialisation of the concentrations 
 Cch            = Cc
 Ccu            = Cc 

!print*, 'VERIF ', ' Ca ', Cc, ' Oa ', Oa  
!print*, 'VERIF ', ' Ta ', Ta, ' pa ', pa  
!print*, 'VERIF ', ' ea ', ea!, ' u ', u  
!print*, 'VERIF ', ' lam ', lam  
!print*, 'VERIF ', ' Rin ', Rin, ' Rli ', Rli  
!print*, 'VERIF ', ' LAI ', LAI, ' Vcmo ', Vcmo  
!print*, 'VERIF ', ' Jmo ', Jmo
!print*, 'VERIF ', ' LAI ', LAI, ' Vcmo ', Vcmo  
!print*, 'VERIF ', ' tts ', tts, ' psi ', psi, ' tto ', tto   
!print*, 'VERIF ', ' z ', z, ' hc ', hc   
!print*, 'VERIF ', ' zo ', zo, ' d ', d   
!print*, 'VERIF ', ' Cab ', Cab, ' Cdm ', Cdm   
!print*, 'VERIF ', ' Cw ', Cw, ' Csm ', Csm   
!print*, 'VERIF ', ' N ',N

! Now we call the routines 
! 1. RTMO 
CALL rtmo(Rin,Rli,Ta,LAI,tts,tto,psi,Ps,Po,Pso,km, Kext, &
        & Esun_,Esky_,P, fEsuno,fEskyo,fEsunt,fEskyt, Eplu_, Emin_, &
        & Lo_, Eout_, Eouto,Eoutt, Rnhs, Rnus, Rnhc, Rnuc, Pnhc, Pnuc,&
        & Pnhc_Cab, Pnuc_Cab, rho, tau, rs, kChlrel, lidf, q)

! Matrix containing values for 1-Ps and Ps of soil
        Fs_mat(1) = 1.-Ps(size(Ps))
        Fs_mat(2) = Ps(size(Ps))

! Matrix containing values for Ps of canopy
         Fc   = (1-Ps(1:size(Ps)-1))/nl     ! Matrix containing values for Ps of canopy

!print*, ' Fc ', minval(Fc), maxval(Fc)


! % 2.3. Biochemical processes
!% photosynthesis (A), fluorescence factor (F), and stomatal resistance (rcw),
!for shaded (h) and sunlit (u) leaves

! Correction of the temperature ... To fix this 
       tempcor = 1

! Shaded leaves 
!print*, ' Pnhc*1E6 ', minval(Pnhc*1E6), maxval(Pnhc*1E6), sum(Pnhc*1E6)
!print*, ' Tch ', minval(Tch), maxval(Tch), sum(Tch)
!print*, ' Cch ', minval(Cch), maxval(Cch), sum(Cch)
!print*, ' Q ', minval(Pnhc*1E6), maxval(Pnhc*1E6), sum(Pnhc*1E6)
!print*, ' Oa ', Oa
!print*, ' ea ', ea
!print*, ' pa ', pa
!print*, ' Vcmo ', Vcmo
!print*, ' Rdparam' ,Rdparam
!print*, ' Tparams ', Tparams
!print*, 'm ', m
!print*, 'tempcor ',tempcor 
!print*, 'stressfactor ' , stressfactor

! Faquahar used for GPP computation
! Not tested for the moment 
IF (faq .eq.1) THEN 
! Computation of the ratio of fluorescences over the shaded leaves using
CALL biochemical_faq(nlh,Ccc,Pnhc*1e6,Tch,ea0(jj),Occ,pres(jj),vcmo,Jmo,EC(jl),EO(jl),&
   EV(jl),ER(jl),EK(jl),kc0(jl)*1e6,ko0(jl)*1e3,option,gcmethod,A0,F0a,rcw0,W0,Ag0,pft,0)


! Farquahar used for GPP computation
! Computation of the ratio of fluorescences over the shaded leaves using
CALL biochemical_faq(nlh,Ccc,Pnuc*1e6,Tcu,ea0(jj),Occ,pres(jj),vcmo,Jmo,EC(jl),EO(jl),&
   EV(jl),ER(jl),EK(jl),kc0(jl)*1e6,ko0(jl)*1e3,option,gcmethod,A1,F1a,rcw1,W1,Ag1,pft,0)


! When using Faquar To come here
! % convert fluorescence factor from units of aPAR to units of indident PAR
!  (iPAR)


!    % convert fluorescence factor from units of aPAR to units of indident PAR
!    (iPAR)
    F0          = F0a*Pnhc/Rnhc
    F1          = F1a*reshape(Pnuc,(/nli*nlazi*nl/))/reshape(Rnuc,(/nli*nlazi*nl/))


 type_integration = 'angles_and_layers'     ! This is working
               np = 1
IF (.NOT.ALLOCATED(Fout)) ALLOCATE(Fout(np))
             Fout = 0.

! Net photosynthesis
CALL integration(1,A1,type_integration,Ps(1:nl),lidf,Fout)
Actot_faq   = LAI*(dot_product(Fc,A0) + Fout(1))

! GPP
CALL integration(1,Ag1,type_integration,Ps(1:nl),lidf,Fout)
Agtot_faq   = LAI*(dot_product(Fc,Ag0) + Fout(1))

ENDIF 


! Collatz used or GPP 
! Shaded leaves 
CALL biochemical(nlh,Pnhc_Cab*1E6,Tch,Cch,ea,Oa,pa,kc0(jl)*1e6,ko0(jl)*1e3,Vcmo,&
               & option,Agh,Ah,Fh,rcwh,Cih,avovc,ardvc,phi_p_h,phi_fs_h,phi_npq_h)

! Sunlit leaves 
CALL  biochemical(nlu,reshape(Pnuc_Cab,(/nlu/))*1E6,Tcu,Ccu,ea,Oa,pa,kc0(jl)*1e6,ko0(jl)*1e3,Vcmo,&
                & option,Agu,Au,Fu,rcwu,Ciu,avovc,ardvc,phi_p_u,phi_fs_u,phi_npq_u)

print*,' Fh:',shape(Fh),minval(Fh),maxval(Fh)
print*,' Fu:',shape(Fu),minval(Fu),maxval(Fu)
print*,'shaded'
print*,' phi_p, (1)=',phi_p_h(1),' (nl)=',phi_p_h(nl)
print*,' phi_npq, (1)=',phi_npq_h(1),' (nl)=',phi_npq_h(nl)

!print*,'   sum'
!print*,'   min=',minval(phi_p_h+phi_fs_h+phi_npq_h)
!print*,'   max=',maxval(phi_p_h+phi_fs_h+phi_npq_h)
print*,' phi_p_h:',shape(phi_p_h),minval(phi_p_h),maxval(phi_p_h)
print*,' phi_fs_h:',shape(phi_fs_h),minval(phi_fs_h),maxval(phi_fs_h)
print*,' phi_npq_h:',shape(phi_npq_h),minval(phi_npq_h),maxval(phi_npq_h)
print*,'sunlit'
print*,'   sum'
print*,'   min=',minval(phi_p_u+phi_fs_u+phi_npq_u)
print*,'   max=',maxval(phi_p_u+phi_fs_u+phi_npq_u)
print*,' phi_p_u:',shape(phi_p_u),minval(phi_p_u),maxval(phi_p_u)
print*,' phi_fs_u:',shape(phi_fs_u),minval(phi_fs_u),maxval(phi_fs_u)
print*,' phi_npq_u:',shape(phi_npq_u),minval(phi_npq_u),maxval(phi_npq_u)
! 3. Calculation of the fluorescence 
CALL rtmf(Esun_, transpose(Emin_), transpose(Eplu_),Fh,reshape(Fu,(/nli,nlazi,nl/)),&
           & LAI,Po,Ps,Pso,tts,tto,psi,LoF,Fhem,Fiprof,MfI,MbI,MfII,MbII,rho,tau,rs,lidf)

!Integration over the layers 
 type_integration = 'angles_and_layers'     ! This is working
               np = 1
IF (.NOT.ALLOCATED(Fout)) ALLOCATE(Fout(np))
             Fout = 0.

! GPP  ??? without substracting the dark respiration 
CALL integration(1,Agu,type_integration,Ps(1:nl),lidf,Fout)
           Agtot  = LAI*(dot_product(Fc,Agh) + Fout(1))

! NPP ???  To ckarify this 
CALL integration(1,Au,type_integration,Ps(1:nl),lidf,Fout)
            Actot = LAI*(dot_product(Fc,Ah) + Fout(1))

!% Net PAR 
CALL integration(1,reshape(Pnuc,(/nlu/)),type_integration,Ps(1:nl),lidf,Fout)
           Pntot  = LAI*(dot_product(Fc,Pnhc) + Fout(1))

!% net PAR_Cab 
CALL integration(1,reshape(Pnuc_Cab,(/nli*nlazi*nl/)),type_integration,Ps(1:nl),lidf,Fout)
        Pntot_Cab = LAI*(dot_product(Fc,Pnhc_Cab) + Fout(1))

! Vertical profile calculations
 type_integration = 'angles'
 np = nl
IF (.NOT.ALLOCATED(Fout_profile)) ALLOCATE(Fout_profile(np))
             Fout_profile = 0.
!PSII fluorescence yields
! Sunlit fluorescence yield per layer
CALL integration(1,reshape(phi_fs_u,(/nli,nlazi,nl/)),type_integration,Ps(1:nl),lidf,Fout_profile)
        fs_yieldu_toc = Fout_profile(1)
! Shaded fluorescence yield at top-of-canopy
        fs_yieldh_toc = phi_fs_h(1)
!PSII photochemical yields
! Sunlit photochemical yield per layer
CALL integration(1,reshape(phi_p_u,(/nli,nlazi,nl/)),type_integration,Ps(1:nl),lidf,Fout_profile)
        p_yieldu_toc = Fout_profile(1)
! Shaded photochemical yield at top-of-canopy
        p_yieldh_toc = phi_p_h(1)
!PSII non-photchemical quenching yields
! Sunlit NPQ yield per layer
CALL integration(1,reshape(phi_npq_u,(/nli,nlazi,nl/)),type_integration,Ps(1:nl),lidf,Fout_profile)
        npq_yieldu_toc = Fout_profile(1)
! Shaded NPQ yield at top-of-canopy
        npq_yieldh_toc = phi_npq_h(1)


! Fluo per jl and for the wavelenght of the studied satellite 
                  LoF_jl = LoF(ifreq_sat)

! We free the local arrays 
DEALLOCATE(Agh,Ah,rcwh,Fh)
DEALLOCATE(A0,Ag0,rcw0,F0a,F0,W0)
DEALLOCATE(Cih,Cch,Tch)
DEALLOCATE(Fout,Fout_profile)
DEALLOCATE(phi_p_h,phi_fs_h,phi_npq_h)
DEALLOCATE(phi_p_u,phi_fs_u,phi_npq_u)

DEALLOCATE(Agu,Au,rcwu,Fu)
DEALLOCATE(A1,Ag1,rcw1,F1a,F1,W1)
DEALLOCATE(Ciu,Ccu,Tcu)

! This activate if we want to calculate the number of layers as function of
! pft... There is memory leaks so far ... 
! Deallocate fields dependent to the layer
!CALL fieldlayer_deallocate

! For the cell jj
! We compute the total of the fluorescence for the selected BETHY grid cell
! by weighting ith the fraction of the PFT in the selected grid cell of BETHY

! Fluorescence

ENDIF      ! Test on LAI if > 0 then calculation made

!    For daily output files
        dayfluo(jl) = LoF_jl
        daygpp(jl) = Agtot

                        LoF_jl = LoF(ifreq_sat)
! IF (isNaN(LoF_jl))     LoF_jl = 0.
        rfluo(iyear,imonth,jj) = rfluo(iyear,imonth,jj) +LoF_jl*frac1
!        print*,'LoF_jl*frac1 equals: ',LoF_jl*frac1
! GPP
     rgppfluo(iyear,imonth,jj) = rgppfluo(iyear,imonth,jj)+Agtot*frac1
!     print*,'Agtot*frac1 equals: ',Agtot*frac1

     zgppfluo(jl) = zgppfluo(jl)+Agtot*frac1     

     ! Convert zgppfluo from umolC/m2/s to gC/m2/s
     zgppfluo(jl) = zgppfluo(jl)*1.E-6*12.0

     ! Output data per vegetation point and per diurnal time-step 
     !it = MOD(t,24.)+1     ! diurnal (sub-day) time index. Allows time (t) 24 to be the first element in the output array
     rfluo_diurnal(iyear,iday0,t,jl) = LoF_jl*frac1
     rgppfluo_diurnal(iyear,iday0,t,jl) = Agtot*frac1
     rlai_diurnal(iyear,iday0,t,jl) = LAI
     rapar_diurnal(iyear,iday0,t,jl) = Pntot*1e6*frac1
     raparcab_diurnal(iyear,iday0,t,jl) = Pntot_Cab*1e6*frac1
     rpar_diurnal(iyear,iday0,t,jl) = P*1e3*frac1
     rswdown_diurnal(iyear,iday0,t,jl) = Rin    ! hourly shortwave radiation down
     rta_diurnal(iyear,iday0,t,jl) = Ta         ! hourly air temperature

     rfsyieldu_toc_diurnal(iyear,iday0,t,jl) = fs_yieldu_toc
     rfsyieldh_toc_diurnal(iyear,iday0,t,jl) = fs_yieldh_toc
     rpyieldu_toc_diurnal(iyear,iday0,t,jl) = p_yieldu_toc
     rpyieldh_toc_diurnal(iyear,iday0,t,jl) = p_yieldh_toc
     rnpqyieldu_toc_diurnal(iyear,iday0,t,jl) = npq_yieldu_toc
     rnpqyieldh_toc_diurnal(iyear,iday0,t,jl) = npq_yieldh_toc

!INCIDENT PAR computed from mo_rtmo for the selected grid cell
   PAR_scope(iyear,imonth,jj)  =  PAR_scope(iyear,imonth,jj) + Pntot*1e6*frac1
PAR_scope_cab(iyear,imonth,jj) =  PAR_scope_cab(iyear,imonth,jj) +Pntot_Cab*1e6*frac1

!Deallocate fluspect output arrays
DEALLOCATE(MfI)
DEALLOCATE(MbI)
DEALLOCATE(MfII)
DEALLOCATE(MbII)
DEALLOCATE(rho)
DEALLOCATE(tau)
DEALLOCATE(rs)
DEALLOCATE(kChlrel)

!Deallocate leafangles output
DEALLOCATE(lidf)

! write(6,'(a4,3(1x,i5),1x,f6.3,8(1x,f7.2),2(1x,f5.2),1x,i4,3(1x,f5.1),2(1x,f8.2))') 'FLUO',imonth,jj,jl,LoF_jl,Agtot,Actot,Long,Lati,&
!& Rin,Ta,ea,pa,LAI,frac1,pft,Vcmo,Cab,tts,Pntot*1e6,Pntot_Cab*1.e6

if (faq.eq.1) then 
 write(6,'(a4,2(1x,i5),1x,f6.3,6(1x,f7.2),1x,i3,2(1x,f6.2))') 'FAQ',jj,jl,LoF_jl,Agtot,Agtot_faq,Actot,Actot_faq,Long,Lati,pft,Vcmo,Cab 
endif 

! Deallocate fields dependent to the layer
!IF (ALLOCATED(xlay)) CALL fieldlayer_deallocate

!END DO  ! Loop on the jd 

END Do     ! loop on jl

!$OMP END PARALLEL DO

ial = 1


END SUBROUTINE  fluorescence


! This routine considers only the fluorescence for which observations exist 
! This then creates a pair of fluorescence data, (obs,mod)
SUBROUTINE fluo_obs(f,fs) 
    IMPLICIT NONE
    REAL, INTENT(in), DIMENSION(:,:,:)     :: f ! simulated fluoresences 
    REAL, INTENT(inout), DIMENSION(:)      :: fs ! simulated fs corresponding to the observed data 
    INTEGER                                :: i

    DO i = 1, nobs_fs
       fs(i)=f(indy_fs(i),indm_fs(i),indng_fs(i))
       !print*, ' i ', indy_fs(i),indm_fs(i),indng_fs(i)
       !print*, ' rgpp ', f(indy_fs(i),indm_fs(i),indng_fs(i)),fs(i)
      ! print*, ' COMP fs_obs  fs_unc fs ', date_fs(i),fs_obs(i), fs_unc(i),fs(i)
!       print*
    END DO

    RETURN

END SUBROUTINE fluo_obs



! Routine pre-computes some quantities of the rtmo and rtmf and this to speed of
! the process of optimisation. To make the pre-computation only 

!SUBROUTINE precom_rtm
!END SUBROUTINE precom_rtm 


!% function [Fi_s,Fi_gs,Fi_g]= calczenithangle(Doy,t,Omega_g,Fi_gm,Long,Lat)
!SUBROUTINE calczenithangle(Doy,t,Omega_g,Fi_gm,Long,Lat,Fi_s,Fi_gs,Fi_g)
SUBROUTINE calczenithangle(Doy,t,tm,Omega_g,Fi_gm,Long,Lat,Fi_s)

USE mo_constants

!%
!% calculates pi/2-the angle of the sun with the slope of the surface.
!%
!% input:
!% Doy       day of the year
!% t         time of the day (hours, local SOLAR time)
!% Omega_g   slope azimuth angle (deg)
!% Fi_gm     slope of the surface (deg)
!% Long      Longitude (decimal)
!% Lat       Latitude (decimal)
!%
!% output:
!% Fi_s      'classic' zenith angle: perpendicular to horizontal plane
!% Fi_gs     solar angle perpendicular to surface slope
!% Fi_g      projected slope of the surface in the plane through the solar beam and the vertical
!% 

IMPLICIT NONE 

! Input variables 
REAL, INTENT(IN)                 ::  Doy,t,Omega_g,Fi_gm,Long,Lat

!Output variables 
REAL, INTENT(OUT)                :: Fi_s 
REAL                             :: Fi_gs,Fi_g 
REAL,INTENT(OUT)                 :: tm

! Local variables
!REAL                             :: GG, d, Et,tm,Omega_s
REAL                             :: GG, d, Et,Omega_s
REAL                             :: Omega_g1,Fi_gm1, Lat1
REAL                             :: Fi_s_max 

!%convert angles into radials
GG               =   (Doy-1)/365*2*pi       !    % converts day of year to radials
Omega_g1         =   Omega_g/180*pi         !    % converts direction of slope to radials
Fi_gm1           =   Fi_gm/180*pi           !    % converts maximum slope to radials
Lat1             =   Lat/180*pi             !    % converts latitude to radials

!%computes the declination of the sun
d               =   0.006918-0.399912*cos(GG  )+ 0.070247*sin(GG  )- &
                   &  0.006758*cos(2*GG)+ 0.000907*sin(2*GG)- &
                   &  0.002697*cos(3*GG)+ 0.00148*sin(3*GG)
                                
!%Equation of time
Et              =   0.017 + .4281 * cos(GG) - 7.351 * sin(GG) - & 
                    & 3.349 * cos(2*GG) - 9.731 * sin(2*GG)

!%computes the time of the day when the sun reaches its highest angle                                
!tm              =   12+(4*(-Long)-Et)/60    !  % de Pury and Farquhar (1997), Iqbal (1983)
tm              =   12+(Et)/60 

!%computes the hour angle of the sun
Omega_s         =   (t-tm)/12*pi

!%computes the zenithangle (equation 3.28 in De Bruin)
Fi_s            =   acos(sin(d)*sin(Lat1)+cos(d)*cos(Lat1)*cos(Omega_s))

! EK 11/02/2013  I do not consider Omega_s for test and this gives the time at
! which the maximum of radiation is obtained at the selected location  
Fi_s_max           =   acos(sin(d)*sin(Lat1)+cos(d)*cos(Lat1))


!%computes the slope of the surface Fi_g in the same plane as the solar beam
Fi_g            =   atan(tan(Fi_gm1)*cos(Omega_s-Omega_g1))

!%computes the angle of the sun with the vector perpendicular to the surface
Fi_gs           =   Fi_s + Fi_g

!print*, ' GEOM Fi_s Fi_g Fi_gs ', Fi_s,Fi_g ,Fi_gs

END SUBROUTINE   calczenithangle


! This routine approximates  the simulated fluorescence data by the GPP   
! Conversion of gpp into fs 
SUBROUTINE conv_gppfs (afs,f,fs,ifs)
  USE mo_constants
  USE constants_pjr

    IMPLICIT NONE 
    REAL, INTENT(in)                       :: afs !  
    INTEGER, INTENT(in)                    :: ifs !  
    REAL, INTENT(in), DIMENSION(:,:,:)     :: f ! simulated bethy fluxes  
    REAL, INTENT(inout), DIMENSION(:)      :: fs !  simulated fs as function of f and afs 
    INTEGER                                :: i 
 
    if (ifs/=0) then 
      DO i = 1, nobs_fs 
       fs(i)=afs*f(indy_fs(i),indm_fs(i),indng_fs(i))
       !print*, ' i ', indy_fs(i),indm_fs(i),indng_fs(i)         
       !print*, ' rgpp ', f(indy_fs(i),indm_fs(i),indng_fs(i)),fs(i)
       !print*, ' fs_obs  fs_unc fs ', date_fs(i),fs_obs(i), fs_unc(i),fs(i)
!       print*
      END DO 
    else 
    fs = 0.
    endif 

    ! stop 

    RETURN

END SUBROUTINE conv_gppfs 


subroutine lower_case(word)
! convert a word to lower case
character (len=*) , intent(inout) :: word
integer :: i,ic,nlen
nlen = len(word)
do i=1,nlen
   ic = ichar(word(i:i))
   if (ic >= 65 .and. ic < 90) word(i:i) = char(ic+32)
end do
end subroutine lower_case



SUBROUTINE modtran_file(month, lon)

USE fluo_param, ONLY : atmos_file,modtran_trop, modtran_sum, modtran_wint 

IMPLICIT NONE 

INTEGER, INTENT(IN)   :: month 
REAL, INTENT(IN)      :: lon 


! Use of tropical atmosphere 
IF ((lon.le.30.).and.(lon.ge.-30.)) atmos_file=trim(modtran_trop)


! NORTH HEMISPHERE 
! Winter in mid-latitude in northern hemisphere 
If (lon.gt.30.) then 
 if ((month.ge.10).or.(month.le.3)) atmos_file=trim(modtran_wint) 
endif 

! Summer  in mid-latitude in northern hemisphere
If (lon.gt.30.) then
 if ((month.ge.4).and.(month.le.9)) atmos_file=trim(modtran_sum)
endif

! SOUTH HEMISPHERE 
! Summer in mid-latitude in southern hemisphere
If (lon.lt.-30.) then
 if ((month.ge.10).or.(month.le.3)) atmos_file=trim(modtran_sum)
endif

! Winter in mid-latitude in northern hemisphere
If (lon.lt.-30.) then
 if ((month.ge.4).and.(month.le.9)) atmos_file=trim(modtran_wint)
endif

END SUBROUTINE modtran_file 

!SUBROUTINE read_modtran_files()

!USE fluo_param, ONLY : read_atm_files,atmos_file,modtran_trop, modtran_sum, modtran_wint

!IMPLICIT NONE

!integer :: iatmosfile

! Use of tropical atmosphere
!iatmosfile = 1
!atmos_file = trim(modtran_trop)
!call read_atm_files(iatmosfile, atmos_file)
!IF ((lon.le.30.).and.(lon.ge.-30.)) atmos_file=trim(modtran_trop)


! NORTH HEMISPHERE 
! Winter in mid-latitude in northern hemisphere 
!iatmosfile = 2
!atmos_file = trim(modtran_wint)
!call read_atm_files(iatmosfile, atmos_file)
!If (lon.gt.30.) then
! if ((month.ge.10).or.(month.le.3)) atmos_file=trim(modtran_wint)
!endif

! Summer  in mid-latitude in northern hemisphere
!iatmosfile = 3
!atmos_file = trim(modtran_sum)
!call read_atm_files(iatmosfile, atmos_file)
!If (lon.gt.30.) then
! if ((month.ge.4).and.(month.le.9)) atmos_file=trim(modtran_sum)
!endif

! SOUTH HEMISPHERE 
! Summer in mid-latitude in southern hemisphere
!If (lon.lt.-30.) then
! if ((month.ge.10).or.(month.le.3)) atmos_file=trim(modtran_sum)
!endif

! Winter in mid-latitude in northern hemisphere
!If (lon.lt.-30.) then
! if ((month.ge.4).and.(month.le.9)) atmos_file=trim(modtran_wint)
!endif

!END SUBROUTINE read_modtran_files

SUBROUTINE modtran_ifile(month, lat, iatmos_file)

USE fluo_param, ONLY : atmos_file,modtran_trop, modtran_sum, modtran_wint

IMPLICIT NONE

INTEGER, INTENT(IN)   :: month
REAL, INTENT(IN)      :: lat
INTEGER, INTENT(OUT)  :: iatmos_file

! Use of standard atmosphere
iatmos_file = 1

! Use of tropical atmosphere 
IF ((lat.le.30.).and.(lat.ge.-30.)) iatmos_file=2


! NORTH HEMISPHERE 
! Winter in mid-latitude in northern hemisphere 
If (lat.gt.30.) then
 if ((month.ge.10).or.(month.le.3)) iatmos_file=3
endif

! Summer  in mid-latitude in northern hemisphere
If (lat.gt.30.) then
 if ((month.ge.4).and.(month.le.9)) iatmos_file=4
endif

! SOUTH HEMISPHERE 
! Summer in mid-latitude in southern hemisphere
If (lat.lt.-30.) then
 if ((month.ge.10).or.(month.le.3)) iatmos_file=4
endif

! Winter in mid-latitude in northern hemisphere
If (lat.lt.-30.) then
 if ((month.ge.4).and.(month.le.9)) iatmos_file=3
endif



END SUBROUTINE modtran_ifile

! This routine is to correct the actual local time of BETHY 
SUBROUTINE pb_hour_bethy(hm)
IMPLICIT NONE 
REAL, INTENT(OUT)     :: hm(24) 

hm(1)=12
hm(2)=13
hm(3)=14
hm(4)=15
hm(5)=16
hm(6)=17
hm(7)=18
hm(8)=19
hm(9)=20
hm(10)=21
hm(11)=22
hm(12)=23
hm(13)=24
hm(14)=1
hm(15)=2
hm(16)=3
hm(17)=4
hm(18)=5
hm(19)=6
hm(20)=7
hm(21)=8
hm(22)=9
hm(23)=10
hm(24)=11



END SUBROUTINE pb_hour_bethy 


END MODULE fluo
  
