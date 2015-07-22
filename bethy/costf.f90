!------------------------------------------------------------------
! cost function module
!------------------------------------------------------------------
MODULE costf

  USE mo_mapping
  USE mo_climate
 
  REAL :: fautleafs, ccosts, lq10ss, lq10fs, ltaufs, laws, lfracss, betas
  REAL :: ECs, EOs, EVs, ERs, EKs, tgams, alphas, alc4s, jmts, kc0s, ko0s
  REAL, ALLOCATABLE :: obf(:,:,:)

  ! include pjr cost module here too
  REAL, ALLOCATABLE, DIMENSION(:)    :: c_obs, c_unc ! data and its unc
  REAL, ALLOCATABLE, DIMENSION(:)    :: flux_obs, flux_unc ! data and its unc
  REAL, ALLOCATABLE, DIMENSION(:,:,:)    :: faparl_obs, faparl_unc ! data and its unc
  REAL, ALLOCATABLE, DIMENSION(:,:,:)    :: faparg_obs, faparg_unc ! data and its unc

! HEW-ADD-03/09/08 for new mapping
  REAL, ALLOCATABLE, DIMENSION(:)    :: px0, p0su ! first guess px0 | uncertainty p0su (physical parameters)

! HEW-DEL-03/09/08 x0s removed ( x0 and x are in units of p0sl rsp. p0su , so x0s=1 )
  REAL, ALLOCATABLE, DIMENSION(:)    :: x0 ! natural initial parameters and their unc 
  INTEGER, ALLOCATABLE, DIMENSION(:)    :: xpf ! flag for parameter transformation
  REAL, ALLOCATABLE, DIMENSION(:)    :: xpfa ! lower boundary if xpf = 2, 3 or 4
  REAL, ALLOCATABLE, DIMENSION(:)    :: xpfb ! uppper bounary if xpf = 4
  REAL, ALLOCATABLE, DIMENSION(:)    :: a ! factor for log transformation

  INTEGER, ALLOCATABLE, DIMENSION(:) :: x_a, x_e, y_a, y_e

CONTAINS

  INTEGER FUNCTION count_params (param_file)
    IMPLICIT NONE 
    CHARACTER(len=80) :: param_file

    OPEN(unit=1, file=param_file, status='old')
    REWIND 1
    READ(1,*) ! skip heading
    READ(1,*) count_params
    RETURN
  END FUNCTION count_params


  SUBROUTINE init_param (np, PARAM_FILE)

    USE mo_trafo

    IMPLICIT NONE

    REAL, ALLOCATABLE, DIMENSION(:) :: px0up, px0dn
    INTEGER, INTENT(in) ::  np ! size of control & param vector
    CHARACTER(len=*) :: param_file ! location of initial values
    ! local vars
    INTEGER :: i, ncheck
    !REAL :: e

    !e = EXP(1.)

    OPEN(1,file=param_file,status='old')
    REWIND 1
    READ(1,*)
    READ(1,*) ncheck

    IF(ncheck > np) THEN
       WRITE(6,*) 'init_param: too many parameters in file ',param_file
       STOP
    ENDIF
    WRITE(6,*) '# number of params to optimise:',ncheck
    ALLOCATE(x0(np), xpf(np), xpfa(np), xpfb(np))
! HEW-ADD-03/09/08 allocate new arrays
!HEW-CHG-041209 : a and b added :
    ALLOCATE(px0(np), px0up(np), px0dn(np), p0su(np), a(np)) 
    a = 0.

! Input of a priori values of parameters
!---------------------------------------

    DO i=1,np
       ! HEW-CHG-03/09/08 read new param file
       READ(1,*) px0(i), p0su(i), xpf(i), xpfa(i), xpfb(i)
       ! HEW-CHG-04/03/03 xpf parameter transformation flag
       IF (abs(xpf(i)) == 3) THEN
          a(i) = LOG(p0su(i)+px0(i)-xpfa(i)) - LOG(px0(i)-xpfa(i))
       ELSE IF (abs(xpf(i)) == 2) THEN
          a(i) = ((p0su(i)-xpfa(i))/(SQRT(p0su(i)+px0(i)-xpfa(i)) + SQRT(px0(i)-xpfa(i))))**2
       ELSE IF (abs(xpf(i)) == 4) THEN
          a(i) = -log((xpfb(i)-xpfa(i))/(p0su(i)+px0(i)-xpfa(i)) - 1) + log((xpfb(i)-xpfa(i))/(px0(i)-xpfa(i)) - 1)
       ELSE
          a(i) = 0.
       ENDIF
    END DO
    CLOSE(1)

    x0 = p2x(px0, xpf, xpfa, xpfb, px0, p0su, a, np)
    px0up = x2p(x0+1,x0,xpf,xpfa,xpfb,px0,p0su,a,np)
    px0dn = x2p(x0-1,x0,xpf,xpfa,xpfb,px0,p0su,a,np)

    OPEN(2,file='params_out.txt')
    WRITE(2,*) 'No.   Mean   Mean-sigma   Mean+sigma   Flag'
    DO i=1,np
!       IF (abs(xpf(i)) == 3) THEN
!          WRITE(2,'(i3,3(e18.10),i3)') i,px0(i),EXP((x0(i)-1)*a(i))+xpfa(i),EXP((x0(i)+1)*a(i))+xpfa(i),xpf(i)
!       ELSE IF (abs(xpf(i)) == 2) THEN
!          WRITE(2,'(i3,3(e18.10),i3)') i,px0(i),a(i)*(x0(i)-1)**2+xpfa(i),a(i)*(x0(i)+1)**2+xpfa(i),xpf(i)
!       ELSE IF (abs(xpf(i)) == 4) then
!          WRITE(2,'(i3,3(e18.10),i3)') i,px0(i),(xpfb(i)-xpfa(i))/(1+exp(-(x0(i)-1)*a(i)))+xpfa(i), &
!               & (xpfb(i)-xpfa(i))/(1+exp(-(x0(i)+1)*a(i)))+xpfa(i),xpf(i)
!       ELSE
!          WRITE(2,'(i3,3(e18.10),i3)') i,px0(i),px0(i)-p0su(i),px0(i)+p0su(i),xpf(i)
!       ENDIF
       WRITE(2,'(i3,3(e18.10),i3)') i,px0(i),px0dn(i),px0up(i),xpf(i)
    ENDDO
    CLOSE(2)

    
    RETURN
  END SUBROUTINE INIT_PARAM
  
  
  INTEGER FUNCTION   init_conc( datfile, firstdat, lastdat, datdir)
  ! gets concentrations and their unc from Globalview as listed in datafile,  and initializes  the jacobiann as well
  USE constants_pjr
  IMPLICIT NONE
  CHARACTER( len=*), INTENT(in) :: datfile ! list of stations we'll use
  INTEGER, INTENT(in) :: firstdat, lastdat  ! first and last years of data
  
  CHARACTER( len=*), OPTIONAL, INTENT(in) :: datdir ! locations for data and jacobians!
  ! local variables
  
  CHARACTER(len=80) :: in_datdir ! internal copies of arguments if present
  INTEGER  i, j, k, l,n ! index variables
  
  CHARACTER(len=80) :: obs_name, unc_name
  INTEGER, PARAMETER :: obs_unit = 10
  INTEGER, PARAMETER :: unc_unit=11
  
  INTEGER, PARAMETER ::  inunit=2
  INTEGER i_stat ! index variable for stations
  REAL t, x, ref, xunc, diff ! things we read in
  
  LOGICAL :: fex

  IF( PRESENT( datdir)) THEN
     in_datdir = datdir
  ELSE
     in_datdir = "./input/gv/"
  ENDIF


  ! count number of stations

  OPEN(unit=inunit,file=TRIM(datfile),status='old')
  REWIND inunit
  n_stats = 0
  DO WHILE (.TRUE.)
     READ(inunit,*,END=1000) ! just skip data
     n_stats = n_stats + 1
  END DO
1000 CONTINUE
  IF (n_stats>0) PRINT*, 'Using ', n_stats, ' flask stations for simulation'
  
  ! handle allocation
  ALLOCATE( c_obs( n_stats *months_per_year *(lastdat - firstdat + 1)))
  c_obs = 0.
  ALLOCATE( c_unc ( n_stats * months_per_year *(lastdat -firstdat + 1)))
  c_unc = 1.
  ALLOCATE( stat_names( n_stats))
  
  REWIND inunit
  DO i_stat = 1, n_stats
     READ(inunit,*) stat_names( i_stat)
     ! create name of data file
     obs_name = TRIM(in_datdir)//TRIM(stat_names(i_stat))//'_01D0_ext.co2' !  only for NOAA stations
     ! and std-dev file
     unc_name = TRIM(in_datdir)//TRIM(stat_names(i_stat))//'_01D0_wts.co2'

     INQUIRE (file=TRIM(obs_name), exist=fex)
     IF (FEX) THEN
        OPEN(unit=obs_unit,file=TRIM(obs_name),status='old')
        REWIND obs_unit 
        OPEN(unit=unc_unit,file=unc_name,status='old')
        REWIND unc_unit 
        ! skip headings
        DO i=1,17 ! should skip the header + line for 1979.00 which we don't want 
           READ(obs_unit,*)
           READ(unc_unit,*) 
        END DO

        DO WHILE(.TRUE.)
           READ(obs_unit,*,END=2000) t,x,ref ,diff
           IF( MOD(t, 1.0) > 0.01 .AND. MOD(t, 1.0) < 0.03)THEN

              READ(unc_unit,*) xunc,xunc ! start of a year, read the unc

           ENDIF

           IF ( in_range( t)) THEN 
              ! we want to ad it to our data but this is very tricky
              ! the Jacobian has each year of data together so that's the way we want it
              n = INT(t -0.001  - firstdat) *months_per_year * n_stats + (i_stat -1)*months_per_year &
                   & + INT( MODULO( t, 1.) * months_per_year) + 1
              c_obs (n) = ref + diff
              c_unc(n) = SQRT(0.25 + xunc*xunc)
              ! there's 4/month in GV, read the other 3
              DO i=1,3
                 READ(obs_unit,*) t,x,ref,diff ! will crash on EOF there but so it should 
                 c_obs (n) = c_obs(n) + ref + diff
              END DO
              c_obs(n) = c_obs (n) /4.
           ENDIF
        END DO
2000    CONTINUE
        CLOSE( obs_unit)
        CLOSE( unc_unit)
     ELSE
        PRINT *, 'No data for station ', TRIM(stat_names(i_stat)), ' using zero instead'
     ENDIF
  END DO
  CLOSE( inunit)
  init_conc = n_stats
  
  RETURN

CONTAINS
  
  LOGICAL FUNCTION in_range(t)
    REAL t
    IF((INT(t -0.001) >= firstdat) .AND. (INT(t-0.001) <= lastdat)) THEN
       in_range = .TRUE.
    ELSE
       in_range = .FALSE.
    ENDIF
    RETURN
  END FUNCTION in_range
  
END FUNCTION init_conc



  INTEGER FUNCTION   init_flux( site_file, firstdat, lastdat, fluxdir)
  ! gets eddy fluxes and their unc from sites as listed in datafile

  USE constants_pjr
  USE mo_namelist, ONLY: year0_site, year1_site

  IMPLICIT NONE
  CHARACTER( len=*), INTENT(in) :: site_file ! list of stations we'll use
  INTEGER, INTENT(in) :: firstdat, lastdat  ! first and last years of data
  
  CHARACTER( len=*), OPTIONAL, INTENT(in) :: fluxdir ! locations for data and jacobians!
  ! local variables
  
  CHARACTER(len=80) :: in_fluxdir ! internal copies of arguments if present
  INTEGER  i, j, k, l,n ! index variables
  
  CHARACTER(len=80) :: obs_name
  INTEGER, PARAMETER :: obs_unit = 10

  CHARACTER (len=4) :: cy
  CHARACTER(len=80) :: header
  
  INTEGER, PARAMETER ::  inunit=2
  INTEGER i_site ! index variable for stations
  REAL t, x, ref, xunc, diff ! things we read in
  INTEGER st_name
  REAL dummyr
  INTEGER dummyi
  
  LOGICAL :: fex

  IF( PRESENT( fluxdir)) THEN
     in_fluxdir = fluxdir
  ELSE
     in_fluxdir = "./input/eddy_sites/"
  ENDIF
          
  ! count number of sites

  OPEN(unit=inunit,file=TRIM(site_file),status='old')
  REWIND inunit
  READ(inunit,*) header
  n_sites = 0
  DO WHILE (.TRUE.)
     READ(inunit,*,END=1000) 
     n_sites = n_sites + 1
  END DO
1000 CONTINUE
  IF (n_sites>0) PRINT*, 'Using ', n_sites, ' eddy flux sites (including sub_sites) for simulation'
  
  ! handle allocation
!!! MAS add 19/10/07
!!! needs to be revised according to temporal resolution of eddy fluxes !!!
  ALLOCATE( flux_obs( n_sites *months_per_year *(lastdat - firstdat + 1)))
  flux_obs = 0.
  ALLOCATE( flux_unc ( n_sites * months_per_year *(lastdat -firstdat + 1)))
  flux_unc = 1.
  ALLOCATE( site_names( n_sites))
  ALLOCATE( site_clim( n_sites))
  ALLOCATE( x_a(n_sites),x_e(n_sites),y_a(n_sites),y_e(n_sites))
  
  WRITE(cy,'(i4)') year0_site
  
  REWIND inunit
  READ(inunit,*) header
  IF (n_sites>0) THEN
     r_sites=1
  ELSE
     r_sites=0
  ENDIF
  DO i_site = 1, n_sites
     READ(inunit,'(2f7.2,i6,3i3,3f5.2,4i4,a)') dummyr,dummyr,dummyi,dummyi,dummyi,dummyi, &
    &        dummyr,dummyr,dummyr,x_a(i_site),x_e(i_site),y_a(i_site),y_e(i_site),header
     
     st_name = INDEX(header,'!')+2
     site_names( i_site) = TRIM(header(st_name:))

     IF (i_site >1) THEN
        IF (site_names(i_site) /= site_names(i_site-1)) r_sites = r_sites + 1
     ENDIF

     ! create name of data file
     IF (site_names(i_site)=='hainich') THEN
        obs_name = TRIM(in_fluxdir)//TRIM(site_names(i_site))//'_'//cy//'.txt'
     ELSE
        obs_name= TRIM(in_fluxdir)//TRIM(site_names(i_site))//'.'//TRIM(cy)//'.flux.hourly.bin'  
     ENDIF

!     IF (site_clim(i_site)== 0) THEN
        INQUIRE (file=TRIM(obs_name), exist=fex)
        IF (FEX) THEN
!!MS$        OPEN(unit=obs_unit,file=TRIM(obs_name),status='old')
!!MS$        REWIND obs_unit 
!!MS$
!!MS$        ! skip headings
!!MS$        DO i=1,17 ! should skip the header + line for 1979.00 which we don't want 
!!MS$           READ(obs_unit,*)
!!MS$           READ(unc_unit,*) 
!!MS$        END DO
!!MS$
!!MS$        DO WHILE(.TRUE.)
!!MS$           READ(obs_unit,*,END=2000) t,x,ref ,diff
!!MS$           IF( MOD(t, 1.0) > 0.01 .AND. MOD(t, 1.0) < 0.03)THEN
!!MS$
!!MS$              READ(unc_unit,*) xunc,xunc ! start of a year, read the unc
!!MS$
!!MS$           ENDIF
!!MS$
!!MS$           IF ( in_range( t)) THEN 
!!MS$              ! we want to ad it to our data but this is very tricky
!!MS$              ! the Jacobian has each year of data together so that's the way we want it
!!MS$              n = INT(t -0.001  - firstdat) *months_per_year * n_stats + (i_stat -1)*months_per_year &
!!MS$                   & + INT( MODULO( t, 1.) * months_per_year) + 1
!!MS$              c_obs (n) = ref + diff
!!MS$              c_unc(n) = SQRT(0.25 + xunc*xunc)
!!MS$              ! there's 4/month in GV, read the other 3
!!MS$              DO i=1,3
!!MS$                 READ(obs_unit,*) t,x,ref,diff ! will crash on EOF there but so it should 
!!MS$                 c_obs (n) = c_obs(n) + ref + diff
!!MS$              END DO
!!MS$              c_obs(n) = c_obs (n) /4.
!!MS$           ENDIF
!!MS$        END DO
2000       CONTINUE
           CLOSE( obs_unit)
        ELSE
!           PRINT *, 'No flux data for site ', TRIM(site_names(i_site))
        ENDIF
!     ENDIF
  END DO
  CLOSE( inunit)

  ALLOCATE(sub_sites(r_sites))
  DO i=1,n_sites
     IF (i==1) THEN
        k=1
        sub_sites(:) = 1
     ELSE
        IF (i>1 .AND. site_names(i) == site_names(i-1)) THEN
           sub_sites(k) = sub_sites(k)+1
        ELSE
           k=k+1
        ENDIF
     ENDIF
  ENDDO


  init_flux = n_sites
  
  RETURN

CONTAINS
  
  LOGICAL FUNCTION in_range(t)
    REAL t
    IF((INT(t -0.001) >= firstdat) .AND. (INT(t-0.001) <= lastdat)) THEN
       in_range = .TRUE.
    ELSE
       in_range = .FALSE.
    ENDIF
    RETURN
  END FUNCTION in_range
  
END FUNCTION init_flux



  SUBROUTINE   init_faparl(n_sites, nrs, outt )
  ! gets FAPAR and unc from sites

  USE constants_pjr
  USE mo_namelist, ONLY: year0_site, year1_site, faparfile
  USE mo_netcdf
  USE mo_calendar, ONLY: sdays

  IMPLICIT NONE
  INTEGER, INTENT(in) :: n_sites  ! number of sites
  INTEGER, INTENT(in) :: nrs  ! number of site simulation years
  INTEGER, INTENT(in) :: outt  ! time resolution
  
  ! local variables
  
  CHARACTER(len=80) :: in_sitedir ! internal copies of arguments if present
  INTEGER  i, j, k, l,n ! index variables

  CHARACTER (len=4) :: cy
  
  INTEGER, PARAMETER ::  inunit=2
  INTEGER i_site, i_year ! index variables
  INTEGER ind, year
!  INTEGER, DIMENSION(15,15) :: help
  
  LOGICAL :: fex

  TYPE (ncfile) :: infile
  TYPE (ncvar) :: nc_doy, nc_year, nc_fapar, nc_errfapar
  INTEGER, DIMENSION(3) :: iload_dims
  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: iload1, iload2
  REAL, ALLOCATABLE, DIMENSION(:) :: iload3, iload4

  in_sitedir = "./input/eddy_sites/"

  ALLOCATE( faparl_obs( nrs, outt, n_sites))
  faparl_obs = 0.
  ALLOCATE( faparl_unc (nrs, outt, n_sites))
  faparl_unc = 1.e9
  
  DO i_site = 1, n_sites
  
     infile%name=TRIM(in_sitedir)//TRIM(faparfile)//TRIM(site_names(i_site))//'.nc'
     nc_doy%name='DOY'
     nc_year%name='year'
     nc_fapar%name='fapar'
     nc_errfapar%name='errfapar'
     

     CALL ncopen(infile)
     
     CALL ncvarinfo(infile,nc_fapar,dimlens=iload_dims)
     
     ALLOCATE(iload1(iload_dims(1),iload_dims(2),iload_dims(3)))
     ALLOCATE(iload2(iload_dims(1),iload_dims(2),iload_dims(3)))
     ALLOCATE(iload3(iload_dims(3))) ! third dimension contains time information
     ALLOCATE(iload4(iload_dims(3))) ! third dimension contains time information

     CALL ncread(infile,nc_doy,iload3)
     CALL ncread(infile,nc_year,iload4)
     
     CALL ncread(infile,nc_fapar,iload1)
     CALL ncread(infile,nc_errfapar,iload2)

     
     CALL ncclose(infile)
     
     DO i_year = 1, nrs
        year = year0_site+i_year-1
        IF (iload4(1) <= year .AND. iload4(iload_dims(3)) >= year) THEN 
           ind = (365*(year-iload4(1))+1)
           i=0
           DO j=ind,ind+364
              i=i+1
              faparl_obs(i_year, i, i_site) = &
                SUM(iload1(x_a(i_site):x_e(i_site),y_a(i_site):y_e(i_site),j) / &
               &    iload2(x_a(i_site):x_e(i_site),y_a(i_site):y_e(i_site),j)**2) / &
               &SUM(1.0/iload2(x_a(i_site):x_e(i_site),y_a(i_site):y_e(i_site),j)**2)
              faparl_unc(i_year, i, i_site) = &
                SUM(1.0/iload2(x_a(i_site):x_e(i_site),y_a(i_site):y_e(i_site),j)) / &
               &SUM(1.0/iload2(x_a(i_site):x_e(i_site),y_a(i_site):y_e(i_site),j)**2)  
              IF (i==365) i=0
           ENDDO
        ENDIF
     ENDDO
     DEALLOCATE(iload1,iload2,iload3,iload4)
  ENDDO
   
  RETURN
  
END SUBROUTINE init_faparl

SUBROUTINE   init_faparg (ns, nrs, outt)
  ! gets FAPAR and unc from sites

  USE constants_pjr
  USE mo_namelist, ONLY: year0_site, year1_site, faparfile, faparfile_global, year0, year1
  USE mo_netcdf
  USE mo_calendar, ONLY: sdays
  
  IMPLICIT NONE
  INTEGER, INTENT(in) :: ns  ! number of sites
  INTEGER, INTENT(in) :: nrs  ! number of site simulation years
  INTEGER, INTENT(in) :: outt  ! time resolution
  
  ! local variables
  
  CHARACTER(len=80) :: in_sitedir ! internal copies of arguments if present
  INTEGER  i, j, k, l, n, t ! index variables

  CHARACTER (len=4) :: cy
  
  INTEGER, PARAMETER ::  inunit=2
  INTEGER i_site, i_year, i_month ! index variables
  INTEGER ind, year
  REAL r_year
!  INTEGER, DIMENSION(15,15) :: help
  
  LOGICAL :: fex

  TYPE (ncfile) :: infile
  TYPE (ncvar) :: nc_doy, nc_year, nc_fapar, nc_errfapar
  INTEGER, DIMENSION(3) :: iload_dims
  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: iload1, iload2
  REAL, ALLOCATABLE, DIMENSION(:) :: iload3, iload4


  ALLOCATE( faparg_obs( nrs, outt, ns))
  faparg_obs = 0.
  ALLOCATE( faparg_unc (nrs, outt, ns))
  faparg_unc = 1.e9
  
  in_sitedir = "./inputs/fapar_global/"
  
  infile%name=TRIM(in_sitedir)//TRIM(faparfile_global)
  nc_year%name='time'
  nc_fapar%name='fapar'
  nc_errfapar%name='errfapar'
  
  CALL ncopen(infile)
  
  CALL ncvarinfo(infile,nc_fapar,dimlens=iload_dims)
  
  ALLOCATE(iload1(iload_dims(1),iload_dims(2),iload_dims(3)))
  ALLOCATE(iload2(iload_dims(1),iload_dims(2),iload_dims(3)))
  ALLOCATE(iload3(iload_dims(3))) ! third dimension contains time information
  
  CALL ncread(infile,nc_fapar,iload1)
  CALL ncread(infile,nc_errfapar,iload2)
  CALL ncread(infile,nc_year,iload3)
  
  CALL ncclose(infile)
  
  t=1
  DO i_year = year0, year1
     DO i_month = 1, months_per_year
        r_year = REAL(i_year) +(i_month -0.5)/months_per_year
        IF (INT(r_year*100) >= INT(iload3(1)*100) .AND. INT(r_year*100) < INT(iload3(iload_dims(3))*100)) THEN
           n = 1
           DO i = iload_dims(2),1,-1
              DO j = 1,iload_dims(1)
                 IF (iload1(j,i,t)>=0.) THEN  
                    faparg_obs(i_year-year0+1,i_month,n) = iload1(j,i,t)
                    faparg_unc(i_year-year0+1,i_month,n) = iload2(j,i,t)
                    n=n+1
                 ENDIF
              ENDDO
           ENDDO
           t = t + 1
        ENDIF
     END DO
  END DO
  DEALLOCATE(iload1,iload2,iload3)
  
END SUBROUTINE init_faparg


! *********   cost_c returns costfunction/mismatch between 
!                  CALCULATED (model) and A PRIORI (observed) CONCENTRATIONs 
SUBROUTINE cost_c(fc, c) ! ccn part of cost function
  REAL :: fc ! returned cost function value
  REAL, DIMENSION(:) :: c ! input ccn vector

  fc = 0.5 * (SUM(((c -c_obs) / c_unc)**2))

  WRITE(6,*) 'cost_c: ',fc

END SUBROUTINE cost_c

! *********   cost_faparl returns costfunction/mismatch between LOCAL
!                  CALCULATED (model) and A PRIORI (observed) FAPAR
SUBROUTINE cost_faparl(fc, fapar, nrs, outt, ng) ! fapar part of cost function
  IMPLICIT NONE
  REAL :: fc ! returned cost function value
  INTEGER :: nrs, outt, ng
  REAL, DIMENSION(nrs,outt,ng) :: fapar ! input ccn vector

  INTEGER :: i,j,k

!!MS$  fc = 0.
!!MS$  do i=1,nrs
!!MS$     do j=1,outt
!!MS$        do k=1,1,ng
!!MS$           if (faparl_unc(i,j,k) .ne. 1.e9) then
!!MS$              fc = fc + 0.5 * ((fapar(i,j,k)-faparl_obs(i,j,k))/faparl_unc(i,j,k))**2
!!MS$           endif
!!MS$        enddo
!!MS$     enddo
!!MS$  enddo

  fc = 0.5 * (SUM(((fapar -faparl_obs) / faparl_unc)**2))

  WRITE(6,*) 'cost_faparl: ',fc

END SUBROUTINE cost_faparl

! *********   cost_faparg returns costfunction/mismatch between GLOBAL
!                  CALCULATED (model) and A PRIORI (observed) FAPAR
SUBROUTINE cost_faparg(fc, fapar, nrs, outt, ng) ! fapar part of cost function
  USE mo_grid, ONLY: rfapar0
  IMPLICIT NONE
  REAL :: fc ! returned cost function value
  INTEGER :: nrs, outt, ng, i, j, k
  REAL, DIMENSION(nrs,outt,ng) :: fapar ! input ccn vector
  REAL :: bias
!  fc = 0.5 * SUM(((fapar-2.92702*faparg_obs) / faparg_unc)**2)
!  fc = 0.5 * SUM(((fapar -faparg_obs) / faparg_unc)**2)
!  write (*,*) 'fc_faparg: ', fc

!  OPEN(unit=1, file='input/fapar_global/bias.txt', status='old')
!  READ(1,*) bias
!  close(1)
  

  fc = 0.
  do i=1,nrs
     do j=1,outt
        do k=1,ng
!           write (9,'(3I4,3E20.12)') i, j, k, fapar(i,j,k), faparg_obs(i,j,k), faparg_unc(i,j,k)
!           biasfac = fbiasa*(faparg_unc(i,j,k)-0.1)
!           biasfac = fbiasa+fbiasb*faparg_unc(i,j,k)
!           bias = (1.-rfapar0(k))*fapar(i,j,k)
           bias = 0.
           fc = fc + 0.5 * ((fapar(i,j,k)-bias-faparg_obs(i,j,k))/(faparg_unc(i,j,k)))**2
!           if (faparg_unc(i,j,k).ge.1e3) biasfac = 0.
!           fc = fc + 0.5 * ((fapar(i,j,k)-(faparg_obs(i,j,k)*(1.+biasfac)))/(faparg_unc(i,j,k)*(1.+biasfac)))**2
!           fc = fc + 0.5 * ((fapar(i,j,k)-(faparg_obs(i,j,k) + bias(i,j,k)))/faparg_unc(i,j,k))**2
        enddo
     enddo
  enddo

  WRITE(6,*) 'cost_faparg: ',fc

END SUBROUTINE cost_faparg


! *********   cost_flux returns costfunction/mismatch between 
!                  CALCULATED (model) and A PRIORI (observed) FLUXes at sites 
SUBROUTINE cost_flux(fc, flux, nrs, outt, ng) ! flux part of cost function
  IMPLICIT NONE 
  REAL :: fc				      ! return cost function value
  INTEGER :: nrs, outt, ng                    ! dimensions
  REAL :: flux(nrs,outt,ng)                   ! input flux vector
  INTEGER, PARAMETER :: uflux_obs = 10        ! flux observations
  LOGICAL :: lvar
  fc = 0.
  


END SUBROUTINE cost_flux


! *********   cost_p returns costfunction of 
!                  OPTIMIZED and A PRIORI PARAMETERs 
SUBROUTINE cost_p( fc, x ) ! cost of prior mismatch
  REAL :: fc ! returned cost function value
  REAL, DIMENSION(:) :: x ! md values in natural units
  INTEGER :: nx, np
! note that x0 and x0s are in physical units while x here is in natural units
! (multiples of unc) so that the calculation looks different from cost_c
!!$  WRITE(6,*) 'cost_p: variables'
!!$  WRITE(6,*) 'x:',x
!!$  WRITE(6,*) 'x0:',x0
!!$  WRITE(6,*) 'x0s:',x0s
!!$  WRITE(6,*) 'cost_p: end variables'

  np = SIZE(x0)
  nx = SIZE(x)
  fc =0
  
  do i=1,np
     if (xpf(i)>0) then
        fc = fc + 0.5 * ((x(i) - x0(i))**2)
     else
     endif
  enddo

!!MS$  if (nx == np) then
!!MS$     fc = 0.5 * (SUM((x - x0)**2))
!!MS$  elseif (nx > np) then
!!MS$     fc = 0.5 * (SUM((x(1:np) - x0)**2) + sum(((x(np+1:nx) - 0.)/frac_u)**2))
!!MS$  ELSE
!!MS$     print*,' wrong parameter numbers, nx = ',nx,' and np = ',np
!!MS$     stop
!!MS$  endif
!!MS$  WRITE(6,*) 'cost_p: variables'
!!MS$  WRITE(6,*) x
!!MS$  WRITE(6,*) x0
!!MS$  WRITE(6,*) x0s
!!MS$  WRITE(6,*) 'cost_p: end variables'

  WRITE(6,*) 'cost_p: ',fc

  RETURN
END SUBROUTINE cost_p

SUBROUTINE cost_pft( fc, x, frac_u) ! ccn part of cost function
  REAL, DIMENSION(:) :: x, frac_u ! input ccn vector
  
  fc = 0.5 * (SUM(((x - 0.) / frac_u)**2))

  WRITE(6,*) 'cost_pft: ',fc

  RETURN
END SUBROUTINE cost_pft

SUBROUTINE cost_wmax( fc, x, pasmmax_u) ! ccn part of cost function
  REAL :: fc ! returned cost function value
  REAL, DIMENSION(:) :: x, pasmmax_u ! input ccn vector
  
  fc = 0.5 * (SUM(((x - 0.) / pasmmax_u)**2))

  WRITE(6,*) 'cost_wmax: ',fc

  RETURN
END SUBROUTINE cost_wmax


! *********   cost_tot returns TOTAL costfunction, sum of 
!                  PARAMETER and CONCENTRATION mismatch
!SUBROUTINE cost_tot( fc, c, x) ! combined cost function
!  REAL :: fc ! returned cost function value
!  REAL, DIMENSION(:) :: x, c ! md values and calculated ccns
!  fc = cost_c(c) + cost_p(x)
!  
!  RETURN
!END SUBROUTINE cost_tot


SUBROUTINE getlsq( n, x, m,y )
  USE tm
  IMPLICIT NONE
  INTEGER :: n,m,nconc
  REAL    :: x(n), y(m)

  nconc = SIZE(conc)
  y(1:nconc) = (conc(:) -c_obs(:)) / c_unc(:)
  y(nconc+1:nconc+n) = x(:) - x0(:)

END SUBROUTINE getlsq

END MODULE costf
