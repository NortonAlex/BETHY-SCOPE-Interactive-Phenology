MODULE tm

  REAL,   POINTER, DIMENSION(:, :)        :: jac     ! mapping from flux 
                                                     ! field to data we're using
  INTEGER                                 :: obs_per_year
  REAL, DIMENSION (:, :, :), ALLOCATABLE  :: f_tm
  REAL, DIMENSION(:), ALLOCATABLE         :: conc
  INTEGER                                 :: totaldat

CONTAINS

  SUBROUTINE   init_tm( datfile, firstdat, lastdat, jacdir, n_stats)
    !  initializes  the jacobiann  for stations listed in datfile
    ! pjr modifies to use netcdf files 02/03/21
    USE constants_pjr
    USE mo_netcdf
    IMPLICIT NONE
    CHARACTER( len=*), INTENT(in) :: datfile ! list of stations we'll use
    INTEGER, INTENT(in) :: firstdat, lastdat, n_stats


    CHARACTER( len=*), OPTIONAL, INTENT(in) :: jacdir ! locations for jacobians


    ! local variables   
    CHARACTER(len=80) :: in_jacdir ! internal copies of arguments if present
    INTEGER  i, j, k, l,n ! index variables
    
    CHARACTER(len=80) :: jac_name, stat_name 
    INTEGER, PARAMETER :: jac_unit = 12
    INTEGER, PARAMETER ::  inunit=2
    INTEGER i_stat ! index variable for stations

    TYPE(ncfile) :: jac_file
    TYPE(ncvar) :: jac_var
    INTEGER, DIMENSION(NF_MAX_DIMS) :: theLens
    REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: in_jac

!HEW-DEL 031211    WRITE(6,*) 'initializing jacobians'
    IF( PRESENT( jacdir)) THEN
       in_jacdir = jacdir
    ELSE
       in_jacdir = "../input/jactd/"
    ENDIF
    
    obs_per_year = months_per_year * n_stats 
    ! handle allocation
    ! jacobian stores sensitivity of one year of ccns at each station to four
    ! years of fluxes from every gridpoint
    
    ALLOCATE( jac( n_lats * n_lons * months_per_year * 4, obs_per_year))
    jac = 0.
    
    OPEN(unit=inunit,file=datfile,status='old')
    REWIND inunit
    DO i_stat = 1, n_stats
       READ(inunit,*) stat_name
       jac_file %name = TRIM(in_jacdir)//"mat_"//TRIM(stat_name)//"_87_cg.nc"

       CALL ncopen( jac_file)
       jac_var %name = 'jacobian'
       CALL ncvarinfo( jac_file, jac_var, dimlens =thelens)
       IF( .NOT. ALLOCATED( in_jac)) ALLOCATE( in_jac( thelens(1), thelens(2), thelens(3), thelens(4)))
       CALL ncread( jac_file, jac_var, in_jac)
       jac(:, (i_stat -1)*months_per_year + 1:  i_stat*months_per_year) =&
            &RESHAPE( in_jac, (/ PRODUCT( thelens( 1:3)), thelens(4) /))
       CALL ncclose( jac_file)
    END DO
    DEALLOCATE( in_jac)
    totaldat = lastdat -firstdat +1
    CALL tm_allocate(n_stats)
    RETURN
  END SUBROUTINE init_tm
  

  SUBROUTINE tm_allocate(n_stats)
    USE constants_pjr
    USE mo_namelist
    USE mo_constants
    USE mo_grid, ONLY: ng
    
    IMPLICIT NONE 
    
    INTEGER, INTENT(in) :: n_stats

    ! allocate f_tm here to make taf happy I hope
    ALLOCATE(f_tm( n_lons, n_lats, months_per_year*totaldat))
    f_tm = 0.0
    ! allocate c here to make TAF happy I hope
    ALLOCATE( conc(n_stats *  months_per_year*totaldat))
    conc = 0.0
    
  END SUBROUTINE tm_allocate


  SUBROUTINE tm_deallocate
    USE constants_pjr
    IMPLICIT NONE 
    DEALLOCATE(f_tm,conc)
  END SUBROUTINE tm_deallocate


  SUBROUTINE run_tm (offset, f, c)
    USE constants_pjr
    IMPLICIT NONE 
    REAL, INTENT(in) :: offset ! base level ccn 
    REAL, INTENT(in), DIMENSION(:,:,:) :: f ! fluxes as x,y,t array
!MAS/TXK: activate for transport test  REAL, DIMENSION(:,:,:) :: f ! fluxes as x,y,t array
    REAL, INTENT(inout), DIMENSION(:) :: c ! ccns as single fector  
    
   
    ! some local variables

    INTEGER i, j, k
    REAL, PARAMETER :: convfac = 1e-12/2.123 ! kilograms to ppmv 
    INTEGER i_year, n_year
    INTEGER n_months, month, m
    INTEGER, PARAMETER :: jdpyear = 365
    INTEGER, DIMENSION(12), PARAMETER :: rdays = &
         & (/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)
    REAL :: comp
    INTEGER :: nlon, nlat, ibase, kbase, l
    
!MAS/TXK: test output
! the purpose of this part of the code is to check the implementation
! of the Jacobian*flux product
! to activate it just remove comments below
! and above where f is declared
!       integer iflux,ig,imonf,nmonf,ic,nfpy
!       i= 20
!       j=15
!       imonf = 6
!       nmonf = 12
!       ic = 2+12
!       ig=i+(j-1)*(n_lons)
!       iflux=ig+(imonf-1)*n_lons * n_lats
!       nfpy = n_lons * n_lats * nmonf
! look at an individual component of the Jacobian
!       print*, 'i,j,imonf,ic', i,j,imonf,ic
!       print*, 'matrix, background flux', jac(iflux,ic),f(i, j, imonf)
!       print*, 'matrix, background flux y2', jac(iflux+nfpy,ic),f(i, j, imonf+12)
!       print*, 'matrix, background flux y3', jac(iflux+2*nfpy,ic),f(i, j, imonf+2*12)
!       print*, 'matrix, background flux y4', jac(iflux+3*nfpy,ic),f(i, j, imonf+3*12)
!       f = 0
! first test constant flux into a single grid cell (1GTC/year)        
!       do k = 1, months_per_year*totaldat
!            f(i,j,k) = 1.e12*2.123
!       enddo
! second test: pulse in January (1GTC in each year)     
 !      do k = 1, months_per_year*totaldat, 12
 !          f(i,j,k) = 1.e12*2.123*jdpyear/rdays(1)
 !      enddo
!MAS/TXK end test output

    c = offset ! start with base level 

   
!$TAF loop = parallel
    DO i = SIZE( c) -obs_per_year, 0, -obs_per_year
       i_year = i / obs_per_year + 1 !counting back the years
! cutoff if not 4 consecutive years of fluxes available, eg 1979-1982
       n_year = MIN( i_year, 4) ! how many years of fluxes are involved 

       ibase = (4 - n_year)*n_lons*n_lats*months_per_year
       kbase = (i_year-n_year)*months_per_year !starting time index for fluxes
       DO j = 1, obs_per_year
          comp = 0.0
          DO k = 1, n_year*months_per_year
             DO nlat = 1, n_lats
                DO nlon = 1, n_lons
                   l = nlon + (nlat-1)*n_lons + (k-1)*n_lons*n_lats
                   comp = comp + jac(ibase+l,j) * f(nlon, nlat, kbase+k)
                END DO
             END DO
          END DO
          c(i+j) = c(i+j) + comp
       END DO

       ! with luck the transpose means we're using the efficient dimension 
       ! of both matrices
       ! now do the bit before the start of the Jacobian
       DO m = 1, (i_year -n_year)*months_per_year
          month = MOD(m -1, months_per_year) + 1
          c (i+1: i+obs_per_year) = c (i+1: i+obs_per_year) +&
               & SUM( f(:, :, m)) * convfac *rdays(month) /jdpyear
       END DO
    END DO

    RETURN
  END SUBROUTINE run_tm
  
END MODULE tm



  
