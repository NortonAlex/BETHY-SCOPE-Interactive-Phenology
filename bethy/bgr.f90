MODULE bgr ! for background fluxes, basically ocean and fossil and land-use
! FastOpt 08/03: made attempt to avoid hard-wired directory names in code
  
  REAL, ALLOCATABLE, DIMENSION(:, :, :), SAVE :: f_bgr !background flux
  REAL, ALLOCATABLE, DIMENSION( :), SAVE :: time_bgr ! time for background flux
  INTEGER, SAVE :: n_direct_flux = 0 !  count direct flux params, to avoid recomp. in run_bgr
  
CONTAINS
  
  SUBROUTINE init_bgr ! read in background land and oc fluxes
    USE constants_pjr
    USE mo_constants
    USE  costf, ONLY :xpf
    USE mo_netcdf
    USE mo_namelist
    
    IMPLICIT NONE 
    
    
    ! locals
    
    INTEGER i
    ! pjr 050520: generalize fossil to a third dimension being the number of patterns
    REAL, DIMENSION(:,:, :), ALLOCATABLE :: fossil_multi
    REAL, DIMENSION (:,:,:), ALLOCATABLE ::  seas_field
    REAL, DIMENSION( :, :), ALLOCATABLE :: ann_field, luc
    REAL, DIMENSION( :), POINTER :: time => NULL()
    REAL, DIMENSION( :, :), POINTER :: fos_amp => NULL() 
    REAL, DIMENSION( :, :), POINTER :: luc_amp => NULL() 
    REAL :: x, t, fos_tot
    REAL, POINTER, DIMENSION(:) :: ocean_time
    REAL, ALLOCATABLE, DIMENSION(:,:,:) :: ocean_fluxes, tm2_fluxes
    INTEGER j,i1,j1, n_years,k
    CHARACTER( len=80) ::  head ! for reading graphz files
    INTEGER :: min_time ! index for start date into fossil amplitude array
    INTEGER :: i_time ! time index variable
    
    
    
    TYPE (ncfile) :: infile
    TYPE (ncvar) :: nc_ocean_lats, nc_ocean_lons, nc_ocean_time, nc_ocean_fluxes
    TYPE (ncvar) :: nc_fossil_fluxes, nc_luc_fluxes, nc_ocean_tak_fluxes
    
    INTEGER, DIMENSION( NF_MAX_DIMS) :: theLens
    TYPE (ncfile) :: flux_file
    TYPE (ncvar) :: flux_var
    REAL, ALLOCATABLE, DIMENSION(:, :, :) :: flux_pattern, flux_total
    INTEGER :: i_flux ! index variable into flux arrays
    REAL :: first_month ! starting date for each basis function
    INTEGER :: start_index, end_index ! array index corresponding to first_month
    INTEGER :: r_length, t_length
    LOGICAL :: fex
    
    !HEW-DEL031211    WRITE(6,*) 'inititializing background fluxes'
    ! FastOpt attempt to avoid hard-wired directory names in code
    
    r_length = (year1 -year0 +1)*months_per_year
    t_length = (yearin1 -yearin0 +1)*months_per_year
    ALLOCATE( f_bgr ( n_lons, n_lats, r_length)); f_bgr = 0.0
    ALLOCATE( time_bgr( (r_length)))
    time_bgr = (/ (REAL(year0) +(i-1.)/12.,i=1,SIZE( time_bgr)) /)
    
    ALLOCATE( luc( n_lons, n_lats)); luc = 0.0
    ALLOCATE( ann_field( n_lons, n_lats)); ann_field = 0.
    ALLOCATE (seas_field( n_lons, n_lats, months_per_year));    seas_field = 0.


!!! FOSSIL BACKGROUND !!!
    ! we now have a file containing multiple patterns and another scaling file for each of these
    infile %name = TRIM(bgrdir)//"fossil_multi.nc"
    CALL ncopen( infile)
    nc_fossil_fluxes %name = 'fossil_multi'
          CALL ncvarinfo( infile, nc_fossil_fluxes, dimlens=thelens)
          ALLOCATE( fossil_multi( thelens(1), thelens(2), thelens(3) ) )
    CALL ncread( infile, nc_fossil_fluxes, fossil_multi)
    CALL ncclose( infile) 
    ! now normalize it to 1 GtC so we can scale by the changing amplitude
    DO i_flux = 1, SIZE( fossil_multi,3)
       fossil_multi( :, :, i_flux) = 1.0e12 * fossil_multi( :, :, i_flux)/SUM( fossil_multi( :, :, i_flux))
    END DO
    

    ! get the time varying amplitude for both patterns
    ! FastOpt OPEN (1, file='../input/background/fossil_amp.txt', status='old') ! shifted to background dir for consistency 
    OPEN (1, file=TRIM(bgrdir)//"fossil_amp_multi.txt", status='old') ! shifted to background dir for consistency 
    REWIND 1
    CALL read_new( 1, head, time, fos_amp,t_length)
    ! now find first time that corresponds with our study period
    min_time = SUM( MINLOC( time, time > REAL( year0)))
    ! check that the fossil amplitudes stretch far enough to cover our study period
    IF( time( SIZE( time)) < REAL( year1)) THEN
       WRITE(0,*) 'init_bgr: not enough fossil amplitudes in ',TRIM(bgrdir)//"fossil_amp_multi.txt"
       WRITE(0,*) 'init_bgr: filling in with last value of ',fos_amp( SIZE( time),2),' for the year ',time(SIZE(time))
       DO i=SIZE(time)+1,SIZE(fos_amp,1)
          fos_amp(i,:)=fos_amp(SIZE(time),:)
       ENDDO
    ENDIF
    
    ! now add in the scaled paterns 
    ! the number of patterns is the minimum of the number of columns in fossil_multi_amp and the number of patterns in 
    !fossil_multi.nc,  note that we can change the number of columns by changing the second line of fossil_multi_amp
    ! note that we increment f_bgr just in case we add something in earlier in the code
    DO i_time = 1,  SIZE( f_bgr, 3)
       DO i_flux = 1, MIN( SIZE( fossil_multi, 3), SIZE( fos_amp, 2))
          f_bgr(:, :, i_time) = f_bgr( :, :,  i_time) +&
               & fos_amp( min_time + i_time -1, i_flux) * fossil_multi(:, :, i_flux)
       END DO
    END DO
    

!!! LAND-USE BACKGROUND !!! using 1990 spatial distribution and an interannual scaling factor from Hougton, 2008
    ! FastOpt    infile %name = '../input/background/luc_flux.nc'
    infile %name = TRIM(bgrdir)//"luc_flux.nc"
    CALL ncopen( infile)
    nc_luc_fluxes %name = 'luc_90'
    CALL ncread( infile, nc_luc_fluxes, luc)
    CALL ncclose( infile)
 ! now normalize it to 1 GtC so we can scale by the changing amplitude
    luc=luc/SUM(luc)

    ! get the time varying amplitude 
    OPEN (1, file=TRIM(bgrdir)//"luc_amp.txt", status='old') ! shifted to background dir for consistency 
    REWIND 1
    CALL read_new( 1, head, time, luc_amp, t_length)
    ! now find first time that corresponds with our study period
    min_time = SUM( MINLOC( time, time > REAL( year0)))
    ! check that the fossil amplitudes stretch far enough to cover our study period
    IF( time( SIZE( time)) < REAL( year1)) THEN
       WRITE(0,*) 'init_bgr: not enough luc amplitudes in ',TRIM(bgrdir)//"luc_amp.txt"
       WRITE(0,*) 'init_bgr: filling in with last value of ',luc_amp( SIZE( time),1),' for the year ',time(SIZE(time))
       DO i=SIZE(time)+1,SIZE(luc_amp,1)
          luc_amp(i,:)=luc_amp(SIZE(time),:)
       ENDDO
    ENDIF
    
    ! now add in the scaled paterns, only one pattern at the moment, therefor 3rd dimension is 1 
    ! note that we increment f_bgr just in case we add something in earlier in the code
    DO i_time = 1,  SIZE( f_bgr, 3)
       f_bgr(:, :, i_time) = f_bgr( :, :,  i_time) +&
            & 1.0e09 * luc_amp( min_time + i_time -1,1) * luc(:, :)
    END DO

!!! OCEAN BACKGROUND !!! climatology, takahashi fluxes
    ! FastOpt    infile %name = '../input/background/ocean_tak.nc'
    infile %name = TRIM(bgrdir)//"ocean_tak.nc"
    CALL ncopen( infile)
    nc_ocean_tak_fluxes %name = 'takahashi_flux'
    CALL ncread( infile, nc_ocean_tak_fluxes, seas_field)
    CALL ncclose( infile)
    DO i = 1, SIZE(f_bgr, 3), months_per_year
       f_bgr( :,:, i:i+months_per_year -1) = f_bgr( :, :, i:i+months_per_year-1) + seas_field
    END DO

    ! now get oc fluxes from lequere et al 2007
    infile %name = TRIM(bgrdir)//"ocean_fluxes-79-07.nc"
    CALL ncopen(infile)
    nc_ocean_fluxes%name = "flux"
    CALL ncvarinfo( infile, nc_ocean_fluxes, dimlens=thelens)
    ALLOCATE( ocean_fluxes( thelens(1), thelens(2), thelens(3) ) )
    ocean_fluxes=0.
    CALL ncread(infile, nc_ocean_fluxes, ocean_fluxes)
    ALLOCATE( ocean_time(thelens(3)))
    ocean_time(:)=0.0
    nc_ocean_time%name = "time"
    CALL ncread(infile, nc_ocean_time, ocean_time)

    ! get anomalies
    ALLOCATE( tm2_fluxes( thelens(1), thelens(2), thelens(3) ) )
    n_years = thelens(3)/12
    DO i= 1, 12       
       ann_field = SUM(ocean_fluxes(:, :, i:i+(n_years -1)*12:12), 3)/n_years ! average
       DO k = 0,n_years -1
          tm2_fluxes(:, :, i+k*12) = ocean_fluxes(:, :, i+k*12) -ann_field
       END DO
    END DO

    ! now scale,  it's in moles/second originally convert to kg/yr 
    tm2_fluxes = 365* 86400* 0.012 * tm2_fluxes 
    
    ! the fluxes *seem* to be from atm to oc , so reverse them
    tm2_fluxes = - tm2_fluxes
    ! and finally select the years we need
    ! first year
    i=1
    DO 
       IF (ocean_time(i)+0.001 > float(year0)) EXIT
       i=i+1
    ENDDO
    ! last year, either end of period or end of tm2_fluxes
    k = MIN( i + (year1 - year0 +1)*12 -1, SIZE( tm2_fluxes, 3))
    f_bgr(:,:,1:k-i+1) = f_bgr(:,:,1:k-i+1) +tm2_fluxes(:, :, i:k)
    
    DEALLOCATE( tm2_fluxes, ocean_time,ann_field, seas_field, fos_amp, luc_amp, & 
         & ocean_fluxes, time, luc, fossil_multi)

!!! WEAK CONSTRAINT !!!    
    ! pjr: count direct fluxes, try to avoid recomputation in run_bgr
    n_direct_flux = COUNT( xpf == direct_flux)
    IF( n_direct_flux > 0) THEN
       INQUIRE (file=flux_temp_file, exist=fex)
       IF (.NOT. fex) THEN
          ALLOCATE( flux_total( SIZE(f_bgr,1), SIZE(f_bgr, 2), SIZE(f_bgr, 3)))
          flux_file%name = TRIM( bgrdir)//TRIM( pattern_file)
          CALL ncopen( flux_file)
          
          ! now temporary file 
          OPEN(unit=1, file=flux_temp_file, form='unformatted', access='direct', &
               &recl=SIZE( flux_total)*KIND( flux_total))
          DO i_flux = 1, n_direct_flux
             flux_total = 0.
             flux_var%name(1:5) = "flux."
             WRITE( flux_var%name(6:10), '(i5.5)') i_flux
             CALL ncvarinfo( flux_file, flux_var, dimlens=thelens)
             ALLOCATE( flux_pattern( thelens(1), thelens(2), thelens(3)))
             CALL ncread( flux_file, flux_var, flux_pattern )
             CALL ncgetatt( flux_file, flux_var, "first_month", first_month) ! get the month corresponding to the start of this pattern
             
             ! now find the position in the array to starte  the patern, it's the closest time match in the time_bgr array 
             start_index = SUM( MINLOC( ABS( time_bgr - first_month)))
             end_index = MIN( SIZE( flux_total, 3), start_index + thelens(3) -1)
             flux_total( :, :, start_index:end_index) = flux_pattern(:, :, 1:end_index -start_index +1)
             WRITE(1, rec=i_flux) flux_total
             DEALLOCATE( flux_pattern)
          END DO
          CLOSE(1)
          CALL ncclose( flux_file)
       ELSE
          ! do nothing
       ENDIF
    ENDIF
    

  END SUBROUTINE init_bgr



  
  SUBROUTINE run_bgr (f, p)
  ! modified pjr 04/04/27 to deal with parameterized background fluxes

    USE mo_namelist
    IMPLICIT NONE 
    REAL, DIMENSION(:, :, :), INTENT(inout) :: f ! biospheric flux modified
    REAL, INTENT( in), DIMENSION (:) :: p ! amplitudes for pattern fluxes


! local variables
    INTEGER :: i_flux
    REAL, DIMENSION( SIZE(f, 1), SIZE(f, 2), SIZE(f, 3)) :: flux_pattern

    ! here to be total
    ! test version 
!!$    f = f_bgr
    ! real version below
    
    f = -f + f_bgr ! note sign change for fluxes *TO* atm

!!! ADD FLUX PATTERN !!!
    IF( n_direct_flux > 0) THEN
       OPEN(unit=1, file=flux_temp_file, form='unformatted', access='direct', &
            &recl=SIZE( f_bgr)*KIND( f_bgr))
       ! we have some parameterized flux patterns, add them
!$TAF loop = parallel
       DO i_flux = 1,n_direct_flux ! loop over all flux params
          READ(1, rec=i_flux) flux_pattern
           f = f + p( i_flux) * flux_pattern
       END DO
       CLOSE(1)
    ENDIF
    
    ! version without background fluxes
    !!$    f = -f
    RETURN

  END SUBROUTINE run_bgr
  
  
  SUBROUTINE read_graphz(n_unit, header, x, y)
    ! subroutine for reading graphz type file,
    ! pjr 04/04/27
    ! arguments
    IMPLICIT NONE 
    INTEGER n_unit ! n_unit number
    CHARACTER( len= *) :: header
    REAL, POINTER, DIMENSION( :) :: x
    REAL, POINTER, DIMENSION (:, :) :: y
    
    ! local variables
    INTEGER :: n_points, n_columns, i_point, i_column
    
    ! at the moment this dissociates the pointers if they're associated, to avoid memory leak
    ! TODO might be to check if things fit then just reuse
    IF( ASSOCIATED( x)) DEALLOCATE(x)
    IF( ASSOCIATED( y)) DEALLOCATE(y)
    READ( n_unit, *) header
    READ( n_unit,*) n_points, n_columns
    ALLOCATE( x( n_points))
    ALLOCATE( y( n_points, n_columns))
    DO i_point = 1, n_points
       READ( n_unit,*) x( i_point), y( i_point, :)
    END DO
    RETURN
    
  END SUBROUTINE read_graphz


  SUBROUTINE read_new(n_unit, header, x, y, r_length)
    ! subroutine for reading graphz type file,
    ! pjr 04/04/27
    ! arguments
    IMPLICIT NONE 
    INTEGER n_unit ! n_unit number
    INTEGER r_length ! n_points for actual run
    CHARACTER( len= *) :: header
    REAL, POINTER, DIMENSION( :) :: x
    REAL, POINTER, DIMENSION (:, :) :: y
    
    ! local variables
    INTEGER :: n_points, n_columns, i_point, i_column
    
    ! at the moment this dissociates the pointers if they're associated, to avoid memory leak
    ! TODO might be to check if things fit then just reuse
    IF( ASSOCIATED( x)) DEALLOCATE(x)
    IF( ASSOCIATED( y)) DEALLOCATE(y)
    READ( n_unit, *) header
    READ( n_unit,*) n_points, n_columns
    ALLOCATE( x( n_points))
    IF (n_points > r_length) THEN
       ALLOCATE( y( n_points, n_columns))
    ELSE
       ALLOCATE( y( r_length, n_columns))
    ENDIF
    DO i_point = 1, n_points
       READ( n_unit,*) x( i_point), y( i_point, :)
    END DO
    RETURN
    
  END SUBROUTINE read_new
  
  
END MODULE bgr

