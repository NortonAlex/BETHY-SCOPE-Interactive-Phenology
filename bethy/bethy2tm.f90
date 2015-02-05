MODULE bethy2tm
! module for transformation from bethy to tm grid
  REAL, ALLOCATABLE, DIMENSION(:,:) :: ag_mat
  REAL, ALLOCATABLE, DIMENSION(:, :, :)    :: netflux

CONTAINS

SUBROUTINE init_bethy2tm(filename)
  ! pjr modifies to use netcdf file 02/03/21
  USE constants_pjr
  USE mo_constants
  USE mo_netcdf
  USE mo_grid, ONLY: ng
  ! read the interpolation matrix
  ! pjr modifies to use netcdf files 02/03/21
  IMPLICIT NONE    
  CHARACTER(len=*) ::  filename

  ! local variables
  TYPE( ncfile) :: bethy2tm_file
  TYPE( ncvar) :: bethy2tm_var
  
  ALLOCATE ( ag_mat (n_lons*n_lats, ng)); ag_mat = 0.

  bethy2tm_file %name = filename
  CALL ncopen( bethy2tm_file)
  bethy2tm_var %name = 'bethy2tm'
  CALL ncread( bethy2tm_file, bethy2tm_var, ag_mat)
  CALL ncclose( bethy2tm_file)
!HEW-DEL031211  PRINT*,'bethy2tm: ',MINVAL( ag_mat), MAXVAL( ag_mat)
  RETURN
END SUBROUTINE init_bethy2tm


SUBROUTINE run_bethy2tm ( f_tm)
  USE constants_pjr
  USE mo_constants
  USE mo_grid, ONLY: ng

  ! takes fluxes on bethy grid defined over a period and maps them to tm grid
  ! note the curious dimensions, bethy has a list of coordinates but the 
  ! tm is a lat,lon array
  ! restricted transformation to non zero elements of ag_mat
  !                                  jan 03, FastOpt
  IMPLICIT NONE 
  REAL, INTENT(out), DIMENSION(:, :, :) :: f_tm

  ! local variables
  REAL, DIMENSION(ng) :: area
  REAL,ALLOCATABLE, DIMENSION(:,:,:) :: f_bethy
  real maxi, mini !FO
  INTEGER n_years, i_year , i, month, m, n, j, lat, lon, n_months, n_cells
  INTEGER, PARAMETER :: jdpy = 365

  n_years = SIZE(netflux,1) ! first dimension 
  n_months = SIZE(netflux,2) ! second dimension 
  n_cells = SIZE(netflux,3) ! third dimension 
  ALLOCATE(f_bethy(n_years,n_months,n_cells))
  f_bethy=netflux

!$taf store n_years = top_tape, rec=1

! .. for diagnostics 
!!MS$  OPEN(2,file='../input/area.lst',status='old')
!!MS$  DO i=1,ng_bethy
!!MS$     READ(2,*) area(i)
!!MS$  END DO
!!MS$  CLOSE(2)

! first transform from bethy time units to kg/m^2/yr
  DO i_year = 1, n_years

! .. diagnostic output
!!     WRITE(6,*) SUM(f_bethy( i_year, :, :)*area)/1e12 
    
     DO n = 1, ng
        f_bethy(i_year, :, n)=f_bethy(i_year, :, n)/float(rdays)*jdpy
     ENDDO

! .. diagnostic output
!!     DO n=1,ng_bethy
!!        WRITE(38,*) SUM( f_bethy(:, :, n))/12.
!!     END DO
  
  END DO

!FO begin for sparsity
  i = 0 ! index into f_tm array
!!  WRITE(6,*) 'run_bethy2tm: checking sums'
!$TAF loop = parallel
  DO i_year = 1, n_years 
!$TAF loop = parallel
     DO m = 1, months_per_year
        i = i + 1
!!MS$        f_tm (:, :, i) = RESHAPE (MATMUL( ag_mat, f_bethy(i_year, m, :)), &
!!MS$             & (/n_lons, n_lats/))
!!MS$ .. diagnostic output
!!MS$        print*,''
!!MS$        WRITE(6,*) 'old: ',SUM( f_tm(:, :, i))/1e12&
!!MS$             &,SUM(f_bethy( i_year, m, :) *area)/1e12 
        j = 0
        DO lat = 1, n_lats
           DO lon = 1, n_lons
              j=j+1
              f_tm(lon,lat,i)=0.
              DO n = 1, ng
                 if(ag_mat(j,n).gt.1.e5) then
                    f_tm(lon,lat,i)=f_tm(lon,lat,i)+ag_mat(j,n)*f_bethy(i_year,m,n)
                 endif
              ENDDO
           ENDDO
        ENDDO     
!!MS$ .. diagnostic output
!!MS$        WRITE(6,*) 'new: ',SUM( f_tm(:, :, i))/1e12&
!!MS$             &,SUM(f_bethy( i_year, m, :) *area)/1e12   
     END DO
  END DO

  DEALLOCATE(f_bethy)
  
  RETURN
END SUBROUTINE run_bethy2tm


END MODULE bethy2tm
