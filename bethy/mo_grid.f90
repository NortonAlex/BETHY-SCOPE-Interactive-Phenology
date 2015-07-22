MODULE mo_grid

!CCCC surface parameters (not modified)
!     declarations /psurface/

  IMPLICIT NONE

! .. Grid declaration
!!MS$  INTEGER, PARAMETER :: ng = 62483  ! 0.5 degree regular grid
!!MS$  INTEGER :: nlon = 720
!!MS$  INTEGER :: nlat = 360
!!MS$  INTEGER, PARAMETER :: ng = 11069  ! 1 deg equal area grid
!!MS$  INTEGER :: nlon = 360
!!MS$  INTEGER :: nlat = 180

!--------------------------------------------
 INTEGER, ALLOCATABLE, DIMENSION(:) :: vg_nv
 INTEGER, ALLOCATABLE, DIMENSION(:,:) :: gridvp
!--------------------------------------------

 INTEGER, PARAMETER :: ng = 170   ! TM2 grid
 INTEGER :: nlon = 36
 INTEGER :: nlat = 24


!cccc ng         number of land grid points !cccc nlon       number of longitudes        >  grid resolution dependant

!cccc nlat       number of latitudes        / 

! .. Local Arrays ..

  INTEGER, DIMENSION (ng) :: elev, vg, vpb, vpe, rfapar0
  INTEGER, ALLOCATABLE, DIMENSION(:) :: gridp
  REAL, ALLOCATABLE, DIMENSION(:) :: frac, frac_p, frac_u, sumfrac
  REAL, DIMENSION (ng) :: lon, lat  

!cccc elev     elevation [m]
!cccc vg       vegetation types per gridcell
!cccc vpb      index to first sub-pixel (1...vp) within grid cell
!cccc vpe      index to last sub-pixel (1...vp) within grid cell

! .. Local Scalars ..
    INTEGER :: vp, sp, nrs

!cccc vp        for grid, total number of vegetation points=sum(pfts per gridcell)
!cccc sp        for sites, total number of vegetation points=sum(pfts per sites)
!cccc sp        for sites, total number simulation years

END MODULE mo_grid
