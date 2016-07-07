module mo_config
! .. Use Statements ..
  USE mo_constants
  USE mo_netcdf
  USE mo_grid
  USE mo_carparams
  USE mo_vegetation
! WOK-ADD-070704
  USE mo_pheno, ONLY: xph

  IMPLICIT NONE
  INTEGER,DIMENSION(ng,nv) :: vtype,help
  REAL,DIMENSION(ng,nv) :: vfrac
  INTEGER :: i1,i2

contains

SUBROUTINE init_config (npoint) 

USE mo_namelist, ONLY: nblocks,iblock 
USE break_jobs
 
! .. Arguments 
  integer :: npoint
! .. Locals
  INTEGER :: CTYPE, j, vl, k

  ! ANorton ...load initialization fields for fluorescence calcs
  ALLOCATE (vg_nv(vp))
  ALLOCATE (gridvp(ng,nv))


  ! Load initialization fields
  
  help=0
  WHERE (vtype > 0) help=1
  vp = SUM(help)

  ! Split veg-points into blocks 
       IF ((iblock==-1) .OR. (nblocks==-1)) THEN
!           print*,'  *  vp block par check 1 * :: vp=',vp
           i1 = 1
           i2 = vp
       else if( (nblocks .ge. 1) .and. (iblock .ge. 1) .and. (iblock .le. nblocks)) then
!           print*,'  *  vp block par check 2 * :: vp=',vp
           call get_splits(vp, nblocks, iblock, i1, i2)
       else
           stop 'iblock and/or nblocks not defined correctly...'
       end if
       write(*,*) 'veg-point block parallelisation: (i1,i2) =',i1,i2
       write(*,*) 'nblocks, iblock =', nblocks, iblock


  call config_allocate(vp,npoint)

  sumfrac = 0.
  vg=0
  vpb=0
  vpe=0
  k=0
  ca=355.e-6  ! mol co2/mol air
  tmade=.FALSE.

  !ANorton, for additional fluorescence fields.   
  gridvp = -1
  DO j=1,ng
     DO vl=1,nv
        CTYPE=vtype(j,vl)
        IF (ctype /= 0) THEN
           k=k+1
           vg_nv(k) = CTYPE
           gridvp(j,vl) = k
           frac(k)=vfrac(j,vl)
        ENDIF
     ENDDO
  ENDDO

k=0

  DO j=1,npoint
     DO vl=1,nv
        CTYPE=vtype(j,vl)  
        IF (ctype /= 0) THEN
           k=k+1
           gridp(k)=j                 
           vg(j)=vg(j)+1
           vpe(j)=k
           frac_p(k)=vfrac(j,vl)
! mas-add-081013 simple algorithm (0.1) to assign PFT frac uncertainties
           frac_u(k)= 0.1
           ph(k)=xph(CTYPE)
           if (ph(k)==0) write (*,*) CTYPE
           c4flg(k)=xc4flg(CTYPE)
           class(k)=xclass(CTYPE)
           hv(k)=xhv(CTYPE)/10.
! WOK-ADD-070629 array of specific leaf area
           sla(k)=slav(CTYPE)
           pft(k)=CTYPE
        ENDIF
     ENDDO
  ENDDO
  vpb=vpe-vg+1
  where(vg==0) vpb=vpe

!txk: MAS to check  IF (vp<npoint .OR. vp/=SUM(vg) .or. vp/=sp) THEN       
  IF (vp<npoint .OR. vp/=SUM(vg)) THEN       
     WRITE(6,*) 'initialize: fewer vegetation points than sites/grid points '
     WRITE(6,*) vp,SUM(vg),npoint
     STOP
  ENDIF

  
  
END SUBROUTINE init_config

SUBROUTINE config_sites (site_file, n_sites) 

  IMPLICIT NONE   

  ! .. Arguments 
  integer :: n_sites
  CHARACTER(len=*), INTENT(in) :: site_file
  ! .. Locals
  INTEGER :: i, j
  INTEGER :: n
  character(len=80) :: header

  OPEN (120,file=trim(site_file),status='unknown',form='formatted')
  READ (120,*) !header

  vtype = 0.
  DO i = 1, n_sites
     READ(120,'(2f7.2,i6,3i3,3f5.2)') lat(i), lon(i), elev(i), vtype(i,:), vfrac(i,:)
  ENDDO
  
  CLOSE(120)

END SUBROUTINE config_sites

SUBROUTINE config_global (grid_file) 

  IMPLICIT NONE   

  ! .. Arguments 
  CHARACTER(len=*), INTENT(in) :: grid_file

  ! .. Local Variables ..
  TYPE(ncfile) :: infile
  TYPE(ncvar) :: nclon 
  TYPE(ncvar) :: nclat
  TYPE(ncvar) :: ncelev
  TYPE(ncvar) :: ncvtype        
  TYPE(ncvar) :: ncfrac        
  TYPE(ncvar) :: ncrfapar0        
  REAL, ALLOCATABLE, DIMENSION (:) :: ilon
  REAL, ALLOCATABLE, DIMENSION (:) :: ilat
  REAL, ALLOCATABLE, DIMENSION(:,:) :: ielev
  INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: ivtype
  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: ifrac
  REAL, ALLOCATABLE, DIMENSION(:,:) :: irfapar0
  INTEGER :: i, j
  INTEGER :: n

  ! .. allocate temporary memory
  ALLOCATE (ilon(nlon),ilat(nlat))
  ALLOCATE (ielev(nlon,nlat))
  ALLOCATE (ivtype(nlon,nlat,nv),ifrac(nlon,nlat,nv))
  ALLOCATE (irfapar0(nlon,nlat))

  ! ..load data from netCDF
  infile%name=grid_file

  CALL ncopen(infile)
  ! get lons
  nclon%name='lon'
  CALL ncread (infile,nclon,ilon)
  ! get lats
  nclat%name='lat'
  CALL ncread (infile,nclat,ilat)
  ! get elevation
  ncelev%name='elev'
  CALL ncread (infile,ncelev,ielev)
  ! get vegetation type
  ncvtype%name='type'
  CALL ncread (infile,ncvtype,ivtype)   
  ! get fractional cover of vegetation type
  ncfrac%name='frac'
  CALL ncread (infile,ncfrac,ifrac)   
  ! get bias scaling factor
!MAS/TXK: allow simulations in case no scaling factor is available
  irfapar0 = 1.
!  ncrfapar0%name='ratio'
!  CALL ncread (infile,ncrfapar0,irfapar0)   

  vtype = 0.
  n = 1
  DO i = nlat,1,-1
     DO j = 1,nlon
        IF (ivtype(j,i,1)/=-9999) THEN  
           lon(n)=ilon(j)
           lat(n)=ilat(i)         
           elev(n)=ielev(j,i)
           vtype(n,:)=ivtype(j,i,:)
!           vfrac(n,:)=ifrac(j,i,:)
!          scale PFT fractions here:
           vfrac(n,:)=ifrac(j,i,:)*irfapar0(j,i)
           rfapar0(n)=irfapar0(j,i)
           n=n+1
        ENDIF
     ENDDO
  ENDDO

  ! .. check if number of grid points matches
  IF (n-1/=ng) THEN
    WRITE (*,*) 'Number of grid cells in input does not match grid information.'
    WRITE (*,*) 'n, ng: ', n, ng
    WRITE (*,*) 'Stopped in config_global.'
    STOP
  ENDIF

  DEALLOCATE (ilon,ilat,ielev,ivtype,ifrac)

END SUBROUTINE config_global

SUBROUTINE config_allocate(vp,ng)

  USE mo_grid, ONLY: gridp, frac, sumfrac, frac_p, frac_u, vg_nv, gridvp
  USE mo_vegetation, ONLY: class,c4flg,hv,ph,sla,pft

  IMPLICIT NONE
  INTEGER, INTENT(in) :: vp,ng
  ALLOCATE (class(vp),c4flg(vp),hv(vp),gridp(vp),frac(vp),frac_p(vp),frac_u(vp))
  ALLOCATE (sumfrac(ng),ph(vp),sla(vp),pft(vp))

END SUBROUTINE config_allocate

!SUBROUTINE config_deallocate
SUBROUTINE config_deallocate(vp,ng) ! note: use these if passing dimension to TAF

  USE mo_grid, ONLY: gridp, frac, frac_p, sumfrac
  USE mo_vegetation, ONLY: class,c4flg,hv,ph,sla,pft

  IMPLICIT NONE
  INTEGER, INTENT(in) :: vp,ng

!$taf next required = vp,ng
  DEALLOCATE (class,c4flg,hv,gridp,frac,frac_p,frac_u,sumfrac,ph,sla,pft)

END SUBROUTINE config_deallocate

end module mo_config
