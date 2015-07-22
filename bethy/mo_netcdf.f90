!$TAF FUNCTION nf_get_vara_real  INPUT = 1,2,3,4
!$TAF FUNCTION nf_get_vara_real OUTPUT =         5

!$TAF MODULE mo_netcdf SUBROUTINE ncopen INPUT = 1,2
!$TAF MODULE mo_netcdf SUBROUTINE ncopen OUTPUT = 1

!$TAF MODULE mo_netcdf SUBROUTINE ncclose INPUT = 1
!$TAF MODULE mo_netcdf SUBROUTINE ncclose OUTPUT = 1

!$TAF MODULE mo_netcdf SUBROUTINE get_var_int1d  INPUT = 1,2  ,4,5,6,7
!$TAF MODULE mo_netcdf SUBROUTINE get_var_int1d OUTPUT =   2,3
!$TAF MODULE mo_netcdf SUBROUTINE get_var_int2d  INPUT = 1,2  ,4,5,6,7
!$TAF MODULE mo_netcdf SUBROUTINE get_var_int2d OUTPUT =   2,3
!$TAF MODULE mo_netcdf SUBROUTINE get_var_int3d  INPUT = 1,2  ,4,5,6,7
!$TAF MODULE mo_netcdf SUBROUTINE get_var_int3d OUTPUT =   2,3
!$TAF MODULE mo_netcdf SUBROUTINE get_var_int4d  INPUT = 1,2  ,4,5,6,7
!$TAF MODULE mo_netcdf SUBROUTINE get_var_int4d OUTPUT =   2,3

!$TAF MODULE mo_netcdf SUBROUTINE get_var_real0d  INPUT = 1,2
!$TAF MODULE mo_netcdf SUBROUTINE get_var_real0d OUTPUT =     3

!$TAF MODULE mo_netcdf SUBROUTINE get_var_real1d  INPUT = 1,2  ,4,5,6,7
!$TAF MODULE mo_netcdf SUBROUTINE get_var_real1d OUTPUT =   2,3
!$TAF MODULE mo_netcdf SUBROUTINE get_var_real2d  INPUT = 1,2  ,4,5,6,7
!$TAF MODULE mo_netcdf SUBROUTINE get_var_real2d OUTPUT =   2,3
!$TAF MODULE mo_netcdf SUBROUTINE get_var_real3d  INPUT = 1,2  ,4,5,6,7
!$TAF MODULE mo_netcdf SUBROUTINE get_var_real3d OUTPUT =   2,3
!$TAF MODULE mo_netcdf SUBROUTINE get_var_real4d  INPUT = 1,2  ,4,5,6,7
!$TAF MODULE mo_netcdf SUBROUTINE get_var_real4d OUTPUT =   2,3

!$TAF MODULE mo_netcdf SUBROUTINE get_var_double0d  INPUT = 1,2
!$TAF MODULE mo_netcdf SUBROUTINE get_var_double0d OUTPUT =     3

!$TAF MODULE mo_netcdf SUBROUTINE get_var_double1d  INPUT = 1,2  ,4,5,6,7
!$TAF MODULE mo_netcdf SUBROUTINE get_var_double1d OUTPUT =   2,3
!$TAF MODULE mo_netcdf SUBROUTINE get_var_double2d  INPUT = 1,2  ,4,5,6,7
!$TAF MODULE mo_netcdf SUBROUTINE get_var_double2d OUTPUT =   2,3
!$TAF MODULE mo_netcdf SUBROUTINE get_var_double3d  INPUT = 1,2  ,4,5,6,7
!$TAF MODULE mo_netcdf SUBROUTINE get_var_double3d OUTPUT =   2,3
!$TAF MODULE mo_netcdf SUBROUTINE get_var_double4d  INPUT = 1,2  ,4,5,6,7
!$TAF MODULE mo_netcdf SUBROUTINE get_var_double4d OUTPUT =   2,3

!$TAF MODULE mo_netcdf SUBROUTINE ncVarInfo  INPUT = 1,2
!$TAF MODULE mo_netcdf SUBROUTINE ncVarInfo OUTPUT =   2,3,4,5
!$TAF FUNCTION nf_get_varm_real  INPUT = 1,2,3,4,5,6
!$TAF FUNCTION nf_get_varm_real OUTPUT =             7

!$TAF FUNCTION nf_get_vars_real  INPUT = 1,2,3,4,5
!$TAF FUNCTION nf_get_vars_real OUTPUT =           6


MODULE mo_netcdf
  ! **********************************************************************
  !
  ! Name
  !   mo_netcdf
  !
  ! Description
  !   Defines general procedures for data read and write to netCDF files.
  !   The variety of original netCDF routines for reading and writing data
  !   are gathered within two (heavily) overloaded procedures, ncread and
  !   ncwrite. Scalars and arrays up to 4 dimensions (type real and integer)
  !   are handled. File opening and closing is done by ncopen and ncclose,
  !   respectively. New variable types for netCDF files (ncfile), netCDF
  !   variables (ncvar) and netCDF dimensions (ncdim) are defined. All
  !   procedures use error handling implicitly. Any error results in program
  !   stop.
  !
  ! Procedure
  !   The original netCDF procedures, defined in libnetcdf.a are collected
  !   withing interfaces. The problem: One procedure for each kind and rank
  !   of data array is needed. count, stride and imap parameters are
  !   handeled as optional arguments. If the netCDF procedures produce an
  !   error the handle_ncerr subroutine is called which stops program
  !   execution and prints an error message.
  !
  !   Each netCDF property (like files and variables) is represented by
  !   a type collecting all important information. 
  !
  ! EXAMPLES
  !
  !   1. Defining a new netCDF file from scratch:
  !      A netCDF file savdat.nc containing the variable fields hammer (3D) and
  !      nail(2D) and the corresponding dimension variables x, y and z. The
  !      fields are of size (nx,ny,nz) and (nx,ny), respectively. Write
  !      some data to the file.
  !
  !  use mo_netcdf
  !
  !  [...]
  !
  !  type(ncFile) :: DataFile
  !  type(ncDim), TARGET :: xDim, yDim, zDim
  !  type(ncVar)  :: xVar, yVar, zVar
  !  type(ncVar)  :: hammerVar, nailVar
  !
  !  real, dimension(nx) :: xDat
  !  real, dimension(ny) :: yDat
  !  real, dimension(nz) :: zDat
  !
  !  real, dimension(nx,ny,nz) :: hammer
  !  real, dimension(nx,ny) :: nail
  !
  !  dataFile %name = "savdat.nc"
  !
  !  xDim %name = "x"
  !  xDim %len  = nx
  !
  !  yDim %name = "y"
  !  yDim %len  = ny
  !
  !  zDim %name = "z"
  !  zDim %len  = nz
  !
  !  xVar %name = "x"
  !  xVar %xType= NF_FLOAT
  !  xVar %nDims = 1
  !  xVar %Dims(1) %p => xDim
  !
  !  yVar %name = "y"
  !  yVar %xType= NF_FLOAT
  !  yVar %nDims = 1
  !  yVar %Dims(1) %p => yDim
  !
  !  zVar %name = "z"
  !  zVar %xType= NF_FLOAT
  !  zVar %nDims = 1
  !  zVar %Dims(1) %p => zDim
  !
  !  hammerVar %name = "hammer"
  !  hammerVar %xType= NF_FLOAT
  !  hammerVar %nDims = 3
  !  hammerVar %Dims(1) %p => xDim
  !  hammerVar %Dims(2) %p => yDim
  !  hammerVar %Dims(3) %p => zDim
  !
  !  nailVar %name = "nail"
  !  nailVar %xType= NF_FLOAT
  !  nailVar %nDims = 2
  !  nailVar %Dims(1) %p => xDim
  !  nailVar %Dims(2) %p => yDim
  !
  !  call ncCreate(dataFile)
  !  call ncPutAtt(dataFile,'creator','Schorsch')
  !
  !  call ncDimDef(dataFile, xDim)
  !  call ncDimDef(dataFile, yDim)
  !  call ncDimDef(dataFile, zDim)
  !
  !  call ncVarDef(dataFile, xVar)
  !  call ncVarDef(dataFile, yVar)
  !  call ncVarDef(dataFile, zVar)
  !
  !  call ncVarDef(dataFile,hammerVar)
  !  call ncPutAtt(dataFile,hammerVar,'long_name','der absolute Hammer')
  !  call ncPutAtt(dataFile,hammerVar,'units','gramm')
  !
  !  call ncVarDef(dataFile,nailVar)
  !  call ncPutAtt(dataFile,nailVar,'tag',12.0)
  !
  !  call ncEndDef(dataFile)
  !
  !  call ncPrint(dataFile, xVar, xDat)
  !  call ncPrint(dataFile, yVar, yDat)
  !  call ncPrint(dataFile, zVar, zDat)
  !  call ncPrint(dataFile, hammerVar, hammer)
  !  call ncPrint(dataFile, nailVar, nail)
  !  
  !  call ncClose(dataFile)
  !
  !  
  !  2. Read data from a netCDF file
  !
  !
  !  use mo_netcdf
  !
  !  [...]
  !
  !  type(ncFile), TARGET :: DataFile
  !  type(ncVar)  :: xVar, yVar, zVar
  !  type(ncVar)  :: hammerVar, nailVar
  !
  !  real, dimension(:), allocatable :: ReadX
  !  real, dimension(:), allocatable :: ReadY
  !  real, dimension(:), allocatable :: ReadZ
  !  real, dimension(:,:,:), allocatable :: ReadHammer
  !  real, dimension(:,:), allocatable :: ReadNail
  !
  !  integer, dimension(NF_MAX_DIMS) :: theLens
  !
  !  dataFile %name = "savdat.nc"
  !
  !  xVar %name = "x"
  !  yVar %name = "y"
  !  zVar %name = "z"
  !  hammerVar %name = "hammer"
  !  nailVar %name = "nail"
  !
  !  call ncOpen(dataFile)
  !
  !  call ncVarInfo(dataFile,xVar,DimLens=theLens)
  !  allocate(ReadX(theLens(1)))
  !  call ncRead(dataFile,xVar,ReadX)
  !     
  !  call ncVarInfo(dataFile,yVar,DimLens=theLens)
  !  allocate(ReadY(theLens(1)))
  !  call ncRead(dataFile,yVar,ReadY)
  !     
  !  call ncVarInfo(dataFile,zVar,DimLens=theLens)
  !  allocate(ReadZ(theLens(1)))
  !  call ncRead(dataFile,zVar,ReadZ)
  !     
  !  call ncVarInfo(dataFile,hammerVar,DimLens=theLens)
  !  allocate(ReadHammer(theLens(1),theLens(2),theLens(3)))
  !  call ncRead(dataFile,hammerVar,ReadHammer)
  !     
  !  call ncVarInfo(dataFile,nailVar,DimLens=theLens)
  !  allocate(ReadNail(theLens(1),theLens(2)))
  !  call ncRead(dataFile,nailVar,ReadNail)
  !
  !  call ncClose(dataFile)
  !
  !  [...]
  !  
  !  deallocate(ReadNail,ReadHammer,ReadZ,ReadY,ReadX)
  !
  !      
  !
  ! AUTHOR
  !   Georg Baeuml, MPI Met, Hamburg <baeuml@dkrz.de>
  ! 
  ! Version control information:
  !   $Date: 2009-02-21 04:29:52 +1100 (Sat, 21 Feb 2009) $
  !   $Revision: 1135 $
  !   $Name$
  !   $Source$
  ! 
  ! **********************************************************************

  IMPLICIT NONE

  INCLUDE 'netcdf.inc'

  ! TODO alle Leseroutinen unterscheiden sich in der Anzahl ihrer Argumente
  ! TODO Annahme: Keine Typ-Konversion erwünscht.

  ! ************************************************************
  !
  ! ncread(file, var, value, start, count, stride, imap)
  !
  ! ************************************************************

  INTERFACE ncread
!!$     MODULE PROCEDURE nf_get_var_text
!!$     MODULE PROCEDURE nf_get_var_int1
!!$     MODULE PROCEDURE nf_get_var_int2
     MODULE PROCEDURE get_var_int1d
     MODULE PROCEDURE get_var_int2d
     MODULE PROCEDURE get_var_int3d
     MODULE PROCEDURE get_var_int4d
     MODULE PROCEDURE get_var_real0d
     MODULE PROCEDURE get_var_real1d
     MODULE PROCEDURE get_var_real2d
     MODULE PROCEDURE get_var_real3d
     MODULE PROCEDURE get_var_real4d
     MODULE PROCEDURE get_var_double0d
     MODULE PROCEDURE get_var_double1d
     MODULE PROCEDURE get_var_double2d
     MODULE PROCEDURE get_var_double3d
     MODULE PROCEDURE get_var_double4d
  END INTERFACE

  INTERFACE ncprint
     MODULE PROCEDURE put_var_real0d
     MODULE PROCEDURE put_var_real1d
     MODULE PROCEDURE put_var_real2d
     MODULE PROCEDURE put_var_real3d
     MODULE PROCEDURE put_var_real4d
     MODULE PROCEDURE put_var_double0d
     MODULE PROCEDURE put_var_double1d
     MODULE PROCEDURE put_var_double2d
     MODULE PROCEDURE put_var_double3d
     MODULE PROCEDURE put_var_double4d
  END INTERFACE

  INTERFACE ncGetAtt
!!$     MODULE PROCEDURE ncGetAtt_text
     MODULE PROCEDURE ncGetAtt_real
     MODULE PROCEDURE ncGetAtt_double
!!$     MODULE PROCEDURE ncGetAtt_dblarr
!!$     MODULE PROCEDURE ncGlobalGetAtt_text
!!$     MODULE PROCEDURE ncGlobalGetAtt_real
!!$     MODULE PROCEDURE ncGlobalGetAtt_double
!!$     MODULE PROCEDURE ncGlobalGetAtt_dblarr
  END INTERFACE

  INTERFACE ncPutAtt
     MODULE PROCEDURE ncPutAtt_text
     MODULE PROCEDURE ncPutAtt_real
     MODULE PROCEDURE ncPutAtt_double
     MODULE PROCEDURE ncPutAtt_dblarr
     MODULE PROCEDURE ncGlobalPutAtt_text
     MODULE PROCEDURE ncGlobalPutAtt_real
     MODULE PROCEDURE ncGlobalPutAtt_double
     MODULE PROCEDURE ncGlobalPutAtt_dblarr
  END INTERFACE

  INTERFACE ncPresent
     MODULE PROCEDURE present_file
     MODULE PROCEDURE present_var
  END INTERFACE
  

  TYPE ncfile
     INTEGER :: id=0              ! ID of opened file
     CHARACTER (len=200) :: name  ! path to file
     LOGICAL :: opened=.FALSE.    ! .true. for opened, .false. for closed
     LOGICAL :: defmode=.FALSE.   ! .true. if in DEFINE mode
  END TYPE ncfile

  TYPE ncDim
     INTEGER :: id=0                     ! ID of dimension
     CHARACTER (len=NF_MAX_NAME) :: name ! name in netCDF file
     INTEGER :: len=0                    ! length of dimension
  END TYPE ncDim

  TYPE ncDimPtr
     TYPE(ncDim), POINTER :: p
  END TYPE ncDimPtr

  TYPE ncVar
     INTEGER :: id=0                     ! ID of variable
     CHARACTER (len=NF_MAX_NAME) :: name ! name in netCDF file
     INTEGER :: xtype                    ! external type
     ! NF_BYTE, NF_CHAR, NF_SHORT, NF_INT, NF_FLOAT, NF_DOUBLE
     INTEGER :: nDims                    ! number of dimensions
     TYPE(ncDimPtr), DIMENSION(NF_MAX_DIMS) :: dims
  END TYPE ncvar

  TYPE(ncdim), DIMENSION(0), SAVE :: ncscalar

  INTEGER, PRIVATE :: ncstatus ! use as status variable for netCDF commands
  INTEGER, PARAMETER, PRIVATE :: dblkind=KIND(1.0D0)
                               ! kind number of double precission
  
CONTAINS

  SUBROUTINE ncVarInfo(file,var,dimIDs,DimLens,IsDef)
    ! This routine gets the information about var. The name of var has to be
    ! known. The id and nDims is filled in by ncVarInfo. If supplied, IsDef
    ! is set to true if the variable has been found in file, false otherwise.
    ! The optional argument dimIDs contains the IDs of the dimensions used
    ! (unused is set to -1). DimLens holds the dimension lenghts.
    ! If the size of DimLens or dimIDs is not large enough (at least nDim)
    ! all entries are set to -1.
    TYPE(ncFile), INTENT(in) :: file
    TYPE(ncVar), INTENT(inout) :: var
    LOGICAL, INTENT(out), OPTIONAL :: IsDef
    INTEGER, DIMENSION(:), INTENT(out), OPTIONAL :: dimIDs
    INTEGER, DIMENSION(:), INTENT(out), OPTIONAL :: DimLens

    LOGICAL :: theStatus
    INTEGER, DIMENSION(NF_MAX_DIMS) :: theDimIDs
    INTEGER :: nAtts

    INTEGER :: jd

    theStatus = ncPresent(file,var)
    IF (.NOT. theStatus) THEN
       var %id = -1
       var %xType = -1
       var %nDims = -1
       IF(PRESENT(IsDef)) THEN
          IsDef = .FALSE.
       END IF
       RETURN
    ELSE
       IF(PRESENT(IsDef)) THEN
          IsDef = .TRUE.
       END IF
    END IF
    
    ncstatus = nf_inq_varid(file %id, var %name, var %id)
    IF (ncstatus .NE. NF_NOERR) THEN
       CALL handle_ncerr('ncVarInfo',ncstatus)
    END IF

    ncstatus = nf_inq_var(file %id, var %id, var %name, var %xType, &
         & var %nDims, theDimIDs, nAtts)
    IF (ncstatus .NE. NF_NOERR) THEN
       CALL handle_ncerr('ncVarInfo',ncstatus)
    END IF

    IF(PRESENT(dimIDs)) THEN
       dimIDs(:) = -1
       IF(SIZE(dimIDs) .GE. var %nDims) THEN
          dimIDs(1:var %nDims) = theDimIDs(1:var %nDims)
       END IF
    END IF
    
    IF(PRESENT(dimLens)) THEN
       dimLens(:) = -1
       IF(SIZE(dimLens) .GE. var %nDims) THEN
          DO jd=1,var %nDims
             ncstatus = nf_inq_dimlen(file %id, theDimIDs(jd), dimLens(jd))
             IF (ncstatus .NE. NF_NOERR) THEN
                CALL handle_ncerr('ncVarInfo',ncstatus)
             END IF
          END DO
       END IF
    END IF
  END SUBROUTINE ncVarInfo
  
  SUBROUTINE ncopen(file, omode)
    TYPE(ncfile), INTENT(inout) :: file
    INTEGER, OPTIONAL, INTENT(in) :: omode
    INTEGER :: theomode

    IF (.NOT. PRESENT(omode)) THEN
       theomode=NF_NOWRITE
    ELSE
       theomode=omode
    END IF
    
    ncstatus = nf_open(file%name, theomode, file%id)
    IF (ncstatus .NE. NF_NOERR) THEN
       file%opened=.FALSE.
       CALL handle_ncerr('ncopen',ncstatus)
    ELSE
       file%opened=.TRUE.
    END IF

  END SUBROUTINE ncopen

  SUBROUTINE ncclose(file)
    TYPE(ncfile), INTENT(inout) :: file

    IF (file%opened .EQV. .FALSE.) THEN
       PRINT *, 'attempt to close already closed file ', file%name
       RETURN
    END IF

    ncstatus = nf_close(file%id)
    IF (ncstatus .NE. NF_NOERR) THEN
       CALL handle_ncerr('ncclose',ncstatus)
    ELSE
       file%opened = .FALSE.
       file%defmode=.FALSE.
    END IF
  END SUBROUTINE ncclose

  SUBROUTINE ncDimLen(file, dim)
    TYPE(ncfile), INTENT(in) :: file
    TYPE(ncdim), INTENT(inout) :: dim

    ncstatus = nf_inq_dimid(file%id, dim%name, dim%id)
    IF (ncstatus .NE. NF_NOERR) THEN
       CALL handle_ncerr('ncdimlen',ncstatus)
    END IF

    ncstatus = nf_inq_dimlen(file%id, dim%id, dim%len)
    IF (ncstatus .NE. NF_NOERR) THEN
       CALL handle_ncerr('ncdimlen',ncstatus)
    END IF
  END SUBROUTINE ncDimLen


  SUBROUTINE handle_ncerr(routine,status)
    CHARACTER(*) :: routine            ! calling routine
    INTEGER, INTENT(in) :: status       ! error status
    IF (status .NE. NF_NOERR) THEN
       PRINT *, 'netCDF error in ', TRIM(routine), ':'
       PRINT *, nf_strerror(status)
!       STOP 'stopped: netCDF error'
    END IF
  END SUBROUTINE handle_ncerr

  ! **************************************************
  ! Subroutines for the ncread interface.
  ! Common procedure: inquire VARID and get specified data
  ! One subroutine for each variable type and 2, 3 or 4 dimensional array.
  ! **************************************************

  ! INTEGER
  
  SUBROUTINE get_var_int1d(file, var, value, start, count, stride, imap)

    IMPLICIT NONE

    INTEGER, PARAMETER :: ncdim=1

    TYPE(ncfile), INTENT(in) :: file
    TYPE(ncvar), INTENT(inout) :: var
    INTEGER, DIMENSION(:), INTENT(out) :: value
    INTEGER, DIMENSION(ncdim), OPTIONAL, INTENT(in) :: start
    INTEGER, DIMENSION(ncdim), OPTIONAL, INTENT(in) :: count
    INTEGER, DIMENSION(ncdim), OPTIONAL, INTENT(in) :: stride
    INTEGER, DIMENSION(ncdim), OPTIONAL, INTENT(in) :: imap

    INTEGER, DIMENSION(ncdim) :: ncstart
    INTEGER, DIMENSION(ncdim) :: nccount
    INTEGER, DIMENSION(ncdim) :: ncstride

    ncstatus = nf_inq_varid(file%id, var%name, var%id) 
    IF (ncstatus .NE. NF_NOERR) THEN
       CALL handle_ncerr('ncread (get_var_int1d)',ncstatus)
    END IF

    ! default for start is origin of array
    IF(PRESENT(start)) THEN
       ncstart=start
    ELSE
       ncstart=1
    END IF

    ! If count is not given, determine the size of value. A subarray
    ! of that size is read from the netCDF file
    IF(PRESENT(count)) THEN
       nccount=count
    ELSE
       nccount=SHAPE(value)
    END IF

    ! TODO: moeglicher Weise ein Problem mit default count, weil beim
    ! subsampling ja effektiv weniger Werte ausgelesen werden. Evtl.
    ! muesste dann nccount=count*stride sein oder so...
    IF(PRESENT(imap)) THEN
       IF(PRESENT(stride)) THEN
          ncstride=stride
       ELSE
          ncstride=1
       END IF

       ncstatus = nf_get_varm_int(file%id, var%id,&
            & ncstart, nccount, ncstride, imap, value)
       
    ELSE IF(PRESENT(stride)) THEN
       ncstatus=nf_get_vars_int(file%id, var%id,&
            & ncstart, nccount, stride, value)
    ELSE
       ncstatus = nf_get_vara_int(file%id, var%id,&
            & ncstart, nccount, value)
    END IF

    IF (ncstatus .NE. NF_NOERR) THEN
       CALL handle_ncerr('ncread (get_var_int1d)',ncstatus)
    END IF
  END SUBROUTINE get_var_int1d

  SUBROUTINE get_var_int2d(file, var, value, start, count, stride, imap)

    IMPLICIT NONE

    INTEGER, PARAMETER :: ncdim=2

    TYPE(ncfile), INTENT(in) :: file
    TYPE(ncvar), INTENT(inout) :: var
    INTEGER, DIMENSION(:,:), INTENT(out) :: value
    INTEGER, DIMENSION(ncdim), OPTIONAL, INTENT(in) :: start
    INTEGER, DIMENSION(ncdim), OPTIONAL, INTENT(in) :: count
    INTEGER, DIMENSION(ncdim), OPTIONAL, INTENT(in) :: stride
    INTEGER, DIMENSION(ncdim), OPTIONAL, INTENT(in) :: imap

    INTEGER, DIMENSION(ncdim) :: ncstart
    INTEGER, DIMENSION(ncdim) :: nccount
    INTEGER, DIMENSION(ncdim) :: ncstride

    ncstatus = nf_inq_varid(file%id, var%name, var%id) 
    IF (ncstatus .NE. NF_NOERR) THEN
       CALL handle_ncerr('ncread (get_var_int2d)', ncstatus)
    END IF

    ! default for start is origin of array
    IF(PRESENT(start)) THEN
       ncstart=start
    ELSE
       ncstart=1
    END IF

    ! If count is not given, determine the size of value. A subarray
    ! of that size is read from the netCDF file
    IF(PRESENT(count)) THEN
       nccount=count
    ELSE
       nccount=SHAPE(value)
    END IF

    ! TODO: moeglicher Weise ein Problem mit default count, weil beim
    ! subsampling ja effektiv weniger Werte ausgelesen werden. Evtl.
    ! muesste dann nccount=count*stride sein oder so...
    IF(PRESENT(imap)) THEN
       IF(PRESENT(stride)) THEN
          ncstride=stride
       ELSE
          ncstride=1
       END IF

       ncstatus = nf_get_varm_int(file%id, var%id,&
            & ncstart, nccount, ncstride, imap, value)
       
    ELSE IF(PRESENT(stride)) THEN
       ncstatus=nf_get_vars_int(file%id, var%id,&
            & ncstart, nccount, stride, value)
    ELSE
       ncstatus = nf_get_vara_int(file%id, var%id,&
            & ncstart, nccount, value)
    END IF

    IF (ncstatus .NE. NF_NOERR) THEN
       CALL handle_ncerr('ncread (get_var_int2d)', ncstatus)
    END IF
  END SUBROUTINE get_var_int2d

  SUBROUTINE get_var_int3d(file, var, value, start, count, stride, imap)

    IMPLICIT NONE

    INTEGER, PARAMETER :: ncdim=3

    TYPE(ncfile), INTENT(in) :: file
    TYPE(ncvar), INTENT(inout) :: var
    INTEGER, DIMENSION(:,:,:), INTENT(out) :: value
    INTEGER, DIMENSION(ncdim), OPTIONAL, INTENT(in) :: start
    INTEGER, DIMENSION(ncdim), OPTIONAL, INTENT(in) :: count
    INTEGER, DIMENSION(ncdim), OPTIONAL, INTENT(in) :: stride
    INTEGER, DIMENSION(ncdim), OPTIONAL, INTENT(in) :: imap

    INTEGER, DIMENSION(ncdim) :: ncstart
    INTEGER, DIMENSION(ncdim) :: nccount
    INTEGER, DIMENSION(ncdim) :: ncstride

    ncstatus = nf_inq_varid(file%id, var%name, var%id) 
    IF (ncstatus .NE. NF_NOERR) THEN
       CALL handle_ncerr('ncread (get_var_int3d)', ncstatus)
    END IF

    ! default for start is origin of array
    IF(PRESENT(start)) THEN
       ncstart=start
    ELSE
       ncstart=1
    END IF

    ! If count is not given, determine the size of value. A subarray
    ! of that size is read from the netCDF file
    IF(PRESENT(count)) THEN
       nccount=count
    ELSE
       nccount=SHAPE(value)
    END IF

    ! TODO: moeglicher Weise ein Problem mit default count, weil beim
    ! subsampling ja effektiv weniger Werte ausgelesen werden. Evtl.
    ! muesste dann nccount=count*stride sein oder so...
    IF(PRESENT(imap)) THEN
       IF(PRESENT(stride)) THEN
          ncstride=stride
       ELSE
          ncstride=1
       END IF

       ncstatus = nf_get_varm_int(file%id, var%id,&
            & ncstart, nccount, ncstride, imap, value)
       
    ELSE IF(PRESENT(stride)) THEN
       ncstatus=nf_get_vars_int(file%id, var%id,&
            & ncstart, nccount, stride, value)
    ELSE
       ncstatus = nf_get_vara_int(file%id, var%id,&
            & ncstart, nccount, value)
    END IF

    IF (ncstatus .NE. NF_NOERR) THEN
       CALL handle_ncerr('ncread (get_var_int3d)', ncstatus)
    END IF
  END SUBROUTINE get_var_int3d

  SUBROUTINE get_var_int4d(file, var, value, start, count, stride, imap)

    IMPLICIT NONE

    INTEGER, PARAMETER :: ncdim=4

    TYPE(ncfile), INTENT(in) :: file
    TYPE(ncvar), INTENT(inout) :: var
    INTEGER, DIMENSION(:,:,:,:), INTENT(out) :: value
    INTEGER, DIMENSION(ncdim), OPTIONAL, INTENT(in) :: start
    INTEGER, DIMENSION(ncdim), OPTIONAL, INTENT(in) :: count
    INTEGER, DIMENSION(ncdim), OPTIONAL, INTENT(in) :: stride
    INTEGER, DIMENSION(ncdim), OPTIONAL, INTENT(in) :: imap

    INTEGER, DIMENSION(ncdim) :: ncstart
    INTEGER, DIMENSION(ncdim) :: nccount
    INTEGER, DIMENSION(ncdim) :: ncstride

    ncstatus = nf_inq_varid(file%id, var%name, var%id) 
    IF (ncstatus .NE. NF_NOERR) THEN
       CALL handle_ncerr('ncread (get_var_int4d)', ncstatus)
    END IF

    ! default for start is origin of array
    IF(PRESENT(start)) THEN
       ncstart=start
    ELSE
       ncstart=1
    END IF

    ! If count is not given, determine the size of value. A subarray
    ! of that size is read from the netCDF file
    IF(PRESENT(count)) THEN
       nccount=count
    ELSE
       nccount=SHAPE(value)
    END IF

    ! TODO: moeglicher Weise ein Problem mit default count, weil beim
    ! subsampling ja effektiv weniger Werte ausgelesen werden. Evtl.
    ! muesste dann nccount=count*stride sein oder so...
    IF(PRESENT(imap)) THEN
       IF(PRESENT(stride)) THEN
          ncstride=stride
       ELSE
          ncstride=1
       END IF

       ncstatus = nf_get_varm_int(file%id, var%id,&
            & ncstart, nccount, ncstride, imap, value)
       
    ELSE IF(PRESENT(stride)) THEN
       ncstatus=nf_get_vars_int(file%id, var%id,&
            & ncstart, nccount, stride, value)
    ELSE
       ncstatus = nf_get_vara_int(file%id, var%id,&
            & ncstart, nccount, value)
    END IF

    IF (ncstatus .NE. NF_NOERR) THEN
       CALL handle_ncerr('ncread (get_var_int4d', ncstatus)
    END IF
  END SUBROUTINE get_var_int4d

  ! **************************************************
  ! REAL
  ! **************************************************

  SUBROUTINE get_var_real0d(file, var, value)

    IMPLICIT NONE

    TYPE(ncfile), INTENT(in) :: file
    TYPE(ncvar), INTENT(inout) :: var
    REAL, INTENT(out) :: value

    ncstatus = nf_inq_varid(file%id, var%name, var%id) 
    IF (ncstatus .NE. NF_NOERR) THEN
       CALL handle_ncerr('ncread (get_var_real0d)', ncstatus)
    END IF

    ncstatus = nf_get_var_real(file%id, var%id, value)
       
    IF (ncstatus .NE. NF_NOERR) THEN
       CALL handle_ncerr('ncread (get_var_real0d)', ncstatus)
    END IF
  END SUBROUTINE get_var_real0d

  SUBROUTINE get_var_real1d(file, var, value, start, count, stride, imap)

    IMPLICIT NONE

    INTEGER, PARAMETER :: ncdim=1

    TYPE(ncfile), INTENT(in) :: file
    TYPE(ncvar), INTENT(inout) :: var
    REAL, DIMENSION(:), INTENT(out) :: value
    INTEGER, DIMENSION(ncdim), OPTIONAL, INTENT(in) :: start
    INTEGER, DIMENSION(ncdim), OPTIONAL, INTENT(in) :: count
    INTEGER, DIMENSION(ncdim), OPTIONAL, INTENT(in) :: stride
    INTEGER, DIMENSION(ncdim), OPTIONAL, INTENT(in) :: imap

    INTEGER, DIMENSION(ncdim) :: ncstart
    INTEGER, DIMENSION(ncdim) :: nccount
    INTEGER, DIMENSION(ncdim) :: ncstride

    ncstatus = nf_inq_varid(file%id, var%name, var%id) 
    IF (ncstatus .NE. NF_NOERR) THEN
       CALL handle_ncerr('ncread (get_var_real1d)', ncstatus)
    END IF

    ! default for start is origin of array
    IF(PRESENT(start)) THEN
       ncstart=start
    ELSE
       ncstart=1
    END IF

    ! If count is not given, determine the size of value. A subarray
    ! of that size is read from the netCDF file
    IF(PRESENT(count)) THEN
       nccount=count
    ELSE
       nccount=SHAPE(value)
    END IF

    ! TODO: moeglicher Weise ein Problem mit default count, weil beim
    ! subsampling ja effektiv weniger Werte ausgelesen werden. Evtl.
    ! muesste dann nccount=count*stride sein oder so...
    IF(PRESENT(imap)) THEN
       IF(PRESENT(stride)) THEN
          ncstride=stride
       ELSE
          ncstride=1
       END IF

       ncstatus = nf_get_varm_real(file%id, var%id,&
            & ncstart, nccount, ncstride, imap, value)
       
    ELSE IF(PRESENT(stride)) THEN
       ncstatus=nf_get_vars_real(file%id, var%id,&
            & ncstart, nccount, stride, value)
    ELSE
       ncstatus = nf_get_vara_real(file%id, var%id,&
            & ncstart, nccount, value)
    END IF

    IF (ncstatus .NE. NF_NOERR) THEN
       CALL handle_ncerr('ncread (get_var_real1d)', ncstatus)
    END IF
  END SUBROUTINE get_var_real1d

  SUBROUTINE get_var_real2d(file, var, value, start, count, stride, imap)

    IMPLICIT NONE

    INTEGER, PARAMETER :: ncdim=2

    TYPE(ncfile), INTENT(in) :: file
    TYPE(ncvar), INTENT(inout) :: var
    REAL, DIMENSION(:,:), INTENT(out) :: value
    INTEGER, DIMENSION(ncdim), OPTIONAL, INTENT(in) :: start
    INTEGER, DIMENSION(ncdim), OPTIONAL, INTENT(in) :: count
    INTEGER, DIMENSION(ncdim), OPTIONAL, INTENT(in) :: stride
    INTEGER, DIMENSION(ncdim), OPTIONAL, INTENT(in) :: imap

    INTEGER, DIMENSION(ncdim) :: ncstart
    INTEGER, DIMENSION(ncdim) :: nccount
    INTEGER, DIMENSION(ncdim) :: ncstride

    ncstatus = nf_inq_varid(file%id, var%name, var%id) 
    IF (ncstatus .NE. NF_NOERR) THEN
       CALL handle_ncerr('ncread (get_var_real2d)', ncstatus)
    END IF

    ! default for start is origin of array
    IF(PRESENT(start)) THEN
       ncstart=start
    ELSE
       ncstart=1
    END IF

    ! If count is not given, determine the size of value. A subarray
    ! of that size is read from the netCDF file
    IF(PRESENT(count)) THEN
       nccount=count
    ELSE
       nccount=SHAPE(value)
    END IF

    ! TODO: moeglicher Weise ein Problem mit default count, weil beim
    ! subsampling ja effektiv weniger Werte ausgelesen werden. Evtl.
    ! muesste dann nccount=count*stride sein oder so...
    IF(PRESENT(imap)) THEN
       IF(PRESENT(stride)) THEN
          ncstride=stride
       ELSE
          ncstride=1
       END IF

       ncstatus = nf_get_varm_real(file%id, var%id,&
            & ncstart, nccount, ncstride, imap, value)
       
    ELSE IF(PRESENT(stride)) THEN
       ncstatus=nf_get_vars_real(file%id, var%id,&
            & ncstart, nccount, stride, value)
    ELSE
       ncstatus = nf_get_vara_real(file%id, var%id,&
            & ncstart, nccount, value)
    END IF

    IF (ncstatus .NE. NF_NOERR) THEN
       CALL handle_ncerr('ncread (get_var_real2d)', ncstatus)
    END IF
  END SUBROUTINE get_var_real2d

  SUBROUTINE get_var_real3d(file, var, value, start, count, stride, imap)

    IMPLICIT NONE

    INTEGER, PARAMETER :: ncdim=3

    TYPE(ncfile), INTENT(in) :: file
    TYPE(ncvar), INTENT(inout) :: var
    REAL, DIMENSION(:,:,:), INTENT(out) :: value
    INTEGER, DIMENSION(ncdim), OPTIONAL, INTENT(in) :: start
    INTEGER, DIMENSION(ncdim), OPTIONAL, INTENT(in) :: count
    INTEGER, DIMENSION(ncdim), OPTIONAL, INTENT(in) :: stride
    INTEGER, DIMENSION(ncdim), OPTIONAL, INTENT(in) :: imap

    INTEGER, DIMENSION(ncdim) :: ncstart
    INTEGER, DIMENSION(ncdim) :: nccount
    INTEGER, DIMENSION(ncdim) :: ncstride

    ncstatus = nf_inq_varid(file%id, var%name, var%id) 
    IF (ncstatus .NE. NF_NOERR) THEN
       CALL handle_ncerr('ncread (get_var_real3d)', ncstatus)
    END IF

    ! default for start is origin of array
    IF(PRESENT(start)) THEN
       ncstart=start
    ELSE
       ncstart=1
    END IF

    ! If count is not given, determine the size of value. A subarray
    ! of that size is read from the netCDF file
    IF(PRESENT(count)) THEN
       nccount=count
    ELSE
       nccount=SHAPE(value)
    END IF

    ! TODO: moeglicher Weise ein Problem mit default count, weil beim
    ! subsampling ja effektiv weniger Werte ausgelesen werden. Evtl.
    ! muesste dann nccount=count*stride sein oder so...
    IF(PRESENT(imap)) THEN
       IF(PRESENT(stride)) THEN
          ncstride=stride
       ELSE
          ncstride=1
       END IF

       ncstatus = nf_get_varm_real(file%id, var%id,&
            & ncstart, nccount, ncstride, imap, value)
       
    ELSE IF(PRESENT(stride)) THEN
       ncstatus=nf_get_vars_real(file%id, var%id,&
            & ncstart, nccount, stride, value)
    ELSE
       ncstatus = nf_get_vara_real(file%id, var%id,&
            & ncstart, nccount, value)
    END IF

    IF (ncstatus .NE. NF_NOERR) THEN
       CALL handle_ncerr('ncread (get_var_real3d)', ncstatus)
    END IF
  END SUBROUTINE get_var_real3d

  SUBROUTINE get_var_real4d(file, var, value, start, count, stride, imap)

    IMPLICIT NONE

    INTEGER, PARAMETER :: ncdim=4

    TYPE(ncfile), INTENT(in) :: file
    TYPE(ncvar), INTENT(inout) :: var
    REAL, DIMENSION(:,:,:,:), INTENT(out) :: value
    INTEGER, DIMENSION(ncdim), OPTIONAL, INTENT(in) :: start
    INTEGER, DIMENSION(ncdim), OPTIONAL, INTENT(in) :: count
    INTEGER, DIMENSION(ncdim), OPTIONAL, INTENT(in) :: stride
    INTEGER, DIMENSION(ncdim), OPTIONAL, INTENT(in) :: imap

    INTEGER, DIMENSION(ncdim) :: ncstart
    INTEGER, DIMENSION(ncdim) :: nccount
    INTEGER, DIMENSION(ncdim) :: ncstride

    ncstatus = nf_inq_varid(file%id, var%name, var%id) 
    IF (ncstatus .NE. NF_NOERR) THEN
       CALL handle_ncerr('ncread (get_var_real4d)', ncstatus)
    END IF

    ! default for start is origin of array
    IF(PRESENT(start)) THEN
       ncstart=start
    ELSE
       ncstart=1
    END IF

    ! If count is not given, determine the size of value. A subarray
    ! of that size is read from the netCDF file
    IF(PRESENT(count)) THEN
       nccount=count
    ELSE
       nccount=SHAPE(value)
    END IF

    ! TODO: moeglicher Weise ein Problem mit default count, weil beim
    ! subsampling ja effektiv weniger Werte ausgelesen werden. Evtl.
    ! muesste dann nccount=count*stride sein oder so...
    IF(PRESENT(imap)) THEN
       IF(PRESENT(stride)) THEN
          ncstride=stride
       ELSE
          ncstride=1
       END IF

       ncstatus = nf_get_varm_real(file%id, var%id,&
            & ncstart, nccount, ncstride, imap, value)
       
    ELSE IF(PRESENT(stride)) THEN
       ncstatus=nf_get_vars_real(file%id, var%id,&
            & ncstart, nccount, stride, value)
    ELSE
       ncstatus = nf_get_vara_real(file%id, var%id,&
            & ncstart, nccount, value)
    END IF

    IF (ncstatus .NE. NF_NOERR) THEN
       CALL handle_ncerr('ncread (get_var_real4d)', ncstatus)
    END IF
  END SUBROUTINE get_var_real4d

  ! **************************************************
  ! REAL
  ! **************************************************

  SUBROUTINE get_var_double0d(file, var, value)

    IMPLICIT NONE

    TYPE(ncfile), INTENT(in) :: file
    TYPE(ncvar), INTENT(inout) :: var
    REAL(kind=dblkind), INTENT(out) :: value

    ncstatus = nf_inq_varid(file%id, var%name, var%id) 
    IF (ncstatus .NE. NF_NOERR) THEN
       CALL handle_ncerr('ncread (get_double0d)', ncstatus)
    END IF

    ncstatus = nf_get_var_double(file%id, var%id, value)
       
    IF (ncstatus .NE. NF_NOERR) THEN
       CALL handle_ncerr('ncread (get_var_double0d)', ncstatus)
    END IF
  END SUBROUTINE get_var_double0d

  SUBROUTINE get_var_double1d(file, var, value, start, count, stride, imap)

    IMPLICIT NONE

    INTEGER, PARAMETER :: ncdim=1

    TYPE(ncfile), INTENT(in) :: file
    TYPE(ncvar), INTENT(inout) :: var
    REAL(kind=dblkind), DIMENSION(:), INTENT(out) :: value
    INTEGER, DIMENSION(ncdim), OPTIONAL, INTENT(in) :: start
    INTEGER, DIMENSION(ncdim), OPTIONAL, INTENT(in) :: count
    INTEGER, DIMENSION(ncdim), OPTIONAL, INTENT(in) :: stride
    INTEGER, DIMENSION(ncdim), OPTIONAL, INTENT(in) :: imap

    INTEGER, DIMENSION(ncdim) :: ncstart
    INTEGER, DIMENSION(ncdim) :: nccount
    INTEGER, DIMENSION(ncdim) :: ncstride

    ncstatus = nf_inq_varid(file%id, var%name, var%id) 
    IF (ncstatus .NE. NF_NOERR) THEN
       CALL handle_ncerr('ncread (get_var_double1d)', ncstatus)
    END IF

    ! default for start is origin of array
    IF(PRESENT(start)) THEN
       ncstart=start
    ELSE
       ncstart=1
    END IF

    ! If count is not given, determine the size of value. A subarray
    ! of that size is read from the netCDF file
    IF(PRESENT(count)) THEN
       nccount=count
    ELSE
       nccount=SHAPE(value)
    END IF

    ! TODO: moeglicher Weise ein Problem mit default count, weil beim
    ! subsampling ja effektiv weniger Werte ausgelesen werden. Evtl.
    ! muesste dann nccount=count*stride sein oder so...
    IF(PRESENT(imap)) THEN
       IF(PRESENT(stride)) THEN
          ncstride=stride
       ELSE
          ncstride=1
       END IF

       ncstatus = nf_get_varm_double(file%id, var%id,&
            & ncstart, nccount, ncstride, imap, value)
       
    ELSE IF(PRESENT(stride)) THEN
       ncstatus=nf_get_vars_double(file%id, var%id,&
            & ncstart, nccount, stride, value)
    ELSE
       ncstatus = nf_get_vara_double(file%id, var%id,&
            & ncstart, nccount, value)
    END IF

    IF (ncstatus .NE. NF_NOERR) THEN
       CALL handle_ncerr('ncread (get_var_double1d)', ncstatus)
    END IF
  END SUBROUTINE get_var_double1d

  SUBROUTINE get_var_double2d(file, var, value, start, count, stride, imap)

    IMPLICIT NONE

    INTEGER, PARAMETER :: ncdim=2

    TYPE(ncfile), INTENT(in) :: file
    TYPE(ncvar), INTENT(inout) :: var
    REAL(kind=dblkind), DIMENSION(:,:), INTENT(out) :: value
    INTEGER, DIMENSION(ncdim), OPTIONAL, INTENT(in) :: start
    INTEGER, DIMENSION(ncdim), OPTIONAL, INTENT(in) :: count
    INTEGER, DIMENSION(ncdim), OPTIONAL, INTENT(in) :: stride
    INTEGER, DIMENSION(ncdim), OPTIONAL, INTENT(in) :: imap

    INTEGER, DIMENSION(ncdim) :: ncstart
    INTEGER, DIMENSION(ncdim) :: nccount
    INTEGER, DIMENSION(ncdim) :: ncstride

    ncstatus = nf_inq_varid(file%id, var%name, var%id) 
    IF (ncstatus .NE. NF_NOERR) THEN
       CALL handle_ncerr('ncread (get_var_double2d)', ncstatus)
    END IF

    ! default for start is origin of array
    IF(PRESENT(start)) THEN
       ncstart=start
    ELSE
       ncstart=1
    END IF

    ! If count is not given, determine the size of value. A subarray
    ! of that size is read from the netCDF file
    IF(PRESENT(count)) THEN
       nccount=count
    ELSE
       nccount=SHAPE(value)
    END IF

    ! TODO: moeglicher Weise ein Problem mit default count, weil beim
    ! subsampling ja effektiv weniger Werte ausgelesen werden. Evtl.
    ! muesste dann nccount=count*stride sein oder so...
    IF(PRESENT(imap)) THEN
       IF(PRESENT(stride)) THEN
          ncstride=stride
       ELSE
          ncstride=1
       END IF

       ncstatus = nf_get_varm_double(file%id, var%id,&
            & ncstart, nccount, ncstride, imap, value)
       
    ELSE IF(PRESENT(stride)) THEN
       ncstatus=nf_get_vars_double(file%id, var%id,&
            & ncstart, nccount, stride, value)
    ELSE
       ncstatus = nf_get_vara_double(file%id, var%id,&
            & ncstart, nccount, value)
    END IF

    IF (ncstatus .NE. NF_NOERR) THEN
       CALL handle_ncerr('ncread (get_var_double2d)', ncstatus)
    END IF
  END SUBROUTINE get_var_double2d

  SUBROUTINE get_var_double3d(file, var, value, start, count, stride, imap)

    IMPLICIT NONE

    INTEGER, PARAMETER :: ncdim=3

    TYPE(ncfile), INTENT(in) :: file
    TYPE(ncvar), INTENT(inout) :: var
    REAL(kind=dblkind), DIMENSION(:,:,:), INTENT(out) :: value
    INTEGER, DIMENSION(ncdim), OPTIONAL, INTENT(in) :: start
    INTEGER, DIMENSION(ncdim), OPTIONAL, INTENT(in) :: count
    INTEGER, DIMENSION(ncdim), OPTIONAL, INTENT(in) :: stride
    INTEGER, DIMENSION(ncdim), OPTIONAL, INTENT(in) :: imap

    INTEGER, DIMENSION(ncdim) :: ncstart
    INTEGER, DIMENSION(ncdim) :: nccount
    INTEGER, DIMENSION(ncdim) :: ncstride

    ncstatus = nf_inq_varid(file%id, var%name, var%id) 
    IF (ncstatus .NE. NF_NOERR) THEN
       CALL handle_ncerr('ncread (get_var_double3d)', ncstatus)
    END IF

    ! default for start is origin of array
    IF(PRESENT(start)) THEN
       ncstart=start
    ELSE
       ncstart=1
    END IF

    ! If count is not given, determine the size of value. A subarray
    ! of that size is read from the netCDF file
    IF(PRESENT(count)) THEN
       nccount=count
    ELSE
       nccount=SHAPE(value)
    END IF

    ! TODO: moeglicher Weise ein Problem mit default count, weil beim
    ! subsampling ja effektiv weniger Werte ausgelesen werden. Evtl.
    ! muesste dann nccount=count*stride sein oder so...
    IF(PRESENT(imap)) THEN
       IF(PRESENT(stride)) THEN
          ncstride=stride
       ELSE
          ncstride=1
       END IF

       ncstatus = nf_get_varm_double(file%id, var%id,&
            & ncstart, nccount, ncstride, imap, value)
       
    ELSE IF(PRESENT(stride)) THEN
       ncstatus=nf_get_vars_double(file%id, var%id,&
            & ncstart, nccount, stride, value)
    ELSE
       ncstatus = nf_get_vara_double(file%id, var%id,&
            & ncstart, nccount, value)
    END IF

    IF (ncstatus .NE. NF_NOERR) THEN
       CALL handle_ncerr('ncread (get_var_double3d)', ncstatus)
    END IF
  END SUBROUTINE get_var_double3d

  SUBROUTINE get_var_double4d(file, var, value, start, count, stride, imap)

    IMPLICIT NONE

    INTEGER, PARAMETER :: ncdim=4

    TYPE(ncfile), INTENT(in) :: file
    TYPE(ncvar), INTENT(inout) :: var
    REAL(kind=dblkind), DIMENSION(:,:,:,:), INTENT(out) :: value
    INTEGER, DIMENSION(ncdim), OPTIONAL, INTENT(in) :: start
    INTEGER, DIMENSION(ncdim), OPTIONAL, INTENT(in) :: count
    INTEGER, DIMENSION(ncdim), OPTIONAL, INTENT(in) :: stride
    INTEGER, DIMENSION(ncdim), OPTIONAL, INTENT(in) :: imap

    INTEGER, DIMENSION(ncdim) :: ncstart
    INTEGER, DIMENSION(ncdim) :: nccount
    INTEGER, DIMENSION(ncdim) :: ncstride

    ncstatus = nf_inq_varid(file%id, var%name, var%id) 
    IF (ncstatus .NE. NF_NOERR) THEN
       CALL handle_ncerr('ncread (get_var_double4d)', ncstatus)
    END IF

    ! default for start is origin of array
    IF(PRESENT(start)) THEN
       ncstart=start
    ELSE
       ncstart=1
    END IF

    ! If count is not given, determine the size of value. A subarray
    ! of that size is read from the netCDF file
    IF(PRESENT(count)) THEN
       nccount=count
    ELSE
       nccount=SHAPE(value)
    END IF

    ! TODO: moeglicher Weise ein Problem mit default count, weil beim
    ! subsampling ja effektiv weniger Werte ausgelesen werden. Evtl.
    ! muesste dann nccount=count*stride sein oder so...
    IF(PRESENT(imap)) THEN
       IF(PRESENT(stride)) THEN
          ncstride=stride
       ELSE
          ncstride=1
       END IF

       ncstatus = nf_get_varm_double(file%id, var%id,&
            & ncstart, nccount, ncstride, imap, value)
       
    ELSE IF(PRESENT(stride)) THEN
       ncstatus=nf_get_vars_double(file%id, var%id,&
            & ncstart, nccount, stride, value)
    ELSE
       ncstatus = nf_get_vara_double(file%id, var%id,&
            & ncstart, nccount, value)
    END IF

    IF (ncstatus .NE. NF_NOERR) THEN
       CALL handle_ncerr('ncread (get_var_double4d)', ncstatus)
    END IF
  END SUBROUTINE get_var_double4d



  SUBROUTINE put_var_real0d(file, var, value)
    
    IMPLICIT NONE

    TYPE(ncfile), INTENT(in) :: file
    TYPE(ncvar), INTENT(inout) :: var
    REAL, INTENT(in) :: value

    ncstatus = nf_inq_varid(file%id,var%name,var%id)
    IF (ncstatus .NE. NF_NOERR) THEN
       CALL handle_ncerr('ncprint (put_var_real0d)', ncstatus)
    END IF

    ncstatus=nf_put_var_real(file%id,var%id,value)
    IF (ncstatus .NE. NF_NOERR) THEN
       CALL handle_ncerr('ncprint (put_var_real0d)', ncstatus)
    END IF
  END SUBROUTINE put_var_real0d
  
  SUBROUTINE put_var_real1d(file, var, value, start, count, stride, imap)

    IMPLICIT NONE

    INTEGER, PARAMETER :: ncdim=1

    TYPE(ncfile), INTENT(in) :: file
    TYPE(ncvar), INTENT(inout) :: var
    REAL, DIMENSION(:), INTENT(in) :: value
    INTEGER, DIMENSION(ncdim), OPTIONAL, INTENT(in) :: start
    INTEGER, DIMENSION(ncdim), OPTIONAL, INTENT(in) :: count
    INTEGER, DIMENSION(ncdim), OPTIONAL, INTENT(in) :: stride
    INTEGER, DIMENSION(ncdim), OPTIONAL, INTENT(in) :: imap

    INTEGER, DIMENSION(ncdim) :: ncstart
    INTEGER, DIMENSION(ncdim) :: nccount
    INTEGER, DIMENSION(ncdim) :: ncstride

    ncstatus = nf_inq_varid(file%id, var%name, var%id) 
    IF (ncstatus .NE. NF_NOERR) THEN
       CALL handle_ncerr('ncprint (put_var_real1d)', ncstatus)
    END IF

    ! default for start is origin of array
    IF(PRESENT(start)) THEN
       ncstart=start
    ELSE
       ncstart=1
    END IF

    ! If count is not given, determine the size of value. A subarray
    ! of that size is read from the netCDF file
    IF(PRESENT(count)) THEN
       nccount=count
    ELSE
       nccount=SHAPE(value)
    END IF

    ! TODO: moeglicher Weise ein Problem mit default count, weil beim
    ! subsampling ja effektiv weniger Werte ausgelesen werden. Evtl.
    ! muesste dann nccount=count*stride sein oder so...
    IF(PRESENT(imap)) THEN
       IF(PRESENT(stride)) THEN
          ncstride=stride
       ELSE
          ncstride=1
       END IF

       ncstatus = nf_put_varm_real(file%id, var%id,&
            & ncstart, nccount, ncstride, imap, value)
       
    ELSE IF(PRESENT(stride)) THEN
       ncstatus=nf_put_vars_real(file%id, var%id,&
            & ncstart, nccount, stride, value)
    ELSE
       ncstatus = nf_put_vara_real(file%id, var%id,&
            & ncstart, nccount, value)
    END IF

    IF (ncstatus .NE. NF_NOERR) THEN
       CALL handle_ncerr('ncprint (put_var_real1d)', ncstatus)
    END IF
  END SUBROUTINE put_var_real1d
  
  SUBROUTINE put_var_real2d(file, var, value, start, count, stride, imap)

    IMPLICIT NONE

    INTEGER, PARAMETER :: ncdim=2

    TYPE(ncfile), INTENT(in) :: file
    TYPE(ncvar), INTENT(inout) :: var
    REAL, DIMENSION(:,:), INTENT(in) :: value
    INTEGER, DIMENSION(ncdim), OPTIONAL, INTENT(in) :: start
    INTEGER, DIMENSION(ncdim), OPTIONAL, INTENT(in) :: count
    INTEGER, DIMENSION(ncdim), OPTIONAL, INTENT(in) :: stride
    INTEGER, DIMENSION(ncdim), OPTIONAL, INTENT(in) :: imap

    INTEGER, DIMENSION(ncdim) :: ncstart
    INTEGER, DIMENSION(ncdim) :: nccount
    INTEGER, DIMENSION(ncdim) :: ncstride

    ncstatus = nf_inq_varid(file%id, var%name, var%id) 
    IF (ncstatus .NE. NF_NOERR) THEN
       CALL handle_ncerr('ncprint (put_var_real2d)', ncstatus)
    END IF

    ! default for start is origin of array
    IF(PRESENT(start)) THEN
       ncstart=start
    ELSE
       ncstart=1
    END IF

    ! If count is not given, determine the size of value. A subarray
    ! of that size is read from the netCDF file
    IF(PRESENT(count)) THEN
       nccount=count
    ELSE
       nccount=SHAPE(value)
    END IF

    ! TODO: moeglicher Weise ein Problem mit default count, weil beim
    ! subsampling ja effektiv weniger Werte ausgelesen werden. Evtl.
    ! muesste dann nccount=count*stride sein oder so...
    IF(PRESENT(imap)) THEN
       IF(PRESENT(stride)) THEN
          ncstride=stride
       ELSE
          ncstride=1
       END IF

       ncstatus = nf_put_varm_real(file%id, var%id,&
            & ncstart, nccount, ncstride, imap, value)
       
    ELSE IF(PRESENT(stride)) THEN
       ncstatus=nf_put_vars_real(file%id, var%id,&
            & ncstart, nccount, stride, value)
    ELSE
       ncstatus = nf_put_vara_real(file%id, var%id,&
            & ncstart, nccount, value)
    END IF

    IF (ncstatus .NE. NF_NOERR) THEN
       CALL handle_ncerr('ncprint (put_var_real2d)', ncstatus)
    END IF
  END SUBROUTINE put_var_real2d

  SUBROUTINE put_var_real3d(file, var, value, start, count, stride, imap)

    IMPLICIT NONE

    INTEGER, PARAMETER :: ncdim=3

    TYPE(ncfile), INTENT(in) :: file
    TYPE(ncvar), INTENT(inout) :: var
    REAL, DIMENSION(:,:,:), INTENT(in) :: value
    INTEGER, DIMENSION(ncdim), OPTIONAL, INTENT(in) :: start
    INTEGER, DIMENSION(ncdim), OPTIONAL, INTENT(in) :: count
    INTEGER, DIMENSION(ncdim), OPTIONAL, INTENT(in) :: stride
    INTEGER, DIMENSION(ncdim), OPTIONAL, INTENT(in) :: imap

    INTEGER, DIMENSION(ncdim) :: ncstart
    INTEGER, DIMENSION(ncdim) :: nccount
    INTEGER, DIMENSION(ncdim) :: ncstride

    ncstatus = nf_inq_varid(file%id, var%name, var%id) 
    IF (ncstatus .NE. NF_NOERR) THEN
       CALL handle_ncerr('ncprint (put_var_real3d)', ncstatus)
    END IF

    ! default for start is origin of array
    IF(PRESENT(start)) THEN
       ncstart=start
    ELSE
       ncstart=1
    END IF

    ! If count is not given, determine the size of value. A subarray
    ! of that size is read from the netCDF file
    IF(PRESENT(count)) THEN
       nccount=count
    ELSE
       nccount=SHAPE(value)
    END IF

    ! TODO: moeglicher Weise ein Problem mit default count, weil beim
    ! subsampling ja effektiv weniger Werte ausgelesen werden. Evtl.
    ! muesste dann nccount=count*stride sein oder so...
    IF(PRESENT(imap)) THEN
       IF(PRESENT(stride)) THEN
          ncstride=stride
       ELSE
          ncstride=1
       END IF

       ncstatus = nf_put_varm_real(file%id, var%id,&
            & ncstart, nccount, ncstride, imap, value)
       
    ELSE IF(PRESENT(stride)) THEN
       ncstatus=nf_put_vars_real(file%id, var%id,&
            & ncstart, nccount, stride, value)
    ELSE
       ncstatus = nf_put_vara_real(file%id, var%id,&
            & ncstart, nccount, value)
    END IF

    IF (ncstatus .NE. NF_NOERR) THEN
       CALL handle_ncerr('ncprint (put_var_real3d)', ncstatus)
    END IF
  END SUBROUTINE put_var_real3d

  SUBROUTINE put_var_real4d(file, var, value, start, count, stride, imap)

    IMPLICIT NONE

    INTEGER, PARAMETER :: ncdim=4

    TYPE(ncfile), INTENT(in) :: file
    TYPE(ncvar), INTENT(inout) :: var
    REAL, DIMENSION(:,:,:,:), INTENT(in) :: value
    INTEGER, DIMENSION(ncdim), OPTIONAL, INTENT(in) :: start
    INTEGER, DIMENSION(ncdim), OPTIONAL, INTENT(in) :: count
    INTEGER, DIMENSION(ncdim), OPTIONAL, INTENT(in) :: stride
    INTEGER, DIMENSION(ncdim), OPTIONAL, INTENT(in) :: imap

    INTEGER, DIMENSION(ncdim) :: ncstart
    INTEGER, DIMENSION(ncdim) :: nccount
    INTEGER, DIMENSION(ncdim) :: ncstride

    ncstatus = nf_inq_varid(file%id, var%name, var%id) 
    IF (ncstatus .NE. NF_NOERR) THEN
       CALL handle_ncerr('ncprint (put_var_real4d)', ncstatus)
    END IF

    ! default for start is origin of array
    IF(PRESENT(start)) THEN
       ncstart=start
    ELSE
       ncstart=1
    END IF

    ! If count is not given, determine the size of value. A subarray
    ! of that size is read from the netCDF file
    IF(PRESENT(count)) THEN
       nccount=count
    ELSE
       nccount=SHAPE(value)
    END IF

    ! TODO: moeglicher Weise ein Problem mit default count, weil beim
    ! subsampling ja effektiv weniger Werte ausgelesen werden. Evtl.
    ! muesste dann nccount=count*stride sein oder so...
    IF(PRESENT(imap)) THEN
       IF(PRESENT(stride)) THEN
          ncstride=stride
       ELSE
          ncstride=1
       END IF

       ncstatus = nf_put_varm_real(file%id, var%id,&
            & ncstart, nccount, ncstride, imap, value)
       
    ELSE IF(PRESENT(stride)) THEN
       ncstatus=nf_put_vars_real(file%id, var%id,&
            & ncstart, nccount, stride, value)
    ELSE
       ncstatus = nf_put_vara_real(file%id, var%id,&
            & ncstart, nccount, value)
    END IF

    IF (ncstatus .NE. NF_NOERR) THEN
       CALL handle_ncerr('ncprint (put_var_real4d)', ncstatus)
    END IF
  END SUBROUTINE put_var_real4d

  SUBROUTINE put_var_double0d(file, var, value)
    
    IMPLICIT NONE

    TYPE(ncfile), INTENT(in) :: file
    TYPE(ncvar), INTENT(inout) :: var
    REAL(kind=dblkind), INTENT(in) :: value

    ncstatus = nf_inq_varid(file%id,var%name,var%id)
    IF (ncstatus .NE. NF_NOERR) THEN
       CALL handle_ncerr('ncprint (put_var_double0d)', ncstatus)
    END IF

    ncstatus=nf_put_var_double(file%id,var%id,value)
    IF (ncstatus .NE. NF_NOERR) THEN
       CALL handle_ncerr('ncprint (put_var_double0d)', ncstatus)
    END IF
  END SUBROUTINE put_var_double0d
  
  SUBROUTINE put_var_double1d(file, var, value, start, count, stride, imap)

    IMPLICIT NONE

    INTEGER, PARAMETER :: ncdim=1

    TYPE(ncfile), INTENT(in) :: file
    TYPE(ncvar), INTENT(inout) :: var
    REAL(kind=dblkind), DIMENSION(:), INTENT(in) :: value
    INTEGER, DIMENSION(ncdim), OPTIONAL, INTENT(in) :: start
    INTEGER, DIMENSION(ncdim), OPTIONAL, INTENT(in) :: count
    INTEGER, DIMENSION(ncdim), OPTIONAL, INTENT(in) :: stride
    INTEGER, DIMENSION(ncdim), OPTIONAL, INTENT(in) :: imap

    INTEGER, DIMENSION(ncdim) :: ncstart
    INTEGER, DIMENSION(ncdim) :: nccount
    INTEGER, DIMENSION(ncdim) :: ncstride

    ncstatus = nf_inq_varid(file%id, var%name, var%id) 
    IF (ncstatus .NE. NF_NOERR) THEN
       CALL handle_ncerr('ncprint (put_var_double1d)', ncstatus)
    END IF

    ! default for start is origin of array
    IF(PRESENT(start)) THEN
       ncstart=start
    ELSE
       ncstart=1
    END IF

    ! If count is not given, determine the size of value. A subarray
    ! of that size is read from the netCDF file
    IF(PRESENT(count)) THEN
       nccount=count
    ELSE
       nccount=SHAPE(value)
    END IF

    ! TODO: moeglicher Weise ein Problem mit default count, weil beim
    ! subsampling ja effektiv weniger Werte ausgelesen werden. Evtl.
    ! muesste dann nccount=count*stride sein oder so...
    IF(PRESENT(imap)) THEN
       IF(PRESENT(stride)) THEN
          ncstride=stride
       ELSE
          ncstride=1
       END IF

       ncstatus = nf_put_varm_double(file%id, var%id,&
            & ncstart, nccount, ncstride, imap, value)
       
    ELSE IF(PRESENT(stride)) THEN
       ncstatus=nf_put_vars_double(file%id, var%id,&
            & ncstart, nccount, stride, value)
    ELSE
       ncstatus = nf_put_vara_double(file%id, var%id,&
            & ncstart, nccount, value)
    END IF

    IF (ncstatus .NE. NF_NOERR) THEN
       CALL handle_ncerr('ncprint (put_var_double1d)', ncstatus)
    END IF
  END SUBROUTINE put_var_double1d
  
  SUBROUTINE put_var_double2d(file, var, value, start, count, stride, imap)

    IMPLICIT NONE

    INTEGER, PARAMETER :: ncdim=2

    TYPE(ncfile), INTENT(in) :: file
    TYPE(ncvar), INTENT(inout) :: var
    REAL(kind=dblkind), DIMENSION(:,:), INTENT(in) :: value
    INTEGER, DIMENSION(ncdim), OPTIONAL, INTENT(in) :: start
    INTEGER, DIMENSION(ncdim), OPTIONAL, INTENT(in) :: count
    INTEGER, DIMENSION(ncdim), OPTIONAL, INTENT(in) :: stride
    INTEGER, DIMENSION(ncdim), OPTIONAL, INTENT(in) :: imap

    INTEGER, DIMENSION(ncdim) :: ncstart
    INTEGER, DIMENSION(ncdim) :: nccount
    INTEGER, DIMENSION(ncdim) :: ncstride

    ncstatus = nf_inq_varid(file%id, var%name, var%id) 
    IF (ncstatus .NE. NF_NOERR) THEN
       CALL handle_ncerr('ncprint (put_var_double2d)', ncstatus)
    END IF

    ! default for start is origin of array
    IF(PRESENT(start)) THEN
       ncstart=start
    ELSE
       ncstart=1
    END IF

    ! If count is not given, determine the size of value. A subarray
    ! of that size is read from the netCDF file
    IF(PRESENT(count)) THEN
       nccount=count
    ELSE
       nccount=SHAPE(value)
    END IF

    ! TODO: moeglicher Weise ein Problem mit default count, weil beim
    ! subsampling ja effektiv weniger Werte ausgelesen werden. Evtl.
    ! muesste dann nccount=count*stride sein oder so...
    IF(PRESENT(imap)) THEN
       IF(PRESENT(stride)) THEN
          ncstride=stride
       ELSE
          ncstride=1
       END IF

       ncstatus = nf_put_varm_double(file%id, var%id,&
            & ncstart, nccount, ncstride, imap, value)
       
    ELSE IF(PRESENT(stride)) THEN
       ncstatus=nf_put_vars_double(file%id, var%id,&
            & ncstart, nccount, stride, value)
    ELSE
       ncstatus = nf_put_vara_double(file%id, var%id,&
            & ncstart, nccount, value)
    END IF

    IF (ncstatus .NE. NF_NOERR) THEN
       CALL handle_ncerr('ncprint (put_var_double2d)', ncstatus)
    END IF
  END SUBROUTINE put_var_double2d

  SUBROUTINE put_var_double3d(file, var, value, start, count, stride, imap)

    IMPLICIT NONE

    INTEGER, PARAMETER :: ncdim=3

    TYPE(ncfile), INTENT(in) :: file
    TYPE(ncvar), INTENT(inout) :: var
    REAL(kind=dblkind), DIMENSION(:,:,:), INTENT(in) :: value
    INTEGER, DIMENSION(ncdim), OPTIONAL, INTENT(in) :: start
    INTEGER, DIMENSION(ncdim), OPTIONAL, INTENT(in) :: count
    INTEGER, DIMENSION(ncdim), OPTIONAL, INTENT(in) :: stride
    INTEGER, DIMENSION(ncdim), OPTIONAL, INTENT(in) :: imap

    INTEGER, DIMENSION(ncdim) :: ncstart
    INTEGER, DIMENSION(ncdim) :: nccount
    INTEGER, DIMENSION(ncdim) :: ncstride

    ncstatus = nf_inq_varid(file%id, var%name, var%id) 
    IF (ncstatus .NE. NF_NOERR) THEN
       CALL handle_ncerr('ncprint (put_var_double3d)', ncstatus)
    END IF

    ! default for start is origin of array
    IF(PRESENT(start)) THEN
       ncstart=start
    ELSE
       ncstart=1
    END IF

    ! If count is not given, determine the size of value. A subarray
    ! of that size is read from the netCDF file
    IF(PRESENT(count)) THEN
       nccount=count
    ELSE
       nccount=SHAPE(value)
    END IF

    ! TODO: moeglicher Weise ein Problem mit default count, weil beim
    ! subsampling ja effektiv weniger Werte ausgelesen werden. Evtl.
    ! muesste dann nccount=count*stride sein oder so...
    IF(PRESENT(imap)) THEN
       IF(PRESENT(stride)) THEN
          ncstride=stride
       ELSE
          ncstride=1
       END IF

       ncstatus = nf_put_varm_double(file%id, var%id,&
            & ncstart, nccount, ncstride, imap, value)
       
    ELSE IF(PRESENT(stride)) THEN
       ncstatus=nf_put_vars_double(file%id, var%id,&
            & ncstart, nccount, stride, value)
    ELSE
       ncstatus = nf_put_vara_double(file%id, var%id,&
            & ncstart, nccount, value)
    END IF

    IF (ncstatus .NE. NF_NOERR) THEN
       CALL handle_ncerr('put_var_double3d', ncstatus)
    END IF
  END SUBROUTINE put_var_double3d

  SUBROUTINE put_var_double4d(file, var, value, start, count, stride, imap)

    IMPLICIT NONE

    INTEGER, PARAMETER :: ncdim=4

    TYPE(ncfile), INTENT(in) :: file
    TYPE(ncvar), INTENT(inout) :: var
    REAL(kind=dblkind), DIMENSION(:,:,:,:), INTENT(in) :: value
    INTEGER, DIMENSION(ncdim), OPTIONAL, INTENT(in) :: start
    INTEGER, DIMENSION(ncdim), OPTIONAL, INTENT(in) :: count
    INTEGER, DIMENSION(ncdim), OPTIONAL, INTENT(in) :: stride
    INTEGER, DIMENSION(ncdim), OPTIONAL, INTENT(in) :: imap

    INTEGER, DIMENSION(ncdim) :: ncstart
    INTEGER, DIMENSION(ncdim) :: nccount
    INTEGER, DIMENSION(ncdim) :: ncstride

    ncstatus = nf_inq_varid(file%id, var%name, var%id) 
    IF (ncstatus .NE. NF_NOERR) THEN
       CALL handle_ncerr('ncprint (put_var_double4d)', ncstatus)
    END IF

    ! default for start is origin of array
    IF(PRESENT(start)) THEN
       ncstart=start
    ELSE
       ncstart=1
    END IF

    ! If count is not given, determine the size of value. A subarray
    ! of that size is read from the netCDF file
    IF(PRESENT(count)) THEN
       nccount=count
    ELSE
       nccount=SHAPE(value)
    END IF

    ! TODO: moeglicher Weise ein Problem mit default count, weil beim
    ! subsampling ja effektiv weniger Werte ausgelesen werden. Evtl.
    ! muesste dann nccount=count*stride sein oder so...
    IF(PRESENT(imap)) THEN
       IF(PRESENT(stride)) THEN
          ncstride=stride
       ELSE
          ncstride=1
       END IF

       ncstatus = nf_put_varm_double(file%id, var%id,&
            & ncstart, nccount, ncstride, imap, value)
       
    ELSE IF(PRESENT(stride)) THEN
       ncstatus=nf_put_vars_double(file%id, var%id,&
            & ncstart, nccount, stride, value)
    ELSE
       ncstatus = nf_put_vara_double(file%id, var%id,&
            & ncstart, nccount, value)
    END IF

    IF (ncstatus .NE. NF_NOERR) THEN
       CALL handle_ncerr('ncprint (put_var_double4d)', ncstatus)
    END IF
  END SUBROUTINE put_var_double4d

  SUBROUTINE ncDimDef(file,dim)

    IMPLICIT NONE

    TYPE(ncfile), INTENT(in) :: file
    TYPE(ncdim), INTENT(inout) :: dim

    ncstatus=nf_def_dim(file%id,dim%name,dim%len,dim%id)
    
  END SUBROUTINE ncDimDef

  SUBROUTINE ncVarDef(file,var)

    IMPLICIT NONE

    TYPE(ncFile), INTENT(in) :: file
    TYPE(ncVar), INTENT(inout) :: var

    TYPE(ncDim), DIMENSION(var%nDims) :: dims
    INTEGER :: jd

    DO jd=1,var%nDims
       dims(jd)=var %dims(jd) %p
    END DO

    IF(var%nDims .EQ. 0) THEN
       ncstatus=nf_def_var(file%id,var%name,var%xType,var%nDims,(/ 0 /), &
            &var%id)
    ELSE
       ncstatus=nf_def_var(file%id,var%name,var%xType,var%nDims,&
            &dims%id,var%id)
    END IF
    IF(ncstatus .NE. NF_NOERR) THEN
       CALL handle_ncerr('ncvardef', ncstatus)
    END IF
  END SUBROUTINE ncVarDef

  SUBROUTINE ncGetAtt_real(file,var,AttName,AttVal)

    IMPLICIT NONE

    TYPE(ncfile), INTENT(in) :: file
    TYPE(ncvar), INTENT(in) :: var
    CHARACTER (len=*), INTENT(in) :: AttName
    REAL, INTENT(out) :: AttVal

    INTEGER :: attlen

    attlen=1


    ncstatus=nf_get_att_real(file%id,var%id,AttName,AttVal)
    IF(ncstatus .NE. NF_NOERR) THEN
       CALL handle_ncerr('ncGetAtt (ncGetAtt_real)', ncstatus)
    END IF
  END SUBROUTINE ncGetAtt_real

  SUBROUTINE ncGetAtt_double(file,var,AttName,AttVal)

    IMPLICIT NONE

    TYPE(ncfile), INTENT(in) :: file
    TYPE(ncvar), INTENT(in) :: var
    CHARACTER (len=*), INTENT(in) :: AttName
    REAL(kind=dblkind), INTENT(out) :: AttVal

    INTEGER :: attlen

    attlen=1


    ncstatus=nf_get_att_double(file%id,var%id,AttName, AttVal)
    IF(ncstatus .NE. NF_NOERR) THEN

       CALL handle_ncerr('ncGetAtt (ncGetAtt_double)', ncstatus)
    END IF
  END SUBROUTINE ncGetAtt_double


  SUBROUTINE ncPutAtt_text(file,var,AttName,AttVal)

    IMPLICIT NONE

    TYPE(ncfile), INTENT(in) :: file
    TYPE(ncvar), INTENT(in) :: var
    CHARACTER (len=*), INTENT(in) :: AttName
    CHARACTER (len=*), INTENT(in) :: AttVal

    INTEGER :: attlen

    attlen=LEN_TRIM(AttVal)

    ncstatus=nf_put_att_text(file%id,var%id,AttName,attlen,AttVal)
    IF(ncstatus .NE. NF_NOERR) THEN
       CALL handle_ncerr('ncPutAtt (ncPutAtt_text)', ncstatus)
    END IF
  END SUBROUTINE ncPutAtt_text

  SUBROUTINE ncPutAtt_real(file,var,AttName,AttVal)

    IMPLICIT NONE

    TYPE(ncfile), INTENT(in) :: file
    TYPE(ncvar), INTENT(in) :: var
    CHARACTER (len=*), INTENT(in) :: AttName
    REAL, INTENT(in) :: AttVal

    INTEGER :: attlen

    attlen=1

    ncstatus=nf_put_att_real(file%id,var%id,AttName,NF_REAL,attlen,AttVal)
    IF(ncstatus .NE. NF_NOERR) THEN
       CALL handle_ncerr('ncPutAtt (ncPutAtt_real)', ncstatus)
    END IF
  END SUBROUTINE ncPutAtt_real

  SUBROUTINE ncPutAtt_double(file,var,AttName,AttVal)

    IMPLICIT NONE

    TYPE(ncfile), INTENT(in) :: file
    TYPE(ncvar), INTENT(in) :: var
    CHARACTER (len=*), INTENT(in) :: AttName
    REAL(kind=dblkind), INTENT(in) :: AttVal

    INTEGER :: attlen

    attlen=1

    ncstatus=nf_put_att_double(file%id,var%id,AttName,NF_REAL,attlen,AttVal)
    IF(ncstatus .NE. NF_NOERR) THEN
       CALL handle_ncerr('ncPutAtt (ncPutAtt_double)', ncstatus)
    END IF
  END SUBROUTINE ncPutAtt_double

  SUBROUTINE ncPutAtt_dblarr(file,var,AttName,AttVal)

    IMPLICIT NONE

    TYPE(ncfile), INTENT(in) :: file
    TYPE(ncvar), INTENT(in) :: var
    CHARACTER (len=*), INTENT(in) :: AttName
    REAL(kind=dblkind), DIMENSION(:), INTENT(in) :: AttVal

    INTEGER :: attlen

    attlen=SIZE(AttVal)

    ncstatus=nf_put_att_double(file%id,var%id,AttName,NF_DOUBLE,attlen,AttVal(1))
    IF(ncstatus .NE. NF_NOERR) THEN
       CALL handle_ncerr('ncPutAtt (ncPutAtt_dblarr)', ncstatus)
    END IF
  END SUBROUTINE ncPutAtt_dblarr

  SUBROUTINE ncGlobalPutAtt_text(file,AttName,AttVal)

    IMPLICIT NONE

    TYPE(ncfile), INTENT(in) :: file
    CHARACTER (len=*), INTENT(in) :: AttName
    CHARACTER (len=*), INTENT(in) :: AttVal

    INTEGER :: attlen

    attlen=LEN_TRIM(AttVal)

    ncstatus=nf_put_att_text(file%id,NF_GLOBAL,AttName,attlen,AttVal)
    IF(ncstatus .NE. NF_NOERR) THEN
       CALL handle_ncerr('ncPutAtt (ncPutAtt_text)', ncstatus)
    END IF
  END SUBROUTINE ncGlobalPutAtt_text

  SUBROUTINE ncGlobalPutAtt_real(file,AttName,AttVal)

    IMPLICIT NONE

    TYPE(ncfile), INTENT(in) :: file
    CHARACTER (len=*), INTENT(in) :: AttName
    REAL, INTENT(in) :: AttVal

    INTEGER :: attlen

    attlen=1

    ncstatus=nf_put_att_real(file%id,NF_GLOBAL,AttName,NF_REAL,attlen,AttVal)
    IF(ncstatus .NE. NF_NOERR) THEN
       CALL handle_ncerr('ncPutAtt (ncPutAtt_real)', ncstatus)
    END IF
  END SUBROUTINE ncGlobalPutAtt_real

  SUBROUTINE ncGlobalPutAtt_double(file,AttName,AttVal)

    IMPLICIT NONE

    TYPE(ncfile), INTENT(in) :: file
    CHARACTER (len=*), INTENT(in) :: AttName
    REAL(kind=dblkind), INTENT(in) :: AttVal

    INTEGER :: attlen

    attlen=1

    ncstatus=nf_put_att_double(file%id,NF_GLOBAL,AttName,NF_REAL,attlen,AttVal)
    IF(ncstatus .NE. NF_NOERR) THEN
       CALL handle_ncerr('ncPutAtt (ncPutAtt_double)', ncstatus)
    END IF
  END SUBROUTINE ncGlobalPutAtt_double

  SUBROUTINE ncGlobalPutAtt_dblarr(file,AttName,AttVal)

    IMPLICIT NONE

    TYPE(ncfile), INTENT(in) :: file
    CHARACTER (len=*), INTENT(in) :: AttName
    REAL(kind=dblkind), DIMENSION(:), INTENT(in) :: AttVal

    INTEGER :: attlen

    attlen=SIZE(AttVal)

    ncstatus=nf_put_att_double(file%id,NF_GLOBAL,AttName,NF_DOUBLE,&
         &attlen,AttVal(1))
    IF(ncstatus .NE. NF_NOERR) THEN
       CALL handle_ncerr('ncPutAtt (ncPutAtt_dblarr)', ncstatus)
    END IF
  END SUBROUTINE ncGlobalPutAtt_dblarr

  
  SUBROUTINE ncCreate(file,incmode)

    IMPLICIT NONE

    TYPE(ncfile), INTENT(inout) :: file
    INTEGER, INTENT(in), OPTIONAL :: incmode

    INTEGER :: cmode

    IF(PRESENT(incmode)) THEN
       cmode=incmode
    ELSE
       cmode=NF_NOCLOBBER
    END IF

    ncstatus=nf_create(file%name, cmode, file%id)
    IF(ncstatus .NE. NF_NOERR) THEN
       CALL handle_ncerr('nccreate', ncstatus)
    END IF
    
    file%opened=.TRUE.
    file%defmode=.TRUE.
  END SUBROUTINE ncCreate
  
  SUBROUTINE ncReDef(file)

    IMPLICIT NONE

    TYPE(ncfile), INTENT(inout) :: file

    ncstatus=nf_redef(file%id)
    IF(ncstatus .NE. NF_NOERR) THEN
       CALL handle_ncerr('ncdef', ncstatus)
    END IF

    file%defmode=.TRUE.
  END SUBROUTINE ncReDef
  

  SUBROUTINE ncEndDef(file)

    IMPLICIT NONE

    TYPE(ncfile), INTENT(inout) :: file

    ncstatus=nf_enddef(file%id)
    IF(ncstatus .NE. NF_NOERR) THEN
       CALL handle_ncerr('ncenddef', ncstatus)
    END IF
    
    file%defmode=.FALSE.
  END SUBROUTINE ncenddef

  ! -------------------------------------------------------------------------
  ! INQUIRY FUNCTIONS
  ! -------------------------------------------------------------------------

  ! *********
  ! ncpresent
  ! *********
  
  FUNCTION present_file(file)
    
    IMPLICIT NONE

    TYPE(ncfile), INTENT(in) :: file
    LOGICAL :: present_file

    ! Inquire ID of file/variable/attribute. 
    ncstatus = nf_open(file%name, NF_NOWRITE, file%id)
    
    ! Test for error flag.
    ! If not error, this means the property defined in the file: return true
    ! otherwise return false.
    IF(ncstatus .EQ. NF_NOERR) THEN
       present_file=.TRUE.
    ELSE
       present_file=.FALSE.
       RETURN
    END IF

    ncstatus = nf_close(file%id)
    IF (ncstatus .NE. NF_NOERR) THEN
       CALL handle_ncerr('present_file', ncstatus)
    END IF
    
  END FUNCTION present_file
  
  FUNCTION present_var(file,var)
    
    IMPLICIT NONE

    TYPE(ncfile), INTENT(in) :: file
    TYPE(ncvar), INTENT(in) :: var
    LOGICAL :: present_var

    ! Inquire ID of file/variable/attribute. 
    ncstatus=nf_inq_varid(file%id,var%name,var%id)
    
    ! Test for error flag.
    ! If not error, this means the property defined in the file: return true
    ! otherwise return false.
    IF(ncstatus .EQ. NF_NOERR) THEN
       present_var=.TRUE.
    ELSE
       present_var=.FALSE.
    END IF
  END FUNCTION present_var
  

  FUNCTION ncinqdims(file,ndims) RESULT(status)
    IMPLICIT NONE
    TYPE(ncfile), INTENT(in) :: file
    INTEGER, INTENT(out) :: ndims
    INTEGER :: status

    status=nf_inq_ndims(file%id,ndims)
       PRINT *, nf_strerror(status)
    

  END FUNCTION ncinqdims
  

    
END MODULE mo_netcdf

 
