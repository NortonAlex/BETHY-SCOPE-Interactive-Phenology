MODULE mo_constants

  IMPLICIT NONE

! bounds for time stepping
!  INTEGER, PARAMETER :: maxyears = 42
  INTEGER, PARAMETER :: maxyears = 29
!  INTEGER, PARAMETER :: maxyears = 6
  INTEGER, PARAMETER :: maxdpy = 366
  INTEGER, PARAMETER :: maxtspd  = 24
  INTEGER, PARAMETER :: maxndays = 1
  INTEGER, PARAMETER :: maxdiu = (1+(maxdpy/10))*maxyears*maxtspd

!CCCC DIMENSION SETTINGS

  INTEGER, PARAMETER :: nv = 3
  INTEGER, PARAMETER :: nprof = 134
  INTEGER, PARAMETER :: maxhor = 4
  INTEGER, PARAMETER :: nplayer = 5
  INTEGER, PARAMETER :: npft = 1
  INTEGER, PARAMETER :: nci = 3

!cccc nv         maximum number of vegetation types at one grid point
!cccc nshmax     maximum number of layers (or shaded layers) for 2-flux scheme
!cccc            (see 'dayphtrans' and 'fracpar')
!cccc maxiayears maximum number of years for interannual calculations
!cccc            should be set to number of years for which interannual
!cccc            data are available (excludes spin-up)
!cccc ngyears    maximum number of years for storing daily weather
!cccc            should be set to total number of simulation years,
!cccc            including spin-up
!cccc            (memomory-critical, check subroutine 'dailyweather')
!cccc nm         maximum number of months
!cccc npft 	 number of pft's used
!cccc nci        number of Carbon-tracers (12C,13C,14C) for fractionation

!CCCC physical CONSTANTS

  REAL, PARAMETER :: pi = 3.1415926535879
  REAL, PARAMETER :: ar = 8.314
  REAL, PARAMETER :: cpd = 1005.46
  REAL, PARAMETER :: g = 9.80665
  REAL, PARAMETER :: stbo = 5.6697e-8
  REAL, PARAMETER :: amd = 28.964

!cccc pi         half circle radiants
!cccc ar         gas constant [J / K mol]
!cccc cpd        specific heat of dry air at constant pressure [J / kg K]
!cccc g          acceleration due to gravity [m / s**2]
!cccc stbo       Stefan-Boltzmann constant [W / m**2 K**4]
!cccc amd        molar mass of air [g / mol]

  REAL, PARAMETER :: rd    = 287.05      ! gas constant for dry air in J/K/kg
  REAL, PARAMETER :: rv    = 461.51      ! gas constant for water vapour in J/K/kg
  REAL, PARAMETER :: vtmpc1= rv/rd-1.    ! dimensionless auxiliary constant
  REAL, PARAMETER :: zeps=rd/rv
  REAL, PARAMETER :: secpd=86400.        ! seconds per day

! Constants used for computation of saturation mixing ratio
! over liquid water (*c_les*) or ice(*c_ies*)
  REAL, PARAMETER :: c1es  = 610.78              ! 
  REAL, PARAMETER :: c2es  = c1es*zeps           ! 
  REAL, PARAMETER :: c3les = 17.269              ! 
  REAL, PARAMETER :: c3ies = 21.875              ! 
  REAL, PARAMETER :: c4les = 35.86               ! 
  REAL, PARAMETER :: c4ies = 7.66                ! 


!CCCC tuning CONSTANTS

  REAL, PARAMETER :: wmaxmin0 = 10.
  REAL, PARAMETER :: relh0 = 0.96
  REAL, PARAMETER :: reld0 = 0.49
  REAL, PARAMETER :: karman = 0.41
  REAL, PARAMETER :: href = 2.
  REAL, PARAMETER :: rz0 = 0.0999995
  REAL, PARAMETER :: az0 = 4.35680
! wmaxmin0  minimum bucket size [mm]
! relh0     humidity parameter 
! reld0     humidity parameter
! karman    karman constant used for computation of aerodynamic conductance
! href      reference height   -"-
! rz0                          -"-     
! az0                          -"-


!CCCC CALENDER CONSTANTS


   INTEGER, PARAMETER :: jan = 1, feb = 1, mar = 1, apr = 1, may = 1, jun = 1
   INTEGER, PARAMETER :: jul = 1, aug = 1, sep = 1, oct = 1, nov = 1, dec = 1
!   INTEGER, PARAMETER :: jan = 31, feb = 28, mar = 31, apr = 30, may = 31, jun = 30
!   INTEGER, PARAMETER :: jul = 31, aug = 31, sep = 30, oct = 31, nov = 30, dec = 31

   INTEGER, PARAMETER :: janf = 1
   INTEGER, PARAMETER :: febf = janf + jan
   INTEGER, PARAMETER :: marf = febf + feb
   INTEGER, PARAMETER :: aprf = marf + mar
   INTEGER, PARAMETER :: mayf = aprf + apr
   INTEGER, PARAMETER :: junf = mayf + may
   INTEGER, PARAMETER :: julf = junf + jun
   INTEGER, PARAMETER :: augf = julf + jul
   INTEGER, PARAMETER :: sepf = augf + aug
   INTEGER, PARAMETER :: octf = sepf + sep
   INTEGER, PARAMETER :: novf = octf + oct
   INTEGER, PARAMETER :: decf = novf + nov

   INTEGER, PARAMETER :: janl = febf - 1
   INTEGER, PARAMETER :: febl = marf - 1
   INTEGER, PARAMETER :: marl = aprf - 1
   INTEGER, PARAMETER :: aprl = mayf - 1
   INTEGER, PARAMETER :: mayl = junf - 1
   INTEGER, PARAMETER :: junl = julf - 1
   INTEGER, PARAMETER :: jull = augf - 1
   INTEGER, PARAMETER :: augl = sepf - 1
   INTEGER, PARAMETER :: sepl = octf - 1
   INTEGER, PARAMETER :: octl = novf - 1
   INTEGER, PARAMETER :: novl = decf - 1
   INTEGER, PARAMETER :: decl = decf + dec - 1

   INTEGER, PARAMETER :: jan2 = janf + jan/2
   INTEGER, PARAMETER :: feb2 = febf + feb/2
   INTEGER, PARAMETER :: mar2 = marf + mar/2
   INTEGER, PARAMETER :: apr2 = aprf + apr/2
   INTEGER, PARAMETER :: may2 = mayf + may/2
   INTEGER, PARAMETER :: jun2 = junf + jun/2
   INTEGER, PARAMETER :: jul2 = julf + jul/2
   INTEGER, PARAMETER :: aug2 = augf + aug/2
   INTEGER, PARAMETER :: sep2 = sepf + sep/2
   INTEGER, PARAMETER :: oct2 = octf + oct/2
   INTEGER, PARAMETER :: nov2 = novf + nov/2
   INTEGER, PARAMETER :: dec2 = decf + dec/2


   INTEGER, PARAMETER :: jan_f = 31, feb_f = 28, mar_f = 31, apr_f = 30, may_f = 31, jun_f = 30
   INTEGER, PARAMETER :: jul_f = 31, aug_f = 31, sep_f = 30, oct_f = 31, nov_f = 30, dec_f = 31

   INTEGER, PARAMETER :: janf_f = 1
   INTEGER, PARAMETER :: febf_f = janf_f + jan_f
   INTEGER, PARAMETER :: marf_f = febf_f + feb_f
   INTEGER, PARAMETER :: aprf_f = marf_f + mar_f
   INTEGER, PARAMETER :: mayf_f = aprf_f + apr_f
   INTEGER, PARAMETER :: junf_f = mayf_f + may_f
   INTEGER, PARAMETER :: julf_f = junf_f + jun_f
   INTEGER, PARAMETER :: augf_f = julf_f + jul_f
   INTEGER, PARAMETER :: sepf_f = augf_f + aug_f
   INTEGER, PARAMETER :: octf_f = sepf_f + sep_f
   INTEGER, PARAMETER :: novf_f = octf_f + oct_f
   INTEGER, PARAMETER :: decf_f = novf_f + nov_f

   INTEGER, PARAMETER :: janl_f = febf_f - 1
   INTEGER, PARAMETER :: febl_f = marf_f - 1
   INTEGER, PARAMETER :: marl_f = aprf_f - 1
   INTEGER, PARAMETER :: aprl_f = mayf_f - 1
   INTEGER, PARAMETER :: mayl_f = junf_f - 1
   INTEGER, PARAMETER :: junl_f = julf_f - 1
   INTEGER, PARAMETER :: jull_f = augf_f - 1
   INTEGER, PARAMETER :: augl_f = sepf_f - 1
   INTEGER, PARAMETER :: sepl_f = octf_f - 1
   INTEGER, PARAMETER :: octl_f = novf_f - 1
   INTEGER, PARAMETER :: novl_f = decf_f - 1
   INTEGER, PARAMETER :: decl_f = decf_f + dec_f - 1

   INTEGER, PARAMETER :: jdpyear = decl
   INTEGER, PARAMETER :: nmonth = 12


! ..!CCCC OPTIONS DECLARATIONS

   LOGICAL :: ragain
   INTEGER :: tspd, tspy, maxts, outt, outyear
   INTEGER :: ryear, rmonth, jday, sday, gday, inho
   INTEGER, DIMENSION(12) :: tspm, atspm
   REAL :: ftspd
   INTEGER :: snflg
   INTEGER :: fflg, stflg
   INTEGER :: nsh, sunflg, nsclflg
   REAL :: zdiagt
   REAL :: frd, frcost
   REAL :: ca     ! ambient air CO2 concentration [mol(CO2)/mol(air)] 
   REAL :: u0=3.  ! wind speed with default value 3.0 m/s
   REAL :: fci0c3, fci0c4
   REAL :: hcpl, bw0, be0, kf0, psiw0
   REAL :: rdpthmin, hovp
   REAL, DIMENSION (nplayer) :: player
   REAL, DIMENSION (5) :: rdpc, rdpsc
   REAL, DIMENSION (20:90) :: frtable
   LOGICAL :: tmade
   INTEGER :: jdlast

!cccc irestart   0: no restart run, 1: restart run
!cccc year0      initial year for simulation
!cccc nrun       number of years in this simulation
!cccc tspd       time steps per day
!cccc tspm       time steps per month
!cccc atspm      accumulated time steps monthly
!cccc ftspd      float of tspd
!cccc maxts      numberof time steps in this run
!cccc outt       number of output values per year
!cccc lonw       westernmost longitude executed (<0: W)
!cccc lone       easternmost longitude executed (<0: W)
!cccc            default is -180 and 180 (if both are set 0 or not specified 
!cccc            in 'bethy.options'
!cccc latn       first and northernmost latitude executed (<0: S)
!cccc lats       last and southernmost latitude executed (<0: S)
!cccc            order of latitudes in input data must be north to south
!cccc fautleaf   ratio of leaf to total maintenance respiration
!cccc ccost      construction cost of biosynthesis [mol(C)/mol(C)]
!cccc fparmax    maximum landscape-scale FPAR (0.97, corr. to LAI = 6)
!cccc frcost     plant respiration costs as fraction of GPP [fraction]
!cccc q10        Q10 for heterotroph respiration related to air temperature 
!cccc dflg       switch for daily or monthly calculation
!cccc wflg       switch for stochastic weather generator

! CALENDER SETTINGS
   
   INTEGER, DIMENSION (12) :: jdfst, jdlst, midday, ndays, midmoday, rdays,jdfst_f,jdlst_f
   
! .. Data Statements ..
   DATA ndays  /jan, feb, mar, apr, may, jun, jul, aug, sep, oct, nov, dec/
   DATA jdfst  /janf, febf, marf, aprf, mayf, junf, julf, augf, sepf, octf, &
             &  novf, decf/
   DATA jdlst  /janl, febl, marl, aprl, mayl, junl, jull, augl, sepl, octl, &
             &  novl, decl/
   DATA jdfst_f  /janf_f, febf_f, marf_f, aprf_f, mayf_f, junf_f, julf_f, augf_f, &
             & sepf_f, octf_f, novf_f, decf_f/
   DATA jdlst_f  /janl_f, febl_f, marl_f, aprl_f, mayl_f, junl_f, jull_f, augl_f, &
             & sepl_f, octl_f, novl_f, decl_f/
   DATA midday /jan2, feb2, mar2, apr2, may2, jun2, jul2, aug2, sep2, oct2, &
            &  nov2, dec2/
   DATA midmoday /16, 44, 75, 105, 135, 166, 197, 228, 258, 289, 319, 350 /
   DATA rdays /31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/

   ! pjr 04/04/28 some flags for different parameter bounding
   integer, parameter ::  bounded_below = 1
   integer, parameter :: bounded_above = -1
   integer, parameter :: bounded = 2
   integer, parameter :: unbounded = 0
   integer, parameter :: direct_flux = 0  ! it's a basis function -> xpf=0.

   ! spin up years in cbalance
   INTEGER, PARAMETER :: aspin = 100

 ! ANorton 12/2014 .. constants that are not defined in this module neither
 ! fluo_param for fluorescence calculation

  REAL, PARAMETER       ::       Av  = 6.02214E23  ! [mol-1]         Constant of Avogadro
  REAL, PARAMETER       ::        h  = 6.6262E-34  ! [J s]           Planck's constant
  REAL, PARAMETER       ::    kappa  =  0.4        ! % []            Von Karman constant
  REAL, PARAMETER       ::     MH2O  = 18          !% [g mol-1]    Molecular mass of water
  REAL, PARAMETER       ::     Mair  = 28.96       !% [g mol-1]    Molecular mass of dry air
  REAL, PARAMETER       ::     MCO2  = 44          !% [g mol-1]    Molecular mass of carbon dioxide
  REAL, PARAMETER       ::  sigmaSB  = 5.67E-8    !% [W m-2 K-4]  Stefan Boltzman constant
  REAL, PARAMETER       ::  deg2rad  = pi/180     !% [rad]         Conversion from deg to rad
  REAL, PARAMETER       ::     rhoa  = 1.2047     !% [kg m-3]      Specific mass of air
                                                  ! This can be calculated from


END MODULE mo_constants
