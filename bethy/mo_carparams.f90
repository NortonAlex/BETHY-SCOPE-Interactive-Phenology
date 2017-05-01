MODULE mo_carparams

  IMPLICIT NONE

! WOK-CHG-2007-07-01 removed tphen0, tphena0, tlmax0, will not be used
  REAL, PARAMETER  :: LAIMAX0       = 6.
  REAL, PARAMETER  :: FCMAX0        = 1.0
  REAL, PARAMETER  :: LAILIM0       = 3.0      
!  REAL, PARAMETER  :: TPHEN0        = 5.0
!  REAL, PARAMETER  :: TPHENA0       = 12.0
!  REAL, PARAMETER  :: TLMAX0        = 15.0
  REAL, PARAMETER  :: AFAUTLEAF     = 0.40
  REAL, PARAMETER  :: ACCOST        = 1.25
  REAL, PARAMETER  :: AQ10F         = 1.5   
  REAL, PARAMETER  :: AQ10S         = 1.5   
  REAL, PARAMETER  :: atauf         = 1.5
  REAL, PARAMETER  :: AAW           = 1.
  REAL, PARAMETER  :: ABETA         = 1.
  REAL, PARAMETER  :: afracs        = 0.2
  REAL, PARAMETER  :: CW0           = 0.5
  REAL, PARAMETER  :: ALBV0         = 0.15
  REAL, PARAMETER  :: Asoil0        = 0.05 
  REAL, PARAMETER  :: ALBS0         = 0.64 
  REAL, PARAMETER  :: FCI1C3        = 0.65 
  REAL, PARAMETER  :: FCI1C4        = 0.37 
  REAL, PARAMETER  :: P0            = 1.01325E5  
  REAL, PARAMETER  :: LAPS          = 6.0E-3    
  REAL, PARAMETER  :: MMC           = 12.E-3   
  REAL, PARAMETER  :: CDRM          = 0.45     
  REAL, PARAMETER  :: EPAR          = 2.2E5  
! WOK-CHG-2007-06-26 removed gddbase (was not used)
!  REAL, PARAMETER  :: gddbase       = 5.0
! gddbase    base temperature for growing degree days
! p0         standard air pressure
! laps       constant lapse rate assumed for calculating actual pressure [K / m]
! mmc        molar mass of carbon [kgC / mol(CO2)]
! cdrm       carbon content of dry matter [gC/g(dry matter)] 
! epar       energy content of PAR [J / mol(photons)]=(4.6 mol/MJ PAR)**-1
  REAL, PARAMETER  :: KARMANN       = 0.41    
  REAL, PARAMETER  :: RLAM0         = 3.1512E6 
  REAL, PARAMETER  :: RLAM1         = -2.38E3 
  REAL, PARAMETER  :: RLAM2         = 3.834E6  
  REAL, PARAMETER  :: ZSPEEDMIN     = 0.5      
  REAL, PARAMETER  :: ZLAIMIN       = 1.0e-12
  REAL, PARAMETER  :: ZRHOSPARMIN   = 0.0      
  REAL, PARAMETER  :: Z0VEGMIN      = 1.0e-2   
  REAL, PARAMETER  :: ZGAMIN        = 3.0e-2    
  REAL, PARAMETER  :: ZFCMIN        = 1.0e-3   
! WOK-CHG-090401: the next parameter is now obsolete:
!  REAL, PARAMETER  :: ZGCMIN        = 1.0e-12     
  REAL, PARAMETER  :: ZWATERMIN     = 1.0e-15
 ! REAL, PARAMETER  :: ZGROWTHRESMIN = 0.0
  REAL, PARAMETER  :: ZZENITMIN     = 1.0e-3    
  REAL, PARAMETER  :: JMAXMIN       = 1.0e-12   
! WOK-GHG-090401: realistic cuticular conductance [m/s]
! Meidner, H. 1986. Cuticular conductance and the humidity response of stomata
! J. Exp. Bot. 37:571-525
  REAL, PARAMETER  :: ZGSMIN        = 4e-5
!  REAL, PARAMETER  :: ZGSMIN        = 0.0     
  REAL, PARAMETER  :: ZMUEMIN       = 0.0174524 
  REAL, PARAMETER  :: ZRIMAX        = 0.137
  REAL, PARAMETER  :: ZBTHICK       = 9.830 
  REAL, PARAMETER  :: ZROOTTHICK    = 1.230 
  REAL, PARAMETER  :: ZFFSCALE      = 5.834   
  REAL, PARAMETER  :: ZBBSCALE      = 3.14
  REAL, PARAMETER  :: ZEPSW         = 1.E-3
  REAL, PARAMETER  :: TMELT         = 273.16

  INTEGER,  PARAMETER :: NL         = 3        ! number of layers within canopy 



!-------------------------------------------------------------------------
! PFT definitions:
! 1: tropical broadleaf evergreen tree
! 2: tropical broadleaf deciduous tree
! 3: temperate broadleaf evergreen tree
! 4: temperate broadleaf deciduous tree
! 5: evergreen coniferous tree
! 6: deciduous coniferous tree
! 7: evergreen shrub
! 8: deciduous shrub
! 9: C3 grass
!10: C4 grass
!11: tundra
!12: swamp
!13: arable crop
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
! specific leaf area, in Diss. only for decidous Plants, values from Schulze
!                     et al. (1994); in m^2 / kg(dry mass).
! Reference: Schulze, E.-D., Kelliher, F.M., Kšrner, Ch., Lloyd, J. and Leuning, R. 1994. 
! Relationships among maximum stomatal conductance, ecosystem surface 
! conductance, carbon assimilation rate, and plant nitrogen nutrition: a 
! global ecology scaling exercise. Ann. Rev. Ecol. Syst. 25, 629-660.
!-------------------------------------------------------------------------
! WOK-CHG-0706 did not conform to new PFTs, updated from Diss., but need data for evergreens
!  REAL, PARAMETER  :: SLAV(0:11)=(/0., 7.8, 13.2, 3.4, 51., 114., 3.5, 11.3, &
!		& 16.9, 16.9, 23.6, 263. /)
!  REAL, PARAMETER  :: SLAV(0:13)=(/0., 9.9, 14.1, 5.7, 11.5, 4.1, 11.3, 6.9, &
!		& 11.5, 16.9, 16.9, 20.0, 16.9, 25.3 /)
!WOK-CHG-070926 set SLA for tundra to more typical value of evergreen shrub
  REAL, PARAMETER  :: SLAV(0:13)=(/0., 9.9, 14.1, 5.7, 11.5, 4.1, 11.3, 6.9, &
		& 11.5, 16.9, 16.9, 6.9, 16.9, 25.3 /)
  REAL, PARAMETER  :: RSOIL(6) =(/0.07, 0.15, 0.10, 0.20, 0.18, 0.35/)
  REAL, PARAMETER  :: ZBSOIL(5)=(/0.065, 0.254, 0.913, 2.902, 5.700/)

  INTEGER, DIMENSION (0:13), PARAMETER :: xc4flg=(/0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, &
                                                & 0, 0, 0/)
  INTEGER, DIMENSION (0:13), PARAMETER :: xclass=(/6, 1, 1, 1, 1, 1, 1, 2, 2, 3, 3, &
                                                 & 5, 4, 4/)
  REAL, DIMENSION(0:13), PARAMETER :: xhv = (/0, 300, 150, 150, 150, 150, 150, 10, &
                                            & 10, 10, 10, 3, 3, 6/)
  REAL, DIMENSION(0:13), PARAMETER :: xvm = (/0, 60, 90, 41, 35, 29, 53, 52, 160, &
                                            & 42, 8, 20, 20, 117/)
  REAL, DIMENSION(0:13), PARAMETER :: xjm = (/0, 118, 179, 82, 70, 52, 95, 102, 266, &
                                           & 80, 140, 37, 37, 220/)


  REAL, PARAMETER :: AEC=59356.    ! EC      ACTIVATION ENERGY FOR KC [J/MOL]
  REAL, PARAMETER :: AEO=35948.    ! EO      ACTIVATION ENERGY FOR KO [J/MOL]
  REAL, PARAMETER :: AEV=58520.    ! EV      ACTIVATION ENERGY FOR VCMAX [J/MOL]
  REAL, PARAMETER :: AER=45000.    ! ER      ACTIVATION ENERGY FOR DARK RES [J/MOL]
  REAL, PARAMETER :: AEK=50967.    !         Q10=2 (Collatz et al. 1992)
  REAL, PARAMETER :: Atgam=1.7e-6
  REAL, PARAMETER :: akc0=460.E-6  ! KC0     MICHAELIS-MENTEN CONSTANT FOR CO2 
                                   !         AT 25C [MOL(CO2) / MOL(AIR)]
  REAL, PARAMETER :: ako0=330.E-3  ! KO0     MICHAELIS-MENTEN CONSTANT FOR O2 
  REAL, PARAMETER :: aALPHA=0.28   ! ALPHA   EFFICIENCY OF OF PHOTON CAPTURE 
                                   !         AT 25C [MOL(O2) / MOL(AIR)]  
  REAL, PARAMETER :: aalc4=0.04    ! ALC4    EFFECTIVE QUANTUM EFFICIENCY
  REAL, PARAMETER :: ajmt=1/25.


 ! ANorton 12/2014 fluorescence calculation
REAL, PARAMETER   :: Hkin   = 220000.
REAL, PARAMETER   :: akp    = 80.E-3   ! KP0  MICHAELIS-MENTEN CONSTANT FOR CO2 FOR C4 PLANT ???
REAL, PARAMETER   :: Rd_f   = 0.5    ! [umol/m2/s] dark respiration rate at 25 oC
REAL, PARAMETER   :: Vpmo   = 220    ! [umol/m2/s] maximum PEP regeneration rate (C4 only)(~50;60;220)
REAL, PARAMETER   :: gbs    = 3E3    ! [umol/m2/s]  bundle sheath conductance to CO2(C4 only)
REAL, PARAMETER   :: gcmin  = 1E-6   ! [m s-1] minimum stomatal conductance (if stomata are closed).
REAL, PARAMETER   :: Ejmax  = 45000  ! [J/mol] for electron transport
REAL, PARAMETER   :: lam    = 750.   ! []  marginal cost of assimilation (~700)

!  parameters (at optimum temperature)
REAL, PARAMETER   :: fipo   = 0.82   !  dark photochemistry fraction (Genty et al., 1989)
REAL, PARAMETER   ::   al   = 1E1    !  smoothing parameter  (choose a value [3,100])
REAL, PARAMETER   :: x_mas  = .6     !  partitioning factor of the electron transport rate (Massad et al., 2007)

 
END MODULE mo_carparams
