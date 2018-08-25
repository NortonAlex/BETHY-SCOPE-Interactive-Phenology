MODULE fluo_param

!function [varargout]=parameters(varargin)
!% 
!% parameters.m
!% authors: Dr. ir. Christiaan van der Tol (tol@itc.nl)
!%          ir. Joris Timmermans (j_timmermans@itc.nl)
!% 
!% date:      7	Dec 2007
!% update:   13	Dec 2007    included fluorescence matrix (CvdT)
!% update:   29  Jan 2008    converted into a function (JT)
!% update    13  Feb 2008        - Derived parameters are all in
!!% subfunctions
!%                               - subfunctions can be called from outside
!%                                 the main function (JT)
!% update    14  Feb 2008        - parameter hc is added =2 (JT). 
!%                               - removed prmave         
!% update    18  Feb 2008   changed htable in hctable (as well as in "parameters)
!%                               changed input_files (seperate files instead
!%                               of single file)... (JT)
!% update    18  Feb 2008   Have put the tables (LAI,z,hc) in Load_Files (JT)
!% update    20  Oct 2008   Changed dx = x(2)-x(1) into dx = abs(x(2)-x(1)) (JT and CvdT)

!% update     7  Nov 2008   Changed layout (CvdT)
!  Fortran version : June 2012 Ernest Koffi, IPSL, Jusssieu, Paris, France 
!%
!%The function ‘parameters’ contains the parameters used in the model ‘ebal’. 
!% parent:       master.m
!% children:     prospect.m and leafangles.m
!%
!% Table of contents of parameters.m
!% 
!% 0. Initialsiations        Not to be changed
!% 1. Simulation options     Specify what you want to calculate (for 
!%                           example: fluorescence), start and end date of 
!%                           the simulations (in case of a time series), and 
!%                           methods to be used
!% 2. Declare file names     Specify the location and file names of model
!%                           input
!% 3. Location               Specify the Longitude and Latitude (needed to 
!%                           calculate solar angles)
!% 4. Numerical parameters	Specify number of layers, number of azimuth 
!%                           angles, maximum number of iterations, maximum 
!%                           allowed energy balance error, and the wavelength resolution.
!% 5. Radiative Transfer     Specify Prospect parameters for reflectance and 
!%                           transmittance of leaves, and Hapke’s parameters 
!%                           for soil, or a table for soil reflectance as a 
!%                           function of wavelength
!% 6. MODTRAN output         Specify a spectrum of incoming light, if this 
!%                           is not available for each time step separately
!% 7.  Fluorescence          Specify the fluorescence yield efficiency and
!%                           the peak ratio of PSI/PSII
!% 8. Energy balance, micrometeorological and physiological parameters
!% 	8.1 Vegetation          Specify type of vegetation (C3 or C4), and bio-
!%                           chemical parameters for photosynthesis
!% 	8.2 Canopy structure	Specify height of vegetation, LAI, roughness 
!%                           length, etc. Either fixed or as a function of time
!% 	8.3 Soil                Thermal properties and surface resistance of the soil
!% 
!% Appendix A-H              Not to be changed
!% Appendix A                Fixed parameters (SAIL)
!% Appendix A                Secondary parameters: these parameters are 
!%                           derived from the above specified parameters
!% Appendix C                Derived Leaf function (not to be changed)
!% Appendix D                Derived Radiation function
!% Appendix E                Derived Leaf function 
!% Appendix F                Soil parameters function
!% Appendix G                Layer parameters function
!% Appendix H                Look-up table for inversion of Planck function

!%% 0. Initialisations
!% Subfunctions call
!if nargin>0    
!    varargout   =   feval(varargin{1},varargin{2:end});   
!    if ~iscell(varargout)
!        varargout={varargout};
!    end
!    return    
!end

USE fluo_func

IMPLICIT NONE 

!% globals
!global wl wlo wlt nwl nwlo nwlt iwlo iwlt litab lazitab nlazi nli x dx nl resolution
!global psi tto tto_ psi_ noa 
!global Mf Mb wlfi wlfo iwlfi iwlfo nwlfi nwlfo
!global Vcmax_file Jmax_file Lambda_file
!global Vcmax Vpmo Jmax lambda Rdopt gcmin Type 
!global gbs  PSIs  hc z zo d rb Cd w rbs rss rwc LAI
!global prm_ prm q GAM fEsun fEsky fEskyt lidf T_LUT
!global CR CD1 Psicor CSSOIL


!global wl wlo wlt nwl nwlo nwlt iwlo iwlt litab lazitab nlazi nli x dx nl
!resolution
!global psi tto tto_ psi_ noa 

INTEGER                                  :: nlazi,nli,noa,nwl,nl
INTEGER                                  :: nwlfi,nwlfo,nwlt 
INTEGER                                  :: nwlP,nwlF,nwlS,nwlO,nwlE,nwlPAR

INTEGER                                  :: spectral_nreg 
REAL, DIMENSION(3)                       :: spectral_start, spectral_end ,spectral_res

REAL, ALLOCATABLE, DIMENSION(:)          :: wlS,wlP,wlT,wlF,wlE,wlPAR,wlO 
REAL, ALLOCATABLE, DIMENSION(:)          :: wl,resolution,xlay 
REAL, ALLOCATABLE, DIMENSION(:)          :: psi_,tto_
REAL, ALLOCATABLE, DIMENSION(:)          :: fEsun, fEsky
REAL, ALLOCATABLE, DIMENSION(:)          :: wlfi,wlfo
REAL, ALLOCATABLE, DIMENSION(:)          :: lidf
REAL, ALLOCATABLE, DIMENSION(:)          :: kChlrel, tran, refl 
REAL, ALLOCATABLE, DIMENSION(:)          :: rho,tau,rs 
REAL, ALLOCATABLE, DIMENSION(:)          :: prm
REAL, ALLOCATABLE, DIMENSION(:,:)        :: MfI, MbI
REAL, ALLOCATABLE, DIMENSION(:,:)        :: MfII, MbII
REAL, ALLOCATABLE, DIMENSION(:,:)        :: prm_

INTEGER, ALLOCATABLE, DIMENSION(:)       :: iwlo,iwlt, litab,lazitab
INTEGER, ALLOCATABLE, DIMENSION(:)       :: iwlfi,iwlfo
INTEGER, ALLOCATABLE, DIMENSION(:)       :: IwlP

REAL                                     :: dx
REAL                                     :: gbs,w 
REAL                                     :: Vcmax_ref ! %[umol/m2/s]        maximum carboxylation capacity
REAL                                     :: Jmax_ref  ! %[umol/m2/s]        maximum electron transport rate (~100;120;270)


INTEGER, PARAMETER                       :: calcfluor = 1  ! %   calculate !chlorophyll !fluorescence
INTEGER, PARAMETER                       :: simulation = 1


INTEGER, PARAMETER                       :: iec = 0    ! 1, we make speaked the model 
INTEGER                                  :: ial  =0    ! The first time the fluo module is clalled, it is put to 1 

INTEGER,PARAMETER                        :: nregion=20   ! Number of regions where the spectrum is computed

! Output variables from rtmo 
REAL                                     :: Rnhs , Rnus
REAL                                     :: Eouto,Eoutt
REAL                                     :: km, Kext 
REAL                                     :: P 
REAL,ALLOCATABLE,DIMENSION(:)            :: Eout_,Lo_, Lout_
REAL,ALLOCATABLE,DIMENSION(:,:,:)        :: Rnuc,Pnuc,Pnuc_Cab
REAL,ALLOCATABLE,DIMENSION(:)            :: Rnhc,Pnhc,Pnhc_Cab
REAL,ALLOCATABLE, DIMENSION(:,:)         :: Eplu_, Emin_
REAL,ALLOCATABLE, DIMENSION(:)           :: Esun_, Esky_
REAL,ALLOCATABLE, DIMENSION(:)           :: fEsuno,fEskyo,fEsunt,fEskyt
REAL,ALLOCATABLE, DIMENSION(:)           :: Ps , Po,  Pso

! Output variables from rtmt
REAL                                     :: Lot, Loutt_t
REAL                                     :: Rnust,Rnhst
REAL,ALLOCATABLE,DIMENSION(:)            :: Rnhct 
REAL,ALLOCATABLE, DIMENSION(:)           :: Eplu, Emin
REAL,ALLOCATABLE,DIMENSION(:,:,:)        :: Rnuct

! Output from rtmf 
REAL,  ALLOCATABLE, DIMENSION(:)         :: LoF,Fhem,Fiprof

REAL,ALLOCATABLE,DIMENSION(:)            :: Fc !% fraction of sunlit leaves per leaf layer

! Output variables from biochemical 
REAL, ALLOCATABLE, DIMENSION(:)          :: An,Fl,rcw,f,Ci
REAL, ALLOCATABLE, DIMENSION(:)          :: Wc,Wj,Vcmax
REAL, ALLOCATABLE, DIMENSION (:)         :: Jmax,zetan,Fe

! Other output 
REAL, ALLOCATABLE,DIMENSION(:)           :: Rnch
REAL, ALLOCATABLE,DIMENSION(:,:,:)       :: Rncu
REAL                                     :: Rnsh,Rnsu


!!!! PARAMETERS 
! 1. PROSPECT 
!%1.1. green leaves
REAL                                     :: Cab    ! %[ug/cm2] Concentration chlorophyl (~60)
REAL                                     :: Cdm    ! %[g/cm2]  Concentration dry material (~0.012;0.07)
REAL                                     :: Cw     ! %[cm]     Concentration water (~0.009)
REAL                                     :: Csm    ! % Concentration senescent material (~0.0)
REAL                                     :: N      ! %  thickness parameter (~1.4;2)
                                    ! We need to take care of this
                                    ! parameter in Fortran , must be nn to
                                    ! avoid with n already set maybe
                                    ! somewhere

!%1.2 leaf parameters
REAL                                     :: f_b    ! %   fraction of brown material
REAL, DIMENSION(9)                       :: leafbio 

!%1.2. brown leaves (currently not used)
REAL                                     :: Cab_b  ! %[ug/cm2] Concentration chlorophyl (~60)
REAL                                     :: Cdm_b  ! %[g/cm2]  Concentration dry material (~0.005)
REAL                                     :: Cw_b   ! %[cm]     Concentration water (~0.009)
REAL                                     :: Cs_b   ! %   Concentration senescent material (~0.0)
REAL                                     :: N_b    ! %  thickness parameter (~2.0)

!% 1.3 Spectral parameter
REAL                                     :: rho_thermal ! % Broadband thermal reflectance
REAL                                     :: tau_thermal !% Broadband thermal translisttance



! 2. Leaf biochemical
REAL                                      :: Vcmo    ! (um m-2 s) Maximum carboxylation capacity (at optimum temperature)
REAL                                      :: m       ! Ball-Berry stomatal conductance  parameter
INTEGER                                   :: option  ! Photoochemicl pathway : 0=C3, 1=C4
REAL                                      :: kV      ! Extinction coefficient for Vcmax in the vertical (maximum at the
                                                     ! top). 0 for uniform Vcmax
REAL                                      :: Rdparam ! Respiration  = Rdparam * Vcmax
REAL, DIMENSION(5)                        :: Tparams ! See PFT.xls file. These are five parameters specifying the temperature response
                                           ! TO CLARIFY THIS .....
INTEGER                                   :: tempcor ! [] boolean (0 or 1) whether or not temperature correction to Vcmax has to be
                                                     ! applied. 
REAL                                      ::  stressfactor !  []   optional input: stress factor to reduce Vcmax (for example
                                                     ! soil moisture, leaf age). Default value = 1.

REAL                                      :: Kcopt   !% [ubar] kinetic coefficient for CO2 (Von Caemmerer and Furbank, 1999)
REAL                                      :: Koopt   !% [mbar] kinetic coeeficient for  O2 (Von Caemmerer and Furbank, 1999)
REAL                                      :: Kf0     !% [] rate constant for fluorescence
REAL                                      :: Kd0     !% [] rate constant for thermal deactivation at Fm
REAL                                      :: Kpc0    !% [] rate constant for photochemisty
REAL                                      :: atheta 
REAL                                      :: lam     ![]    Cowan's stomatal parameter (not used in this version of SCOPE) 
REAL                                      :: Jmo     ! Jmax : maximum electron transport rate table  (not used in this version of scope)



INTEGER                                   ::  model_choice    ! Fluorescence model_choice integer with
                                       !  0: with sustained quenching after
                                       !  drought, as in Lee et al. (2013)
                                       !  1: calibrated for cotton data set: no
                                       !  drought
!Leaf biochemical (Magnani model) 
REAL                                       :: Tyear  ! % [oC] mean annual temperature
REAL                                       :: beta   ! % [] fraction of photons partitioned to PSII !(0.507 for C3, 0.4 for C4; Yin et al. 2006; Yin and Struik 2012)
REAL                                       :: NPQs   ! !% NPQs   % [s-1] rate constant of sustained thermal !dissipation, normalized to (kf+kD) (=kNPQs'; Porcar-Castell 2011)  
REAL                                       :: qLs   !  !% qLs  % [] fraction of functional reaction centres !(Porcar-Castell 2011)
 

!3. Fluorescence
REAL                                      :: fqe1     ! Fluorescence quantum yield efficiency at photosystem level
REAL                                      :: fqe2     ! Fluorescence quantum yield efficiency at photosystem level
REAL                                      :: prat    ! %  PSI/PSII peak ratio
REAL                                      :: gcmin   ! %[m s-1] minimum stomatal conductance (if stomata are closed).

REAL                                      :: freq_Fsmin                         ! Minimum frequeency to compute Fs
REAL                                      :: freq_Fsmax                         ! Maximum frequency to compute Fs



!4. Soil
INTEGER                                   :: spectrum   ! Spectrum number (column in the data base soil_file)
REAL                                      :: rss        ! (s m-1) Soil resistane for evaporation from the pore space
REAL                                      :: rs_thermal ! Broadband soil reflectance in the thermal range (1-emissivity)
REAL                                      :: cs         ! (J m-2 K-1) Volumetric heat capacity of the soil
REAL                                      :: rhos       ! (kg m-3) Specific mass of the soil
REAL                                      :: lambdas    ! (J m-1 K-1) Heat conductivity of the soil
REAL                                      :: SMC        ! Volumetric soil moisture content in the root zone
REAL                                      :: GAM

!5. Meteo (values in data files. Thes values are per default, otherwise files are provided)
REAL                                      :: z        ! (m) measurement height of meteorological data
REAL                                      :: Rin      ! (W m-2) brodband incoming shortwave radiation (0.4 - 2.5 um)
REAL                                      :: Ta       !  Air temperature
REAL                                      :: Rli      ! (W m-2) broadland incoming longwave radiation (2.5 - 50 um)
REAL                                      :: pa       ! (hPa) Air pressure
REAL                                      :: ea       ! (hPa) Atmospheric vapour pressure
REAL                                      :: u        ! (m s-1) Wind speed at height z
REAL                                      :: Ca       ! (ppm) Atmospheric concentration
REAL                                      :: Oa       ! (ppm) Atmospheric O2 concentration

!6. Canopy
REAL                                      :: LAI        ! (m2 m-2) Leaf area index
REAL                                      :: hc         ! (m)Leaf vegetation height
REAL                                      :: zo         ! %[m] roughness height for momentum transport (~0.12*hc)
REAL                                      :: LIDFa      ! Leaf inclination
REAL                                      :: LIDFb      ! Variation in leaf inclination
REAL                                      :: leafwidth  ! (m) Leaf width
REAL                                      :: q          ! %   hot spot parameter (width leaf/height vegetation)

! 7. Aerodynamic
REAL                                      :: z_0        ! (m) Roughness length for momentum of the canopy
REAL                                      :: d          ! (m) Displacement height
REAL                                      :: Cd         ! Leaf drag coefficient
REAL                                      :: rb         ! Leaf boundary resistance
REAL                                      :: CR         ! Verhoef et al. (1997) drag coefficient for isolated tree
REAL                                      :: CD1        ! Verhoef et al. (1997) fitting parameter
REAL                                      :: Psicor     ! Verhoef et al. (1997) roughness layer correction
REAL                                      :: CSSOIL     ! Verhoef et al. (1997) drag coefficient for soil
REAL                                      :: rbs        ! (s m-1) Soil boundary layer resistance
REAL                                      :: rwc        ! (s m-1) Within canopy layer resitance
REAL                                      :: PSIs       ! %[J kg-1] soil water potential
INTEGER                                   :: SoilHeatMethod  ! Flag for the Method to compute the soil heat 
REAL, DIMENSION(12,2)                     :: Tsold0

! 8. Angles
REAL                                      :: tts        ! (deg) Solar zenith angle
REAL                                      :: tto        ! (deg) Observation zenith angle
REAL                                      :: psi        ! (deg) Azimuthal difference between solar and observations angles

!
!%9. Energy balance
! 9.1 Iteration stops
INTEGER                                   :: maxit       ! %        maximum number of iterations
REAL                                      :: maxEBer     ! %[W m-2] maximum accepted error in energy bal.

!%9.2 Weigting factors
REAL                                      :: Wc_Tc        ! %   Weight coefficient for iterative calculation of Tc


! SATELLITE frequency used for fluorescence data retrieval 
REAL                                      :: freq_sat     ! in nm 
INTEGER                                   :: ifreq_sat    ! indice of this frequency in the wave length wlF 

REAL, ALLOCATABLE, DIMENSION(:)           :: list_freq_sat 
INTEGER, ALLOCATABLE, DIMENSION(:)        :: ilist_freq_sat 

INTEGER                                   :: nfreq       ! number of frequency
REAL                                      :: dfreq       ! step of frequency
REAL                                      :: freq0     ! Start of the calcualtion


CHARACTER(len=80)                        :: line

! 10. Files and related fields 
! 10.1 Leaf file
CHARACTER(len=80)                        :: leaf_file = '../input/fluspect_parameters/optipar_fluspect.txt'  !
INTEGER                                  :: nopti ! lines in the optipar file with
REAL,  ALLOCATABLE,DIMENSION(:,:)        :: opticoef


! 10.2  Soil file
CHARACTER(len=80)                        :: soil_file = '../input/soil_spectrum/soilnew.txt' ! soil file file
INTEGER                                  :: nsoil 
REAL, ALLOCATABLE, DIMENSION(:,:)        :: rsfile

! 10.3 MODTRAN output spectrum (used as input for SCOPE)
CHARACTER(len=80)                        :: path_atmos_file = '../input/radiationdata/'
CHARACTER(len=80)                        :: atmos_file
CHARACTER(len=80)                        :: modtran_trop   ! The atmosphere files, sum for summer, wint for winter and trop for tropics
CHARACTER(len=80)                        :: modtran_wint   ! The atmosphere files, sum for summer, wint for winter and trop for tropics
CHARACTER(len=80)                        :: modtran_sum   ! The atmosphere files, sum for summer, wint for winter and trop for tropics

REAL, ALLOCATABLE,DIMENSION(:,:)         :: atmoM


! 10.4 observation angles in case of BRDF calculation
CHARACTER(len=80)                        :: brdf_file = '../input/directional/brdf_angles2.dat'
INTEGER                                  :: nangles      ! lines brdf_file file
REAL, ALLOCATABLE, DIMENSION(:,:)        :: angles





! ------------------------------------------------------------------------------------------------


CONTAINS 

SUBROUTINE initparam 

IMPLICIT NONE 

! ----------------------------------------------------------------
! Declarations of the parameters 

!%4.2 Spectrum
!REAL                               :: VISresolution    ! %[um]  Resolution of spectrum in VISible
!REAL                               :: NIRresolution    ! %[um]  Resolution of spectrum in Near InfraRed
!REAL                               :: MWIRresolution   ! %[um]  Resolution of spectrum in Mid Wave InfraRed
!REAL                               :: TIRresolution    ! %[um]  Resolution of spectrum in Thermal InfraRed

LOGICAL,DIMENSION(:),ALLOCATABLE     :: mask
!INTEGER, ALLOCATABLE, DIMENSION(:)   :: wle,wlf

INTEGER                              :: i,j
INTEGER, PARAMETER                   :: inunit = 2

INTEGER                              :: ireg 

! ---------------------------------------------------------------------
! Now we assign the default values of these parameters 

!1. PROSPECT 
!%1.1 green leaves
Cab          = 80.00                       ! %[ug/cm2] Chorophyl AB content (~60)
Cdm          =  0.012                      ! %[g/cm2]  Dry matter content (~0.012;0.07)
Cw           =  0.009                      ! %[cm]     Leaf water equivalent laye r(~0.009)
Csm          =  0.0                        ! %         Senescent material fraction (~0.0)
N            =  1.4                        ! %         Leaf thickness parameter (~1.4;2)

!%1.2. brown leaves (currently not used)
Cab_b        = 0.000                       ! %[ug/cm2]  Concentration chlorophyl (~60)
Cdm_b        =  0.01                       ! %[g/cm2]   Concentration dry material (~0.005)
Cw_b         =  0.000                      ! %[cm]      Concentration water (~0.009)
Cs_b         =  0.7                        ! %          Concentration senescent material (~0.0)
N_b          =  2.0                        ! %          thickness parameter (~2.0)
f_b          = .0                       ! %  fraction of brown material

!% 1.3 Spectral parameter
rho_thermal  = 0.01                        ! % Broadband thermal reflectance 
tau_thermal  = 0.01                        ! % Broadband thermal translisttance

! 2. Leaf biochemical 
Vcmo         = 60                          ! (um m-2 s) Maximum carboxylation capacity (at optimum temperature)  
Jmo          = 2.*Vcmo                     ! Jmax : maximum electron transport rate table  (not used in this version of scope) 
m            = 8                           ! Ball-Berry stomatal conductance parameter 
option       = 0                           ! Photoochemicl pathway : 0=C3, 1=C4  
kV           = 0                           ! Extinction coefficient for Vcmax in the vertical (maximum at the
                                           ! top). 0 for uniform Vcmax 
Rdparam      = 0.015                       ! Respiration  = Rdparam * Vcmax 
Tparams      = [0.2, 0.3,281., 308., 328.] ! See PFT.xls file. These are five parameters specifying the temperature response 
                                           ! TO CLARIFY THIS .....

tempcor      = 0                           ! [] boolean (0 or 1) whether or not temperature correction to Vcmax has to be
                                           ! applied.
stressfactor = 1.                          !  []   optional input: stress factor to reduce Vcmax (for example
                                           !  soil moisture, leaf age). Default value =  1.
Kcopt        = 350                         !% [ubar] kinetic coefficient for CO2 (Von Caemmerer and Furbank, 1999)
Koopt        = 450                         !% [mbar] kinetic coeeficient for  O2 (Von Caemmerer and Furbank, 1999)
Kf0          = 0.05                        !% [] rate constant for fluorescence
Kd0          = 0.95                        !% [] rate constant for thermal deactivation at Fm
Kpc0         = 4.0                         !% [] rate constant for photochemisty
atheta       = 0.8
lam          = 750.                        ![]    Cowan's stomatal parameter (not used in this version of SCOPE)

model_choice = 0                       ! Fluorescence model_choice integer with
                                       !  0: with sustained quenching after
                                       !  drought, as in Lee et al. (2013)
                                       !  1: calibrated for cotton data set: no
                                       !  drought
!Leaf biochemical (Magnani model)
 Tyear    =   15.                      ! % [oC] mean annual temperature
 beta     = 0.507                      ! % [] fraction of photons partitioned to PSII !(0.507 for C3, 0.4 for C4; Yin et al. 2006; Yin and Struik 2012)
 NPQs     = 0.                         ! !% NPQs   % [s-1] rate constant of sustained thermal !dissipation, normalized to (kf+kD) (=kNPQs'; Porcar-Castell 2011)
 qLs      = 1.                         !  !% qLs  % [] fraction of functional reaction centres !(Porcar-Castell 2011)

                                           ! TO CLARIFY THIS .....
!3. Fluorescence 
fqe2          = 0.02                     ! Fluorescence quantum yield efficiency at photosystem level  
fqe1          = fqe2/5.                  ! Fluorescence quantum yield efficiency at photosystem level  
prat         = 1.                        ! %  PSI/PSII peak ratio

! Here we can choose the range of frequency we want to process  --> 640-850  
!freq_Fsmin   = 755.                      ! Minimum frequeency to compute Fs
!freq_Fsmax   = 755.                      ! Maximum frequency to compute Fs 

freq_Fsmin   = 755.                      ! Minimum frequeency to compute Fs
freq_Fsmax   = 757.                      ! Maximum frequency to compute Fs


!4. Soil 
spectrum     = 1                           ! Spectrum number (column in the data base soil_file)
rss          = 500.                        ! (s m-1) Soil resistane for evaporation from the pore space 
rs_thermal   = 0.06                        ! Broadband soil reflectance in the thermal range (1-emissivity)
cs           = 1.18E+03                    ! (J m-2 K-1) Volumetric heat capacity of the soil
rhos         = 1.80E+3                     ! (kg m-3) Specific mass of the soil  
lambdas      = 1.55                        ! (J m-1 K-1) Heat conductivity of the soil 
SMC          = 0.25                        ! Volumetric soil moisture content in the root zone 

!5. Meteo (values in data files. Thes values are per default, otherwise files are provided)
z            = 10.                         ! (m) measurement height of meteorological data  
Rin          = 600.                        ! (W m-2) brodband incoming shortwave radiation (0.4 - 2.5 um)
Ta           = 20.                         ! (oC) Air temperature 
Rli          = 350.                        ! (W m-2) broadland incoming longwave radiation (2.5 - 50 um) 
pa           = 970.                        ! (hPa) Air pressure 
ea           = 15.                         ! (hPa) Atmospheric vapour pressure 
u            = 2.                          ! (m s-1) Wind speed at height z 
Ca           = 380.                        ! (ppm) Atmospheric concentration
Oa           = 209.                        ! (ppm) Atmospheric O2 concentration

!6. Canopy 
LAI          = 2.                          ! (m2 m-2) Leaf area index  
hc           = 1.                          ! (m)Leaf vegetation height 
LIDFa        = -0.35                       ! Leaf inclination 
LIDFb        = -0.15                       ! Variation in leaf inclination 
leafwidth    = 0.1                         ! (m) Leaf width 
q            =  leafwidth/hc               ! %   hot spot parameter (width leaf/ height vegetation)

! 7. Aerodynamic 
zo           = 0.123                       ! (m) Roughness length for momentum of the canopy 
d            = 0.67                        ! (m) Displacement height 
Cd           = 0.30                        ! Leaf drag coefficient 
rb           = 10.                         ! Leaf boundary resistance 
CR           = 0.35                        ! Verhoef et al. (1997) drag coefficient for isolated tree 
CD1          = 20.6                        ! Verhoef et al. (1997) fitting parameter  
Psicor       = 0.2                         ! Verhoef et al. (1997) roughness layer correction  
CSSOIL       = 0.01                        ! Verhoef et al. (1997) drag coefficient for soil 
rbs          = 10.                         ! (s m-1) Soil boundary layer resistance 
rwc          = 0.                          ! (s m-1) Within canopy layer resitance 
PSIs         = 0.                          ! %[J kg-1]      soil water potential
SoilHeatMethod =1

! 8. Angles 
tts          = 30.                         ! (deg) Solar zenith angle  
tto          = 0.                          ! (deg) Observation zenith angle 
psi          = 90.                          ! (deg) Azimuthal difference between solar and observations angles 

! 
!%9. Energy balance 
! 9.1 Iteration stops
maxit           = 200                      ! %        maximum number of iterations
maxEBer         = 1.                       ! %[W m-2] maximum accepted error in energy bal.

!%9.2 Weigting factors
Wc_Tc           = 0.35                     ! %   Weight coefficient for iterative calculation of Tc



!%% 10. Numerical parameters
nl = 60           ! number of leaf layers


!%% 10. Other Energy balance, micrometeorological and physiological parameters
!% 10.1 vegetation 
Vcmax_ref        = 75.                     ! %[umol/m2/s]        maximum carboxylation capacity
gbs             = 3E3                     ! %[umol/m2/s]        bundle sheath conductance to CO2(C4 only)
Jmax_ref        = 2.5*Vcmax_ref           ! %[umol/m2/s]        maximum electron transport rate (~100;120;270)
gcmin           = 1E-6                    ! %[m s-1]            minimum stomatal conductance (if stomata are closed). 
                                          ! %                   ... this parameter is used to avoid instability
                                                                    
! File MODTRAN for default  standard atmosphere
atmos_file = trim(path_atmos_file)//'FLEX-S3_std.atm'

! File for tropics
modtran_trop = trim(path_atmos_file)//'FLEX-S3_Trop.atm'

! Winter
modtran_wint = trim(path_atmos_file)//'FLEX-S3_Wint.atm'

! Summer
modtran_sum = trim(path_atmos_file)//'FLEX-S3_Mar.atm'


! Here select SIF frequency (nm) to write to output (e.g. SIF retrieval wavelength)
freq_sat =  757.                      

! To compute the sensibility of  Fs to the frequency 
nfreq=31       ! number of frequency 
dfreq =5       ! step of frequency 
freq0=650.     ! freq from which we save the fluorescence

! ============================================================================
!%% Appendix A. Fixed parameters (not to be changed)
!% leaf angle distribution classes (as in SAIL)
! Creation of litab table 
!litab           = ([linspace(5,75,8) linspace(81,89,5)])'  ! %        leaf inclination table
nli  = 13  
IF (.NOT. ALLOCATED(litab)) ALLOCATE(litab(nli))
j=1
do i=5,75,10 
 litab(j)=i 
 j=j+1
end do 
do i=81,89,2 
 litab(j)=i 
 j=j+1
end do 

! Creation of lazitab table 
!lazitab         = linspace(180/36,360-180/36,36);%         !         leaf azimuth table
nlazi =  36 
IF (.NOT. ALLOCATED(lazitab)) ALLOCATE(lazitab(nlazi))
j=1
do i=5,355,10 
 lazitab(j)=i
 j=j+1
end do 

! ==================================================================================
! Define bands 
CALL define_bands 

!print*, ' size wlE ', size(wlE), ' size wlF ', size(wlF)
! ---------------

!!! Read some input files relevant for some parameters/data  .....
!Read the MODTRAN rad and then re_arrange for the model  
! We need to open nregion files if necessary **** here 25/03/2013 E Koffi
!DO ireg = 1, nregion  

!print*,' before aggreg' 
!print*, 'spectral_nreg ', spectral_nreg
!print*, 'spectral_start ', spectral_start
!print*,'spectral_end ', spectral_end
!print*, ' spectral_res ', spectral_res

CALL aggreg (atmos_file,spectral_nreg,spectral_start,spectral_end,spectral_res)
!print*,minval(atmoM), maxval(atmoM), sum(atmoM)

!END DO 

! ========   The leaf file ===========================================
! Read leaf file 
OPEN(unit=inunit,file=leaf_file,status='old')
REWIND inunit
 nopti = 0
    DO WHILE (.TRUE.)
       READ(inunit,*,END=1000) line ! just skip data
       nopti =  nopti +1
    END DO
1000 CONTINUE
CLOSE(inunit)
!print*, ' nopti ', nopti

IF (.NOT.ALLOCATED(opticoef)) ALLOCATE(opticoef(8,nopti))

OPEN(unit=inunit,file=leaf_file,status='old')
REWIND inunit
READ(inunit,*)opticoef
CLOSE(inunit)
! ----------------------------------------------------------------------

!print*, ' opticoef 1 ', opticoef(1,:)

!  =========  The soil file ==============
OPEN(unit=inunit,file=soil_file,status='old')
REWIND inunit
 nsoil = 0
    DO WHILE (.TRUE.)
       READ(inunit,*,END=1001) line ! just skip data
       nsoil =  nsoil +1
    END DO
1001 CONTINUE
CLOSE(inunit)
!print*, ' nsoil ', nsoil

IF (.NOT.ALLOCATED(rsfile)) ALLOCATE(rsfile(nsoil,4))


OPEN(unit=inunit,file=soil_file,status='old')
REWIND inunit
DO i=1, nsoil
READ(inunit,*)rsfile(i,:)
END DO
CLOSE(inunit)

!print*, ' rsfile 1 ', rsfile(:,2)
!--------------------------------------------------------------------------

! ===========  BRDF file  =======================================
!angles          = load([path_input,'directional/',brf_file]); %     Multiple
!observation angles in case of BRDF calculation

OPEN(unit=inunit,file = brdf_file,status='old')
REWIND inunit
nangles = 0
    DO WHILE (.TRUE.)
       READ(inunit,*,END=1002) line ! just skip data
       nangles =  nangles +1
    END DO
1002 CONTINUE
CLOSE(inunit)
!print*, ' nangles ', nangles

IF (.NOT.ALLOCATED(angles)) ALLOCATE(angles(nangles,2))

OPEN(unit=inunit,file = brdf_file,status='old')
REWIND inunit
DO i=1, nangles
READ(inunit,*)angles(i,:)
END DO
CLOSE(inunit)
!print*, ' angles 1, ', angles(:,1)

! -----------------------------------------------------------------------

! Call fluspect 
!print*, ' nwlS ', nwlS 
!print*, ' nwlP ', nwlP 
!print*, ' nopti ', nopti 
nwl = nwlS 
IF (.NOT.ALLOCATED(wl)) ALLOCATE(wl(nwlS))
wl =  wlS

IF (.NOT.ALLOCATED(tau)) ALLOCATE(rho(nwl))
IF (.NOT.ALLOCATED(tau)) ALLOCATE(tau(nwl))
IF (.NOT.ALLOCATED(rs)) ALLOCATE(rs(nwl))

IF (.NOT.ALLOCATED(tran)) ALLOCATE(tran(nwlP))
IF (.NOT.ALLOCATED(refl)) ALLOCATE(refl(nwlP))
IF (.NOT.ALLOCATED(kChlrel)) ALLOCATE(kChlrel(nwlP))


leafbio(1)  = Cab
leafbio(2)  = Cdm
leafbio(3)  = Cw
leafbio(4)  = Csm
leafbio(5)  = N
leafbio(6)  = fqe1                           
leafbio(7)  = fqe2                            
leafbio(8)  = rho_thermal
leafbio(9)  = tau_thermal

CALL  fluspect(leafbio)

!print*, 'size wlT ', size(wlT)  
!print*, 'size wlP ', size(wlP)  
!print*, 'size wlS ', size(wlS)  

!IF (.NOT. ALLOCATED(mask)) ALLOCATE(mask(size(wlS)))
!mask = ((wlS>=minval(wlT)).and.(wlS<=maxval(wlT)))
!IF (.NOT.ALLOCATED(IwlT)) ALLOCATE(IwlT(count(mask)))
!CALL mask_ind(size(mask),mask,IwlT)


!IF (.NOT. ALLOCATED(mask)) ALLOCATE(mask(size(wlS)))
!mask = ((wlS>=minval(wlP)).and.(wlS<=maxval(wlP)))
!IF (.NOT.ALLOCATED(IwlP)) ALLOCATE(IwlP(count(mask)))
!CALL mask_ind(size(mask),mask,IwlP)

!print*, 'min max IwlP ', minval(IwlP),maxval(IwlP)  
!print*, 'min max IwlT ', minval(IwlT),maxval(IwlT)  
!print*, 'size rsfile ', size(rsfile(:,2)), 'IwlP ', size(IwlP)

!rho(IwlP) = refl  
!tau(IwlP) = tran  
! rs(IwlP) = rsfile(:,2)  

!rho(IwlT) = rho_thermal
!tau(IwlT) = tau_thermal
!rs(IwlT)  = rs_thermal

!print*, ' rho(iwlfo) ', minval(rho(iwlfo)), maxval(rho(iwlfo)),sum(rho(iwlfo))
!print*, ' tau(iwlfo) ', minval(tau(iwlfo)), maxval(tau(iwlfo)),sum(tau(iwlfo))


!% B.1 Output directory

!% B.2 number of classes, leaf angle distribution, observation angles
IF (.NOT.ALLOCATED(lidf)) ALLOCATE(lidf(nli))
CALL leafangles(LIDFa,LIDFb)

! Define layers
!IF (.NOT.ALLOCATED(xlay)) ALLOCATE(xlay(nl))
!CALL layers
!dx     = abs(xlay(2)-xlay(1))                          !  %

tto_   = angles(:,1)      ! %[deg]  Observation zenith Angles for calcbrdf
psi_   = angles(:,2)      !  %[deg]  Observation zenith Angles for calcbrdf
noa    = size(tto_,1)     !  %       Number of Observation Angles 

!print*, ' noa ', noa 


!% B.4 Soil thermal inertia (may be not to be used)
![GAM]           = Soil_Intertia(cs,rhos,lambdas)
CALL  Soil_Intertia(cs,rhos,lambdas,GAM)

!print*, ' GAM ', GAM 

!% B.5 Inversion of planck by Lookup Table  (not used in this version)
! We do not consider this for the moment  
![T_LUT]    =   Temperature_Lookup_Table;

! Define layers 
CALL Nlayers(LAI)
CALL layers

!print*, ' Nlayers ', nl 
! ALLOCATE fields 
CALL field_allocate       
!print*, ' apres field_allocate'
CALL fieldlayer_allocate 

END SUBROUTINE initparam 

SUBROUTINE field_allocate
IMPLICIT NONE

! Output variables from rtmo
ALLOCATE(Eout_(nwl))
ALLOCATE(Lout_(nwl))
ALLOCATE(Lo_(nwl))

! Output from rtmf 
ALLOCATE(LoF(nwlfo))
ALLOCATE(Fhem(nwlfo))
END SUBROUTINE field_allocate


SUBROUTINE field_deallocate
IMPLICIT NONE
DEALLOCATE(Lout_,Eout_,Lo_,LoF,Fhem)
DEALLOCATE(list_freq_sat)
DEALLOCATE(ilist_freq_sat)
END SUBROUTINE field_deallocate

SUBROUTINE fieldlayer_allocate
IMPLICIT NONE

! Output from rtmo 
ALLOCATE(Rnuc(nli,nlazi,nl))
ALLOCATE(Pnuc(nli,nlazi,nl))
ALLOCATE(Pnuc_Cab(nli,nlazi,nl))
ALLOCATE(Rnhc(nl),Pnhc(nl),Pnhc_Cab(nl))
ALLOCATE(Eplu_(nl+1,nwl))
ALLOCATE(Emin_(nl+1,nwl))
!print*, ' 0'
ALLOCATE(Ps(nl+1))
ALLOCATE(Po(nl+1))
ALLOCATE(Pso(nl+1))
!print*, ' 1 '
ALLOCATE(Esun_(nwl))
ALLOCATE(Esky_(nwl))
ALLOCATE(fEsuno(nwl))
ALLOCATE(fEskyo(nwl))
ALLOCATE(fEsunt(nwl))
ALLOCATE(fEskyt(nwl))

 ! print*, ' 2'
! Output variables from rtmt 
ALLOCATE(Eplu(nl+1))
ALLOCATE(Emin(nl+1))
ALLOCATE(Rnhct(nl))
ALLOCATE(Rnuct(nli,nlazi,nl))

 ! print*, ' 3'
! Output from rtmf 
ALLOCATE(Fiprof(nl+1))

  !print*, ' 4'
! Other
ALLOCATE(Rnch(nl))
ALLOCATE(Rncu(nli,nlazi,nl))

ALLOCATE(Fc(nl))
END SUBROUTINE fieldlayer_allocate


SUBROUTINE fieldlayer_deallocate
IMPLICIT NONE

DEALLOCATE(Rnuc,Pnuc,Pnuc_Cab)
DEALLOCATE(Rnhc,Pnhc,Pnhc_Cab)
DEALLOCATE(Eplu_,Emin_)
DEALLOCATE(Esun_,Esky_)
DEALLOCATE(fEsuno,fEskyo)
DEALLOCATE(fEsunt,fEskyt)
DEALLOCATE(Ps,Po,Pso)
DEALLOCATE(xlay,Fc)

DEALLOCATE(Eplu,Emin,Rnhct,Rnuct)
DEALLOCATE(Rnch,Rncu)

DEALLOCATE(Fiprof)

END SUBROUTINE fieldlayer_deallocate


!function [spectral] = define_bands
SUBROUTINE define_bands 

!USE fluo_func, ONLY : mask_ind 

IMPLICIT NONE 

! Output variables 
!REAL, ALLOCATABLE,DIMENSION(:)                          :: wlS,wlP,wlF,wlT
!REAL, ALLOCATABLE,DIMENSION(:)                          :: wlE,wlO, wlPAR 

! Local variables 
REAL, ALLOCATABLE,DIMENSION(:)                          :: reg1,reg2,reg3 
!INTEGER, DIMENSION(3)                                   :: spectral_start,spectral_end,spectral_res 
INTEGER                                                 :: n1,n2,n3, ntot,k 
INTEGER                                                 :: nE, nF, nT 
!INTEGER                                                 :: spectral_nreg 
LOGICAL, ALLOCATABLE, DIMENSION(:)                      :: mask
INTEGER, ALLOCATABLE,DIMENSION(:)                       :: j
REAL, ALLOCATABLE,DIMENSION(:)                          :: vi 

REAL                                                    :: minwlS, maxwlS
REAL                                                    :: minwlE, maxwlE
REAL                                                    :: minwlF, maxwlF

INTEGER                                                 :: i

INTEGER                                                 :: kf


!    % Define spectral regions for SCOPE v_1.40
!    % All spectral regions are defined here as row vectors
!    % WV Jan. 2013
    
!    % 3 spectral regions for SCOPE
    
   ! reg1 =   400 :    1 :  2400;
    n1 = int((2400 -400)/1)+1
    ntot = 1
 IF (.NOT.ALLOCATED(reg1)) ALLOCATE(reg1(n1))
  do k=1,n1
    reg1(k) = 400.+(k-1)*1.
     ntot   = ntot +1
  end do
!  print*, ' n1 ', n1, ' ntot ', ntot 


!    reg2 =  2500 :  100 : 15000;
    n2 = int((15000 -2500)/100)+1
 IF (.NOT.ALLOCATED(reg2)) ALLOCATE(reg2(n2))
   do k=1,n2
   reg2(k) = 2500.+(k-1)*100.
      ntot = ntot +1
  end do
!  print*, ' n2 ', n2, ' ntot ', ntot 
 

  !reg3 = 16000 : 1000 : 50000;
  n3 = int((50000 -16000)/1000)+1
 IF (.NOT.ALLOCATED(reg3)) ALLOCATE(reg3(n3))
   do k=1,n3
   reg3(k) = 16000.+(k-1)*1000.
      ntot = ntot +1
  end do
!  print*, ' n3 ', n3, ' ntot ', ntot

    
!    spectral.wlS  = [reg1 reg2 reg3];
IF (.NOT.ALLOCATED(wlS)) ALLOCATE(wlS(ntot-1))
IF (.NOT.ALLOCATED(wl)) ALLOCATE(wl(ntot-1))
    wlS = [reg1, reg2, reg3]
     wl =  wlS     
!print*, ' min max wlS ', minval(wlS), maxval(wlS), sum(wlS)
!print*, ' wlS ', wlS 

!    % Other spectral (sub)regions
IF (.NOT.ALLOCATED(wlP)) ALLOCATE(wlP(n1))
   ! spectral.wlP   = reg1;                            % PROSPECT data range
    wlP    =  reg1

!print*, ' nwlP ',size(wlP)

!print*, ' min max wlP ', minval(wlP), maxval(wlP),sum(wlP) 
!print*, ' wlP ', wlP 

 !   spectral.wlE   = 400:1:750;                       % excitation in E-F matrix
  nE = int((750 -400)/1)+1
  ntot =  1
 IF (.NOT.ALLOCATED(wlE)) ALLOCATE(wlE(nE))
   do k=1,nE
   wlE(k) = 400.+(k-1)*1.
      ntot = ntot +1
  end do
!  print*, ' nE ', nE, ' ntot ', ntot

!print*, ' min max wlE ', minval(wlE), maxval(wlE),sum(wlE) 
!print*, ' wlE ', wlE 

!    spectral.wlF   = 640:1:850;                       % chlorophyll fluorescence in E-F matrix
  !nF = int((850 -640)/1)+1
  nF = int((freq_Fsmax -freq_Fsmin)/1)+1
  ntot =  1
 IF (.NOT.ALLOCATED(wlF)) ALLOCATE(wlF(nF))



IF (.NOT. ALLOCATED(list_freq_sat)) ALLOCATE(list_freq_sat(nfreq))
IF (.NOT. ALLOCATED(ilist_freq_sat)) ALLOCATE(ilist_freq_sat(nfreq))

  do kf=1,nfreq 
   list_freq_sat(kf) = freq0+(kf-1)*dfreq 

   do k=1,nF
   wlF(k) = freq_Fsmin+(k-1)*1.
      ntot = ntot +1
      
! Find the indice for the frequency freq_gosat (here 755nm) used for fluo data retrieval
      IF (wlF(k) .EQ. freq_sat) ifreq_sat = k 
      IF (wlF(k) .EQ.  list_freq_sat(kf)) ilist_freq_sat(kf) = k 
     
    end do

 end do 





!print*, ' nF ', nF, ' ntot ', ntot, ' notot-1 ', ntot-1
!print*, ' n2+n3 ', n2+n3 , ' n1+n2+n3 ', n1+n2+n3

!print*, ' min max wlF ', minval(wlF), maxval(wlF),sum(wlF) 
!print*, ' wlF ', wlF 

!    spectral.wlO   = reg1;                            % optical part
    IF (.NOT.ALLOCATED(wlO)) ALLOCATE(wlO(n1))
        wlO   =  reg1

!print*, ' min max wlO ', minval(wlO), maxval(wlO),sum(wlO) 
!print*, ' wlO ', wlO 

!    spectral.wlT   = [reg2 reg3];                     % thermal part
  IF (.NOT.ALLOCATED(wlT)) ALLOCATE(wlT(n2+n3))
       wlT   = [reg2, reg3]                    
    
!print*, ' nwT ',size(wlT) 
!print*, ' min max wlT ', minval(wlT), maxval(wlT),sum(wlT) 
!print*, ' wlT ', wlT 
! Already assigned above 
   ! wlS            = spectral.wlS;

!    spectral.wlPAR = wlS(wlS>=400 & wlS<=700);  % PAR range
IF (.NOT. ALLOCATED(mask)) ALLOCATE(mask(size(wlS)))
mask =(( wlS>=400).AND.( wlS<=700))
IF (.NOT.ALLOCATED(j)) ALLOCATE(j(count(mask)))
CALL mask_ind(size(mask),mask,j)
!print*, ' size(j )', size(j)
IF (.NOT.ALLOCATED(wlPAR)) ALLOCATE(wlPAR(size(j)))
 wlPAR  = wlS(j)



!print*, 'nwlP ', size(wlP)
!print*, 'nwlT ', size(wlT)
!print*, 'nwlPAR ', size(wlPAR)
!open(3,file='bands.dat',status='unknown')
!write(3,'(f18.10)')wlP
!write(3,'(f18.10)')wlT
!write(3,'(f18.10)')wlPAR
!close(3)

!print*, ' min max wlPAR ', minval(wlPAR), maxval(wlPAR),sum(wlPAR) 
!print*, ' wlPAR ', wlPAR 
    
  !  % Data used by aggreg routine to read in MODTRAN data
 !   spectral.SCOPEspec.nreg = 3;
    spectral_nreg = 3

!    spectral.SCOPEspec.start = [ 400  2500  16000];
   spectral_start = [ 400,  2500,  16000]

!    spectral.SCOPEspec.end   = [2400 15000  50000];
    spectral_end   = [2400, 15000,  50000]

 !   spectral.SCOPEspec.res   = [   1   100   1000];
    spectral_res   = [   1,   100,   1000]
!  print*, '  spectral_nreg ',  spectral_nreg
!  print*, '  spectral_start  ',  spectral_start 
!  print*, ' spectral_end  ',   spectral_end 
!  print*, ' spectral_res ',  spectral_res


    nwlP = size(wlP)
    nwlF = size(wlF)
    nwlS = size(wlS)
    nwlO = size(wlO)
    nwlE = size(wlE)
 nwlPAR  = size(wlPAR)
    nwl  = size(wlS)


! Define intersection 
![B,iwlfi]    = intersect(wlS,wlE);
![C,iwlfo]    = intersect(wlS,wlF);

minwlS    = minval(wlS)
maxwlS    = maxval(wlS)

minwlE    = minval(wlE)
maxwlE    = maxval(wlE)

minwlF    = minval(wlF)
maxwlF    = maxval(wlF)


!print*, ' minwlS ', minwlS, ' maxwlS ', maxwlS
!print*, ' minwlE ', minwlE, ' maxwlE ', maxwlE
!print*, ' minwlF ', minwlF, ' maxwlF ', maxwlF


![B,iwlfi]    = intersect(wlS,wlE);
IF (.NOT. ALLOCATED(mask)) ALLOCATE(mask(size(wlS)))
mask = ((wlS>=minwlE).and.(wlS<=maxwlE))
IF (.NOT.ALLOCATED(iwlfi)) ALLOCATE(iwlfi(count(mask)))
CALL mask_ind(size(mask),mask,iwlfi)

nwlfi    = size(iwlfi,1)

!print*, ' iwlfi ', minval(iwlfi), maxval(iwlfi)
!print*, ' wlS(iwlfi) ', minval(wlS(iwlfi)), maxval(wlS(iwlfi))

![C,iwlfo]    = intersect(wlS,wlF);
IF (.NOT. ALLOCATED(mask)) ALLOCATE(mask(size(wlS)))
mask = ((wlS>=minwlF).and.(wlS<=maxwlF))
IF (.NOT.ALLOCATED(iwlfo)) ALLOCATE(iwlfo(count(mask)))
CALL mask_ind(size(mask),mask,iwlfo)

nwlfo    = size(iwlfo,1)

!print*, ' iwlfo ', minval(iwlfo), maxval(iwlfo)
!print*, ' wlS(iwlfo) ', minval(wlS(iwlfo)), maxval(wlS(iwlfo))

DEALLOCATE(mask)

END SUBROUTINE define_bands


!%% Appendix F. Soil parameters function (not to be changed)
!function [GAM]  = Soil_Intertia(cs,rhos,lambdas)
SUBROUTINE  Soil_Intertia(cs,rhos,lambdas,GAM)
IMPLICIT NONE

! Input variables
REAL, INTENT(IN)         :: cs, rhos, lambdas

! Input Output
REAL, INTENT(OUT)        :: GAM

!!% soil thermal intertia
GAM  = sqrt(cs*rhos*lambdas)           ! % soil thermal intertia

END SUBROUTINE Soil_Intertia


!%% Appendix G. Layer parameters function (not to be changed)
!% vector describing the relative depth in the canopy (0: top, -1: bottom)
!function [x]    =  layers()


SUBROUTINE layers

IMPLICIT NONE

! Local variables
INTEGER                          :: i
REAL                             :: x0

!x  = (-1/nl:-1/nl:-1)   !%create layering of canopy (linear in this case)
ALLOCATE(xlay(nl))
x0 = -1./nl
DO i=1,nl
xlay(i) = i*x0
END DO

dx     = abs(xlay(2)-xlay(1))

END SUBROUTINE layers


!function [refl,tran,Mb,Mf] = fluspect(leafpar, optipar)
SUBROUTINE fluspect(leafpar)

!% calculates reflection and transmission of a leaf using FLUSPECT
!%
!% Authors: Wout Verhoef, Christiaan van der Tol (tol@itc.nl), Joris Timmermans, 
!% Date: 2007
!% Update from PROSPECT to FLUSPECT: January 2011 (CvdT)
!%
!% usage:
!% [refl, tran] = fluspect(leafpar, optipar, wle, wlf)
!% 
!% input:
!% leafpar: a vector with 5 elements:
!% Cab         = leafpar(1);
!% Cw          = leafpar(2);
!% Cdm         = leafpar(3);
!% Cs          = leafpar(4);
!% N           = leafpar(5);
!% fqe         = leafpar(6); fluorescence quantum yield
!% prat        = leafpar(7); PSI/PSII peak ratio
!%
!% optipar:
!% a matrix with the following columns (rows refer to the wavelength):
!% nr          = optipar(:,1);
!% Kdm         = optipar(:,2);
!% Kab         = optipar(:,3);
!% Kw          = optipar(:,4);
!% Ks          = optipar(:,5);
!% phiI        = optipar(:,6);
!% phiII       = optipar(:,7);
!%
!% wle           index for the excitation wavelengths (corresponding to optipar)
!% wlf           index for the fluorescence wavelengths (corresponding to optipar)
!%
!% output:
!% relf          reflectance
!% tran          transmittance
!% Mb            backward scattering fluorescence matrix
!% Mf            forward scattering fluorescence matrix

IMPLICIT NONE 

! Input variables 
REAL, INTENT(IN), DIMENSION(:)                   :: leafpar 
!REAL, INTENT(IN), DIMENSION(:,:)                 :: optipar 

! Output variables 
!REAL, INTENT(OUT),  DIMENSION(nwl)                :: tran,refl,kChrel

! Local variables 
REAL                                                :: Cab,Csm,Cw,Cdm,N,fqe,prat
REAL, DIMENSION(nwlP)                               :: nr,Kdm,Kab,Kw, Ks 
REAL, DIMENSION(nwlP)                               :: phiI, phiII, Kall,t1,t2,taut
REAL, DIMENSION(nwlP)                               :: talf,ralf
REAL, DIMENSION(nwlP)                               :: t12,r12,t21,r21
REAL, DIMENSION(nwlP)                               :: denom,Taf,Ra
REAL, DIMENSION(nwlP)                               :: t,r,D,rq,tq,a,b 
REAL, DIMENSION(nwlP)                               :: bNm1, bN2, a2
REAL, DIMENSION(nwlP)                               :: Tsub,Rsub,s,k,kChl

!REAL, DIMENSION(nwl)                               :: nr,Kdm,Kab,Kw, Ks
!REAL, DIMENSION(nwl)                               :: phiI, phiII, Kall,t1,t2,taut
!REAL, DIMENSION(nwl)                               :: talf,ralf
!REAL, DIMENSION(nwl)                               :: t12,r12,t21,r21
!REAL, DIMENSION(nwl)                               :: denom,Taf,Ra
!REAL, DIMENSION(nwl)                               :: t,r,D,rq,tq,a,b
!REAL, DIMENSION(nwl)                               :: bNm1, bN2, a2
!REAL, DIMENSION(nwl)                               :: Tsub,Rsub,s,k,kChl



!REAL, DIMENSION(nopti)                               :: nr,Kdm,Kab,Kw, Ks
!REAL, DIMENSION(nopti)                               :: phiI, phiII, Kall,t1,t2,taut
!REAL, DIMENSION(nopti)                               :: talf,ralf
!REAL, DIMENSION(nopti)                               :: t12,r12,t21,r21
!REAL, DIMENSION(nopti)                               :: denom,Taf,Ra
!REAL, DIMENSION(nopti)                               :: t,r,D,rq,tq,a,b
!REAL, DIMENSION(nopti)                               :: bNm1, bN2, a2
!REAL, DIMENSION(nopti)                               :: Tsub,Rsub,s,k,kChl



REAL, ALLOCATABLE,DIMENSION(:)                     :: te,tf,re,rf
REAL, ALLOCATABLE,DIMENSION(:)                     :: xe,xf,ten,tfn
REAL, ALLOCATABLE,DIMENSION(:)                     :: ren,rfn
REAL, ALLOCATABLE,DIMENSION(:,:)                   :: sigmoid,siglelf
REAL, ALLOCATABLE,DIMENSION(:,:)                   :: MfnI,MbnI
REAL, ALLOCATABLE,DIMENSION(:,:)                   :: MfnII,MbnII

REAL                                               :: wII,wI, eps 

LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: mask
INTEGER, ALLOCATABLE, DIMENSION(:)                 :: j,jj
INTEGER, ALLOCATABLE, DIMENSION(:)                 :: I_rt, I_a, I_an
INTEGER, ALLOCATABLE, DIMENSION(:)                 :: Iwle, Iwlf 
REAL, ALLOCATABLE,DIMENSION(:,:)                   :: Ih,Iv,we,wf
REAL, ALLOCATABLE,DIMENSION(:,:)                   :: wxe,wxf,wre,wrf

REAL                                               :: minwle, maxwle 
REAL                                               :: minwlf, maxwlf 

INTEGER                                            :: ndub,i
INTEGER                                            :: l,lm

REAL, DIMENSION(nwl)                               :: intexp 

!%% parameters
!% fixed parameters for the fluorescence module
ndub        = 15

!print*, 'leafpar ', leafpar 


!% Fluspect parameters
!Bio
Cab         = leafpar(1)
Cdm         = leafpar(2)
Cw          = leafpar(3)
Csm         = leafpar(4)
N           = leafpar(5)
fqe1        = leafpar(6)
fqe2        = leafpar(7)
rho_thermal = leafpar(8)
tau_thermal = leafpar(9)


!print*, ' Cab ', Cab 
!print*, ' Cdm ', Cdm 
!print*, ' Cw ', Cw 
!print*, ' Csm ', Csm 
!print*, ' N ', N 
!print*, ' fqe1 ', fqe1 , ' fqe2 ', fqe2
!print*, 'rho_thermal ', rho_thermal 
!print*, 'tau_thermal ', tau_thermal 

!print*, ' size ', size(nr)
!print*, ' opticoef ',minval(opticoef(2,:)), maxval(opticoef(2,:))

nr          = opticoef(2,:)
Kdm         = opticoef(3,:)
Kab         = opticoef(4,:)
Kw          = opticoef(5,:)
Ks          = opticoef(6,:)
phiI        = opticoef(7,:)
phiII       = opticoef(8,:)

!print*, ' nr ', minval(nr), maxval(nr) 
!print*, ' Kdm ',minval(Kdm), maxval(Kdm) 
!print*, ' Kab ', minval(Kab), maxval(Kab) 
!print*, ' Kw ', minval(Kw), maxval(Kw) 
!print*, ' Ks ', minval(Ks), maxval(Ks) 
!print*, ' phiI ', minval(phiI), maxval(phiI) 
!print*, ' phiII ', minval(phiII),maxval(phiII) 


!print*, ' hdjdjdjd ', size(wle), size(wlf), size(wlp)
!print*, wle 

minwle      = minval(wle)
maxwle      = maxval(wle)
!print*, ' minwle ', minwle 
!print*, ' maxwle ', maxwle 
minwlf      = minval(wlf)
maxwlf      = maxval(wlf)

!print*, minwle, maxwle, minwlf,maxwlf 
!% indices of wle and wlf within wlp
!Iwle        = find(wlp>=minwle & wlp<=maxwle);
IF (.NOT. ALLOCATED(mask)) ALLOCATE(mask(size(wlp)))
!print*, ' apres allocate '
mask = ((wlp>=minwle).and.(wlp<=maxwle))
IF (.NOT.ALLOCATED(Iwle)) ALLOCATE(Iwle(count(mask)))
CALL mask_ind(size(mask),mask,Iwle)

!print*, ' dhdhdhdhdh'

!Iwlf        = find(wlp>=minwlf & wlp<=maxwlf);
IF (.NOT. ALLOCATED(mask)) ALLOCATE(mask(size(wlp)))
mask = ((wlp>=minwlf).and.(wlp<=maxwlf))
IF (.NOT.ALLOCATED(Iwlf)) ALLOCATE(Iwlf(count(mask)))
CALL mask_ind(size(mask),mask,Iwlf)

!print*, ' size Iwle ', size(Iwle), ' size Iwlf ', size(Iwlf)


!print*,'nr ', minval(nr), maxval(nr),sum(nr)
!print*,'Kdm ', minval(Kdm), maxval(Kdm), sum(Kdm)
!print*,'Kab ', minval(Kab), maxval(Kab),sum(Kab)
!print*,'Kw ', minval(Kw), maxval(Kw),sum(Kw)
!print*,'Ks ', minval(Ks), maxval(Ks), sum(Ks)
!print*,'phiI', minval(phiI), maxval(phiI), sum(phiI)
!print*,'phiII', minval(phiII), maxval(phiII), sum(phiII)



!%% Prospect calculations
Kall        = (Cab*Kab+Cdm*Kdm+Cw*Kw+Csm*Ks)/N
!print*, ' in fluspect Kall   ', minval(Kall), maxval(Kall),sum(Kall)

!j           = find(Kall>0)
IF (.NOT. ALLOCATED(mask)) ALLOCATE(mask(size(Kall)))
mask = (Kall>0.)
IF (.NOT.ALLOCATED(j)) ALLOCATE(j(count(mask)))
CALL mask_ind(size(mask),mask,j)

t1          = (1-Kall)*exp(-Kall)
CALL expint(1,size(Kall),Kall,intexp)
t2           = Kall**2.*intexp 

!print*, ' in fluspect Kall   ', minval(Kall), maxval(Kall),sum(Kall)
!print*, ' in fluspect intexp   ', minval(intexp), maxval(intexp),sum(intexp)
!print*, ' in fluspect t1   ', minval(t1), maxval(t1),sum(t1)
!print*, ' in fluspect t2   ', minval(t2), maxval(t2),sum(t2)

taut         = 1. 
taut(j)      = t1(j)+t2(j)

kChlrel     = 0.
kChlrel(j)   = Cab*Kab(j)/(Kall(j)*N)

CALL  calctav(59,nr,talf)
ralf        = 1.-talf

CALL  calctav(90,nr,t12)
r12         = 1.-t12
t21         = t12/(nr**2)
r21         = 1.-t21

!% top layer
denom       = 1.-r21*r21*taut**2
!print*, ' in fluspect r21   ', minval(r21), maxval(r21),sum(r21)
!print*, ' in fluspect taut   ', minval(taut), maxval(taut),sum(taut)

Taf          = talf*taut*t21/denom
Ra          = ralf+r21*taut*Taf

!% deeper layers
t           = t12*taut*t21/denom
r           = r12+r21*taut*t

!print*, ' before t ', minval(t), maxval(t), sum(t)
!print*, ' before r ', minval(r), maxval(r), sum(r)

!% Stokes equations to compute properties of next N-1 layers (N real) Normal case
D           = sqrt((1.+r+t)*(1.+r-t)*(1.-r+t)*(1.-r-t))
rq          = r**2
tq          = t**2
a           = (1.+rq-tq+D)/(2*r)
b           = (1.-rq+tq+D)/(2*t)
!print*, ' D = ', minval(D), maxval(D),sum(D)
!print*, ' rq = ', minval(rq), maxval(rq),sum(rq)
!print*, ' tq = ', minval(tq), maxval(tq), sum(tq)
!print*, ' a = ', minval(a), maxval(a), sum(a)
!print*, ' b = ', minval(b), maxval(b), sum(b)

bNm1        = b**(N-1.)                 ! %** CHECK
bN2         = bNm1**2.
a2          = a**2.
denom       = a2*bN2-1.
Rsub        = a*(bN2-1.)/denom
Tsub        = bNm1*(a2-1.)/denom

!print*, ' a2 ', minval(a2), maxval(a2),sum(a2)
!print*, ' bNm1 ', minval(bNm1), maxval(bNm1), sum(bNm1)
!print*, ' denom ', minval(denom), maxval(denom), sum(denom)
!print*, ' BEFORE Rsub ', minval(Rsub), maxval(Rsub), sum(Rsub)
!print*, ' Tsub ', minval(Tsub), maxval(Tsub), sum(Tsub)



s           = r/t    ! Conservative scattering (CS)

!print*, ' s old ', minval(s), maxval(s), sum(s)


! Version from V1.53   ---> EK Sept 2014
!---------------------------------------------

s(j)        = 2.*a(j)/(a(j)**2.-1)*log(b(j))   ! Normal case overwrites CS case 

!print*, ' s new ', minval(s), maxval(s), sum(s)

!print*, 'r+t ', minval(r+t), maxval(r+t), sum(r+t)


! Case of zero absorption 
IF (.NOT. ALLOCATED(mask)) ALLOCATE(mask(size(r+t)))
mask = ((r+t) >=1.)
IF (.NOT.ALLOCATED(jj)) ALLOCATE(jj(count(mask)))
CALL mask_ind(size(mask),mask,jj)


Tsub(jj)     = t(jj)/(t(jj)+(1.-t(jj))*(N-1.))
Rsub(jj)     = 1.-Tsub(jj)


!print*, 'AFTER  Rsub ', minval(Rsub), maxval(Rsub), sum(Rsub)
!print*, ' Tsub ', minval(Tsub), maxval(Tsub), sum(Tsub)

!print*, 'after t ', minval(t), maxval(t),sum(t)
!print*, 'after r ', minval(r), maxval(r),sum(r)


! Reflectance and transmisttance of the leaf: combine top layer with next N-1
! layers 
denom      =  1.-Rsub*r 
 tran      =  Taf*Tsub/denom  
 refl      =  Ra+Taf*Rsub*t/denom 

 t        =   tran
 r        =   refl

!print*, ' after 2 t ', minval(t), maxval(t),sum(t)
!print*, ' after 2 r ', minval(r), maxval(r),sum(r)
!print*, ' r+t ', minval(r+t), maxval(r+t), sum(r+t)

!I_rt     =   (r+t)<1
IF (.NOT. ALLOCATED(mask)) ALLOCATE(mask(size(r)))
mask = ((r+t)<1.)
IF (.NOT.ALLOCATED(I_rt)) ALLOCATE(I_rt(count(mask)))
CALL mask_ind(size(mask),mask,I_rt)


D(I_rt)  =   sqrt((1 + r(I_rt) + t(I_rt)) * (1 + r(I_rt) - t(I_rt)) * &
                  (1 - r(I_rt) + t(I_rt)) * (1 - r(I_rt) - t(I_rt)))
a        = 1. 
b        = 1. 

a(I_rt)  =   (1 + r(I_rt)**2. - t(I_rt)**2. + D(I_rt)) / (2.*r(I_rt))
b(I_rt)  =   (1 - r(I_rt)**2. + t(I_rt)**2. + D(I_rt)) / (2.*t(I_rt))

!a(~I_rt) =   1.
!b(~I_rt) =   1.

!print*, ' D I_rt ', minval(D(I_rt)), maxval(D(I_rt)), sum(D(I_rt)) 
!print*, ' a ', minval(a(I_rt)), maxval(a(I_rt)),sum(a(I_rt))
!print*, ' b ', minval(b(I_rt)), maxval(b(I_rt)), sum(b(I_rt))


!I_a      =   a>1;
IF (.NOT. ALLOCATED(mask)) ALLOCATE(mask(size(a)))
mask = (a>1.)
IF (.NOT.ALLOCATED(I_a)) ALLOCATE(I_a(count(mask)))
CALL mask_ind(size(mask),mask,I_a)

s(I_a)   =   2.*a(I_a) / (a(I_a)**2. - 1.) * log(b(I_a))

!I_an    = a <=1.
IF (.NOT. ALLOCATED(mask)) ALLOCATE(mask(size(a)))
mask = (a<=1.)
IF (.NOT.ALLOCATED(I_an)) ALLOCATE(I_an(count(mask)))
CALL mask_ind(size(mask),mask,I_an)

!print*, ' I_an ', I_an 

s(I_an)  =   r(I_an) / t(I_an)

k        =   (a-1) / (a+1) * log(b)
kChl      =   kChlrel * k
k(jj)     = 0.
kChl(jj)   = 0.


!print*, 's I_a ', minval(s(I_a)), maxval(s(I_a)), sum(s(I_a))
!print*, 's I_an ', minval(s(I_an)), maxval(s(I_an)), sum(s(I_an))

!print*, 'BEFORE TEST  fqe1 ', fqe1, 'fqe2 ', fqe2


! if fqe >0
if ((fqe1.GT.0) .and. (fqe2.gt.0.))  then 

!print*, ' fqe1 ', fqe1, 'fqe2 ', fqe2

!% Reflectance and transmittance of the leaf: combine top layer with next N-1 layers
!denom = 1.-Rsub*r
!tran  = Taf*Tsub/denom
!refl  = Ra+Taf*Rsub*t/denom

!print*, ' in fluspect Taf ', minval(Taf), maxval(Taf),sum(Taf)
!print*, ' in fluspect r ', minval(r), maxval(r),sum(r)
!print*, ' in fluspect Resub ', minval(Rsub), maxval(Rsub),sum(Rsub)
!print*, ' in fluspect denom ', minval(denom), maxval(denom),sum(denom)
!print*, ' in fluspect refl ', minval(refl), maxval(refl),sum(refl)
!print*, ' in fluspect tran ', minval(tran), maxval(tran), sum(tran)
!print*, ' in fluspect kChlrel ', minval(kChlrel), maxval(kChlrel),sum(kChlrel)

!print*, ' size wle ', size(wle), ' size wlE ', size(wlE)
!print*, ' size wlf ', size(wlf), ' size wlF ', size(wlF)

IF (.NOT.ALLOCATED(MbI))  ALLOCATE(MbI(size(wlf),size(wle)))
IF (.NOT.ALLOCATED(MfI))  ALLOCATE(MfI(size(wlf),size(wle)))
IF (.NOT.ALLOCATED(MbII)) ALLOCATE(MbII(size(wlf),size(wle)))
IF (.NOT.ALLOCATED(MfII)) ALLOCATE(MfII(size(wlf),size(wle)))

print*, ' here 1 '

wII         = prat/(1.+prat)
wI          = 1.-wII
eps         = 2.**(-ndub)

!% initialisations
IF (.NOT.ALLOCATED(te))  ALLOCATE(te(size(wle)))
IF (.NOT.ALLOCATED(re))  ALLOCATE(re(size(wle)))
IF (.NOT.ALLOCATED(xe))  ALLOCATE(xe(size(wle)))
IF (.NOT.ALLOCATED(ten)) ALLOCATE(ten(size(wle)))
IF (.NOT.ALLOCATED(ren)) ALLOCATE(ren(size(wle)))
IF (.NOT.ALLOCATED(Ih))  ALLOCATE(Ih(1,size(wle)))
IF (.NOT.ALLOCATED(we))  ALLOCATE(we(1,size(wle)))

!print*, ' here 11 '

IF (.NOT.ALLOCATED(tf))  ALLOCATE(tf(size(wlf)))
IF (.NOT.ALLOCATED(rf))  ALLOCATE(rf(size(wlf)))
IF (.NOT.ALLOCATED(xf))  ALLOCATE(xf(size(wlf)))
IF (.NOT.ALLOCATED(tfn)) ALLOCATE(tfn(size(wlf)))
IF (.NOT.ALLOCATED(rfn)) ALLOCATE(rfn(size(wlf)))
IF (.NOT.ALLOCATED(Iv))  ALLOCATE(Iv(size(wlf),1))
IF (.NOT.ALLOCATED(wf))  ALLOCATE(wf(size(wlf),1))

!print*, ' here 12 '

!print*, ' s ', size(s) , ' k ', size(k), ' te ', size(te)
!print*, ' eps ', eps
te          = 1.-(k(Iwle)+s(Iwle))*eps
!print*, ' here entre deux 0 '
!print*, ' size s  ', size(s), ' size k ', size(k), size(Iwlf), minval(Iwlf), maxval(Iwlf)
tf          = 1.-(k(Iwlf)+s(Iwlf))*eps
!print*, ' entre deux '
re          = s(Iwle)*eps
rf          = s(Iwlf)*eps

Ih(1,:)   = 1.
Iv(:,1)   = 1.

!print*, ' entre deux 2 '
!we(1,:) = wl(Iwle)
we(1,:) = wle
!wf(:,1) = wl(Iwlf)
wf(:,1) = wlf


!print*, ' here iiiii ' 
IF (.NOT.ALLOCATED(siglelf)) ALLOCATE(siglelf(size(wlf),size(wle)))
!siglelf     = 1./(1+exp((Iv*wl(wle)' - wl(wlf)*Ih)*1E2));
siglelf   = 1./(1.+exp((matmul(Iv,we)- matmul(wf,Ih))*1E2))

!print*, ' here 13 '
!sigmoid     = 1./(1+exp(-wlf/10)*exp(wle'/10));  % matrix computed as an outproduct
IF (.NOT.ALLOCATED(sigmoid)) ALLOCATE(sigmoid(size(wlf),size(wle)))
sigmoid     = 1./(1+matmul(exp(-wf/10.),exp(we/10.)))   !  % matrix computed as an outproduct
!print*, ' in fluspect siglelf ', minval(siglelf), maxval(siglelf), sum(siglelf)
!print*, ' in fluspect sigmoid ', minval(sigmoid), maxval(sigmoid), sum(sigmoid)
!print*, ' here 20' 

![MfI,  MbI]  = deal(.5 * fqe * ((.5*phiI( Iwlf))*eps) * kChl(Iwle)'.*sigmoid);
![MfII, MbII] = deal(.5 * fqe * ((.5*phiII(Iwlf))*eps) * kChl(Iwle)'.*sigmoid);


we(1,:) = kChl(Iwle)
!wf(:,1) = ((wI*phiI(Iwlf) + wII*phiII(Iwlf))*eps)
wf(:,1) = 0.5*phiI(Iwlf)*eps


!MfI = .5*fqe* (matmul(wf,we)*sigmoid)
MfI = fqe1* (matmul(wf,we)*sigmoid)
MbI = MfI


!wf(:,1) = (wII*phiII(Iwlf))*eps
wf(:,1) = 0.5*phiII(Iwlf)*eps
!MfII = .5*fqe* (matmul(wf,we)*sigmoid)
MfII = fqe2* (matmul(wf,we)*sigmoid)
MbII = MfII

!print*, ' in fluspect MbI ', minval(MbI), maxval(MbI), sum(MbI)
!print*, ' in fluspect MbII ', minval(MbII), maxval(MbII), sum(MbII)
!print*, ' in fluspect MfI ', minval(MfI), maxval(MfI), sum(MfI)
!print*, ' in fluspect MfII ', minval(MfII), maxval(MfII), sum(MfII)


!print*, ' here 21 '
IF (.NOT.ALLOCATED(MbnI))  ALLOCATE(MbnI(size(wlf),size(wle)))
IF (.NOT.ALLOCATED(MfnI))  ALLOCATE(MfnI(size(wlf),size(wle)))
IF (.NOT.ALLOCATED(MbnII)) ALLOCATE(MbnII(size(wlf),size(wle)))
IF (.NOT.ALLOCATED(MfnII)) ALLOCATE(MfnII(size(wlf),size(wle)))

!print*, ' here 22 '
IF (.NOT.ALLOCATED(wxe)) ALLOCATE(wxe(1,size(wle)))
IF (.NOT.ALLOCATED(wre)) ALLOCATE(wre(1,size(wle)))
IF (.NOT.ALLOCATED(wxf)) ALLOCATE(wxf(size(wlf),1))
IF (.NOT.ALLOCATED(wrf)) ALLOCATE(wrf(size(wlf),1))

!print*, ' here 3'
 
!% doubling routine
DO i = 1,ndub
    xe      = te/(1-re*re)
    xf      = tf/(1-rf*rf)
    ten     = te*xe
    tfn     = tf*xf
    ren     = re*(1+ten)
    rfn     = rf*(1+tfn)

    wxe(1,:) = xe
    wxf(:,1) = xf

    wre(1,:) = re
    wrf(:,1) = rf

 !   MfnI    = MfI  .*(xf*Ih + Iv*xe')         + MbI  .*(xf*xe').*(rf*Ih + Iv*re');
 !   MbnI    = MbI  .*(1+(xf*xe').*(1+rf*re')) + MfI.*((xf.*rf)*Ih+Iv*(xe.*re)');
 !   MfnII   = MfII .*(xf*Ih + Iv*xe')         + MbII .*(xf*xe').*(rf*Ih +Iv*re');
 !   MbnII   = MbII .*(1+(xf*xe').*(1+rf*re')) + MfII .*((xf.*rf)*Ih+Iv*(xe.*re)');


   MfnI  =  MfI*(matmul(wxf,Ih) + matmul(Iv,wxe)) + & 
          &  MbI*matmul(wxf,wxe)*(matmul(wrf,Ih) + matmul(Iv,wre))
   MbnI = MbI*(1.+matmul(wxf,wxe)*(1.+matmul(wrf,wre))) +  &
         & MfI*(matmul((wxf*wrf),Ih)+matmul(Iv,(wxe*wre))) 

   MfnII  =  MfII*(matmul(wxf,Ih) + matmul(Iv,wxe)) + &
          &  MbII*matmul(wxf,wxe)*(matmul(wrf,Ih) + matmul(Iv,wre))
   MbnII = MbII*(1.+matmul(wxf,wxe)*(1.+matmul(wrf,wre))) +  &
         & MfII*(matmul((wxf*wrf),Ih)+matmul(Iv,(wxe*wre)))

    te      = ten
    re      = ren
    tf      = tfn
    rf      = rfn

   MfI      = MfnI
   MbI      = MbnI  
  MfII      = MfnII
  MbII      = MbnII  


END DO 

DEALLOCATE(te,re,xe,ten,ren,Ih,we)
DEALLOCATE(mask,j)
DEALLOCATE(Iwle, Iwlf)
DEALLOCATE(MfnI,MbnI, MfnII,MbnII)

!print*, ' in fluspect refl ', minval(refl), maxval(refl),sum(refl)
!print*, ' in fluspect tran ', minval(tran), maxval(tran), sum(tran)
!print*, ' in fluspect kChlrel ', minval(kChlrel), maxval(kChlrel),sum(kChlrel)

!print*, ' in fluspect MbI ', minval(MbI), maxval(MbI), sum(MbI)
!print*, ' in fluspect MbII ', minval(MbII), maxval(MbII), sum(MbII)
!print*, ' in fluspect MfI ', minval(MfI), maxval(MfI), sum(MfI)
!print*, ' in fluspect MfII ', minval(MfII), maxval(MfII), sum(MfII)

!print*, 'size wlT ', size(wlT)
!print*, 'size wlP ', size(wlP)
!print*, 'size wlS ', size(wlS)

IF (.NOT. ALLOCATED(mask)) ALLOCATE(mask(size(wlS)))
mask = ((wlS>=minval(wlT)).and.(wlS<=maxval(wlT)))
IF (.NOT.ALLOCATED(IwlT)) ALLOCATE(IwlT(count(mask)))
CALL mask_ind(size(mask),mask,IwlT)


IF (.NOT. ALLOCATED(mask)) ALLOCATE(mask(size(wlS)))
mask = ((wlS>=minval(wlP)).and.(wlS<=maxval(wlP)))
IF (.NOT.ALLOCATED(IwlP)) ALLOCATE(IwlP(count(mask)))
CALL mask_ind(size(mask),mask,IwlP)

!print*, 'min max IwlP ', minval(IwlP),maxval(IwlP)
!print*, 'min max IwlT ', minval(IwlT),maxval(IwlT)
!print*, 'size rsfile ', size(rsfile(:,2)), 'IwlP ', size(IwlP)

ENDIF 


rho(IwlP) = refl
tau(IwlP) = tran
 rs(IwlP) = rsfile(:,2)

rho(IwlT) = rho_thermal
tau(IwlT) = tau_thermal
rs(IwlT)  = rs_thermal

!print*, ' rho(iwlfo) ', minval(rho(iwlfo)), maxval(rho(iwlfo)),sum(rho(iwlfo))
!print*, ' tau(iwlfo) ', minval(tau(iwlfo)), maxval(tau(iwlfo)),sum(tau(iwlfo))


END SUBROUTINE fluspect 


!function tav = calctav(alfa,nr)
SUBROUTINE calctav(alfa,nr,tav)

IMPLICIT NONE 

! Input variables 
INTEGER, INTENT(IN)                 :: alfa
REAL, INTENT(IN), DIMENSION(nwl)    :: nr 

! Ouput variables 
REAL, INTENT(OUT), DIMENSION(nwl)   :: tav 

! Local variables  
REAL, PARAMETER                     :: pi = 3.1415926535879 
REAL, DIMENSION(nwl)                :: np,n2,nm,a,k
REAL, DIMENSION(nwl)                :: b1,b2,b,b3,a3,ts
REAL, DIMENSION(nwl)                :: tp1,tp2,tp3,tp4,tp5,tp
REAL                                :: sa, rd
 
INTEGER                             :: i

rd          = pi/180
n2          = nr**2
np          = n2+1
nm          = n2-1
a           = (nr+1)*(nr+1)/2
k           = -(n2-1)*(n2-1)/4
sa          = sin(alfa*rd)



!b1          = (alfa~=90)*sqrt((sa**2-np/2)*(sa**2-np/2)+k)
b1 = 0. 
if (alfa.ne.90)  b1  = sqrt((sa**2-np/2)*(sa**2-np/2)+k)

b2          = sa**2-np/2
b           = b1-b2
b3          = b**3
a3          = a**3

ts          = (k**2./(6*b3)+k/b-b/2)-(k**2./(6*a3)+k/a-a/2)
tp1         = -2*n2*(b-a)/(np**2.)
tp2         = -2*n2*np*log(b/a)/(nm**2)

tp3         = n2*(1./b-1./a)/2
tp4         = 16*n2**2.*(n2**2+1)*log((2*np*b-nm**2)/(2*np*a-nm**2))/(np**3.*nm**2)
tp5         = 16*n2**3.*(1./(2*np*b-nm**2)-1./(2*np*a-nm**2))/(np**3)

tp          = tp1+tp2+tp3+tp4+tp5
tav         = (ts+tp)/(2*sa**2)

END SUBROUTINE calctav


!function [lidf]=  leafangles(a,b)
!% Subroutine FluorSail_dladgen

SUBROUTINE leafangles(a,b)                                     

!function [lidf]=  leafangles(a,b)                                     
!% Subroutine FluorSail_dladgen
!% Version 2.3 
!% For more information look to page 128 of "theory of radiative transfer models applied in optical remote sensing of
!% vegetation canopies"
!%
!% FluorSail for Matlab
!% FluorSail is created by Wout Verhoef, 
!% National Aerospace Laboratory (NLR)
!% Email: verhoef@nlr.nl
!%
!% This code was created by Joris Timmermans, 
!% International institute for Geo-Information Science and Earth Observation. (ITC)
!% Email: j_timmermans@itc.nl
!%
!%% main function

IMPLICIT NONE 

!Input variables 
REAL, INTENT(IN)                    :: a,b

! Local variables 
 REAL, DIMENSION(1,13)                :: F 
 REAL                                 :: theta 
 INTEGER                              :: i

!F           =   zeros(1,13);
F = 0.

DO i = 1, 8                                                               
    theta  =  i*10.             !% theta_l =  10:80
    CALL dcum(a,b,theta,F(1,i))     !% F(theta)
END DO 

DO i=9,12                                                              
    theta   =   80. + (i-8)*2.   !% theta_l = 82:88
    CALL dcum(a,b,theta,F(1,i))  !% F(theta)
END DO 

F(1,13) = 1.                     !%  theta_l = 90:90

lidf = 0. 
DO i=13,2,-1                                                           
    lidf(i) =  F(1,i) -   F(1,i-1)  !% lidf   =   dF/dtheta;
    !lidf(i,1) =  F(1,i) -   F(1,i-1)  !% lidf   =   dF/dtheta;

END DO 

!lidf(1,1) =   F(1,1)                  !% Boundary condition
lidf(1) =   F(1,1)                  !% Boundary condition

END SUBROUTINE leafangles  

!%
!function [F]   =  dcum(a,b,theta)
SUBROUTINE  dcum(a,b,theta,F)

IMPLICIT NONE 

REAL, PARAMETER             :: pi = 3.1415926535879 

! Input variables 
REAL, INTENT(IN)            :: a,b,theta 

! Output variables 
REAL, INTENT(OUT)           :: F 

! Local variables 
REAL                        :: rd,eps,delx,x1
REAL                        :: y,dx1,theta2


rd  =   pi/180.                ! %   Geometrical constant

if ( a>1. ) then  
    F  = 1 - cos(theta*rd)
else
    eps     =   1e-8
    delx    =   1.
    
    x1      =   2*rd *theta
    theta2  =   x1                                                                           
  do  while (delx > eps)
         y   =   a*sin(x1) + 0.5*b*sin(2*x1)
        dx1  =   0.5*(y - x1 + theta2)
        x1   =   x1 + dx1
        delx =   abs(dx1)
  end do 
    F    =   (2*y + theta2)/pi    ! %   Cumulative leaf inclination density function

    !%pag 139 thesis says: F = 2*(y+p)/pi. 
    !%Since theta2=theta*2 (in rad), this makes F=(2*y + 2*theta)/pi    
endif 

END SUBROUTINE dcum 

! Computation of the number of layers according to the LAI
  SUBROUTINE Nlayers(LAI)

   REAL, INTENT(IN)         :: LAI
   nl =anint(LAI) *10
   if (nl ==0) nl =10

END SUBROUTINE Nlayers


 !function [M] = aggreg(atmfile,SCOPEspec)

 SUBROUTINE aggreg (atmosfile,nreg,streg,enreg,width) 

 IMPLICIT NONE 

!% Aggregate MODTRAN data over SCOPE bands by averaging (over rectangular band
!% passes)
! atmosphere file
CHARACTER(len=120)                           :: header 

! Input variables 
INTEGER, INTENT(IN)                          :: nreg 
CHARACTER(len=80), INTENT(IN)                :: atmosfile 
REAL, DIMENSION(3),INTENT(IN)                :: streg,enreg,width 

!Output variables 
!REAL, ALLOCATABLE,DIMENSION(:,:),INTENT(OUT) :: atmoM 

! Local variables 
INTEGER                                      :: nwM 
REAL,   ALLOCATABLE,DIMENSION(:,:)           :: xdata
INTEGER,PARAMETER                            :: inunit = 2
REAL,  ALLOCATABLE,DIMENSION(:)              :: wlM
REAL,  ALLOCATABLE,DIMENSION(:,:)            :: U 
INTEGER                                      :: iwl,r,i,k
INTEGER                                      :: nwS
INTEGER, DIMENSION(3)                        :: nwreg,off  
REAL, ALLOCATABLE,DIMENSION(:,:)             :: S 
INTEGER, ALLOCATABLE,DIMENSION(:)            :: n 
INTEGER, ALLOCATABLE,DIMENSION(:)            :: j 
REAL                                         :: w 
 
!% Aggregate MODTRAN data over SCOPE bands by averaging (over rectangular band
!% passes)

print*, ' nreg ', nreg 
print*,  'streg ', streg 
print*, ' enreg ', enreg
print*, 'width ', width

!% Read .atm file with MODTRAN data
!s   = importdata(atmfile);
OPEN(unit=inunit,file=atmosfile,status='old')
REWIND inunit
READ(inunit,*) header
READ(inunit,*) header
 nwM = 0
    DO WHILE (.TRUE.)
       READ(inunit,*,END=1000) header ! just skip data
       nwM =  nwM+1
    END DO
1000 CONTINUE
CLOSE(inunit)
print*, ' nwM ', nwM 

IF (.NOT.ALLOCATED(xdata)) ALLOCATE(xdata(nwM,20))

OPEN(unit=inunit,file=atmosfile,status='old')
REWIND inunit
READ(inunit,*) header
READ(inunit,*) header
!print*, header

DO i=1, nwM
READ(inunit,*)xdata(i,:)
END DO
CLOSE(inunit)

!print*, xdata(1,:)
!print*
!print*, xdata(i-1,:)

IF (.NOT.ALLOCATED(wlM)) ALLOCATE(wlM(nwM))

wlM = xdata(:,2)
!T  = s.data(:,3:20);

!% Extract 6 relevant columns from T
!%  1: <Eso*costts/pi>
!%  3: <rdd>
!%  4: <tss>
!%  5: <tsd>
!% 12: <tssrdd>
!% 16: <La(b)>

!U     = [T(:,1) T(:,3) T(:,4) T(:,5) T(:,12) T(:,16)];
IF (.NOT.ALLOCATED(U)) ALLOCATE(U(nwM,6))
!U     =  [xdata(:,3),xdata(:,5),xdata(:,6),xdata(:,7),xdata(:,14),xdata(:,18)]
U(:,1) = xdata(:,3)
U(:,2) = xdata(:,5)
U(:,3) = xdata(:,6)
U(:,4) = xdata(:,7)
U(:,5) = xdata(:,14)
U(:,6) = xdata(:,18)

!DO  i=1,6
!print*, ' U ',i, minval(U(:,i)), maxval(U(:,i)), sum(U(:,i))
!print*, U(:,i)
!END DO
!stop 
!nwM   = length(wlM);

!SCOPEspec.nreg;
!nreg  = 3 

!streg = SCOPEspec.start;
!streg = [ 400,  2500,  16000]

!enreg = SCOPEspec.end;
!enreg = [2400, 15000,  50000]

!width = SCOPEspec.res;
!width =  [   1,   100,   1000]


!% Nr. of bands in each region
!nwreg = int32((enreg-streg)./width)+1;
nwreg = int((enreg-streg)/width)+1

print*, ' nwreg', nwreg 

!off   = int32(zeros(nreg,1));
off   =  0

!for i=2:nreg
!    off(i) = off(i-1)+nwreg(i-1);
!end

DO i =2, nreg 
  off(i)= off(i-1) + nwreg(i-1)
END DO 

!print*, 'off ', minval(off), maxval(off), sum(off)

!nwS = sum(nwreg);
nwS = sum(nwreg)
print* , ' nwS ', nwS 

!n   = zeros(nwS,1);    !% Count of MODTRAN data contributing to a band
IF (.NOT.ALLOCATED(n)) ALLOCATE(n(nwS))
n = 0


!S   = zeros(nwS,6);    ! Intialize sums
IF (.NOT.ALLOCATED(S)) ALLOCATE(S(nwS,6))
!S   = zeros(nwS,6);    ! Intialize sums
S   = 0.

!%k   = int32(0);
!j   = int32(zeros(nreg,1))  !% Band index within regions
IF (.NOT.ALLOCATED(j)) ALLOCATE(j(nreg))
j = 0 


!for iwl = 1:nwM
!    w   = wlM(iwl);    % MODTRAN wavelength
!    for r = 1:nreg
!        j(r) = int32(round(w-streg(r))./(width(r)))+1;
!        if j(r)>0 && j(r)<=nwreg(r)                 % test if index is in valid range
!            k      = j(r)+off(r);                   % SCOPE band index
!            S(k,:) = S(k,:)+U(iwl,:);               % Accumulate from contributing MODTRAN data
!            n(k)   = n(k)+1;                        % Increment count
!        end
!    end
!end


DO iwl = 1,nwM
    w   = wlM(iwl)     !    % MODTRAN wavelength
    DO r = 1,nreg
        j(r) = nint(nint(w-streg(r))/(width(r)))+1

            !print*, ' w streg width ', w,streg(r),width(r)
            !print*, ' nint w-streg ratio ', nint(w-streg(r)), nint(w-streg(r))/(width(r))
            !print*, ' iwl nwreg jr ', iwl,nwreg(r),j(r)

        if( (j(r)>0).and.(j(r)<=nwreg(r))) then ! % test if index is in valid  range
            k      = j(r)+off(r)               !    % SCOPE band index
            !print*, ' nwreg ', nwreg(r)
            !print*, ' k ', k, ' j r ', j(r), ' off ', off(r), ' n ', n(k) , iwl 
            !print*, 'BEFORE S ', S(k,:)

            S(k,:) = S(k,:)+U(iwl,:)           !    % Accumulate from contributing MODTRAN data
            n(k)   = n(k)+1                    !    % Increment count 
            !print*, '  U ', iwl,U(k,:) 
            !print*, 'AFTER S ', S(k,:) 

        endif
   END DO
END DO 


IF (.NOT.ALLOCATED(atmoM)) ALLOCATE(atmoM(nwS,6))
atmoM = 0.
DO  i = 1,6
!print*, 'S ',  i,minval(S(:,i) ), maxval(S(:,i)), sum(S(:,i)) 
    atmoM(:,i) = S(:,i)/n      !% Calculate averages per SCOPE band
!print*, 'atmoM ',  i,minval(atmoM(:,i) ), maxval(atmoM(:,i)), sum(atmoM(:,i)) 
!print*, 'n ',  i,minval(n), maxval(n), sum(n) 
END DO 

print*,minval(atmoM), maxval(atmoM), sum(atmoM)

DEALLOCATE(xdata,wlM,S,U,j,n)

END SUBROUTINE aggreg


END MODULE fluo_param 

