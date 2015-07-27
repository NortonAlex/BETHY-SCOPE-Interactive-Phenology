MODULE mo_rtmo  


CONTAINS 


SUBROUTINE rtmo(Rin,Rli,Ta,LAI,tts,tto,psi, Ps, Po, Pso, km, Kext, &
        & Esun_,Esky_,P, fEsuno,fEskyo,fEsunt,fEskyt, Eplu_, Emin_, &
        & Lo_, Eout_, Eouto, Eoutt, Rnhs, Rnus, Rnhc, Rnuc,&
        & Pnhc, Pnuc, Pnhc_Cab, Pnuc_Cab, rho, tau, rs, kClrel)

USE fluo_func

!%% function RTMo

!% calculates the spectra of hemisperical and directional observed visible
!% and thermal radiation (fluxes E and radiances L), as well as the single
!% and bi-directional gap probabilities
!%
!% the function does not require any non-standard Matlab functions. No
!% changes to the code have to be made to operate the function for a
!% particular canopy. All necessary parameters and variables are input or
!% global and need to be specified elsewhere.
!%
!% Authors:      Wout Verhoef            (verhoef@nlr.nl)
!%               Christiaan van der Tol  (tol@itc.nl)
!%               Joris Timmermans        (j_timmermans@itc.nl)
!%
!% updates:      10 Sep 2007 (CvdT)      - calculation of Rn
!%                5 Nov 2007             - included observation direction
!%               12 Nov 2007             - included abs. PAR spectrum output
!%                                       - improved calculation efficiency
!%               13 Nov 2007             - written readme lines
!%               11 Feb 2008 (WV&JT)     - changed Volscat
!%                           (JT)        - small change in calculation Po,Ps,Pso
!%                                        - introduced parameter 'lazitab'
!%                                       - changed nomenclature
!%                                       - Appendix IV: cosine rule
!%               04 Aug 2008 (JT)        - Corrections for Hotspot effect in the probabilities
!%               05 Nov 2008 (CvdT)      - Changed layout
!%               04 Jan 2011 (JT & CvdT) - Included Pso function (Appendix IV)
!%                                       - removed the analytical function (for checking)
!%               02 Oct 2012 (CvdT)      - included incident PAR in output
!%
!%               Jan/Feb 2013 (WV)       - Major revision towards SCOPE version 1.40:
!%                                       - Parameters passed using structures
!%                                       - Improved interface with MODTRAN atmospheric data
!%                                       - Now also calculates 4-stream
!%                                         reflectances rso, rdo, rsd and rdd
!%                                         analytically
!%               Apri 2013 (CvT)         - improvements in variable names
!%                                           and descriptions

!%
!% Table of contents of the function
!%
!%   0.      Preparations
!%       0.1     parameters
!%       0.2     initialisations
!%   1.      Geometric quantities
!%       1.1     general geometric quantities
!%       1.2     geometric factors associated with extinction and scattering
!%       1.3     geometric factors to be used later with rho and tau
!%       1.4     solar irradiance factor for all leaf orientations
!%       1.5     probabilities Ps, Po, Pso
!%   2.      Calculation of upward and downward fluxes
!%   3.      Outgoing fluxes, hemispherical and in viewing direction, spectrum
!%   4.      Net fluxes, spectral and total, and incoming fluxes
!%   A1      functions J1 and J2 (introduced for stable solutions)
!%   A2      function volscat
!%   A3      function e2phot
!%   A4      function Pso
!%
!% references:
!%{1} Verhoef (1998), 'Theory of radiative transfer models applied in
!%    optical remote sensing of vegetation canopies'. PhD Thesis Univ. Wageninegn
!%{2} Verhoef, W., Jia, L., Xiao, Q. and Su, Z. (2007) Unified optical -
!%    thermal four - stream radiative transfer theory for homogeneous
!%    vegetation canopies. IEEE Transactions on geoscience and remote
!%    sensing, 45,6.
!%{3} Verhoef (1985), 'Earth Observation Modeling based on Layer Scattering
!%    Matrices', Remote sensing of Environment, 17:167-175
!%
!% Usage:
!% function [rad,gap] = RTMo(spectral,atmo,soil,leafopt,canopy,angles)
!%
!% The input and output are structures. These structures are further
!% specified in a readme file.
!%
!%
!% Input:
!%   spectral    information about wavelengths and resolutions
!%   atmo        MODTRAN atmospheric parameters
!%   soil        soil properties
!%   leafopt     leaf optical properties
!%   canopy      canopy properties (such as LAI and height)
!%   angles      viewing and observation angles
!%
!% Output: (see at the end of the subroutine )
!%   gap         probabilities of direct light penetration and viewing
!%   rad         a large number of radiative fluxes: spectrally distributed
!%               and integrated, and canopy radiative transfer coefficients.

!% Input:
!%   Symbol  Description                 Unit            Dimension
!%   ------  -----------                 ----            ---------
!%   Esun_   direct solar radiation      (W m-2 um-1)    [nwl]
!%   Esky_   diffuse sky radiation       (W m-1 um-1)    [nwl]
!%   LAI     leaf area index             (m2 m-2)        [1]
!%   tts     solar zenith angle          (degrees)       [1]
!%   tto     viewing angle               (degrees)       [1]
!%   psi     azimuth angle difference between solar and viewing position
!%           (degrees)                   (degrees)       [1]
!%
!% Output:
!%   Symbol  Description                                 Dimension
!%   ------  -----------                                 ---------
!%   Rnhs    Net radiation of shaded soil (W m-2)        [1]
!%   Rnus    Net radiation of sunlit soil                [1]
!%   Rnhc    Net radiation of shaded leaves              [nl]
!%   Rnuc    Net radiation of sunlit leaves              [13,36,nl]
!%   Pnh     Net PAR by shaded leaves (mol m-2 s-1)      [nl]
!%   Pnu     Net PAR by sunlit leaves                    [13,36,nl]
!%   Eplu_   Total upward diffuse radiation (W m-2 um-1) [nl+1]
!%   Emin_   Total downward diffuse radiation            [nl+1]
!%   Louto   Total outgoing hemispherical rad (<=2.5 um) [1]
!%   Loutt   Total outgoing hemishperical rad (>2.5 um)  [1]
!%   Lout_   Spectrum of outgoing hemispherical rad      [nwl]
!%   Lo_     Spectrum of outgoing rad in viewing dir     [nwl]
!%   PARhc   Incident PAR on shaded leaves               [nl]
!%   PARuc   Incident PAR on sunlit leaves               [13,36,nl]
!%   Ps      probability of sunlit leaves                [nl+1]
!%   Po      probability of viewing a leaf or soil       [nl+1]
!%   Pso     probability of viewing a sunlit leaf/soil   [nl+1]
!%   Kext       extinction coefficient in viewing dir       [1]
!%
!% Notes:
!%   nl      number of layers
!%   nwl     number of wavelengths
!% '_'means: a flux at different wavelengths (a vertically oriented vector)
!% page numbers in the explanation refer to Verhoef (1998)
!%
!%% 0. Preparations
!global nwl nl
!global x dx prm_ lidf q
!global litab lazitab nlazi
!global wl resolution 
!global deg2rad

! ---------------------------------------------------------------------------
!  E Koffi, IPSL, Jussieu, Paris, July 2012 
! Translation of the RTMo.m matlab code of Christian Van der Tol in Fortran 
! ----------------------------------------------------------------------------

USE fluo_param, ONLY:  nwl, nl, xlay, dx, prm_, lidf, q, litab, &
                        & lazitab, nli,nlazi, &
                        & wl,wlT,wlP,wlS,wlF,wlF,wlPAR,&
                        & nwlfi, nwlfo, nwlT, nwlP,nwlPAR, &
                        & atmoM

IMPLICIT NONE

REAL, PARAMETER                             :: pi  = 3.1415926535879
REAL, PARAMETER                             :: deg2rad = pi/180.

! Input variables 
REAL, INTENT(IN)                            :: LAI,tts,tto,psi
REAL, INTENT(IN)                            :: Rin,Rli,Ta 
REAL, DIMENSION(:), INTENT(IN)            :: rho,tau,rs,kClrel    ! fluspect output

! Output variables 
REAL, DIMENSION(nl+1,nwl), INTENT(OUT)      :: Eplu_,Emin_
REAL, DIMENSION(nl), INTENT(OUT)            :: Rnhc,Pnhc,Pnhc_Cab
REAL, DIMENSION(nli,nlazi,nl),INTENT(OUT)   :: Rnuc,Pnuc,Pnuc_Cab
REAL, DIMENSION(nwl), INTENT(OUT)           :: Esun_,Esky_
REAL, DIMENSION(nwl), INTENT(OUT)           :: fEsuno,fEskyo,fEsunt,fEskyt
REAL, DIMENSION(nwl), INTENT(OUT)           :: Lo_, Eout_
REAL, DIMENSION(nl+1),INTENT(OUT)           :: Ps , Po,  Pso
REAL, INTENT(OUT)                           :: Rnhs,Rnus
REAL, INTENT(OUT)                           :: Eouto,Eoutt
REAL, INTENT(OUT)                           :: km, Kext 
REAL, INTENT(OUT)                           :: P

! Variables with nwl dimension 
REAL, DIMENSION(nwl)                        :: epsc, epss
REAL, DIMENSION(nwl)                        :: sigb,sigf,sb
REAL, DIMENSION(nwl)                        :: sf,vb,vf,w
REAL, DIMENSION(nwl)                        :: s1,s2,v1,v2
REAL, DIMENSION(nwl)                        :: a,m,rinf,rinf2

REAL, DIMENSION(nwl)                        :: J1Km,J2km
REAL, DIMENSION(nwl)                        :: J1ks,J2ks
REAL, DIMENSION(nwl)                        :: e1,e2
REAL, DIMENSION(nwl)                        :: re,denom
REAL, DIMENSION(nwl)                        :: Pss,Qss
REAL, DIMENSION(nwl)                        :: Poo,Qoo
REAL, DIMENSION(nwl)                        :: tau_dd,rho_dd
REAL, DIMENSION(nwl)                        :: tau_sd,rho_sd
REAL, DIMENSION(nwl)                        :: tau_do,rho_do
REAL, DIMENSION(nwl)                        :: T1,T2
REAL, DIMENSION(nwl)                        :: rso
REAL, DIMENSION(nwl)                        :: rho_sod,rho_sos,rho_so
REAL, DIMENSION(nwl)                        :: t1m,t3m,t4m,t5m,t12m,t16m
REAL, DIMENSION(nwl)                        :: rsd,rdd,rdo
REAL, DIMENSION(nwl)                        :: Eplu_1,Eplu0,Emin_1
REAL, DIMENSION(nwl)                        :: delta1,delta2

REAL, DIMENSION(nwl)                        :: temp,res,Ls
REAL, DIMENSION(nwl)                        :: Fd

REAL, DIMENSION(nwl)                        :: Psun_
REAL, DIMENSION(nwl)                        :: piLoc_ ,piLos_
REAL, DIMENSION(nwl)                        :: piLo_
REAL, DIMENSION(nwl)                        :: E_

! Variables with nl 
REAL, DIMENSION(nl)                         :: Pnh
REAL, DIMENSION(nl)                         :: Pdif
REAL, DIMENSION(nl)                         :: Pndif_Cab 
REAL, DIMENSION(nl)                         :: Pndif,Rndif 
REAL, DIMENSION(nl)                         :: x
REAL, DIMENSION(nl+1)                       :: J1kx,J2kx
REAL, DIMENSION(nl+1)                       :: F1,F2
REAL, DIMENSION(nl+1)                       :: Psoi
REAL, DIMENSION(nl+1)                       :: xl

! Variables with nli 
REAL, DIMENSION(nli)                        :: frho,ftau
REAL, DIMENSION(nli)                        :: chi_s,chi_o
REAL, DIMENSION(nli)                        :: cos_ttli
REAL, DIMENSION(nli)                        :: sin_ttli, ksli
REAL, DIMENSION(nli)                        :: koli, sobli, sofli , bfli
REAL, DIMENSION(nli)                        :: Cs,Ss

! Variables with nlazi
REAL, DIMENSION(nlazi)                      :: cos_ttlo

! Variables with nl, nwl 
REAL, DIMENSION(nl,nwl)                     :: Rndif_
REAL, DIMENSION(nl,nwl)                     :: Pndif_Cab_,Pndif_

! Variables with nli, nlazi 
REAL, DIMENSION(nli,nlazi)                  :: cos_deltas,fs
REAL, DIMENSION(nli,nlazi)                  :: Pdir,Rndir
REAL, DIMENSION(nli,nlazi)                  :: Pndir,Pndir_Cab

! Variables with nli, nl, nlazi 
REAL, DIMENSION(nli,nlazi,nl)               :: Puc, Pnu 

! Variables with nwlP
INTEGER,  DIMENSION(nwlP)                   :: IwlP

!Variables with nwlT
INTEGER,  ALLOCATABLE,DIMENSION(:)          :: IwlT

! Constants
REAL                                        :: cos_tts, tan_tto
REAL                                        :: cos_tto ,sin_tts
REAL                                        :: tan_tts, dso
REAL                                        :: psi1
REAL                                        :: iLAI
REAL, DIMENSION(1)                          :: tx,tm
double precision                            :: integral
!external                                    :: func

REAL                                        :: bf,sob,sof,sdb,sdf
REAL                                        :: ddb,ddf,dob,dof

REAL                                        :: tau_ss,tau_oo, Z
REAL                                        :: Pso2w
REAL                                        :: Esunto, Eskyto,Etoto
REAL                                        :: Esuntt, Eskytt,Etott
REAL                                        :: ap
REAL                                        :: Psun ,Asun, Pnsun, Pnsun_Cab
!REAL                                        :: Louto, Loutt
REAL                                        :: Rndirsoil,Rndifsoil


! Other ALLOCATABL arrays 
LOGICAL, ALLOCATABLE, DIMENSION(:)          :: mask 
INTEGER, ALLOCATABLE, DIMENSION(:)          :: indices 
INTEGER, ALLOCATABLE, DIMENSION(:)          :: J_o, J_t 
REAL, ALLOCATABLE, DIMENSION(:)             :: P_,tempIpar 
INTEGER,ALLOCATABLE,DIMENSION(:)            :: Ipar

INTEGER                                     :: i,j,n
INTEGER                                     :: ike,ke 
INTEGER                                     :: il

REAL                                        :: minPAR, maxPAR


!%% 0. Preparations
minPAR  = minval(wlPAR);
maxPAR  = maxval(wlPAR);

!Ipar    = find(wl>=minPAR & wl<=maxPAR); % Indices for PAR wavelenghts within wl
IF (.NOT. ALLOCATED(mask)) ALLOCATE(mask(size(wlPAR)))
mask = ((wl>=minPAR).AND.(wl<=maxPAR))
IF (.NOT.ALLOCATED(Ipar)) ALLOCATE(Ipar(count(mask)))
CALL mask_ind(size(mask),mask,Ipar)


!% 0.1 parameters
epsc        = 1-rho-tau
epss        = 1-rs                         !% [nwl]     emissivity  (soil)
iLAI        = LAI/nl                       !% [1]       LAI of a layer

!xl          = [0, x]                       !% [nl+1]    leaf layers + soil
xl = 0.
xl(2:nl+1) = xlay

!% 0.2 initialisations (allocation of memory)
!   Rndif         = zeros(nl+1,1)                !% [nl+1]    abs. diffuse rad soil+veg
![PARdif,Pndif]   = deal(zeros(nl,1))          !% [nl]      net PAR veg
![Emin_,Eplu_]    = deal(zeros(nl+1,nwl))       !% [nl,nwl]  up and down diff. rad.
![Rndif_,Pndif_]  = deal(zeros(nl,nwl))       !% [nl,nwl]  abs diff and PAR veg.
![PARuc,Rnuc,Pnu] = deal(zeros(length(litab),nlazi,nl)) ! [nli,nlazi,nl] inc and net rad and PAR sunlit

!print*, ' rho ', minval(rho), maxval(rho),sum(rho)
!print*, ' tau ', minval(tau), maxval(tau), sum(tau)
!print*, ' rs ', minval(rs), maxval(rs), sum(tau)
!print*, ' kClrel ', minval(kClrel), maxval(kClrel)


Rndif  = 0.
Pndif  = 0.
Emin_  = 0.
Eplu_  = 0.
Rndif_ = 0.
Pndif_ = 0.
Rnuc   = 0.
Pnu    = 0.

!%% 1.0 Geometric quantities
!% 1.1 general geometric quantities
!% these variables are scalars
cos_tts     = cos(tts*deg2rad)             !%           cos solar       angle   
tan_tto     = tan(tto*deg2rad)             !%           tan observation angle

cos_tto     = cos(tto*deg2rad)             !%           cos observation angle   
sin_tts     = sin(tts*deg2rad)             !%           sin solar       angle
tan_tts     = tan(tts*deg2rad)             !%           tan observation angle


!psi         = abs(psi-360*round(psi/360)) !%           (to ensure that volscatt is symmetric for psi=90 and psi=270)
psi1         = abs(psi-360*NINT(psi/360))   !%           (to ensure that volscatt is symmetric for psi=90 and psi=270)
!dso         = sqrt(tan_tts.^2 + tan_tto.^2 - 2*tan_tts.*tan_tto.*cos(psi*deg2rad));
dso         = sqrt(tan_tts**2 + tan_tto**2 - 2*tan_tts*tan_tto*cos(psi*deg2rad))

!print*, ' tts =', tts, ' tto ', tto, ' psi1 =', psi1
!print*, ' cos_tts ', cos_tts, ' tan_tto ', tan_tto, ' cos_tto ', cos_tto 
!print*, ' sin_tts ', sin_tts, ' tan_tts ', tan_tts, ' psi ', psi1, ' dso ', dso 

! 1.2 geometric factors associated with extinction and scattering
![chi_s,chi_o,frho,ftau]=volscat(tts,tto,psi,litab)     !% volume scattering
call volscat(tts,tto,psi1,litab,chi_s,chi_o,frho,ftau)     !% volume scattering

cos_ttlo    = cos(lazitab*deg2rad)       ! %    cos leaf orientation angles

cos_ttli    = cos(litab*deg2rad)           !% [13]      cos leaf angles
sin_ttli    = sin(litab*deg2rad)           !% [13]      sinus leaf angles


ksli        = chi_s/cos_tts               !% [13]  p306{1} extinction coefficient in direction of sun per leaf angle
koli        = chi_o/cos_tto               !% [13]   p307{1} extinction coefficient in direction of observer per leaf angle

sobli       = frho*pi/(cos_tts*cos_tto)   !% [13]      pag 309{1} area scattering coefficient fractions
sofli       = ftau*pi/(cos_tts*cos_tto)   !% [13]      pag 309{1}
bfli        = cos_ttli**2.                !% [13]

!print*, 'cos_ttlo ', minval(cos_ttlo), maxval(cos_ttlo), sum(cos_ttlo)
!print*, 'cos_ttli ', minval(cos_ttli), maxval(cos_ttli), sum(cos_ttli)
!print*, 'sin_ttli ', minval(sin_ttli), maxval(sin_ttli), sum(sin_ttli)
!print*, 'ksli ', minval(ksli), maxval(ksli), sum(ksli)
!print*, 'koli ', minval(koli), maxval(koli), sum(koli)
!print*, 'sobli ', minval(sobli), maxval(sobli), sum(sobli)
!print*, 'bfli ', minval(bfli), maxval(bfli), sum(bfli)

!%integration over angles (using a vector inproduct) -> scalars
!k          = ksli'*lidf;                  ! %  pag 306{1}    extinction coefficient in direction of sun.
km          = DOT_PRODUCT(ksli,lidf)  ! % pag 306{1}    extinction coefficient in direction of sun.
Kext        = DOT_PRODUCT(koli,lidf) ! %  pag 307{1}    extinction coefficient in direction of observer
bf          = DOT_PRODUCT(bfli,lidf) ! % '
sob         = DOT_PRODUCT(sobli,lidf) ! % '         weight of specular2directional back    scatter coefficient
sof         = DOT_PRODUCT(sofli,lidf) !  % '          weight of specular2directional forward scatter coefficient

!print*, 'km ', km, ' Kext ', Kext, ' bf ', bf, ' sob ', sob, ' sof ', sof 


!% 1.3 geometric factors to be used later with rho and tau, f1 f2 of pag 304:
!% these variables are scalars 
!sdb         = 0.5*(k+bf)                   !% fs*f1
!sdf         = 0.5*(k-bf)                   !% fs*f2     weight of specular2diffuse     foward  scatter coefficient 
sdb         = 0.5*(km+bf)                   !% fs*f1
sdf         = 0.5*(km-bf)                   !% fs*f2     weight of specular2diffuse     foward  scatter coefficient 

ddb         = 0.5*(1+bf)                   !% f1^2+f2^2 weight of diffuse2diffuse      back    scatter coefficient 
ddf         = 0.5*(1-bf)                   !% 2*f1*f2   weight of diffuse2diffuse      forward scatter coefficient 
dob         = 0.5*(Kext+bf)                   !% fo*f1     weight of diffuse2directional  back    scatter coefficient 
dof         = 0.5*(Kext-bf)                   !% fo*f2     weight of diffuse2directional  forward scatter coefficient 

!print*, ' sdb ', sdb , ' sdf ', sdf, ' ddb ', ddb 
!print*, ' ddf ', ddf , ' dob ', dob, ' dof ', dof 


! 1.4 solar irradiance factor for all leaf orientations
Cs          = cos_ttli*cos_tts             !% [nli]     pag 305 modified by Joris
Ss          = sin_ttli*sin_tts             !% [nli]     pag 305 modified by Joris

!print*, 'Cs ', minval(Cs), maxval(Cs), sum(Cs)
!print*, 'Ss ', minval(Ss), maxval(Ss), sum(Ss)


!cos_deltas  = Cs*ones(1,nlazi) + Ss*cos_ttlo !%[nli,nlazi]
do j=1, nlazi
  do i=1, nli
   cos_deltas(i,j)=Cs(i) + Ss(i)*cos_ttlo(j)
   end do
end do

fs          = abs(cos_deltas/cos_tts)       !% [nli,nlazi] pag 305

!print*, ' cos_deltas ',  minval(cos_deltas), maxval(cos_deltas), sum(cos_deltas) 
!print*, ' fs ',  minval(fs), maxval(fs), sum(fs) 
!print*, ' xl  ', minval(xl), maxval(xl), sum(xl)
!print*, ' dx  ', dx

! 1.5 probabilities Ps, Po, Pso
Ps          =   exp(km*xl*LAI)                            ! % [nl+1]    p154{1} probability of viewing a leaf in solar dir
Po          =   exp(Kext*xl*LAI)                            ! % [nl+1]    p154{1} probability of viewing a leaf in observation dir
Ps(1:nl)    =   Ps(1:nl) *(1-exp(-km*LAI*dx))/(km*LAI*dx)  ! % Correct Ps/Po for finite dx
Po(1:nl)    =   Po(1:nl) *(1-exp(-Kext*LAI*dx))/(Kext*LAI*dx)  ! % Correct Ps/Po for finite dx

!print*, ' q ', q , ' LAI ', LAI, ' dso ', dso
!print *, 'Shape of xl is ', shape(xl)    !ANorton 
do j=1,size(xl)
!    Psoi(j,:)=   quad(@(y)Psofunction(K,k,LAI,q,dso,y),xl(j)-dx,xl(j))/dx
   !print *,'Entered loop over xl for simpson call'      !ANorton
    n= 4
   call simpson(Kext,km,LAI,q,dso,func,xl(j)-dx,xl(j),integral,n)
   Pso(j) = integral/dx
   
end do
!print *, 'Passed the loop over xl for simpson call'     !ANorton

!print*, ' Ps 1 ', minval(Ps), maxval(Ps), sum(Ps)
!print*, ' Po 1 ', minval(Po), maxval(Po), sum(Po)
!print*, ' Pso 1', minval(Pso), maxval(Pso), sum(Pso)

!Pso(Pso >Po) Pso= min([Po(Pso>Po),Ps(Pso>Po)],[],2)    !%takes care of rounding error
!Pso(Pso>Ps)= min([Po(Pso>Ps),Ps(Pso>Ps)],[],2)    !%takes care of rounding error
WHERE (Pso > Po) 
Pso = min(Po,Ps) 
ENDWHERE 

WHERE (Pso > Ps)
Pso = min(Po,Ps)
ENDWHERE

!print*, ' Ps 2 ', minval(Ps), maxval(Ps), sum(Ps)
!print*, ' Po 2 ', minval(Po), maxval(Po), sum(Po)
!print*, ' Pso 2', minval(Pso), maxval(Pso), sum(Pso)

!print*, ' size tau and rho ',  size(tau), size(rho)
!print*, ' tau ',  minval(tau), maxval(tau), sum(tau)
!print*, ' rho ',  minval(rho), maxval(rho), sum(rho)


!%% 2. Calculation of upward and downward fluxes
!% the following are vectors with lenght nwl
sigb        = ddb*rho+ddf*tau              !% [nwl]     sigmab, p305{1} diffuse     backscatter scattering coefficient for diffuse  incidence 
sigf        = ddf*rho+ddb*tau              !% [nwl]     sigmaf, p305{1} diffuse     forward     scattering coefficient for forward  incidence 
sb          = sdb*rho + sdf*tau            !% [nwl]     sb,     p305{1} diffuse     backscatter scattering coefficient for specular incidence 
sf          = sdf*rho + sdb*tau            !% [nwl]     sf,     p305{1} diffuse     foward      scattering coefficient for specular incidence 
vb          = dob*rho + dof*tau            !% [nwl]     vb,     p305{1} directional backscatter scattering coefficient for diffuse  incidence 
vf          = dof*rho + dob*tau            !% [nwl]     vf,     p305{1} directional forward     scattering coefficient for diffuse  incidence 
w           = sob*rho + sof*tau            !% [nwl]     w,      p309{1} bidirectional scattering coefficent (directional-directional)         
a           = 1-sigf                       !% [nwl]     attenuation
!m           = sqrt(a.^2-sigb.^2);          !% [nwl]
m           = sqrt(a**2-sigb**2)           !% [nwl]
rinf        = (a-m)/sigb                  !% [nwl]
rinf2       = rinf*rinf                   !% [nwl]

!print*, ' sigb ', minval(sigb), maxval(sigb), sum(sigb)
!print*, ' sigf ', minval(sigf), maxval(sigf), sum(sigf)
!print*, ' sb ', minval(sb), maxval(sb), sum(sb)
!print*, ' sf ', minval(sf), maxval(sf), sum(sf)
!print*, ' vb ', minval(vb), maxval(vb), sum(vb)
!print*, ' vf ', minval(vf), maxval(vf), sum(vf)
!print*, ' w ', minval(w), maxval(w), sum(w)
!print*, ' a ', minval(a), maxval(a), sum(a)
!print*, ' m ', minval(m), maxval(m), sum(m)
!print*, ' rinf ', minval(rinf), maxval(rinf), sum(rinf)
!print*, ' rinf2 ', minval(rinf2), maxval(rinf2), sum(rinf2)


!% direct solar radiation
!J1ks        = calcJ1(-1, m,km,LAI)         !% [nwl]
CALL calcJ1(nwl,-1.,m,km,LAI,J1ks)         !% [nwl]
!J2ks        = calcJ2( 0, m,km,LAI)         !% [nwl]
CALL  calcJ2(nwl,0.,m,km,LAI,J2ks)         !% [nwl]

!J1K        = calcJ1(-1, m,K,LAI);          % [nwl]   % added for calculation of rdo
CALL calcJ1(nwl,-1.,m,Kext,LAI,J1Km)       !% [nwl]
!J2K        = calcJ2( 0, m,K,LAI);          % [nwl]   % added for calculation of rdo
CALL calcJ2(nwl,0.,m,Kext,LAI,J2Km)        !% [nwl]

!print*, ' J1ks ', minval(J1ks), maxval(J1ks), sum(J1ks)
!print*, ' J2ks ', minval(J2ks), maxval(J2ks), sum(J2ks)
!print*, ' J1Km ', minval(J1Km), maxval(J1Km), sum(J1Km)
!print*, ' J2Km ', minval(J2Km), maxval(J2Km), sum(J2Km)

e1          = exp(-m*LAI)                !% [nwl]
e2          = e1**2                      !% [nwl]
re          = rinf*e1                   !% [nwl]
denom       = 1-rinf2*e2                !% [nwl]

!print*, ' e1 ', minval(e1), maxval(e1), sum(e1)
!print*, ' e2 ', minval(e2), maxval(e2), sum(e2)
!print*, ' re ', minval(re), maxval(re), sum(re)
!print*, ' denom ', minval(denom), maxval(denom), sum(denom)


s1          = sf+rinf*sb
s2          = sf*rinf+sb
v1          = vf+rinf*vb
v2          = vf*rinf+vb

!print*, ' s1 ', minval(s1), maxval(s1), sum(s1)
!print*, ' s2 ', minval(s2), maxval(s2), sum(s2)
!print*, ' v1 ', minval(v1), maxval(v1), sum(v1)
!print*, ' v2 ', minval(v2), maxval(v2), sum(v2)


!Pss         = (sf+rinf*sb)*J1ks;         ! % [nwl]
!Qss         = (sf*rinf+sb)*J2ks;         ! % [nwl]

Pss         = s1*J1ks         ! % [nwl]
Qss         = s2*J2ks         ! % [nwl]

Poo         = v1*J1Km    !    % (nwl)   % added for calculation of rdo
Qoo         = v2*J2Km    !  % [nwl]   % added for calculation of rdo

!print*, ' Pss ', minval(Pss), maxval(Pss), sum(Pss)
!print*, ' Qss ', minval(Qss), maxval(Qss), sum(Qss)
!print*, ' Poo ', minval(Poo), maxval(Poo), sum(Poo)
!print*, ' Qoo ', minval(Qoo), maxval(Qoo), sum(Qoo)

tau_ss      = exp(-km*LAI)               !   % [1]
tau_oo      = exp(-Kext*LAI)            !      % [1]

Z           = (1 - tau_ss * tau_oo)/(Kext + km)  !;  % needed for analytic rso


!print*, ' tau_ss ', tau_ss
!print*, ' tau_oo ', tau_oo
!print*, ' Z ', Z

tau_dd      = (1-rinf2)*e1 /denom     !   % [nwl]
rho_dd      = rinf*(1-e2)  /denom     !   % [nwl]
tau_sd      = (Pss-re*Qss) /denom   !     % [nwl]
tau_do      = (Poo-re*Qoo) /denom   !     % [nwl]
rho_sd      = (Qss-re*Pss) /denom   !     % [nwl]
rho_do      = (Qoo-re*Poo) /denom   !     % (nwl)

!print*, ' tau_dd ', minval(tau_dd), maxval(tau_dd), sum(tau_dd)
!print*, ' rho_dd ', minval(rho_dd), maxval(rho_dd), sum(rho_dd)
!print*, ' tau_sd ', minval(tau_sd), maxval(tau_sd), sum(tau_sd)
!print*, ' tau_do ', minval(tau_do), maxval(tau_do), sum(tau_do)
!print*, ' rho_sd ', minval(rho_sd), maxval(rho_sd), sum(rho_sd)
!print*, ' rho_do ', minval(rho_do), maxval(rho_do), sum(rho_do)


T1          = v2*s1*(Z-J1ks*tau_oo)/(Kext+m)+v1*s2*(Z-J1Km*tau_ss)/(km+m)
T2          = -(Qoo*rho_sd+Poo*tau_sd)*rinf
rho_sod     = (T1+T2)/(1-rinf2)

!print*, ' T1 ', minval(T1), maxval(T1), sum(T1)
!print*, ' T2 ', minval(T2), maxval(T2), sum(T2)
!print*, ' rho_sod ', minval(rho_sod), maxval(rho_sod), sum(rho_sod)


rho_sos     = w * sum(Pso(1:nl))*iLAI
rho_so      = rho_sod + rho_sos

Pso2w       = Pso(nl+1)


!print*, ' rho_sos ', minval(rho_sos), maxval(rho_sos), sum(rho_sos)
!print*, ' rho_so ', minval(rho_so), maxval(rho_so), sum(rho_so)
!print*, ' Pso2w ', Pso2w 

!% Analytical rso following SAIL
rso         = rho_so + rs * Pso2w + &
             &  ((tau_sd+tau_ss*rs*rho_dd)*tau_oo+(tau_sd+tau_ss)*tau_do)*rs/denom

!print*, ' rso ', minval(rso), maxval(rso), sum(rso)

!% Extract MODTRAN atmosphere parameters at the SCOPE wavelengths
t1m  = atmoM(:,1)
t3m  = atmoM(:,2)
t4m  = atmoM(:,3)
t5m  = atmoM(:,4)
t12m = atmoM(:,5)
t16m = atmoM(:,6)

!print*, ' t1m ', minval(t1m), maxval(t1m), sum(t1m)
!print*, ' t3m ', minval(t3m), maxval(t3m), sum(t3m)
!print*, ' t4m ', minval(t4m), maxval(t4m), sum(t4m)
!print*, ' t5m ', minval(t5m), maxval(t5m), sum(t5m)
!print*, ' t12m ', minval(t12m), maxval(t12m), sum(t12m)
!print*, ' t16m ', minval(t16m), maxval(t16m), sum(t16m)

!% radiation fluxes, downward and upward (these all have dimenstion [nwl]
!% first calculate hemispherical reflectances rsd and rdd according to SAIL
!% these are assumed for the reflectance of the surroundings
!% rdo is computed with SAIL as well

denom       = 1.-rs*rho_dd

!print*, ' denom ', minval(denom), maxval(denom), sum(denom)

!% SAIL analytical reflectances

rsd     = rho_sd + (tau_ss + tau_sd)*rs*tau_dd/denom
rdd     = rho_dd + tau_dd*rs*tau_dd/denom

rdo     = rho_do + (tau_oo + tau_do)*rs*tau_dd/denom

!print*, 'rsd ', minval(rsd), maxval(rsd), sum(rsd)
!print*, 'rdd ', minval(rdd), maxval(rdd), sum(rdd)
!print*, 'rdo ', minval(rdo), maxval(rdo), sum(rdo)


!% assume Fd of surroundings = 0 for the momemnt
!% initial guess of temperature of surroundings from Ta;

!Fd      = zeros(nw ', minval(bfli), maxval(bfli), sum(bfli)l,1);
Fd       = 0.
!Ls      = Planck(wl,atmo.Ta+273.15);
CALL  Planck(nwl,wl,Ta+273.15,Ls)

!print*, 'Ta ', Ta
!print*, 'Ls ', minval(Ls), maxval(Ls), sum(Ls)

! Solar and sky irradiance using 6 atmosperic functions
Esun_   = pi*t1m*t4m
Esky_   = pi/(1-t3m*rdd)*(t1m*(t5m+t12m*rsd)+Fd+(1-rdd)*Ls*t3m+t16m)

!print*, 'Esun_ ', minval(Esun_), maxval(Esun_), sum(Esun_)
!print*, 'Esky_ ', minval(Esky_), maxval(Esky_), sum(Esky_)


! ----------------------------------------------------------------------

!% fractional contributions of Esun and Esky to total incident radiation in
!% optical and thermal parts of the spectrum
![fEsuno,fEskyo,fEsunt,fEskyt]          = deal(0*Esun_);   %initialization
fEsuno  = 0.
fEskyo  = 0.
fEsunt  = 0.
fEskyt  = 0.

!J_o             = wl<3000                          %find optical spectrum
!print*, ' size wl ', size(wl) , ' min max wl ', minval(wl), maxval(wl)

DEALLOCATE (mask)
IF (.NOT. ALLOCATED(mask)) ALLOCATE(mask(size(wl)))
mask = (wl<3000.)
IF (.NOT.ALLOCATED(J_o)) ALLOCATE(J_o(count(mask)))
CALL mask_ind(size(mask),mask,J_o)

!print*, ' size J_o ', size(J_o), ' size mask ', size(mask)

!Esunto          = 0.001 * Sint(Esun_(J_o),wl(J_o)) ; %Calculate optical sun fluxes (by Integration), including conversion mW >> W 
CALL Sint(size(J_o),Esun_(J_o),wl(J_o),ap)
Esunto       =  0.001 * ap 

!Eskyto          = 0.001 * Sint(Esky_(J_o),wl(J_o)); %Calculate optical sun fluxes (by Integration)
CALL Sint(size(J_o),Esky_(J_o),wl(J_o),ap)
Eskyto          = 0.001 *ap 

Etoto           = Esunto + Eskyto     !;  %Calculate total fluxes
fEsuno(J_o)     = Esun_(J_o)/Etoto    !;  %fraction of contribution of Sun fluxes to total light
fEskyo(J_o)     = Esky_(J_o)/Etoto    !;  %fraction of contribution of Sky fluxes to total light

DEALLOCATE(mask)

!print*, 'J_o ', minval(J_o), maxval(J_o), sum(J_o)
!print*, 'Esunto ', Esunto
!print*, 'Eskyto ', Eskyto
!print*, 'Etoto ', Etoto
!print*, 'fEsunoJ_o ', minval(fEsuno(J_o)), maxval(fEsuno(J_o)), sum(fEsuno(J_o))
!print*, 'fEskyo J_o', minval(fEskyo(J_o)), maxval(fEskyo(J_o)), sum(fEskyo(J_o)) 

!J_t             = wl>=3000;                         %find thermal spectrum
!print*, ' wl in rtmo ', minval(wl), maxval(wl)
IF (.NOT. ALLOCATED(mask)) ALLOCATE(mask(size(wl)))
mask = (wl>=3000.)
IF (.NOT.ALLOCATED(J_t)) ALLOCATE(J_t(count(mask)))
CALL mask_ind(size(mask),mask,J_t)
DEALLOCATE(mask)

!Esuntt          = 0.001 * Sint(Esun_(J_t),wl(J_t)); %Themal solar fluxes
CALL Sint(size(J_t),Esun_(J_t),wl(J_t),ap)
Esuntt          = 0.001 * ap 

!Eskytt          = 0.001 * Sint(Esky_(J_t),wl(J_t)); %Thermal Sky fluxes
CALL Sint(size(J_t),Esky_(J_t),wl(J_t),ap)
Eskytt          = 0.001 *ap

Etott           = Eskytt + Esuntt  !;                  %Total
fEsunt(J_t)     = Esun_(J_t)/Etott !;                 %fraction from Esun 
fEskyt(J_t)     = Esky_(J_t)/Etott !;                 %fraction from Esky

!print*, 'J_t ', minval(J_t), maxval(J_t), sum(J_t)
!print*, 'Esuntt ', Esuntt
!print*, 'Eskytt ', Eskytt
!print*, 'Etott ', Etott
!print*, 'fEsunt J_t ', minval(fEsunt(J_t)), maxval(fEsunt(J_t)), sum(fEsunt(J_t))
!print*, 'fEskyt J_t', minval(fEskyt(J_t)), maxval(fEskyt(J_t)), sum(fEskyt(J_t))


if (Rin.ne. -999) then 
    Esun_(J_o) = fEsuno(J_o)*Rin
    Esky_(J_o) = fEskyo(J_o)*Rin
    Esun_(J_t) = fEsunt(J_t)*Rli
    Esky_(J_t) = fEskyt(J_t)*Rli


!print*, 'Esun_ J_o ', minval(Esun_(J_o)), maxval(Esun_(J_o)), sum(Esun_(J_o))
!print*, 'Esky_ J_o ', minval(Esky_(J_o)), maxval(Esky_(J_o)), sum(Esky_(J_o))
!print*, 'Esun_ J_t ', minval(Esun_(J_t)), maxval(Esun_(J_t)), sum(Esun_(J_t))
!print*, 'Esky_ J_t ', minval(Esky_(J_t)), maxval(Esky_(J_t)), sum(Esky_(J_t))

endif 

Eplu_1      = rs*((tau_ss+tau_sd)*Esun_+tau_dd*Esky_)/denom;
Eplu0       = rho_sd*Esun_ + rho_dd*Esky_ + tau_dd*Eplu_1
Emin_1      = tau_sd*Esun_ + tau_dd*Esky_ + rho_dd*Eplu_1
delta1      = Esky_  - rinf*Eplu0
delta2      = Eplu_1 - rinf*Emin_1

!print*, 'Eplu_1 ', minval(Eplu_1), maxval(Eplu_1), sum(Eplu_1)
!print*, 'Eplu0 ', minval(Eplu0), maxval(Eplu0), sum(Eplu0)
!print*, 'Emin_1 ', minval(Emin_1), maxval(Emin_1), sum(Emin_1)
!print*, 'delta1 ', minval(delta1), maxval(delta1), sum(delta1)
!print*, 'delta2 ', minval(delta2), maxval(delta2), sum(delta2)

!% calculation of the fluxes in the canopy
!for i = 1:nwl
do i = 1,nwl
   ! J1kx        = calcJ1(xl,m(i),km,LAI)    !%           [nl]
    CALL calcJ1b(nl+1,xl,m(i),km,LAI,J1kx)
    !J2kx        = calcJ2(xl,m(i),km,LAI)    !%           [nl]
    CALL calcJ2b(nl+1,xl,m(i),km,LAI,J2kx)    !%           [nl]

    F1          = Esun_(i)*J1kx*(sf(i)+rinf(i)*sb(i)) + delta1(i)*exp( m(i)*LAI*xl) !      %[nl]
    F2          = Esun_(i)*J2kx*(sb(i)+rinf(i)*sf(i)) + delta2(i)*exp(-m(i)*LAI*(xl+1)) !  %[nl]
    Emin_(:,i)  = (F1+rinf(i)*F2)/(1-rinf2(i))      !%        [nl,nwl]
    Eplu_(:,i)  = (F2+rinf(i)*F1)/(1-rinf2(i))      !%        [nl,nwl]
end do 

!print*, 'J1kx ', minval(J1kx), maxval(J1kx), sum(J1kx)
!print*, 'J2kx ', minval(J2kx), maxval(J2kx), sum(J2kx)
!print*, 'F1 ', minval(F1), maxval(F1), sum(F1)
!print*, 'F2 ', minval(F2), maxval(F2), sum(F2)
!print*, 'Emin_ ', minval(Emin_), maxval(Emin_), sum(Emin_)
!print*, 'Eplu_ ', minval(Eplu_), maxval(Eplu_), sum(Eplu_)

!% Incident and absorbed solar radiation
!Psun        = 0.001 * Sint(e2phot(wlPAR*1E-9,Esun_(Ipar)),wlPAR)  !;   % Incident solar PAR in PAR units
IF (.NOT. ALLOCATED(tempIpar)) ALLOCATE(tempIpar(size(Ipar)))
CALL e2phot(size(Ipar),wlPAR*1E-9,Esun_(Ipar),tempIpar) 
CALL Sint(size(Ipar),tempIpar,wlPAR,ap)
Psun        = 0.001 * ap 

!Asun        = 0.001 * Sint(Esun_.*epsc,wl);  
CALL Sint(nwl,Esun_*epsc,wl,ap)
Asun = ap *0.001

!Pnsun       = 0.001 * Sint(e2phot(wlPAR*1E-9,Esun_(Ipar).*epsc(Ipar)),wlPAR) !;  %Absorbed solar radiation 
CALL e2phot(size(Ipar),wlPAR*1E-9,Esun_(Ipar)*epsc(Ipar),tempIpar)
CALL Sint(size(Ipar),tempIpar,wlPAR,ap)
Pnsun       = 0.001 *ap
                                                                             !  in PAR range in moles m-2 s-1
!Pnsun_Cab   = 0.001 * Sint(e2phot(wlPAR*1E-9,kClrel(Ipar).*Esun_(Ipar).*epsc(Ipar)),wlPAR) ! % Absorbed solar 
CALL e2phot(size(Ipar),wlPAR*1E-9,kClrel(Ipar)*Esun_(Ipar)*epsc(Ipar),tempIpar)
CALL Sint(size(Ipar),tempIpar,wlPAR,ap)
Pnsun_Cab       = 0.001 *ap

!print*, 'Psun ', Psun 
!print*, 'Asun ', Asun 
!print*, 'Pnsun ', Pnsun 
!print*, 'kClrel ', minval(kClrel), maxval(kClrel)
!print*, 'Ipar ', minval(Ipar), maxval(Ipar)
!print*, 'Pnsun_Cab ', Pnsun_Cab 

!%% 3. outgoing fluxes, hemispherical and in viewing direction, spectrum
!% (CT071105: compared with analytical solution: is ok)
!% hemispherical, spectral
!Lout_       = Eplu_(1,:)';              !           [nwl]

!print*, ' size Eplu ', size(Eplu_,1), size(Eplu_,2),' size Eout_ ', size(Eout_)
Eout_   = Eplu_(1,:)                     !           [nwl]

!print*, 'Eout_ ', minval(Eout_), maxval(Eout_), sum(Eout_)

!% in viewing direction, spectral
!piLoc_      = (vb.*(Emin_(1:nl,:)'*Po(1:nl)) +...
!               vf.*(Eplu_(1:nl,:)'*Po(1:nl)) +...
!               w.*Esun_*sum(Pso(1:nl)))*iLAI;
piLoc_ = ( vb*(MATMUL (transpose(Emin_(1:nl,:)),Po(1:nl))) + & 
        &  vf*(MATMUL (transpose(Eplu_(1:nl,:)),Po(1:nl)))  +& 
        &   w*Esun_*sum(Pso(1:nl)))*iLAI 


!piLos_     = rs.*Emin_(nl+1,:)'*Po(nl+1) + rs.*Esun_*Pso(nl+1);
piLos_  =   rs*reshape(Emin_(nl+1,:),(/nwl/))*Po(nl+1) + rs*Esun_*Pso(nl+1)

 piLo_  = piLoc_ + piLos_         ! %   [nwl]

   Lo_  =  piLo_/pi 

!print*, 'piLoc_', minval(piLoc_), maxval(piLoc_), sum(piLoc_)
!print*, 'piLos_', minval(piLos_), maxval(piLos_), sum(piLos_)
!print*, 'piLo_', minval(piLo_), maxval(piLo_), sum(piLo_)
!print*, 'Lo_', minval(Lo_), maxval(Lo_), sum(Lo_)


!% up and down and hemispherical out, cumulative over wavelenght
!IwlP        = spectral.IwlP;
!IF (.NOT. ALLOCATED(IwlP)) ALLOCATE(IwlP(nwlP))
!IwlT        = spectral.IwlT;

IF (.NOT. ALLOCATED(mask)) ALLOCATE(mask(size(wl)))
mask = ((wlS>=minval(wlT)).and.(wlS<=maxval(wlT)))
IF (.NOT.ALLOCATED(IwlT)) ALLOCATE(IwlT(count(mask)))
CALL mask_ind(size(mask),mask,IwlT)


DO il=1, size(wlP)
 IwlP(il) = il 
END DO 


!Eouto       = 0.001 * Sint(Eout_(IwlP),wlP) !;%  [1] hemispherical out, in optical range (W m-2)
CALL Sint(size(IwlP),Eout_(IwlP),wlP,ap )
Eouto       = 0.001 * ap 

!Eoutt       = 0.001 * Sint(Eout_(IwlT),wlT) !%   [1] hemispherical out, in thermal range (W m-2)
CALL Sint(size(IwlT),Eout_(IwlT),wlT,ap )
Eoutt       = 0.001 *ap 

!print*, 'IwlP', minval(IwlP), maxval(IwlP), sum(IwlP)
!print*, 'IwlT', minval(IwlT), maxval(IwlT), sum(IwlT)
!print*, 'Eouto', Eouto
!print*, 'Eoutt', Eoutt

!%% 4. net fluxes, spectral and total, and incoming fluxes
!% incident PAR at the top of canopy, spectral and spectrally integrated
!P_          = e2phot(wl(Ipar)*1E-9,(Esun_(Ipar)+Esky_(Ipar)));
IF (.NOT. ALLOCATED(P_)) ALLOCATE(P_(size(Ipar)))
CALL e2phot(size(Ipar),wl(Ipar)*1E-6,(Esun_(Ipar)+Esky_(Ipar)),P_)
!P           = .001 * Sint(P_,wlPAR);
CALL Sint(size(Ipar),P_,wlPAR,ap )
P           = .001 * ap

!% total direct radiation (incident and net) per leaf area (W m-2 leaf)
Pdir        = fs * Psun                      !  % [13 x 36]   incident
Rndir       = fs * Asun                      !  % [13 x 36]   net
Pndir       = fs * Pnsun                     !  % [13 x 36]   net PAR
Pndir_Cab   = fs * Pnsun_Cab                 !  % [13 x 36]   net PAR Cab

!print*, 'P_', minval(P_), maxval(P_), sum(P_)
!print*, 'P', P
!print*, 'Pdir', minval(Pdir), maxval(Pdir), sum(Pdir)
!print*, 'Rndir', minval(Rndir), maxval(Rndir), sum(Rndir)
!print*, 'Pndir', minval(Pndir), maxval(Pndir), sum(Pndir)
!print*, 'Pndir_Cab', minval(Pndir_Cab), maxval(Pndir_Cab), sum(Pndir_Cab)


! ---------------------------------------------------

!% canopy layers, diffuse radiation
DO  j = 1,nl
    !% diffuse incident radiation for the present layer 'j' (mW m-2 um-1)
    E_         = .5*(Emin_(j,:) + Emin_(j+1,:)+ Eplu_(j,:)+ Eplu_(j+1,:))

    !% incident PAR flux, integrated over all wavelengths (moles m-2 s-1)
    !Pdif(j)    = .001 * Sint(e2phot(wlPAR*1E-9,E_(Ipar)'),wlPAR) ! % [nl]  including conversion mW >> W
    CALL  e2phot(size(wlPAR),wlPAR*1E-9,E_(Ipar),tempIpar)
    CALL Sint(size(wlPAR),tempIpar,wlPAR,ap )
    Pdif(j)    = .001 * ap

    !% net radiation (mW m-2 um-1) and net PAR (moles m-2 s-1 um-1), per wavelength
    !Rndif_(j,:)         = E_.*epsc'   ! % [nl,nwl]  Net (absorbed) radiation by leaves
    Rndif_(j,:)    = E_*epsc   ! % [nl,nwl]  Net (absorbed) radiation by leaves

   ! Pndif_(j,:)         = .001 *(e2phot(wlPAR*1E-9, Rndif_(j,Ipar)'))' !% [nl,nwl]  Net (absorbed) as PAR photons
    CALL e2phot(size(Ipar),wlPAR*1E-9, Rndif_(j,Ipar),Pndif_(j,:))
    Pndif_(j,:)    =  0.001*Pndif_(j,:)
    

    !Pndif_Cab_(j,:)     = .001 *(e2phot(wlPAR*1E-9,kClrel(Ipar).*Rndif_(j,Ipar)'))' !  % [nl,nwl]  Net (absorbed)
    CALL e2phot(size(Ipar),wlPAR*1E-9,kClrel(Ipar)*Rndif_(j,Ipar),tempIpar)
    Pndif_Cab_(j,:) = .001 *tempIpar
                                                                                   !  as PAR photons by Cab

   ! % net radiation (W m-2) and net PAR (moles m-2 s-1), integrated over all wavelengths
   ! Rndif(j)            = .001 * Sint(Rndif_(j,:),wl)           !;% [nl]  Full spectrum net diffuse flux
    CALL Sint(size(wl),Rndif_(j,:),wl,ap )
    Rndif(j)       = .001 * ap 

   ! Pndif(j)            =        Sint(Pndif_(j,Ipar),wlPAR)     !% [nl] Absorbed PAR
    CALL Sint(size(wlPAR),Pndif_(j,Ipar),wlPAR,Pndif(j))

    !Pndif_Cab(j)        =        Sint(Pndif_Cab_(j,Ipar),wlPAR) !;% [nl] Absorbed PAR by Cab integrated
     CALL Sint(size(wlPAR),Pndif_Cab_(j,Ipar),wlPAR,Pndif_Cab(j))

END DO 

!print*, 'E_', minval(E_), maxval(E_), sum(E_)
!print*, 'Pdif', minval(Pdif), maxval(Pdif), sum(Pdif)
!print*, 'Rndif_', minval(Rndif_), maxval(Rndif_), sum(Rndif_)
!print*, 'Pndif_', minval(Pndif_(:,Ipar)), maxval(Pndif_(:,Ipar)), sum(Pndif_(:,Ipar))
!print*, 'Pndif_Cab', minval(Pndif_Cab_(:,Ipar)), maxval(Pndif_Cab_(:,Ipar)), sum(Pndif_Cab_(5:,Ipar))
!print*, 'Rndif', minval(Rndif), maxval(Rndif), sum(Rndif)
!print*, 'Pndif', minval(Pndif), maxval(Pndif), sum(Pndif)
!print*, 'Pndif_Cab', minval(Pndif_Cab), maxval(Pndif_Cab), sum(Pndif_Cab)

!% soil layer, direct and diffuse radiation
!Rndirsoil   = .001 * Sint(Esun_.*epss,wl)       !   % [1] Absorbed solar flux by the soil
CALL Sint(size(wl),Esun_*epss,wl,ap)
Rndirsoil   = .001 * ap 

!Rndifsoil   = .001 * Sint(Emin_(nl+1,:).*epss',wl) !; % [1] Absorbed diffuse downward flux by the soil (W m-2)
CALL Sint(size(wl),Emin_(nl+1,:)*epss,wl,ap)
Rndifsoil   = .001 * ap 

!print*, ' Rndirsoil ', Rndirsoil
!print*, 'Rndifsoil ', Rndifsoil


 ! ---------------------------------------------------

!% net (n) radiation R and net PAR P per component: sunlit (u), shaded (h) soil(s) and canopy (c),
!% [W m-2 leaf or soil surface um-1]
Rnhc        = Rndif       !  % [nl] shaded leaves or needles
Pnhc        = Pndif       !  % [nl] shaded leaves or needles
Pnhc_Cab    = Pndif_Cab   ! % [nl] shaded leaves or needles

!print*, 'Rnhc', minval(Rnhc), maxval(Rnhc), sum(Rnhc)
!print*, 'Pnhc', minval(Pnhc), maxval(Pnhc), sum(Pnhc)
!print*, 'Pnhc_Cab', minval(Pnhc_Cab), maxval(Pnhc_Cab), sum(Pnhc_Cab)



DO j = 1,nl
      Puc(:,:,j)  = Pdir      + Pdif(j)  !% [13,36,nl] Total fluxes on sunlit leaves or needles
      Rnuc(:,:,j) = Rndir     + Rndif(j) ! % [13,36,nl] Total fluxes on sunlit leaves or needles
      Pnuc(:,:,j) = Pndir     + Pndif(j) ! % [13,36,nl] Total fluxes on sunlit leaves or needles
 Pnuc_Cab(:,:,j)  = Pndir_Cab + Pndif_Cab(j)  !;% [13,36,nl] Total fluxes on sunlit leaves or needles
END DO 


!print*, 'Puc', minval(Puc), maxval(Puc), sum(Puc)
!print*, 'Rnuc', minval(Rnuc), maxval(Rnuc), sum(Rnuc)
!print*, 'Pnuc', minval(Pnuc), maxval(Pnuc), sum(Pnuc)
!print*, 'Pnuc_Cab', minval(Pnuc_Cab), maxval(Pnuc_Cab), sum(Pnuc_Cab)

Rnhs        = Rndifsoil      !       % [1] shaded soil
Rnus        = Rndifsoil + Rndirsoil !  % [1] sunlit soil


!print*, 'Rnhs', Rnhs
!print*, 'Rnus', Rnus

!%% place output in structure rad ... this for help 
!gap.k       = k;
!gap.K       = K;
!gap.Ps      = Ps;
!gap.Po      = Po;

!rad.rsd     = rsd;
!rad.rdd     = rdd;
!rad.rdo     = rdo;
!rad.rso     = rso;

!rad.Esun_   = Esun_ !        % [2162x1 double]   incident solar spectrum (mW m-2um-1)
!rad.Esky_   = Esky_;        % [2162x1 double]   incident sky spectrum (mW m-2 um-1)
!rad.inPAR   = P;            % [1 double]        incident spectrally integrated PAR (moles m-2 s-1)

!rad.fEsuno  = fEsuno;       % [2162x1 double]   normalized spectrum of direct light (optical)
!rad.fEskyo  = fEskyo;       % [2162x1 double]   normalized spectrum of diffuse light (optical)
!rad.fEsunt  = fEsunt;       % [2162x1 double]   normalized spectrum of direct light (thermal)
!rad.fEskyt  = fEskyt;       % [2162x1 double]   normalized spectrum of diffuse light (thermal)

!rad.Eplu_   = Eplu_;        % [61x2162 double]  upward diffuse radiation in the canopy (mW m-2 um-1)
!rad.Emin_   = Emin_;        % [61x2162 double]  downward diffuse radiation in the canopy (mW m-2 um-1)

!rad.Lo_     = Lo_;          % [2162x1 double]   TOC radiance in observation direction (mW m-2 um-1 sr-1)
!rad.Eout_   = Eout_;        % [2162x1 double]   TOC upward radiation (mW m-2 um-1)
!rad.Eouto   = Eouto;        % [1 double]        TOC spectrally integrated upward optical radiation (W m-2)
!rad.Eoutt   = Eoutt;        % [1 double]        TOC spectrally integrated upward thermal radiation (W m-2)

!rad.Rnhs    = Rnhs;         % [1 double]        net radiation (W m-2) of shaded soil 
!rad.Rnus    = Rnus;         % [1 double]        net radiation (W m-2) of sunlit soil
!rad.Rnhc    = Rnhc;         % [60x1 double]     net radiation (W m-2) of shaded leaves
!rad.Rnuc    = Rnuc;         % [13x36x60 double] net radiation (W m-2) of sunlit leaves
!rad.Pnh     = Pnhc;         % [60x1 double]     net PAR (moles m-2 s-1) of shaded leaves
!rad.Pnu     = Pnuc;         % [13x36x60 double] net PAR (moles m-2 s-1) of sunlit leaves
!rad.Pnh_Cab = Pnhc_Cab;     % [60x1 double]     net PAR absorbed by Cab (moles m-2 s-1) of shaded leaves
!rad.Pnu_Cab = Pnuc_Cab;     % [13x36x60 double] net PAR absorbed by Cab (moles m-2 s-1) of sunlit leaves

END SUBROUTINE rtmo 


!%% APPENDIX I functions J1 and J2 (introduced for stable solutions)

!function J1 = calcJ1(x,m,k,LAI)
SUBROUTINE  calcJ1(n,x,m,k,LAI,J1)
IMPLICIT NONE 

! Input variables 
INTEGER, INTENT(IN)              :: n
REAL, INTENT(IN)                 :: k,LAI
REAL, INTENT(IN),DIMENSION(n)    :: m
REAL, INTENT(IN)                 :: x

! Ouput variables 
REAL, INTENT(OUT),DIMENSION(n)   :: J1

WHERE(abs(m-k)>1E-3)  
 !J1 = (exp(m*LAI*x)-exp(k*LAI*x))./(k-m)
 J1 = (exp(m*LAI*x)-exp(k*LAI*x))/(k-m)
ELSEWHERE
! J1 = -.5*(exp(m*LAI*x)+exp(k*LAI*x))*LAI.*x.*(1-1/12*(k-m).^2*LAI^2.*x.^2);
 J1 = -.5*(exp(m*LAI*x)+exp(k*LAI*x))*LAI*x*(1.-1./12.*(k-m)**2.*LAI**2.*x**2.)

ENDWHERE 
return
!END FUNCTION calcJ1
END SUBROUTINE calcJ1


SUBROUTINE  calcJ1b(n,x,m,k,LAI,J1)
IMPLICIT NONE

! Input variables
INTEGER, INTENT(IN)              :: n
REAL, INTENT(IN)                 :: k,LAI
REAL, INTENT(IN),DIMENSION(n)    :: x
REAL, INTENT(IN)                 :: m

! Ouput variables
REAL, INTENT(OUT),DIMENSION(n)   :: J1

IF (abs(m-k)>1E-3) THEN 
 !J1 = (exp(m*LAI*x)-exp(k*LAI*x))./(k-m)
 J1 = (exp(m*LAI*x)-exp(k*LAI*x))/(k-m)
ELSE 
! J1 = -.5*(exp(m*LAI*x)+exp(k*LAI*x))*LAI.*x.*(1-1/12*(k-m).^2*LAI^2.*x.^2);
 J1 = -.5*(exp(m*LAI*x)+exp(k*LAI*x))*LAI*x*(1.-1./12.*(k-m)**2.*LAI**2.*x**2.)
ENDIF 
return
!END FUNCTION calcJ1
END SUBROUTINE calcJ1b


!function J2 = calcJ2(x,m,k,LAI)
SUBROUTINE calcJ2(n,x,m,k,LAI,J2)
IMPLICIT NONE
INTEGER, INTENT(IN)                        :: n
REAL, INTENT(IN)                           :: k,LAI
REAL, INTENT(IN), DIMENSION(n)             :: m
REAL, INTENT(IN)                           :: x
REAL, INTENT(OUT), DIMENSION(n)            :: J2

!J2 = (exp(k*LAI*x)-exp(-k*LAI)*exp(-m*LAI*(1+x)))./(k+m);
J2 = (exp(k*LAI*x)-exp(-k*LAI)*exp(-m*LAI*(1.+x)))/(k+m)

return
!END FUNCTION 
END SUBROUTINE  calcJ2


!function J2 = calcJ2(x,m,k,LAI)
SUBROUTINE calcJ2b(n,x,m,k,LAI,J2)
IMPLICIT NONE
INTEGER, INTENT(IN)                        :: n
REAL, INTENT(IN)                           :: k,LAI
REAL, INTENT(IN), DIMENSION(n)             :: x
REAL, INTENT(IN)                           :: m
REAL, INTENT(OUT), DIMENSION(n)            :: J2

!J2 = (exp(k*LAI*x)-exp(-k*LAI)*exp(-m*LAI*(1+x)))./(k+m);
J2 = (exp(k*LAI*x)-exp(-k*LAI)*exp(-m*LAI*(1.+x)))/(k+m)
return
!END FUNCTION
END SUBROUTINE calcJ2b


!%% APPENDIX II function volscat

!function [chi_s,chi_o,frho,ftau]    =   volscat(tts,tto,psi,ttli)
SUBROUTINE volscat(tts,tto,psi,ttli,chi_s,chi_o,frho,ftau)

IMPLICIT NONE 
! Parmaters 
!global deg2rad nli
REAL, PARAMETER                    ::     pi  = 3.1415926535879
REAL, PARAMETER                    :: deg2rad = pi/180.
REAL, PARAMETER                    :: nli =13

!Input variables 
REAL, INTENT(IN)                    :: tts,tto,psi
INTEGER, INTENT(IN),DIMENSION(nli)  :: ttli

! Output variables 
REAL, INTENT(out), DIMENSION(nli)  :: chi_s,chi_o
REAL, INTENT(OUT), DIMENSION(nli)  :: frho,ftau 

! Local variables 
REAL, DIMENSION(nli)               :: psi_rad, Totl,zero   
REAL, DIMENSION(nli)               :: cos_psi, cos_ttli, sin_ttli  
REAL                               :: cos_tts,sin_tts, cos_tto
REAL                               :: sin_tto
REAL, DIMENSION(nli)               :: Cs, Ss,Co,So 
REAL, DIMENSION(nli)               :: As, Ao, bts, bto 
REAL, DIMENSION (nli)              :: delta1, delta2
REAL, DIMENSION(nli)               :: bt1, bt3 ,bt2
REAL, DIMENSION(nli)               :: T1,T2, Tot
REAL, DIMENSION(nli)               :: Jmin, Jplus 

!%Volscatt version 2.
!%created by W. Verhoef
!%edited by Joris Timmermans to matlab nomenclature.
!% date: 11 February 2008
!%tts    [1]         Sun            zenith angle in degrees
!%tto    [1]         Observation    zenith angle in degrees
!%psi    [1]         Difference of  azimuth angle between solar and viewing position
!%ttli   [ttli]      normal of upperside of leaf

!psi_rad         = psi*deg2rad*ones(nli,1)
psi_rad= 1. 
psi_rad         = psi*deg2rad*psi_rad 

cos_psi         = cos(psi*deg2rad)                 !%   cosinus of relative azimuth   angle

cos_ttli        = cos(ttli*deg2rad)                !%   cosinus of normal of upperside of leaf
sin_ttli        = sin(ttli*deg2rad)                !%   sinus   of normal of upperside of leaf

cos_tts         = cos(tts*deg2rad)                 !%   cosinus of sun      zenith    angle
sin_tts         = sin(tts*deg2rad)                 !%   sinus   of sun      zenith    angle


cos_tto         = cos(tto*deg2rad)                 !%   cosinus of observer zenith    angle
sin_tto         = sin(tto*deg2rad)                 !%   sinus   of observer zenith    angle



Cs              = cos_ttli*cos_tts                 !%   p305{1}
Ss              = sin_ttli*sin_tts                 !%   p305{1}                          

Co              = cos_ttli*cos_tto                 !%   p305{1}
So              = sin_ttli*sin_tto                 !%   p305{1}


As              = max(Ss,Cs)
Ao              = max(So,Co)


bts             = acos(-Cs/As)                    !%   p305{1}
bto             = acos(-Co/Ao)                    !%   p305{2}


chi_o           = 2/pi*((bto-pi/2)*Co + sin(bto)*So)
chi_s           = 2/pi*((bts-pi/2)*Cs + sin(bts)*Ss)
 
delta1          = abs(bts-bto)                    ! %   p308{1}
delta2          = pi-abs(bts + bto - pi)          ! %   p308{1}

Tot             = psi_rad + delta1 + delta2        !%   pag 130{1}

bt1             = min(psi_rad,delta1)
bt3             = max(psi_rad,delta2)
bt2             = Tot - bt1 - bt3

T1              = 2.*Cs*Co + Ss*So*cos_psi
T2              = sin(bt2)*(2*As*Ao + Ss*So*cos(bt1)*cos(bt3))

Jmin            = ( bt2)*T1 - T2
Jplus           = (pi-bt2)*T1 + T2

frho            =  Jplus/(2*pi**2)
ftau            = -Jmin /(2*pi**2)

!% pag.309 wl-> pag 135{1}
zero=0.
frho            = max(zero,frho)
ftau            = max(zero,ftau)

return

END SUBROUTINE 

!%% APPENDIX III function e2phot

SUBROUTINE  e2phot(n,lambda,E,sortie)
!%molphotons = e2phot(lambda,E) calculates the number of moles of photons
!%corresponding to E Joules of energy of wavelength lambda (m)
IMPLICIT NONE

! Paramters
!global A
REAL, PARAMETER                  :: A = 6.02214E23       ! % [mol-1] Constant of Avogadro

! Input variables
INTEGER , INTENT(IN)             :: n
REAL, INTENT(IN), DIMENSION(n)   :: lambda,E

! Ouput varaibales
REAL, INTENT(OUT),DIMENSION(n)   :: sortie

! Local variables
REAL, DIMENSION(n)               :: em, photons

CALL ephoton(n,lambda,em)
photons  = E/em
sortie   = photons/A

RETURN
END SUBROUTINE e2phot


SUBROUTINE  ephoton(n,lambda,res)
!%E = phot2e(lambda) calculates the energy content (J) of 1 photon of
!%wavelength lambda (m)
IMPLICIT NONE

! Input variables
INTEGER, INTENT(IN)             :: n
REAL, INTENT(IN),DIMENSION(n)   :: lambda

! Ouput variables
REAL, INTENT(OUT),DIMENSION(n)  :: res

! Local variables
REAL                            :: h,c

h        = 6.6262E-34                   !% [J s]         Planck's constant'
c        = 3E8                         ! % [m s-1]       speed of light
res = h*c/lambda                 ! % [J]           energy of 1 photon

RETURN 
END SUBROUTINE  ephoton


SUBROUTINE Planck(npt,wl,Tb,Lb)

!function Lb = Planck (wl,Tb,em)
IMPLICIT NONE

! Input variables 
INTEGER, INTENT(IN)                    :: npt 
REAL, INTENT(IN)                       :: Tb
REAL, DIMENSION(npt), INTENT(IN)       :: wl

! Output variables 
REAL, DIMENSION(npt), INTENT(OUT)      :: Lb

!Local 
REAL, DIMENSION(npt)                  :: em
REAL                                  :: c1, c2

    c1 = 1.191066e-22
    c2 = 14388.33

    em = 1.  
    Lb = em* c1*(wl*1e-9)**(-5.)/(exp(c2/(wl*1e-3*Tb))-1.)

END SUBROUTINE Planck 

END MODULE  mo_rtmo 


