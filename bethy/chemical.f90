MODULE chemical 

USE fluo_func
USE mo_constants

IMPLICIT NONE

! function [biochem_out] = biochemical(biochem_in)
!%
!% Date: 	21 Sep 2012
!% Update:   20 Feb 2013
!%   
!% Authors: 	Joe Berry and Christiaan van der Tol, contributions of others.
!% Sources: 	
!%           Farquhar et al. 1980, Collatz et al (1991, 1992).
!%
!% This function calculates:
!%    - stomatal resistance of a leaf or needle (s m-1)
!%    - photosynthesis of a leaf or needle (umol m-2 s-1)
!%    - fluorescence of a leaf or needle (fraction of fluor. in the dark)
!%
!% Usage:
!% function [A,Cs,eb,f,rcw] = biochemical(C,Cs,Q,T,ea,eb,O,p,Vcmo,gcparam,Type,tempcor,ra,Tparams,Rdparam,stressfactor)
!% the function was tested for Matlab 7.2.0.232 (R2006a)
!%
!% Calculates net assimilation rate A, fluorescence F using !biochemical model
!%
!% Input (units are important):
!% C         % [umol m-3]            conc. of CO2 in the ambient air
!% Cs        % [umol m-3]            initial estimate of conc. of CO2 in the
!%                                   ...bounary layer of the leaf
!% Q         % [umol photons m-2 s-1]net radiation, PAR
!% T         % [oC or K]             leaf temperature
!% eb        % [hPa]                 intial estimate of the vapour pressure in leaf boundary layer
!% O         % [mmol m-3]            concentration of O2 (in the boundary
!%                                   ...layer, but no problem to use ambient)
!% p         % [hPa]                 air pressure
!% Vcmo      % [umol/m2/s]           maximum carboxylation capacity
!% m         % []                    Ball-Berry coefficient 'm' for stomatal regulation
!% Type      % []                    text parameter, either 'C3' !for C3 or any
!%                                   ...other text for C4
!% tempcor   % []                    boolean (0 or 1) whether or not
!%                                   ...temperature correction to !Vcmax has to be applied.
!% ra        % [s m-1]               aerodynamic + boundary layer !resitance for the
!%                                   ...calculation of boundary !layer concentrations
!% Tsparams  % [],[],[K],[K],[K]     vector of 5 temperature !correction
!%                                   parameters, look in !!spreadsheet of PFTs. Only if tempcor=1, otherwise use
!%                                   dummy values
!% Rdparam   % []                    respiration as fraction of !Vcmax
!% stressfactor []                   optional input: stress !factor to reduce Vcmax (for
!%                                   example soil moisture, leaf age). Default value = 1.
!%
!% Note: always use the prescribed units. Temperature can be either oC or K
!% Note: input can be single numbers, vectors, or n-dimensional
!% matrices
!%
!% Output:
!% A         % [umol/m2/s]           net assimilation rate of the !leaves
!% Cs        % [umol/m3]             CO2 concentration in the !boundary layer!
!% eb        % [hPa]                 vapour pressure in the !boundary layer
!% fs_fo0    % []                    fluorescence as fraction of dark
!%                                   ...adapted
!% rcw       % [s m-1]               stomatal resistance
!% qE        % []                    non photochemical quenching
!% fs        % []                    fluorescence as fraction of !PAR
!% Ci        % [umol/m3]             internal CO2 concentration
!% Kn        % []                    rate constant for excess !heat!
!%                                   ...dissipation in case of NPQ
!% fo0       % []                    unstressed, dark adapted fluorescence (fraction of aPAR)
!% fm0       % []                    unstressed, light saturated fluorescence (fraction of aPAR)
!% fo        % []                    dark adapted fluorescence (fraction of aPAR)
!% fm        % []                    light saturated fluorescence (fraction of aPAR)
!% qQ        % []                    photochemical quenching

CONTAINS 

SUBROUTINE biochemical(npts,Q,T1,Csi,ea,Oa,p,akc,ako,Vcmo,option,Ag,A,eta,rcw,Ci)

USE fluo_param, ONLY : Rdparam,Tparams,m,tempcor,stressfactor
USE fluo_param, ONLY :  Kcopt, Koopt,Kf,Kd,Kpc,atheta

IMPLICIT NONE 

! Input variables
INTEGER,INTENT(IN)                   :: npts 
REAL                                 :: Vcmo,Oa,p,ea
REAL, DIMENSION(npts),INTENT(IN)     :: Q,T1
REAL, DIMENSION(npts), INTENT(IN)    :: Csi
INTEGER ,INTENT(IN)                  :: option
REAL, INTENT(IN)                     :: akc,ako


! Output variables 
REAL, DIMENSION(npts),INTENT(OUT)    :: A,eta,rcw, Ag
REAL, DIMENSION(npts),INTENT(OUT)    :: Ci
REAL, DIMENSION(npts)                :: qE,Vcmax

! Local variables 
REAL, DIMENSION(npts)                :: Cs, O,eb
REAL, DIMENSION(npts)                :: fs
REAL, DIMENSION(npts)                :: T,qt,TH,TL
REAL, DIMENSION(npts)                :: Kc,Ko,Rd,kpopt,kp
REAL, DIMENSION(npts)                :: spfy,gam,Vc,Vs
REAL, DIMENSION(npts)                :: RH,es,s,Je,Ve,V,Ja
REAL, DIMENSION(npts)                :: V0,V1
REAL, DIMENSION(npts)                :: ps,qQ,fm,Kn 
REAL, DIMENSION(npts)                :: effcon,a1,a2,fo
REAL, DIMENSION(npts)                :: temp
REAL                                 :: po0 
REAL                                 :: Kcopt1,Koopt1 
REAL                                 :: fo0,fm0
REAL                                 :: Tref,slti,shti,Thl,Thh,Trdm 
INTEGER                              :: C4
LOGICAL, SAVE 	     		     :: testa=.TRUE.

! Type of plant             C4: 0 (C3 plant), 1 (C4 plant)
C4 = option 

!%% input
!Rdparam       = biochem_in.Rdparam  ! %[umol/m2/s]        dark respiration rate at 25 oC as fraction of Vcmax
!Tparams       = biochem_in.Tparams !  %[]    temperature sensitivities of Vcmax, etc (dummy values of applTcorr was selected above)
!Cs            = biochem_in.Cs
!Q             = biochem_in.Q
!T             = biochem_in.T
!eb            = biochem_in.eb
!O             = biochem_in.O
!p             = biochem_in.p
!Vcmo          = biochem_in.Vcmo
!m             = biochem_in.m
!Type          = biochem_in.Type
!tempcor       = biochem_in.tempcor
!stressfactor  = biochem_in.stressfactor

! Reading the input data
!open(3,file='bio_in.dat',status='old')
!read(3,'(f18.10)')Rdparam 
!read(3,'(f18.10)')Tparams 
!read(3,'(f18.10)')Cs
!read(3,'(f18.10)')Q
!read(3,'(f18.10)')T
!read(3,'(f18.10)')eb
!read(3,'(f18.10)')O
!read(3,'(f18.10)')p
!read(3,'(f18.10)')Vcmo
!read(3,'(f18.10)')m
!read(3,'(a18)')type1
!read(3,'(f18.10)')tempcor
!read(3,'(f18.10)')stressfactor
!close(3)
!print*, ' tempcor in biocchemical ', tempcor 
!print*, ' T1 ', MINVAL(T1), maxval(T1), sum(T1)

WHERE (T1 <100.)
T = T1 + 273.15  ! % convert temperatures to K if not already
ENDWHERE

!rhoa        = 1.2047          ! % [kg m-3]       specific mass of air
!Mair        = 28.96           ! % [g mol-1]      molecular mass of dry air
 
!Kcopt       = 350            !  % [ubar]     kinetic coefficient for CO2 (Von Caemmerer and Furbank, 1999)
!Koopt       = 450 !% [mbar]        kinetic coeeficient for  O2 (Von Caemmerer and Furbank, 1999)
!Kf          = 0.05 ! % []            rate constant for fluorescence
!Kd          = 0.95 !             % []            rate constant for thermal deactivation at Fm
!Kpc         = 4.0 !  % []            rate constant for photochemisty

kpopt       = Vcmo/56*1E6 ! % []     PEPcase rate constant for C02, 
                          ! used here: Collatz et al: Vcmo = 39 umol m-1 s-1; kp = 0.7 mol m-1 s-1.
!print*, ' option ', option 
!print*, ' Kcopt ', Kcopt 
!print*, ' Koopt ', Koopt 
!print*, ' Kf ', Kf, ' Kd', Kd, ' Kpc ', Kpc 
!print*, ' atheta ', atheta
!print*, ' Vcmo ', Vcmo
!print*, ' Kpopt ', minval(Kpopt), maxval(Kpopt), sum(Kpopt)
!print*, ' T ', minval(T), maxval(T), sum(T)
!print*, ' tempcor ', tempcor 
!print*, ' Tparams ', Tparams
!print*, ' Oa ', Oa
!print*, ' p ', p

!%% temperature definitions
Tref        = 25.+273.15 !    % [K]           absolute temperature at 25 oC
slti        = Tparams(1)
shti        = Tparams(2)
Thl         = Tparams(3)
Thh         = Tparams(4)
Trdm        = Tparams(5)

!%% convert all to bar
Cs          = Csi * 1e-6 * p *1E-3

eb          = ea 

!print*, 'Cs ', minval(Cs), maxval(Cs), sum(Cs) 
!print*, 'eb ', minval(eb), maxval(eb), sum(eb) 

!O           = O1 * 1e-3 * p *1E-3 * ~strcmp('C4',Type) !;       % forced to be zero for C4 vegetation (this is a trick to prevent oxygenase)
!print*, ' Oa ', Oa
O           = Oa * 1e-3 * p *1E-3    ! For C3 plant 
!print*, ' O ', minval(O), maxval(O), sum(O)

IF (C4) O   = 0.            ! C4 plant 

Kcopt1       = akc * 1e-6
Koopt1       = ako * 1e-3

!print*, ' CHEM Kcopt ', Kcopt, ' akc ', akc 
!print*, ' CHEM Koopt ', Koopt, ' ako ', ako 


!%% temperature corrections
qt          = 0.1 * (T-Tref) * tempcor
TH          = 1 + tempcor* exp(shti * (T   -Thh))
TL          = 1 + tempcor* exp(slti * (Thl -T))

!print*, 'qt ', minval(qt), maxval(qt), sum(qt) 
!print*, 'TH ', minval(TH), maxval(TH), sum(TH) 
!print*, 'TL ', minval(TL), maxval(TL), sum(TL) 

Kc          = Kcopt1 * 2.1**qt
Ko          = Koopt1 * 1.2**qt
kp          = kpopt*1.8**qt

Rd          = Rdparam * Vcmo * 1.8**qt/(1+exp(1.3*(T-Trdm)))

!print*, ' Kc ', minval(Kc), maxval(Kc), sum(Kc) 
!print*, ' Ko ', minval(Ko), minval(Ko), maxval(Ko) 
!print*, ' Kp ', minval(Kp), minval(Kp), maxval(Kp) 
!print*, ' Rd ', minval(Rd), minval(Rd), maxval(Rd) 


! C4 plant 
 IF (C4) THEN 
 Vcmax   =   Vcmo * 2.1**qt/(TL*TH) * 0.5
 ELSE 
! C3 plant
 Vcmax   =   Vcmo * 2.1**qt/TH * 0.5
 ENDIF

!print*, ' Vcmax ', minval(Vcmax), minval(Vcmax), maxval(Vcmax) 

spfy        = 2600 * 0.75 **qt ! % This is, in theory, Vcmax/Vomax.*Ko./Kc, but used as a separate parameter

!print*, ' spfy ', minval(spfy), minval(spfy), maxval(spfy) 

!%% calculation of potential electron transport rate
po0         = Kpc/(Kf+Kd+Kpc)  !  % dark photochemistry fraction (Genty et al., 1989)
!print*,'Kpc ', Kpc 
!print*,'Kf ', Kf 
!print*,'Kd ', Kd 
!print*, 'po0 ', po0
!print*, ' Q ', minval(Q), maxval(Q), sum(Q) 

Je          = 0.5*po0 * Q   !         % electron transport rate

!print*, 'Je ', minval(Je), maxval(Je), sum(Je) 

!%% calculation of the intersection of enzyme and light limited curves
!% this is the original Farquhar model
gam         = 0.5 /spfy *O  !; %[bar]       compensation point [bar]

!print*, 'gam ', minval(gam), maxval(gam), sum(gam) 

!%% calculation of internal CO2 concentration, photosynthesis
!RH          = eb./satvap(T-273.15)  
CALL satvap(npts,(T-273.15),es,s)
RH = eb/es 

!print*, ' es ' , minval(es), maxval(es), sum(es)
!print*, ' s ' , minval(s), maxval(s), sum(s)
!print*, ' RH ', minval(RH), maxval(RH), sum(RH) 

!Ci          = max(.2*Cs,Cs.*(1-1.6./(m.*RH)))
Ci          = max(.2*Cs,Cs*(1.-1.6/(m*RH)))

!print*, ' O ', minval(O), maxval(O),sum(O)
!print*, ' Ci ', minval(Ci), maxval(Ci), sum(Ci) 
!print*, ' m ', m


IF (C4) THEN                
! C4 plant 
        Vc          = Vcmax
        Vs          = kp*Ci
        effcon      = 1/6.                   ! % Berry and Farquhar (1978): 1/0.167
ELSE 
! C3 plant
        Vc          = Vcmax*(Ci-gam)/((Kc * (1.+O/Ko)) + Ci)
        Vs          = Vcmo/2.* 1.8**qt    
        effcon      = 1./(4.5+10.5*gam/Ci)! ; % Von Caemmerer and Farquhar (1981)

!print*, ' C3 ' 
ENDIF 

!print*, ' gam ', minval(gam ), maxval(gam), sum(gam)
!print*, ' Kc ', minval(Kc), maxval(Kc), sum(Kc) 
!print*, ' Ko ', minval(Ko), maxval(Ko), sum(Ko) 
!print*, ' Je ', minval(Je), maxval(Je), sum(Je) 
!print*, ' effcon ', minval(effcon), maxval(effcon), sum(effcon) 


Ve          = Je*(Ci-gam)/(Ci+2*gam) * effcon

!print*, ' Ve ', minval(Ve), maxval(Ve), sum(Ve) 
!print*, ' Vc ', minval(Vc), maxval(Vc), sum(Vc) 
!print*, ' Vs ', minval(Vs), maxval(Vs), sum(Vs) 


! To see 
![a1,a2]     = abc(atheta,-(Vc+Ve),Vc.*Ve)
temp = atheta
CALL abc(npts,temp,-(Vc+Ve),Vc*Ve,a1,a2)

!print*, ' a1 ', minval(a1), maxval(a1), sum(a1)
!print*, ' a2 ', minval(a2), maxval(a2), sum(a2)


!V           = min(a1,a2)*(Ci>gam) + max(a1,a2)*(Ci<=gam)
!V0 = 0.
!V1 = 0.
WHERE (Ci >gam)
   V  = min(a1,a2)
ELSEWHERE
   V  = max(a1,a2) 
ENDWHERE 


!print*, ' V ', minval(V), maxval(V), sum(V) 
!print*, ' min a1 a2', min(a1,a2)
!print*, ' min ', min(a1,a2)*Ci
!print*, ' max ', max(a1,a2)*Ci

! To see 
![a1,a2]     = abc(0.98,-(V+Vs),V.*Vs)
temp = 0.98
CALL abc(npts,temp,-(V+Vs),V*Vs,a1,a2)

Ag      = min(a1,a2)
A       = Ag - Rd
Ja      = Ag /((Ci-gam)/(Ci+2*gam))/effcon  !;  % actual electron transport rate

!print*,' a1 ', minval(a1), maxval(a1), sum(a1)
!print*,' a2 ', minval(a2), maxval(a2),sum(a2)
!print*, ' Ag ', minval(Ag), maxval(Ag), sum(Ag) 
!print*, ' A ', minval(A), maxval(A), sum(A) 
!print*, ' Ja ', minval(Ja), maxval(Ja), sum(Ja) 



!%% calculation of stomatal resistance
rcw  = 0.625*(Cs-Ci)/A *rhoa/Mair*1E3    * 1e6 / p * 1E3
!print*, ' rcw ', minval(rcw), maxval(rcw), sum(rcw) 
WHERE (A <=0)
rcw  = 0.625*1E6
ENDWHERE

!print*, ' rcw ', minval(rcw), maxval(rcw), sum(rcw) 


!%% fluorescence (Replace this part by Magnani or other model if needed
!print*,' po0 in biochemical = ',po0
!print*,' Ja in biochemical = ',Ja
!print*,' Je in biochemical = ',Je
!print*,' po0 shape in biochemical = ',shape(po0)
!print*,' Ja shape in biochemical = ',shape(Ja)
!print*,' Je shape in biochemical = ',shape(Je)
 WHERE (Ja == 0.) 
 ps = 0.
 ELSEWHERE
 ps   =   po0*Ja/Je         !          % this is the photochemical yield
 ENDWHERE


!ps          = po0*Ja/Je     !          % this is the photochemical yield

!print*, ' ps ', minval(ps), maxval(ps), sum(ps) 

! To see 
![eta,qE,qQ,fs,fo,fm,fo0,fm0,Kn]    = TB12(ps,Kp,Kf,Kd)
CALL TB12 (npts,ps,Kpc,Kf,Kd,eta,qE,qQ,fs,fo,fm,fo0,fm0,Kn)

!print*, ' eta ', minval(eta), maxval(eta), sum(eta)
!print*, 'qE', minval(qE), maxval(qE), sum(qE) 
!print*, 'qQ ', minval(qQ), maxval(qQ), sum(qQ)
!print*, ' fs ', minval(fs), maxval(fs), sum(fs) 
!print*, ' fo ', minval(fo), maxval(fo), sum(fo) 
!print*, ' fm ', minval(fm), maxval(fm), sum(fm) 
!print*, 'fo0', fo0
!print*, 'fm0 ', fm0
!print*, ' Kn  ', minval(Kn), maxval(Kn), sum(Kn)



!%% convert back to ppm
Ci          = Ci*1e6/ p * 1E5

!print*, ' Ci ', minval(Ci), maxval(Ci), sum(Ci) 

!%% Collect outputs
!biochem_out.A       = A
!biochem_out.Ci      = Ci
!biochem_out.eta     = eta
!biochem_out.rcw     = rcw
!biochem_out.qE      = qE
!biochem_out.fs      = fs
!biochem_out.Kn      = Kn
!biochem_out.fo0     = fo0
!biochem_out.fm0     = fm0
!biochem_out.fo      = fo
!biochem_out.fm      = fm
!biochem_out.qQ      = qQ
!biochem_out.Vcmax   = Vcmax

!IF (testa) THEN
!	print *, 'biochemical was called'
!	testa = .FALSE.
!END IF 
!print *, ' biochemical was called '

END SUBROUTINE biochemical  
!%%% end of function biochemical

!%% abc formula
!function [x2,x1] = abc(a,b,c)
!if a == 0
!    x1      = -c./b;
!    x2      = x1;
!else
!    x1      = (-b+sqrt(b.^2-4.*a.*c))./(2.*a);
!    x2      = (-b-sqrt(b.^2-4.*a.*c))./(2.*a);
!end
!return;
!%%% end of abc formula

SUBROUTINE abc(npts,a,b,c,x1,x2)
!%% abc formula
!function [x2,x1] = abc(a,b,c)
IMPLICIT NONE

! Inputs
INTEGER, INTENT(IN)                  :: npts
REAL, DIMENSION(npts), INTENT(IN)    :: a,b,c

! Outputs
REAL, INTENT(OUT),DIMENSION(npts)    :: x1,x2

! Local
REAL, DIMENSION(npts)                :: gc

WHERE (a == 0)
  x1  = x1
  x2  = -c/b
ELSE WHERE
  x1  = (-b-sqrt(b**2-4.*a*c))/(2.*a)
  x2  = (-b+sqrt(b**2-4.*a*c))/(2.*a)
END WHERE


RETURN 
END SUBROUTINE abc 


!SUBROUTINE satvap(npts,T,es,s)
!function [es,s] = satvap(T)
!%% function [es,s]= satvap(T)
!% Author: Dr. ir. Christiaan van der Tol
!% Date: 2003
!%
!% calculates the saturated vapour pressure at
!% temperature T (degrees C)
!% and the derivative of es to temperature s (kPa/C)
!% the output is in mbar or hPa. The approximation formula that is used is:
!% es(T) = es(0)*10^(aT/(b+T));
!% where es(0) = 6.107 mb, a = 7.5 and b = 237.3 degrees C
!% and s(T) = es(T)*ln(10)*a*b/(b+T)^2

! Inputs
!INTEGER, INTENT(IN)               :: npts
!REAL, DIMENSION(npts), INTENT(IN) :: T

! Outputs
!REAL, INTENT(OUT),DIMENSION(npts) :: es,s

! Local
!REAL                              :: a,b

!%% constants
!a  = 7.5
!b  = 237.3              !; %degrees C

!%% calculations
!es  = 6.107*10.**(7.5*T/(b+T))
!s   = es*log(10.)*a*b/(b+T)**2

!RETURN
!END SUBROUTINE satvap

SUBROUTINE TB12 (npt,ps,Kp,Kf,Kd, eta,qE,qQ,fs,fo,fm,fo0,fm0,Kn)
!function [eta,qE,qQ,fs,fo,fm,fo0,fm0,Kn] = TB12(ps,Kp,Kf,Kd)

IMPLICIT NONE 

! Input variables 
INTEGER, INTENT(IN)                :: npt
REAL, DIMENSION(npt), INTENT(IN)   :: ps
REAL, INTENT(IN)                   :: Kp,Kf,Kd


! Out put variables 
REAL, INTENT(OUT)                  :: fo0,fm0
REAL, DIMENSION(npt),INTENT(OUT)   :: fo,eta,qE,qQ,fs,fm,Kn

! Local variables 
REAL, DIMENSION(npt)               :: x
REAL                               :: po0 


po0         = Kp/(Kf+Kd+Kp)  !     % dark photochemistry fraction (Genty et al., 1989)
x           = 1.-ps/po0  ! ;                % degree of light saturation
!%Kn          = (3.9867 * x - 1.0589).*x;  % empirical fit to Flexas, Daumard, Rascher, Berry data
!%p = [4.5531;8.5595;1.8510];
!%Kn1   = p(1)./(p(3)+exp(-p(2)*(x-.5)));
Kn          = (6.2473 * x - 0.5944)*x   !; % empirical fit to Flexas' data

fo0         = Kf/(Kf+Kp+Kd)    !   % dark adapted fluorescence yield Fo
fo          = Kf/(Kf+Kp+Kd+Kn) !;  % dark adapted fluorescence yield Fo
fm          = Kf/(Kf+Kd+Kn) !;  % light adapted fluorescence yield Fm
fm0         = Kf/(Kf+Kd)   !    % light adapted fluorescence yield Fm
fs          = fm*(1.-ps)
eta         = fs/fo0
qQ          = 1.-(fs-fo)/(fm-fo)   !;       % photochemical quenching
qE          = 1.-(fm-fo)/(fm0-fo0) !;      % non-photochemical quenching
return

END SUBROUTINE TB12



!function [A,F,rcw,f,Ci,zetan,Wc,Wj,Vcmax,Jmax] =
!biochemical(C,Q,T,eb,O,p,Vcmo,Jmo,lam,option,gcmethod)
! This is simply a translation of the matlab code of Ver dan Tol model et al.

   
SUBROUTINE biochemical_faq(npts,C1,Q,T1,eb,O1,p, Vcmo,Jmo,Ekc,Eko,Evcmax,Erd, Ekp,&
        &  akc,ako,option,gcmethod,An,Fl,rcw,W,Ag,pft,iflag)


! Note the following parameters are maybe intended to be optimized later ... 
!USE mo_carparams, ONLY : Ejmax, Vpmo, lam, Hkin, akp, Rd_f, &
USE mo_carparams, ONLY : Ejmax, lam, Hkin, akp, Rd_f, &
                     &   gbs, gcmin, fipo, al,  x_mas


!%
!% Date: 	10 Jun 2007
!%   Update: 24 Jul 2007: removed some external functions for unit 
!%                conversions. Only syntax changes, operation unchanged.
!%            3 Aug 2007: removes possible negative values of vapour
!%            pressure deficit, which cause imaginary Ci
!%				27 May 2008, included Leunings model as an alternative for Cowans model (CvdT)
!% Authors: 	Christiaan van der Tol, Wout Verhoef, Andries Rosema, Michiel van der Molen and Reinder Ronda
!% Sources: 	Farquhar et al. 1980, Cowan (1977), Arneth et al. (2002), 
!%           Rosema et al. (1999), and Van der Tol et al. (2009),
!%           Genty et al. (1989)
!   Adapted to FORTRAN code, E koffi, IPSL, Jussieu, Paris, 2012
!%
!% This function calculates:
!%    - stomatal resistance of a leaf or needle (s m-1)
!%    - photosynthesis of a leaf or needle (umol m-2 s-1) !%    - fluorescence of a leaf or needle (fraction of fluor. in the dark)
!%
!% Usage:
!% function [A,F,rcw] = biochemical(C,Q,T,eb,O,p,Vcmo,Jmo,lam,option,gcmethod)
!% the function was tested for Matlab versions 5.1 and 7.2.0.232 (R2006a)
!%
!% Calculates net assimilation rate A, fluorescence F using biochemical model
!%
!% Input (units are important): 
!% npts                              % the number of numerically selected
!calculations of fluo data depending on the type of the light on the leaves 
! we have 60 layers for shaded leaves and 13*36*60 for sunlit leaves with 
! 13 the number of leaf angles, 36 the number of leaf orientations and 60 the
! number of leaf layers 
!% C         % [umol m-3]            conc. of CO2 in leaf boundary layer
!% Q         % [umol photons m-2 s-1]net radiation, PAR
!% T         % [oC or K]             leaf temperature
!% eb        % [hPa]                 vapour pressure in leaf boundary layer
!% O         % [mmol m-3]            concentration of O2
!% p         % [Pa]                  air pressure
!% Vcmo      % [umol/m2/s]           maximum carboxylation capacity
!% Jmo       % [umol/m2/s]           maximum electron transport rate
!% lam		% []                    marginal cost of assimilation
!% option    % []                    1: C4, ~1: C3
!% gcmethod  % []                    default: Cowan (1977), 0: Leuning (1995)
!%
!% globals:
!% rhoa      % [kg m-3]              air mass density
!% Mair      % [g mol-1]             molecular mass of dry air
!% R         % [J mol-1K-1]          molar gas constant
!% gbs       % [umol/m2/s]           bundle sheath conductance to CO2(C4 only)
!% gcmin     % [m s-1]               minimum stomatal conductance (if stomata are closed). 
!% Rdopt     % [umol/m2/s]           dark respiration rate at 25 oC
!%
!% Note: always use the prescribed units. Temperature can be either oC or K
!% Note: input can be single numbers or vertical vectors
!%
!% Output:
!% An         % [umol/m2/s]           net assimilation rate
!% Fl         % []                    phi_f/phi_f0 ratio of fluorescence over fluorescence at low light
!% rcw       % [s m-1]               stomatal resistance for water
!% extra, optional output (in this order)
!% f         % []                    which is limiting: 0= enzyme, 1= light
!% Ci        % [umol m-3]            internal CO2 concentration
!%

! %% parameters (at optimum temperature)
!  fipo  = 0.82     !%  dark photochemistry fraction (Genty et al., 1989)
!    al  = 1E1      !; % smoothing parameter  (choose a value [3,100])
!     x  = .6       !;% partitioning factor of the electron transport rate
!     (Massad et al., 2007)
! x = x_mas in mo_carparams

! Ekc % [J/mol]  for CO2 (Farquhar, Q10=2.24)
! Eko % [J/mol]  for O2 (Farquhar, Q10=1.63)
! Ekp % [J/mol]  for PEP
! Erd % J/mol for dark respiration eg Q10=2 (Farquhar:66405@25C,Q10=2.46)
! Ejmax   % [J/mol] for electron transport
!Evcmax  % [J/mol]  for carboxylation (Farquhar, Q10=2.21)

!Kcopt   !  MICHAELIS-MENTEN CONSTANT FOR CO2
!Koopt  ! MICHAELIS-MENTEN CONSTANT FOR O2
!Kpopt   ! MICHAELIS-MENTEN CONSTANT FOR CO2 FOR C4 PLANT ???

!% 
!%  
!%
!% Parameters Vcmo and Jmo are corrected for temperature effects...
!%  with Boltzman functions
!% Parameters Kc, Ko and Rd are corrected for temperature effects ...
!%  with Arrhenius functions
!% check out the temperature definitions if you wish to change them.
!% CHANGE LINE 119 if TEMPERATURE CORRECTIONS OF PARAMETERS ARE REQUIRED


IMPLICIT NONE 

! Input variables
INTEGER, INTENT(IN)                  :: npts 

REAL, DIMENSION(npts), INTENT(IN)    :: Q,T1 

REAL, INTENT(IN)                     :: C1, O1,eb, p,Vcmo, Jmo
REAL, INTENT(IN)                     :: Ekc,Eko,Evcmax,Erd,Ekp
REAL, INTENT(IN)                     :: akc,ako

INTEGER, INTENT(IN)                  :: option,gcmethod
INTEGER, INTENT(IN)                  :: pft,iflag


! Output variables
REAL, DIMENSION(npts), INTENT(OUT)   :: An,Fl,rcw
REAL, DIMENSION(npts), INTENT(OUT)   :: W,Ag 

! Local variables
INTEGER                              :: i
INTEGER                              :: C4

REAL                                 :: C,O
REAL                                 :: Cs
REAL                                 :: D0
REAL                                 :: aleuning 

REAL, DIMENSION(npts)                :: f,Ci,Wc,Wj,Vcmax
REAL, DIMENSION(npts)                :: Jmax,zetan
REAL, DIMENSION(npts)                :: T 
REAL, DIMENSION(npts)                :: dT, RTT,fip
REAL, DIMENSION(npts)                :: J, Vomax,Vpmax 
REAL, DIMENSION(npts)                :: gam,a1,b1,a2,b2 
REAL, DIMENSION(npts)                :: Fe
!REAL, DIMENSION(npts)                :: Fe,W
REAL, DIMENSION(npts)                :: Ce, Cij,Cic 
REAL, DIMENSION(npts)                :: Kc, Ko, Kp,Rd 
REAL, DIMENSION(npts)                :: D 
REAL,DIMENSION(npts)                 :: es,s,val 

LOGICAL,DIMENSION(:),ALLOCATABLE     :: mask
INTEGER, ALLOCATABLE, DIMENSION(:)   :: indices 
REAL, ALLOCATABLE, DIMENSION(:)      :: wjn, wjd
LOGICAL, SAVE 	     		     :: testb=.TRUE.

REAL                                 ::  x       !;% partitioning factor of the electron transport rate (Massad et al., 2007)
REAL                                 ::  Rdopt   !   %[umol/m2/s]     dark respiration rate at 25 oC
REAL                                 ::  Kcopt   !%460;  % [ubar] kinetic coefficient  for CO2 (Von Caemmerer and Furbank, 1999)
REAL                                 ::  Koopt   !;%330;  % [mbar] kinetic coeeficient for  O2 (Von Caemmerer and Furbank, 1999)
REAL                                 ::  Kpopt   !; % [ubar]    kinetic coefficient for PEPcase for CO2

REAL                                 ::  Tref      ![K]absolute temperature at 25 oC
REAL                                 ::  Toptjmax  ![K] optimum temperature for max.electron transp. (25 oC)
REAL                                 ::  Toptvcmax ![K]  opt. T. for max. carbox.(25 oC)


! x  = .6       !;% partitioning factor of the electron transport rate (Massad et al., 2007)
x = x_mas 


!Rdopt   %[umol/m2/s]     dark respiration rate at 25 oC
Rdopt = Rd_f


! For C3 plant 
Tref = 25.+273.15 
Toptjmax = 298.
Toptvcmax = 298.   

Cs        = 2.*C1         ! %  CO2 concentration in bundle sheath !(only for C4)

!%% input
T = T1
WHERE (T1 <100.) 
T = T1 + 273.15  ! % convert temperatures to K if not already
ENDWHERE 

C4 = option           !% We consider a C3 by default, option = 0 and C4 option 1 

IF (C4) THEN  
Tref      = 35.+273.15
Toptjmax  = 308
Toptvcmax = 308

if (iflag==0)  then  
 open(10,file='verif_c4.dat',status='unknown')
 open(11,file='Qfield.dat',status='unknown')
 open(12,file='Tfield.dat',status='unknown')

write(10,*)'pft ', pft
write(10,*) ' C1 ',C1
write(10,*) ' O1 ',O1
write(10,*) ' eb ',eb
write(10,*) ' p  ',p
write(10,*) ' Vcmo ',Vcmo
write(10,*) ' Jmo ',Jmo
write(10,*) ' Ekc ', Ekc
write(10,*) ' Eko ', Eko
write(10,*) 'Evcmax ', Evcmax 
write(10,*) ' Erd ', Erd
write(10,*) ' Ekp ', Ekp
write(10,*) ' akc ', akc
write(10,*) ' ako ', ako
write(10,*) ' akp ', akp 

write(11,*) Q 
write(12,*) T1 

close(11)
close(12)

endif 
ENDIF    

!%% convert all to bar
C      = C1 * Mair/rhoa*1E-3 * 1e-6 * p *1E-5
Cs     = Cs* Mair/rhoa*1E-3 * 1e-6 * p *1E-5
O      = O1 * Mair/rhoa*1E-3 * 1e-3 * p *1E-5

Kcopt  = akc * 1e-6
Koopt  = ako * 1e-3
Kpopt  = akp * 1e-3   


!%% compute arrhenius functions for Kc, Ko, tau
dT = T - Tref
!dT = 0.
CALL arrhenius(npts,Kcopt,  Ekc,  dT, Tref, T, ar,Kc)

CALL arrhenius(npts,Koopt,  Eko,  dT, Tref, T, ar,Ko)

IF (C4)  CALL arrhenius(npts,Kpopt,  Ekp,  dT, Tref, T, ar,Kp)

CALL arrhenius(npts,Rdopt,  Erd,  dT, Tref, T, ar,Rd) 


!% Baldocchi says Rdopt = Rdz = Vcmo .* 0.004657 = 0.34 and ...
!% ... in light Rdopt = Rdz = 0.4 * Vcmo .* 0.004657

!%% Boltzmann temperature corrections to Jmax and Vcmax
CALL boltzmann(npts,Jmo,  Ejmax,  Toptjmax,T,ar,Hkin,Jmax)

CALL boltzmann(npts,Vcmo, Evcmax, Toptvcmax,T,ar,Hkin,Vcmax)

IF (C4)  CALL boltzmann(npts,Vcmo, Evcmax, Toptvcmax,T,ar,Hkin,Vpmax)

!%% calculation of potential electron transport rate
fip  = fipo* 2* Jmax/ (2*Jmax+fipo*Q)   !;%       photochemical quenching    

J    = 0.5*fip*Q                         ! %       light limitation

!%% calculation of the intersection of enzyme and light limited curves
!% this is the original Farquhar model
Vomax  = 0.21*Vcmax             !  % [umol/m2/s] max. rate oxygenation of Rubisco
gam    = 0.5*Vomax/Vcmax*Kc/Ko*O  !; %[bar]       compensation point [bar]
a1     = Vcmax
b1     = Kc* (1+O/Ko)
a2     = J/4
b2     = 2*gam

IF (C4) THEN 
    a1  = Vpmax
    b1  = Kp
    a2  = x*J/2
    b2  = 0
ENDIF 


!!%concentration at which enzyme and electron limited photosynthesis are equal
Ce = (a2*b1-a1*b2)/(a1-a2)

!%% calculation of internal CO2 concentration with Cowan's model for stomata
!% see Appendix in Arneth et al. (2002), Global Biochem. Cycles 16 (5/1  5/13), doi = "10.1029/2000GB001374"}
!% or Appendix in Van der Tol et al. (2007) Biogeosciences, 4(1), pp.137-154
CALL satvap(npts,(T-273.15),es,s)
val=(es-eb)*1e-3
D  = max(1E-5,val)

!% Cowan (1977)
IF (gcmethod) THEN 
   CALL Cowan(npts,b1,D,gam,lam,C,Cic)      !  %       Vmax case (RuBP lim)
   CALL Cowan(npts,b2,D,gam,lam,C,Cij)      !  %       Jmax case (electron lim)
    Ci  = Ce

WHERE( (Ce > Cic).AND.(Ce >Cij))
    Ci = Cic
ENDWHERE
WHERE( (Ce <= Cic).AND.(Ce <=Cij))
    Ci = Cij
END WHERE

!% f smooths the transition between enzyme and electron limited conditions
CALL smoothing(npts,Ce,Ci,C,al,f)

  Wc   = a1*Ci /(b1+Ci)
  Wj   = a2*Ci /(b2+Ci)
  W    = (1-f)*Wc + f*Wj          !%  gross photosynthesis rate
  Ag   = (1 - gam/Ci) *W          !%  GPP withouT Rd 
  An   = (1 - gam/Ci) *W - Rd     !%  net photosynthesis rate

ELSE
!%Leuning (1995) - Not used here since I do not know the value of aleuning EK
!July 2012 
  Ci=Cs

    DO i = 1,3
        CALL smoothing(npts,Ce,Ci,C,al,f)   !% smooths the transition between enzyme and electron limited conditions
        Wc  = a1*Ci/(b1+Ci)   
        Wj  = a2*Ci/(b2+Ci)
        W   = (1.-f)*Wc + f*Wj        !%  gross photosynthesis rate
        Ag  = (1 - gam/Ci) *W      !;% GPP without Rd
        An  = (1 - gam/Ci) *W - Rd   !;% net photosynthesis rate
     
        CALL leuning(npts,aleuning,An,D,D0,gcmin,C,Ci)
    END DO 
ENDIF 

!%% calculation of photosynthesis
!% this is the original Farquhar model

!% f smooths the transition between enzyme and electron limited conditions
CALL smoothing(npts,Ce,Ci,C,al,f)

Wc   = a1*Ci /(b1+Ci)
Wj   = a2*Ci /(b2+Ci)
W    = (1-f)*Wc + f*Wj          !%  gross photosynthesis rate
An   = (1 - gam/Ci) *W - Rd     !%  net photosynthesis rate

IF (C4) An  = An-gbs*(Cs-Ci)     !  %       leakage

!%% fluorescence
!% see Genty et al (1989), Rosema et al (1998) and Van der Tol et al.
!% (2009, AFM)
Fe       = (1-fip)/(1-fipo)
zetan    = 1.

IF (.NOT. ALLOCATED(mask)) ALLOCATE(mask(size(Wj)))
mask = (Wj>0)
IF (.NOT.ALLOCATED(indices)) ALLOCATE(indices(count(mask)))
IF (.NOT.ALLOCATED(wjn)) ALLOCATE(wjn(count(mask)))
IF (.NOT.ALLOCATED(wjd)) ALLOCATE(wjd(count(mask)))

CALL mask_ind(size(mask),mask,indices)

wjn=W(indices)
wjd=Wj(indices)

!zetan(j) = W(j)./Wj(j)
zetan(indices) = W(indices)/Wj(indices)
Fl       = Fe*zetan                 !% fluorescence (relative factor)

!%% convert back to umol m-3
Ci       = Ci*rhoa/Mair*1E3    * 1e6 / p * 1E5
C        = C *rhoa/Mair*1E3    * 1e6 / p * 1E5
!% Ce     = Ce*rhoa/Mair*1E3    * 1e6 ./ p .* 1E5;
!% Cic    = Cic*rhoa/Mair*1E3    * 1e6 ./ p .* 1E5;
!% Cij    = Cij*rhoa/Mair*1E3    * 1e6 ./ p .* 1E5;


!%% calculation of stomatal resistance
rcw   = 0.625*(C-Ci)/An
WHERE (An <=0) 
 rcw = 0.625/gcmin
ENDWHERE   

if (C4) then 
if (iflag==0) then 
write(10,*) 'rcw ', rcw 
write(10,*) ' W ', W
write(10,*)  ' An ', An
close(10)
! stop 
endif 
endif 

!IF (testb) THEN
!	print *, 'biochemical_faq was called'
!	testb = .FALSE.
!END IF
!print *, ' biochemical_faq was called '

!%%% end of function biochemical
END SUBROUTINE biochemical_faq 


SUBROUTINE  arrhenius(npts,Kopt, E, dT, Tref, T, R,K)  
!%% Arrhenius function
!function K = arrhenius(Kopt, E, dT, Tref, T, R)
IMPLICIT NONE 

! Inputs 
INTEGER, INTENT(IN)               :: npts 
REAL, DIMENSION(npts), INTENT(IN) :: dT,T 
REAL,    INTENT(IN)               :: Kopt,E,Tref 
REAL, INTENT(IN)                  :: R

!Outputs
REAL, INTENT(OUT),DIMENSION(npts) :: K

K = Kopt * exp(dT * E / (Tref* R* T) )


RETURN 
!%%% end of arrhenius function
END SUBROUTINE arrhenius


SUBROUTINE  boltzmann(npts,Vmaxopt, Eakin, Topt, T, R, Hkin,Vmax)  
!%% Boltzmann function
!function Vmax = boltzmann(Vmaxopt, Eakin, Topt, T, R, Hkin)
IMPLICIT NONE 

!Inputs
INTEGER, INTENT(IN)               :: npts 
REAL, DIMENSION(npts), INTENT(IN) :: T
REAL,  INTENT(IN)                 :: Vmaxopt,Eakin,Topt
REAL, INTENT(IN)                  :: R, Hkin

! Outputs 
REAL, INTENT(OUT),DIMENSION(npts) :: Vmax

! Local 
REAL, DIMENSION(npts)             :: dT,RTT,denom, numm

dT     = T - Topt
! We do not correct 
dT = 0.
RTT    = R * Topt * T
numm   = Vmaxopt * Hkin * exp(Eakin* dT/ RTT)
denom  = Hkin - Eakin * (1 - exp(Hkin * dT/ RTT))
Vmax   =  numm / denom

!print*, ' end  boltzmann Vmax  ', minval(Vmax), maxval(Vmax) 
RETURN
!%%% end of Boltzmann function
END SUBROUTINE  boltzmann


SUBROUTINE  smoothing(npts,Ce,Ci,C,a,z)    
!%% smoothing function
!function z = smoothing(Ce,Ci,C,a)
IMPLICIT NONE 

!Inputs 
INTEGER, INTENT(IN)                  :: npts
REAL, DIMENSION(npts), INTENT(IN)    :: Ci,Ce 
REAL, INTENT(IN)                     :: C,a

! Outputs 
REAL, INTENT(OUT),DIMENSION(npts)    :: z

z  = 1./(1+exp(-a*(Ci-Ce)/C))

RETURN 
!%%% end of smoothing function
END SUBROUTINE smoothing 

SUBROUTINE Cowan(npts,b,D,gam,lam,C,Ci)   
!%% Cowan's model, solution of Arneth et al. (2002), GBC 16(1) 5-1 - 5-13
!function Ci = Cowan(b,D,gam,lam,C)
IMPLICIT NONE 

! Inputs 
INTEGER, INTENT(IN)               :: npts
REAL, DIMENSION(npts), INTENT(IN) :: D,b,gam
REAL, INTENT(IN)                  :: C,lam

! Outputs 
REAL, INTENT(OUT),DIMENSION(npts) :: Ci

! Local 
REAL, DIMENSION(npts)             :: x,y,z,Ci0

x    = lam - 1.6*D/(b + gam)         
y    = 1.6*D- 2*C*lam  + 1.6*D*(gam-b)/(b+gam)
z    = (lam*C - 1.6*D)*C  + 1.6*D*gam*b/(b+gam)

CALL abc(npts,x,y,z,Ci0,Ci)  

RETURN 
!%%% end of Cowan's function
END SUBROUTINE  Cowan

SUBROUTINE  leuning(npts,aleuning,A,D,D0,gcmin,C,Ci)  
!%% Leunings model (Wang and Leuning, 1998)
!function Ci = leuning(aleuning,A,D,D0,gcmin,C)
IMPLICIT NONE 

! Inputs
INTEGER, INTENT(IN)                  :: npts
REAL, DIMENSION(npts), INTENT(IN) :: D,A
REAL, INTENT(IN)                  :: C,gcmin,aleuning,D0

! Outputs
REAL, INTENT(OUT),DIMENSION(npts) :: Ci

! Local
REAL, DIMENSION(npts)             :: gc 

gc = gcmin + aleuning*A/(1+D/D0)
Ci = C-A/gc

RETURN 
!%%% end of Leuning's function
END SUBROUTINE leuning 




END MODULE chemical 

