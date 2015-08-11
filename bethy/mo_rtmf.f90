MODULE mo_rtmf 


!function Ltot_ = RTMf(Esun_,    Emin_,  Eplu_   ,...
!                      F0,       F1              ,... 
!                      LAI                       ,... 
!                                Po,     Pso     ,... 
!                      tts,      tto,    psi     )
!
!% function 'RTMf' calculates the spectrum of outgoing fluorescence
!% in viewing direction.
!%
!% Authors: Wout Verhoef and Christiaan van der Tol (tol@itc.nl)
!% Date:     12 Dec 2007
!% Update:   26 Aug 2008 CvdT        small correction to matrices)
!%           07 Nov 2008 CvdT        changed layout
!% Update:   19 Mar 2009 CvdT        major corrections: lines 95-96,
!%                                   101-107, and 119-120.
!% Update:    7 Apr 2009 WV & CvdT   major correction: lines 89-90, azimuth
!%                                   dependence was not there in previous verions (implicit assumption of
!%                                   azimuth(solar-viewing) = 0). This has been corrected
!%
!% Fortran version : June 2012 Ernest Koffi 
!%
!% Table of contents of the function:
!%   0       preparations
!%       0.0     globals
!%       0.1     initialisations
!%       0.2     geometric quantities
!%       0.3     solar irradiance factor and ext. in obs dir for all leaf angle/azumith classes
!%   1.0     calculation of fluorescence flux in observation direction   
!%
!% Usage: function Ltot_ = ...
!%   RTMf2(Esun_,Emin_,Eplu_,F0,F1,LAI,Po,Pso,tts, tto)
!% 
!% Input:
!%   Symbol  Description                                 Unit            Dimension
!%   ------  -----------                                 ----            ---------
!%   Esun_   direct inc. PAR                             (W m-2 um-1)    [nwl]
!%   Emin_   downward fluxes, spectral                   (W m-2 um-1)    [nwl,nl+1] 
!%   Eplu_   upward fluxes, spectral                     (W m-2 um-1)    [nwl,nl+1]
!%   F0      fluorescence factor for shaded leaves                       [nl]
!%   F1      fluorescence factor for sunlit leaves                       [nli,nwlfo,nl]
!%   LAI     leaf area index                             (m2 m-2)        [1]
!%   Po      probability of viewing a leaf or soil                       [nl+1]
!%   Pso     probability of viewing a sunlit leaf/soil                   [nl+1]
!%   tts     solar zenith angle                          (degrees)       [1]
!%   tto     viewing angle                               (degrees)       [1]
!%   psi     azimuth angle                               (degrees)       [1]
!%
!% Output
!%   Symbol  Description                                 Unit            Dimension
!%   ------  -----------                                 ----            ---------
!%   Ltot_   total outgoing fluor in viewing dir         (W m-2 s1-1)    [1]   
!%
!% Notes:
!%   nl      number of layers
!%   nwfli_i   number of wavelengths of input (net PAR)
!%   nwfli_o   number of wavelengths of input (net PAR)
!% '_'means: a flux at different wavelengths (a vertically oriented vector)
!%


IMPLICIT NONE 

CONTAINS 


SUBROUTINE rtmf(Esun_, Emin_, Eplu_,etahi,etaui,LAI,Po,Ps,Pso,tts,tto,psi,LoF,Fhem,Fiprof,& 
                & MfI,MbI,MfII,MbII,rho,tau,rs)


USE fluo_func,  ONLY : Sint 
USE fluo_param, ONLY : nl,nwl,nli,nwlfi,nwlfo,iwlfi,iwlfo,nlazi,nwlF,wlF,&
               &  litab,lazitab,lidf


REAL, PARAMETER                                 :: pi      = 3.1415926535879
REAL, PARAMETER                                 :: deg2rad = pi/180.

!Input variables
REAL,  INTENT(IN)                               :: LAI,tts,tto,psi
REAL,  DIMENSION(nwl),INTENT(IN)                :: Esun_
REAL,  DIMENSION(nwl,nl+1), INTENT(IN)          :: Emin_, Eplu_
REAL,  DIMENSION(nl), INTENT(IN)                :: etahi
REAL,  DIMENSION(nli,nlazi,nl),INTENT(IN)       :: etaui
REAL,  DIMENSION(nl+1),INTENT(IN)               :: Po,Ps,Pso
REAL,  DIMENSION(:,:), INTENT(IN)               :: MfI,MbI,MfII,MbII
REAL,  DIMENSION(:), INTENT(IN)                 :: rho,tau,rs

! Output variables 
REAL,  DIMENSION(nwlfo), INTENT(OUT)            :: LoF,Fhem 
REAL,  DIMENSION(nl+1), INTENT(OUT)             :: Fiprof 

! Local variables 
REAL,  DIMENSION(nwl)                           :: Esun_b
REAL,  DIMENSION(nwl,nl+1)                      :: Emin_b, Eplu_b
REAL,  DIMENSION(nwlfo,nwlfi)                   :: Mb, Mf
REAL,  DIMENSION(nwlfo,nl+1)                    :: MpluEmin, MpluEplu
REAL,  DIMENSION(nwlfo,nl+1)                    :: MminEmin, MminEplu 
REAL,  DIMENSION(nwlfo,nl+1)                    :: Fmin,Fplu 
REAL,  DIMENSION(nl)                            :: piLs, piLu,piLd 
!REAL,  DIMENSION(nl+1)                          :: piLs, piLu,piLd 
REAL,  DIMENSION(nl+1)                          :: G1, G2
REAL,  DIMENSION(nwlfi)                         :: Esunf_
REAL,  DIMENSION(nwlfi,nl+1)                    :: Eminf_,Epluf_
REAL                                            :: Fdmin, Fsmin,Fdplu,Fsplu

REAL,  DIMENSION(nli,nlazi)                     :: cds,cdo,fs,fo,fsfo    
REAL,  DIMENSION(nli,nlazi)                     :: absfo,absfsfo,foctl 
REAL,  DIMENSION(nli,nlazi)                     :: absfs,fsctl 
REAL                                            :: cos_tto, sin_tto
REAL                                            :: cos_tts, sin_tts
REAL, DIMENSION(nli)                            :: cos_ttli, sin_ttli
REAL, DIMENSION(nli,nlazi)                      :: ctl2
REAL, DIMENSION(nlazi)                          :: cos_phils, cos_philo

REAL, DIMENSION(nwlfo)                          :: MpluEsun,MminEsun 
REAL, DIMENSION(nwlfo,nwlfi)                    :: Mplu,Mmin 

REAL, DIMENSION(nli,nlazi)                      :: wEsuni,vbEmini, vfEplui
REAL, DIMENSION(nli,nlazi)                      :: sigfEmini,sigfEplui
REAL, DIMENSION(nli,nlazi)                      :: sigbEmini,sigbEplui
REAL, DIMENSION(nli,nlazi)                      :: sfEsuni,sbEsuni 
REAL, DIMENSION(nli,nlazi)                      :: tpLu,  tpLs
REAL, DIMENSION(nlazi,nli)                      :: ttpLu,  ttpLs
REAL, DIMENSION(nlazi)                          :: cLu, cLs
REAL, DIMENSION(nwl)                            :: resolutionb

REAL, DIMENSION(nl)                             :: etah
REAL, DIMENSION(nli,nlazi,nl)                   :: etau 
REAL, DIMENSION(nwlfo,2)                        :: LoF_, Fhem_
REAL, DIMENSION(nl+1,2)                         :: Fiprofile 
REAL, DIMENSION(nwlfo)                          :: rho1, tau1 
!REAL, DIMENSION(nwl)                            :: rs 

REAL                                            :: xdd2 
REAL                                            :: ap 
REAL                                            :: piLtoti,piLtoti0,iLAI
REAL                                            :: t2,att,sig,m
REAL                                            :: fac, facs,rinf 
REAL                                            :: Gnew
INTEGER                                         :: i,j,k,ips
INTEGER                                         :: ilayer,ip
INTEGER                                         :: cont

INTEGER, PARAMETER                              :: fidf =3
REAL                                            :: Ltot_v


!%% 0.1 initialisations
 MpluEmin(:,:) = 0.
 MpluEplu(:,:) = 0.
 MminEmin(:,:) = 0.
 MminEplu(:,:) = 0.

![piLs,piLu]         = deal(zeros(nl+1,1));
      piLs(:)  = 0. 
      piLu(:)  = 0. 

!LoF_               = zeros(nwlfo,1);
        LoF(:) = 0.
     LoF_(:,:) = 0.
    Fhem_(:,:) = 0.
Fiprofile(:,:) = 0.

       rho1(:) = 0.
       tau1(:) = 0.

!%for speed-up the calculation only uses wavelengthi and wavelengtho part of the spectrum
!Esunf_              = Esun_(iwlfi);
!Eminf_              = Emin_(iwlfi,:);
!Epluf_              = Eplu_(iwlfi,:);
!iLAI                = LAI/nl;                   % LAI of a layer        [1]

Esunf_              = Esun_(iwlfi)
!print*, ' Esunf_ ', minval(Esunf_), maxval(Esunf_), sum(Esunf_)
Eminf_              = (Emin_(:,iwlfi))
!print*, ' Eminf_ ', minval(Eminf_), maxval(Eminf_), sum(Eminf_)
Epluf_              = (Eplu_(:,iwlfi))
!print*, ' Epluf_ ', minval(Epluf_), maxval(Epluf_), sum(Epluf_)

iLAI                = 1.*LAI/nl                     ! LAI of a layer        [1]



!fs    = cds/cos_tts                              !% [13,36]
!fo    = cdo/cos_tto                         !% [13,36]
       fsfo  =    fs*fo 

   absfs     = abs(fs)                                                   !% [13,36]
   absfo     = abs(fo)                                                   !% [13,36]
  absfsfo    = abs(fsfo)                                                 !% [13,36]


!print*, ' fsfo ', minval(fsfo), maxval(fsfo), sum(fsfo)
!print*, ' absfs ', minval(absfs), maxval(absfs), sum(absfs)
!print*, ' absfo ', minval(absfo), maxval(absfo), sum(absfo)

!%% 0.2 geometric quantities
!print*, ' size(rho(iwlfo),1) ', size(rho(iwlfo),1)
!print*, ' size(rho) ', size(rho,1)
!print*, ' size rho1 ', size(rho1,1)

rho1        = rho(iwlfo)  !; % [nwlfo]     leaf/needle reflection
tau1        = tau(iwlfo)  !; % [nwlfo]     leaf/needle transmission

!DO i =1, nwlfo 
!  rho1(i)        = rho(iwlfo(i))  !; % [nwlfo]     leaf/needle reflection
!  tau1(i)        = tau(iwlfo(i))  !; % [nwlfo]     leaf/needle transmission
!print*, i, iwlfo(i), rho(iwlfo(i))
!END DO 

!print*, 'iwlfi ', minval(iwlfi), maxval(iwlfi), sum(iwlfi)
!print*, 'iwlfo ', minval(iwlfo), maxval(iwlfo), sum(iwlfo)
!print*, 'IwlP ', minval(IwlP), maxval(IwlP), sum(IwlP)
!print*, 'LAI ', LAI
!print*, 'litab', minval(litab), maxval(litab), sum(litab)
!print*, 'lazitab', minval(lazitab), maxval(lazitab), sum(lazitab)
!print*, 'lidf', minval(lidf), maxval(lidf), sum(lidf)
!print*, 'Ps', minval(Ps), maxval(Ps), sum(Ps)
!print*, 'Po', minval(Po), maxval(Po), sum(Po)
!print*, 'Pso', minval(Pso), maxval(Pso), sum(Pso)
!print*, 'Esun_', minval(Esun_), maxval(Esun_), sum(Esun_)
!print*, 'Emin_', minval(Emin_), maxval(Emin_), sum(Emin_)
!print*, 'Eplu_', minval(Eplu_), maxval(Eplu_), sum(Eplu_)
!print*, 'rho', minval(rho), maxval(rho), sum(rho)
!print*, 'rho1', minval(rho1), maxval(rho1), sum(rho1)
!print*, 'tau', minval(tau), maxval(tau), sum(tau)
!print*, 'tau1', minval(tau1), maxval(tau1), sum(tau1)
!print*, 'rtmf: MbI ', sum(MbI) !minval(MbI), maxval(MbI), sum(MbI)
!print*, 'rtmf: MbII ', sum(MbII) !minval(MbII), maxval(MbII), sum(MbII)
!print*, 'rtmf: MfI ', sum(MfI) !minval(MfI), maxval(MfI), sum(MfI)
!print*, 'rtmf: MfII ', sum(MfII) !minval(MfII), maxval(MfII), sum(MfII)
!print*, 'rs ', minval(rs), maxval(rs), sum(rs)
!print*,  ' tto ', tto, ' tts ', tts, ' psi ', psi
!print*, 'etahi ', minval(etahi), maxval(etahi), sum(etahi)
!print*, 'etaui ', minval(etaui), maxval(etaui), sum(etaui)



cos_tto             = cos(tto*deg2rad)           ! % cos observation angle
sin_tto             = sin(tto*deg2rad)           ! % sin observation angle

cos_tts             = cos(tts*deg2rad)           ! % sinus solar angle
sin_tts             = sin(tts*deg2rad)           ! % sinus solar angle

cos_ttli            = cos(litab*deg2rad)         ! % cos leaf inclinaation angles
sin_ttli            = sin(litab*deg2rad)         ! % sin leaf inclinaation angles

cos_phils           = cos(lazitab*deg2rad)       ! % cos leaf orientation  angles
cos_philo           = cos((lazitab-psi)*deg2rad) !    % cos leaf orientation  angles

!%% 0.3 solar irradiance factor and ext. in obs dir for all leaf angle/azumith classes

DO j=1,nlazi  
 DO i=1,nli
  cds(i,j)   = cos_ttli(i)*cos_tts + sin_ttli(i)*sin_tts*cos_phils(j)
  cdo(i,j)   = cos_ttli(i)*cos_tto + sin_ttli(i)*sin_tto*cos_philo(j)
  fo(i,j)    = cdo(i,j)/cos_tto                         !% [13,36]
  fs(i,j)    = cds(i,j)/cos_tts                         !% [13,36]
  foctl(i,j) = fo(i,j)*(cos_ttli(i)) 
  fsctl(i,j) = fs(i,j)*(cos_ttli(i)) 
  ctl2(i,j)  = cos_ttli(i)**2.                         
 END DO 
END DO   

!print*, ' Esunf_ ', minval(Esunf_), maxval(Esunf_), sum(Esunf_)
!print*, ' Eminf_ ', minval(Eminf_), maxval(Eminf_), sum(Eminf_)
!print*, ' Epluf_ ', minval(Epluf_), maxval(Epluf_), sum(Epluf_)
!print*, ' cds ', minval(cds), maxval(cds), sum(cds)
!print*, ' cdo ', minval(cdo), maxval(cdo), sum(cdo)
!print*, ' fs ', minval(fs), maxval(fs), sum(fs)
!print*, ' foctl ', minval(foctl), maxval(foctl), sum(foctl)
!print*, ' fsctl ', minval(fsctl), maxval(fsctl), sum(fsctl)
!print*, ' ctl2 ', minval(ctl2), maxval(ctl2), sum(ctl2)


!fs    = cds/cos_tts                              !% [13,36]
!fo    = cdo/cos_tto                         !% [13,36]
       fsfo  =    fs*fo 

   absfs     = abs(fs)                                                   !% [13,36]
   absfo     = abs(fo)                                                   !% [13,36]
  absfsfo    = abs(fsfo)                                                 !% [13,36]


!print*, ' fsfo ', minval(fsfo), maxval(fsfo), sum(fsfo)
!print*, ' absfs ', minval(absfs), maxval(absfs), sum(absfs)
!print*, ' absfo ', minval(absfo), maxval(absfo), sum(absfo)
!print*, ' absfsfo ', minval(absfsfo), maxval(absfsfo), sum(absfsfo)

!2. calculation of fluorescence flux in observation direction
!% fluorescence matrices

DO ips = 2,1,-1 
  IF (ips==1) THEN 
      Mb = MbI
      Mf = MfI 
    etah = 1.
    etau = 1. 
 ENDIF 

 IF (ips==2) THEN
     Mb  = MbII
     Mf  = MfII
    etah = etahi 
    etau = etaui
 ENDIF

Mplu                = 0.5*(Mb+Mf)              !% [nwlfo,nwl]
Mmin                = 0.5*(Mb-Mf)              !% [nwlfo,nwl]

!print*, ' Mb ', minval(Mb), maxval(Mb), sum(Mb)
!print*, ' Mf ', minval(Mf), maxval(Mf), sum(Mf)
!print*, ' Mplu ', minval(Mplu), maxval(Mplu), sum(Mplu)
!print*, ' Mmin ', minval(Mmin), maxval(Mmin), sum(Mmin)

DO  j = 1,nl+1
!     MpluEmin(:,j)   = matmul(Mplu,(Eminf_(:,j)*resolution(iwlfi)))   ! %[nwlfo,nl+1]
!     MpluEplu(:,j)   = matmul(Mplu,(Epluf_(:,j)*resolution(iwlfi)))   !  % [nwlfo,nl+1]
!     MminEmin(:,j)   = matmul(Mmin,(Eminf_(:,j)*resolution(iwlfi)))   ! %[nwlfo,nl+1]
!     MminEplu(:,j)   = matmul(Mmin,(Epluf_(:,j)*resolution(iwlfi)))   !  % [nwlfo,nl+1]


     MpluEmin(:,j)   = matmul(Mplu,Eminf_(:,j))   ! %[nwlfo,nl+1]
     MpluEplu(:,j)   = matmul(Mplu,Epluf_(:,j))   !  %[nwlfo,nl+1]
     MminEmin(:,j)   = matmul(Mmin,Eminf_(:,j))   ! %[nwlfo,nl+1]
     MminEplu(:,j)   = matmul(Mmin,Epluf_(:,j))   !  % [nwlfo,nl+1]

!print*, '  MpluEmin(:,j) ', j,minval( MpluEmin(:,j)), maxval( MpluEmin(:,j)), sum( MpluEmin(:,j))
!print*, '  MpluEplu(:,j) ', j,minval( MpluEplu(:,j)), maxval( MpluEplu(:,j)), sum( MpluEplu(:,j))
!print*, '  MminEmin(:,j) ', j,minval( MminEmin(:,j)), maxval( MminEmin(:,j)), sum( MminEmin(:,j))
!print*, '  MminEplu(:,j) ', j,minval( MminEplu(:,j)), maxval( MminEplu(:,j)), sum( MminEplu(:,j))

END DO

!    MpluEsun   = matmul(Mplu,(Esunf_*resolution(iwlfi))) !      % integration by inproduct
!    MminEsun   = matmul(Mmin,(Esunf_*resolution(iwlfi)))   !    % integration by inproduct


    MpluEsun   = matmul(Mplu,Esunf_) !      % integration by inproduct
    MminEsun   = matmul(Mmin,Esunf_)   !    % integration by inproduct

!print*, '  MpluEsun ', minval( MpluEsun), maxval( MpluEsun), sum(MpluEsun)
!print*, '  MminEsun ', minval( MminEsun), maxval( MminEsun), sum(MminEsun)

         xdd2  =  sum(matmul(transpose(ctl2),lidf))/nlazi 

!print*, ' xdd2 ', xdd2



! we calculate the spectrum for all individual leaves, sunlit and shaded. 
DO i = 1,nwlfo

       wEsuni       = absfsfo * MpluEsun(i) + fsfo  * MminEsun(i) !  % [nli,nlazi] 
       sfEsuni      = absfs   * MpluEsun(i) - fsctl * MminEsun(i) !;  % [nli,nlazi]
       sbEsuni      = absfs   * MpluEsun(i) + fsctl * MminEsun(i) !; %  [nli,nlazi]

!    print*,' wEsuni ' , minval(  wEsuni), maxval( wEsuni), sum( wEsuni) 
!    print*,' sfEsuni ' , minval(  sfEsuni), maxval( sfEsuni), sum( sfEsuni) 
!    print*,' sbEsuni ' , minval(  sbEsuni), maxval( sbEsuni), sum( sbEsuni) 

    DO  j =1,nl

        sigfEmini   = MpluEmin(i,j)   - ctl2 * MminEmin(i,j)  ! [nli,nlazi]
        sigfEplui   = MpluEplu(i,j+1) - ctl2 * MminEplu(i,j+1)  !  % [nli,nlazi]
        sigbEmini   = MpluEmin(i,j)   + ctl2 * MminEmin(i,j)  ! % [nli,nlazi]
        sigbEplui   = MpluEplu(i,j+1) + ctl2 * MminEplu(i,j+1) !  %  [nli,nlazi]

        vbEmini     = absfo     *MpluEmin(i,j)  + foctl *MminEmin(i,j)         ! % [nli,nlazi]
        vfEplui     = absfo     *MpluEplu(i,j+1)  - foctl *MminEplu(i,j+1)        ! % [nli,nlazi]
         
         piLs(j)    = sum(matmul(transpose(etau(:,:,j)*(wEsuni + vbEmini + vfEplui)),lidf))/nlazi
         piLd(j)    = etah(j) * (sum(matmul(transpose(vbEmini + vfEplui),lidf))/nlazi)   ! 
            
         Fsmin      =  sum(matmul(transpose(etau(:,:,j)*(sfEsuni + sigfEmini + sigbEplui)),lidf))/nlazi
         Fsplu      =  sum(matmul(transpose(etau(:,:,j)*(sbEsuni + sigbEmini + sigfEplui)),lidf))/nlazi

         Fdmin      =  etah(j) *(sum(matmul(transpose(sigfEmini + sigbEplui),lidf))/nlazi)
         Fdplu      =  etah(j) * (sum(matmul(transpose(sigbEmini + sigfEplui),lidf))/nlazi)

        Fmin(i,j+1) = Ps(j+1) * Fsmin + (1.-Ps(j+1)) * Fdmin
         Fplu(i,j)  = Ps(j)   * Fsplu + (1.-Ps(j))   * Fdplu


 !   print*,'  sigfEmini ' , minval( sigfEmini), maxval( sigfEmini), sum(sigfEmini) 
 !   print*,'  sigfEplui ' , minval( sigfEplui), maxval( sigfEplui), sum(sigfEplui)
 !   print*,'  sigbEmini ' , minval( sigbEmini), maxval( sigbEmini), sum(sigbEmini)
 !   print*,'  sigbEplui ' , minval( sigbEplui), maxval( sigbEplui), sum(sigbEplui)

 !   print*,' vbEmini ' , minval( vbEmini), maxval( vbEmini), sum(vbEmini)
 !   print*,' vfEplui ' , minval( vfEplui), maxval( vfEplui), sum(vfEplui)

  !  print*,' PiLs ' ,PiLs(j) ,' PiLd ', PiLd(j), ' Fsmin ', Fsmin, ' Fsplu ', Fsplu 
 !   print*,' PiLu ' ,PiLu(j) 

  !   print*,' Fdmin ' , Fdmin
  !   print*,' Fdplu ' , Fdplu
  !   print*,' Fmin ' , Fmin(i,j+1)
  !   print*,' Fplu ' , Fmin(i,j)


    END DO 


    !print*, ' PiLd ', PiLd 

   ! print*,'  PiLs ' , i,minval( PiLs), maxval( PiLs), sum(PiLs)
   ! print*,'  PiLd ' , i,minval( PiLd), maxval( PiLd), sum(PiLd)
   ! print*,'  sigfEmini ' , i,minval( sigfEmini), maxval( sigfEmini), sum(sigfEmini)
   ! print*,'  sigfEmini ' , i,minval( sigfEmini), maxval( sigfEmini), sum(sigfEmini)
   ! print*,'  sigfEplui ' , i,minval( sigfEplui), maxval( sigfEplui), sum(sigfEplui)
   ! print*,'  sigbEmini ' , i,minval( sigbEmini), maxval( sigbEmini), sum(sigbEmini)
   ! print*,'  sigbEplui ' , i,minval( sigbEplui), maxval( sigbEplui), sum(sigbEplui)

   ! print*,' vbEmini ' , i,minval( vbEmini), maxval( vbEmini), sum(vbEmini)
   ! print*,' vfEplui ' , i,minval( vfEplui), maxval( vfEplui), sum(vfEplui)


        Fmin(i,1) = 0.     ! No fluorescence from sky and from the soil 
    Fplu(i,nl+1)  = 0. 

    piLtoti         = iLAI*sum(Pso(1:nl)*piLs + (Po(1:nl)-Pso(1:nl))*piLd)                     !  % total canopy pi*F in obs dir for this ps 

    IF (piLtoti.gt.500) THEN
        print*,'**** piLtoti **** ', piLtoti
        print*,' rtmf:     MbI      ', sum(MbI) !minval(MbI), maxval(MbI), sum(MbI)
        print*,' rtmf:     MbII      ', sum(MbII) !minval(MbII), maxval(MbII), sum(MbII)
        print*,' rtmf:     MfI      ', sum(MfI) !minval(MfI), maxval(MfI), sum(MfI)
        print*,' rtmf:     MfII      ', sum(MfII) !minval(MfII), maxval(MfII), sum(MfII)
!        print*,'    iLAI ', iLAI
!        print*,'    piLs ', minval(piLs), maxval(piLs)
!        print*,'    piLd ', minval(piLd), maxval(piLd)
    ENDIF

    LoF_(i,ips)      = piLtoti/pi



 !print*,' ips ', ips, i,'  piLtoti ' ,  piLtoti
 !print*,' ips ', ips, i,'  Pso*PiLs ' ,  Pso*PiLs
 !print*,' ips ', ips, i,' (Po-Pso)*PiLd ' ,  (Po-Pso)*PiLd
 !print*,' ips ', ips, i,' sum  ' ,  sum(Pso*PiLs +(Po-Pso)*PiLd)
 !print*,' ips ', ips, i,' ips  LoF_ ' , ips, LoF_(i,ips)


              t2    = xdd2*(rho1(i)-tau1(i))/2.
             att    = 1.-(rho1(i)+tau1(i))/2.+t2
             sig    = (rho1(i) + tau1(i))/2.+t2 
               m    = sqrt(att**2.-sig**2.)
            rinf    = (att - m)/sig 
             fac    = 1. -m *iLAI 
            facs    = (rs(i)-rinf)/(1.-rs(i)*rinf)


              G1(1) = 2.
               Gnew = 0.

       ! print*, i,' t2 ', t2 , ' att ', att , ' sig  ', sig 
       ! print*, i,' m ', m , ' rinf  ', rinf , ' fac ', fac 
       ! print*, i,' facs  ', facs 

    DO WHILE ( abs(Gnew-G1(1)) > 1.e-3)

           G1(1) = Gnew 
    
           DO ip = 2,nl+1 
            G1(ip) = fac * G1(ip-1) + (Fmin(i,ip)+rinf*Fplu(i,ip-1)) * iLAI;
           END DO 
     
          G2(nl+1)= G1(nl+1)*facs
    
          DO ip = nl,1,-1
             G2(ip) = fac*G2(ip+1) + (rinf*Fmin(i,ip+1) + Fplu(i,ip))*iLAI
          END DO
  
       Gnew  = -rinf * G2(1)
        
   END DO    ! END WHILE   

   Fhem_(i,ips) = (rinf*G1(1)+G2(1))/(1.-rinf**2.)

 !  print*,i, ' Fhem_ ', Fhem_(i,ips)
   

 END DO 

  DO ilayer =1, nl+1 
  CALL Sint(nwlF,Fplu(:,ilayer),wlF,ap)
  Fiprofile(ilayer,ips) = ap*1.e-3
!  print*, ' Fplu ', minval(Fplu(:,ilayer)), maxval(Fplu(:,ilayer)),sum(Fplu(:,ilayer))
!  print*, ' Fiprofile ', ap*1.e-3

  END DO 

END DO    ! Loop on  ips 

!print*, ' nwlF ',  nwlF
!print*, ' Fplu ', minval(Fplu), maxval(Fplu), sum(Fplu)
!print*, ' wlF ', minval(wlF), maxval(wlF), sum(wlF)


LoF    =  LoF_(:,1) + LoF_(:,2) 
Fhem   =  Fhem_(:,1) + Fhem_(:,2)
Fiprof =  Fiprofile(:,1) + Fiprofile(:,2)



!DO i=1,nwlfo 
!print*,  ' LoFL ', i,LoF (i)
!END DO  

!DO i=1,nwlfo 
!print*,  ' FHEM ', i,Fhem(i)
!END DO  

!DO i=1,nl+1
!print*, ' PROF ',  i, Fiprof(i)
!END DO
!print*, 'IN RMTF  min max sum  ', minval(Fiprof), maxval(Fiprof), sum(Fiprof) 


END SUBROUTINE rtmf 

!SUBROUTINE Sint(nx,y,x,yi)

!    % Simpson integration
!    % x and y must be any vectors (rows, columns), but of the same length
!    % x must be a monotonically increasing series
    
!    % WV Jan. 2013, for SCOPE 1.40
! Translate in Fortran 2013, April, E Koffi 

!IMPLICIT NONE 
 
! Input variables 
!INTEGER, INTENT(IN)                :: nx 
!REAL, DIMENSION(nx), INTENT(IN)    :: x,y

! Ouput variables 
!REAL, INTENT(OUT)                  :: yi 

! Local variables 
!REAL,DIMENSION(nx-1)               :: step, mean 

!    nx   = length(x);
!    if (size(x,1) == 1) THEN 
!        x = transpose(x)
!    endif 

!    if (size(y,1).ne.1) then 
!        y = transpose(y)
!    endif 

!    step = x(2:nx) - x(1:nx-1)
!    mean = .5 * (y(1:nx-1) + y(2:nx))
!    yi  = dot_product(mean, step)

!END SUBROUTINE Sint 

END MODULE mo_rtmf
