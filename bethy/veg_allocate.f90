!*********************************************************
!*  SUBROUTINE veg_allocate
!*********************************************************
SUBROUTINE veg_allocate (vp)
  USE mo_constants
  USE mo_carparams, ONLY: nl
  USE mo_vegetation
  USE mo_surface
  USE mo_beta
  USE mo_carvar
  USE mo_pheno
  USE mo_hydro
  !USE mo_grid, ONLY: ng, vp
 
  IMPLICIT NONE
  INTEGER, INTENT(in) :: vp
  
!WOK-ADD-070612 ztrans
!WOK-ADD-070626 zpcevp, zpsevp
!WOK-CHG-070626 removed zctrans
!WOK-CHG-070711 moved 'lai' here from mo_pheno
  ALLOCATE( rbe(vp) )
  ALLOCATE( vm(vp), jmf(vp), fautleaf(vp), ccost(vp) )
  ALLOCATE( aw(vp), q10f(vp), q10s(vp), fracs(vp), tauf(vp) )
  ALLOCATE( EC(vp), EO(vp), EV(vp), ER(vp), EK(vp), tgam(vp) )
  ALLOCATE( alpha(vp), alc4(vp), kc0(vp), ko0(vp) )
  ALLOCATE( zfpar(vp), zgc(vp), ptv(vp), zfc(vp) )
  ALLOCATE( ztrans(vp), zptrans(vp), zpcevp(vp), zpsevp(vp) )
  ALLOCATE( zassc(vp), zraut(vp), zrphc(vp),zgrowth(vp),zmaint(vp) )
  ALLOCATE( pgs(vp,nl) )
  ALLOCATE( beta(vp), cs(vp), fe (vp) )
  ALLOCATE( pci(vp,nl), zci(vp) )
  ALLOCATE (psoilst(vp))
  ALLOCATE (zlai(vp), lai(vp))
  ALLOCATE( Chl(vp), Cdm_arr(vp), Csm_arr(vp), hc_arr(vp), LIDFa_arr(vp), LIDFb_arr(vp) )
  ALLOCATE( leafwidth_arr(vp) )
  ALLOCATE( vms(vp), kc0s(vp), ko0s(vp), vomf(vp), rdf(vp) )
  zlai = 0.
  zfpar = 0.
  q10s = 0.
  fracs = 0.
  beta = 0.
  aw = 0.
  zgc = 0.
END SUBROUTINE veg_allocate
