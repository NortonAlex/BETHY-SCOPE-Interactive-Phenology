!*********************************************************
!*  SUBROUTINE diagnostics_allocate
!*********************************************************
SUBROUTINE veg_deallocate (vp)
  USE mo_constants
  USE mo_vegetation
  USE mo_surface
  USE mo_beta
  USE mo_carvar
  USE mo_carparams, ONLY: nl
  use mo_pheno
  USE mo_hydro
  !USE mo_grid, ONLY: ng, vp

  IMPLICIT NONE
  INTEGER, INTENT(in) :: vp
  
!WOK-ADD-070612 ztrans
!WOK-ADD-070626 zpcevp, zpsevp
!WOK-CHG-070626 removed zctrans
!WOK-CHG-070711 moved 'lai' here from mo_pheno
!$taf next required = vp
  DEALLOCATE( rbe )
!$taf next required = vp
  DEALLOCATE( vm, jmf, fautleaf, ccost )
!$taf next required = vp
  DEALLOCATE( aw, q10f, q10s, fracs, tauf )
!$taf next required = vp
  DEALLOCATE( EC, EO, EV, ER, EK, tgam )
!$taf next required = vp
  DEALLOCATE( alpha, alc4, kc0, ko0 )
!$taf next required = vp
  DEALLOCATE( zfpar, zgc, ptv, zfc )
!$taf next required = vp
  DEALLOCATE( zassc, zraut, zrphc,zgrowth,zmaint )
!$taf next required = vp
  DEALLOCATE( ztrans, zptrans, zpcevp, zpsevp )
!$taf next required = vp, nl
  DEALLOCATE( pgs )
!$taf next required = vp
  DEALLOCATE( beta, cs, fe )
!$taf next required = vp
  DEALLOCATE( pci, zci )
!$taf next required = vp
  DEALLOCATE( psoilst )
!$taf next required = vp
  DEALLOCATE( zlai, lai )
!$taf next required = vp
  DEALLOCATE( Chl )
!$taf next required = vp
  DEALLOCATE( Cdm_arr )
!$taf next required = vp
  DEALLOCATE( Csm_arr )
!$taf next required = vp
  DEALLOCATE( hc_arr )
!$taf next required = vp
  DEALLOCATE( LIDFa_arr )
!$taf next required = vp
  DEALLOCATE( LIDFb_arr )
!$taf next required = vp
  DEALLOCATE( leafwidth_arr )
!$taf next required = vp
  DEALLOCATE( dapar )

END SUBROUTINE veg_deallocate
