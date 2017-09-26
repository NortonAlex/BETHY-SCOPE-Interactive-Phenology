MODULE mo_mapping 
!------------------------------------------------------------------
! map vector of control parameters to physical model variables
!------------------------------------------------------------------
    
  IMPLICIT NONE
  
  INTEGER, DIMENSION( :, :), ALLOCATABLE :: paramap_global, paramap_site, paramap
  REAL, PUBLIC                           :: glob_offset, fbiasa, fbiasb
  real, parameter                        :: eps = 1.e-4     ! a small number


CONTAINS

  SUBROUTINE init_global_mapping (mapping_file_global)
    USE mo_ctrl
    USE mo_constants
    USE mo_grid, ONLY: vp

    CHARACTER(len=*) :: mapping_file_global
    ! local variables
    INTEGER i,j


!!    WRITE(6,*) '  - init parameter mapping according to PFT grid'
    ALLOCATE( paramap_global (npar, vp))
    OPEN(1,file=mapping_file_global, status='old')
    REWIND 1
    READ(1,*) i,j
    IF((i /= npar) .OR. (j /= vp)) THEN       
       WRITE(6,*) 'init_mapping: inconsistent dimensions in file ',mapping_file_global
       WRITE(6,*) i,j,npar,vp
       STOP
    ENDIF
    DO i=1,vp
       READ(1,'(9999i6)') paramap_global(:, i)
    END DO
    RETURN
  END SUBROUTINE init_global_mapping


  SUBROUTINE init_site_mapping (mapping_file_site)
    USE mo_ctrl
    USE mo_constants
    USE mo_grid, ONLY: sp

    CHARACTER(len=*) :: mapping_file_site 
    ! local variables
    INTEGER i,j

!!    WRITE(6,*) '  - init parameter mapping according PFT grid'
    ALLOCATE( paramap_site (npar, sp))
    OPEN(1,file=mapping_file_site, status='old')
    REWIND 1
    READ(1,*) i,j
    IF((i /= npar) .OR. (j /= sp)) THEN       
       WRITE(6,*) 'init_mapping: inconsistent dimensions in file ',mapping_file_site
       WRITE(6,*) i,j,npar,sp
!       STOP
    ENDIF
    DO i=1,sp
       READ(1,'(9999i6)') paramap_site(:, i)
    END DO
    RETURN
  END SUBROUTINE init_site_mapping
  
!------------------------------------------------------------------
! ********  mapping maps natural params x onto physical params p = (vm, ...) 
!------------------------------------------------------------------
! HEW-ADD-03/09/08 new arguments (phys. params, bounds) for parmapping
SUBROUTINE mapping(par,nvar,nxp, scale)

  USE mo_constants
  USE mo_vegetation
  USE mo_pheno
  USE mo_beta
  USE mo_ctrl
  USE mo_grid
  USE mo_hydro, ONLY: pasmmax, pasmmax0
  USE mo_namelist, ONLY: optpftg, optpftl, optbsmg, optbsml

  IMPLICIT NONE

  INTEGER, INTENT(in) :: nvar, nxp
  REAL, DIMENSION(nvar), INTENT(in) :: par
  INTEGER, INTENT(in) :: scale

  INTEGER i, j, ipar, k, npoints, n
  REAL pasmmax_tree, pasmmax_herb, f_herb, frac_tree, frac_herb, mult_tree, mult_herb

  IF (scale==1) THEN
     ALLOCATE( paramap(npar, vp))
     paramap = paramap_global
     npoints = vp        ! loop over all vege points, map according PFT (vp) dependance
  ELSE ! (scale==2)
     ALLOCATE( paramap(npar, sp))
     paramap = paramap_site
     npoints = sp        ! loop over all vege points, map according PFT (vp) dependance
  END IF

  DO i = 1, npoints
     ! vm
     ipar = paramap(ivm,i)
     vm (i) = par(ipar)

     ! jmf
     ipar = paramap(ijmf,i)
     jmf (i) = par(ipar)

     ! fautleaf
     ipar = paramap(ifautleaf,i)
     fautleaf (i) = par(ipar)

     !  ccost
     ipar = paramap(iccost,i)
     ccost  (i) = par(ipar)

     ! q10s
     ipar = paramap(iq10s,i)
     q10s  (i) = par(ipar)

     ! q10f
     ipar = paramap(iq10f,i)
     q10f  (i) = par(ipar)

     ! tauf
     ipar = paramap(itauf,i)
! MAS20100520:changed parameter for tauf to 1./tauf
     tauf  (i) = par(ipar)

     ! aw
     ipar = paramap(iaw,i)
     aw  (i) = par(ipar)

     ! fracs
     ipar = paramap(ifracs,i)
     fracs  (i) = par(ipar)

     ! er
     ipar = paramap(ier,i)
     er  (i) = par(ipar)

     ! ev
     ipar = paramap(iev,i)
     ev  (i) = par(ipar)

     ! eo
     ipar = paramap(ieo,i)
     eo  (i) = par(ipar)

     ! ec
     ipar = paramap(iec,i)
     ec  (i) = par(ipar)

     ! ek
     ipar = paramap(iek,i)
     ek  (i) = par(ipar)

     ! alpha
     ipar = paramap(ialpha,i)
     alpha  (i) = par(ipar)

     ! alc4
     ipar = paramap(ialc4,i)
     alc4  (i) = par(ipar)

     ! kc0
     ipar = paramap(ikc0,i)
     kc0  (i) = par(ipar)

     ! ko0
     ipar = paramap(iko0,i)
     ko0  (i) = par(ipar)

     ! tgam
     ipar = paramap(itgam,i)
     tgam  (i) = par(ipar)

     ! beta
     ipar = paramap(ibeta,i)
     beta(i) = par(ipar)

     ! plaimax
     ipar = paramap(iplaimax,i)
     plaimax(i) =  par(ipar)

     ! ptphen
     ipar = paramap(iptphen,i)
     ptphen(i) = par(ipar)

     ! ptphenr
     ipar = paramap(iptphenr,i)
     ptphenr(i) = par(ipar)

     ! ptphen
     ipar = paramap(ipdphen,i)
     pdphen(i) = par(ipar)

     ! pdphenr
     ipar = paramap(ipdphenr,i)
     pdphenr(i) = par(ipar)

     ! ptshd
 !    ipar = paramap(iptshd,i)
 !    ptshd(i) = par(ipar)

     ! ptshds
!     ipar = paramap(iptshds,i)
!     ptshds(i) = par(ipar)

     ! plgr
     ipar = paramap(iplgr,i)
     plgr(i) = par(ipar)

     ! pkl
     ipar = paramap(ipkl,i)
     pkl(i) = par(ipar)

     ! ptauw
     ipar = paramap(iptauw,i)
     ptauw(i) = par(ipar)

     ! Chl
     ipar = paramap(iChl,i)
     Chl (i) = par(ipar)

     ! Cdm_arr (vp array)
     ipar = paramap(iCdm,i)
     Cdm_arr (i) = par(ipar)

     ! Csm_arr  (vp array)
     ipar = paramap(iCsm,i)
     Csm_arr (i) = par(ipar)

     ! hc_arr  (vp array)
     ipar = paramap(ihc,i)
     hc_arr (i) = par(ipar)

     ! LIDFa_arr  (vp array)
     ipar = paramap(iLIDFa,i)
     LIDFa_arr (i) = par(ipar)

     ! LIDFb_arr  (vp array)
     ipar = paramap(iLIDFb,i)
     LIDFb_arr (i) = par(ipar)

     ! leafwidth_arr  (vp array)
     ipar = paramap(ileafwidth,i)
     leafwidth_arr (i) = par(ipar)

     ! vms (Vcmax for SCOPE only)
     ipar = paramap(ivms,i)
     vms (i) = par(ipar)

     ! kc0s (kc0 for SCOPE only)
     ipar = paramap(ikc0s,i)
     kc0s (i) = par(ipar)

     ! ko0s (ko0 for SCOPE only)
     ipar = paramap(iko0s,i)
     ko0s (i) = par(ipar)

     ! vomf (Vomax as fraction of Vcmax, for SCOPE only)
     ipar = paramap(ivomf,i)
     vomf (i) = par(ipar)

     ! rdf (Rd as fraction of Vcmax, for SCOPE only)
     ipar = paramap(irdf,i)
     rdf (i) = par(ipar)

     ! pks
!     ipar = paramap(ipks,i)
!     pks(i) = par(ipar)
!     pks(i) = 1./30.

     ! pkm
!     ipar = paramap(ipkm,i)
!     pkm(i) = par(ipar)
!     pkm(i) = 1./30.

     ! pkml
!     ipar = paramap(ipkml,i)
!     pkml(i) = par(ipar)
!     pkml(i) = 1./180.

     n = gridp(i)
!      IF (scale==1) THEN
!         IF (optpftg) THEN
!            frac(i) = frac_p(i) + par(nxp+i)
!            IF (optbsmg) THEN
!               pasmmax(i)= pasmmax0(n) + par(nxp+vp+i)
!            ELSE
!               pasmmax(i)= pasmmax0(n)
!            ENDIF
!         ELSE
!            frac(i) = frac_p(i)
!            IF (optbsmg) THEN
!               pasmmax(i)= pasmmax0(n) + par(nxp+i)
!            ELSE
!               pasmmax(i)= pasmmax0(n)
!            ENDIF
!         ENDIF
!      ELSEIF (scale==2) THEN 
!         IF (optpftl) THEN
!            IF (optpftg .AND. optbsmg) THEN
!               frac(i) = frac_p(i) + par(nxp+vp+vp+i)
!               IF (optbsml) THEN
!                  pasmmax(i) = pasmmax0(n) + par(nxp+vp+vp+sp+i)
!               ELSE
!                  pasmmax(i) = pasmmax0(n)
!               ENDIF
!            ELSEIF ((optpftg .AND. .NOT. optbsmg) .OR. (.NOT. optpftg .AND. optbsmg)) THEN
!               frac(i) = frac_p(i) + par(nxp+vp+i)
!               IF (optbsml) THEN
!                  pasmmax(i) = pasmmax0(n) + par(nxp+vp+sp+i)
!               ELSE
!                  pasmmax(i) = pasmmax0(n)
!               ENDIF
!            ELSE
!               frac(i) = frac_p(i) + par(nxp+i)
!               IF (optbsml) THEN
!                  pasmmax(i) = pasmmax0(n) + par(nxp+sp+i)
!               ELSE
!                  pasmmax(i) = pasmmax0(n)
!               ENDIF
!            ENDIF
!         ELSE
            frac(i) = frac_p(i)
!            IF (optbsml) THEN
!               IF (optpftg .AND. optbsmg) THEN                      
!                  pasmmax(i) = pasmmax0(n) + par(nxp+vp+vp+i)
!               ELSEIF ((optpftg .AND. .NOT. optbsmg) .OR. (.NOT. optpftg .AND. optbsmg)) THEN
!                  pasmmax(i) = pasmmax0(n) + par(nxp+vp+i)
!               ELSE
!                  pasmmax(i) = pasmmax0(n)+ par(nxp+i)
!               ENDIF
!            ELSE
               pasmmax(i) = pasmmax0(n)
!            ENDIF
!         ENDIF
!      ENDIF

  END DO
! ... adjust bucket size for grass or crop sub-pixels if co-occurring with trees
  DO j = 1, ng
     IF (vg(j).GT.0) THEN
        pasmmax_tree = 0.
        pasmmax_herb = 0.
        frac_tree = 0.
        frac_herb = 0.
        DO k = vpb(j), vpe(j)
!           write(*,*) 'j k: ', j, k
           IF(pft(k).LE.6) THEN
              pasmmax_tree = pasmmax_tree + pasmmax(k)*frac(k)
              frac_tree = frac_tree + frac(k)
           ENDIF
           IF(pft(k)==9.OR.pft(k)==10.OR.pft(k)==13) THEN
              pasmmax_herb = pasmmax_herb + pasmmax(k)*frac(k)
              frac_herb = frac_herb + frac(k)
           ENDIF
        ENDDO
!        write(*,*) 'pasmmax_tree: ', pasmmax_tree/(frac_tree+1e-20)
!        write(*,*) 'pasmmax_herb: ', pasmmax_herb/(frac_herb+1e-20)
        IF (pasmmax_tree.GT.1E-20.AND.pasmmax_herb.GT.1E-20) THEN
           mult_tree = frac_tree * (pasmmax_tree+pasmmax_herb) &
                       / pasmmax_tree / (frac_tree+frtdh*frac_herb)
!           mult_herb = mult_herb * frtdh * frac_herb * pasmmax_tree &
!                       / pasmmax_herb / (frac_tree+1e-20)
           mult_herb = frtdh * frac_herb * (pasmmax_tree+pasmmax_herb) &
                       / pasmmax_herb / (frac_tree+frtdh*frac_herb)
!           write(*,*) 'mult_tree: ', mult_tree
!           write(*,*) 'mult_herb: ', mult_herb
                      DO k = vpb(j), vpe(j)
              IF(pft(k).LE.6) pasmmax(k) = pasmmax(k) * mult_tree
              IF(pft(k)==9.OR.pft(k)==10.OR.pft(k)==13) pasmmax(k) = pasmmax(k) * mult_herb
           ENDDO
        ENDIF
!        pasmmax_tree = 0.
!        pasmmax_herb = 0.
!        frac_tree = 0.
!        frac_herb = 0.
!        DO k = vpb(j), vpe(j)
!           IF(pft(k).LE.6) THEN
!              pasmmax_tree = pasmmax_tree + pasmmax(k)*frac(k)
!              frac_tree = frac_tree + frac(k)
!           ENDIF
!           IF(pft(k)==9.OR.pft(k)==10.OR.pft(k)==13) THEN
!              pasmmax_herb = pasmmax_herb + pasmmax(k)*frac(k)
!              frac_herb = frac_herb + frac(k)
!           ENDIF
!        ENDDO
!        write(*,*) 'AFTER: pasmmax_tree: ', pasmmax_tree/(frac_tree+1e-20)
!        write(*,*) 'AFTER: pasmmax_herb: ', pasmmax_herb/(frac_herb+1e-20)
     ENDIF
  ENDDO

  fbiasa = par(nxp-2)
  fbiasb = par(nxp-1)
  glob_offset = par(nxp)

  DEALLOCATE (paramap)

END SUBROUTINE mapping

END MODULE mo_mapping
