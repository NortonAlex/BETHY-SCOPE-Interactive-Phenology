PROGRAM inv_hess
  ! inverts a symmetric matrix and prints out roots of diagonals as well as saving matrix
  ! really a wrapper for symmetric inverse routines from lapack
  ! author Peter Rayner
  ! date:  Sep 2001
  ! modified, Aug 2004, Marko Scholze
  ! modified, Oct 2003, HEW
  
USE mo_namelist
USE mo_trafo

  IMPLICIT NONE 

  REAL, ALLOCATABLE, DIMENSION( :, :) :: hessinverse, hess, hessinverse_denoised, hess_denoised, e, einverse, aux, evec
  REAL, DIMENSION(:), ALLOCATABLE :: evals
  REAL, DIMENSION(1000000) ::  work
  REAL, ALLOCATABLE, DIMENSION(:) :: x, p0, p, p0su, a, pps, pms, xs, xpfa, xpfb, x0
  Integer, ALLOCATABLE, DIMENSION(:) :: xpf
  real, parameter :: eps = 1.e-12
  REAL :: fc
  INTEGER, DIMENSION(1000000) :: lwork
  INTEGER n, i, j, info, k ,l, cutoff
  INTEGER :: ifail
  CHARACTER :: header, w
  CHARACTER (LEN= 20) :: filename, progName

! .. read namelist variables
  CALL get_namelist

! .. read prior parameter file
  OPEN(1,file=param_file,status='old')
  READ(1,*)
  READ(1,*) n

  ALLOCATE(x(n), p0(n), p(n), p0su(n), xpf(n), a(n), pps(n), &
       &      pms(n), xs(n), xpfa(n), xpfb(n), x0(n))
  a = 0.
  DO i=1,n
     READ(1,*) p0(i),p0su(i),xpf(i),xpfa(i),xpfb(i)
     IF (xpf(i)<0) THEN
        a(i) = LOG(p0su(i)+p0(i)) - LOG(p0(i))
     ELSE IF (xpf(i) == 2) THEN
        a(i) = ((p0su(i))/(SQRT(p0su(i)+p0(i)) + SQRT(p0(i))))**2
     ELSE IF (xpf(i) == 4) THEN
        a(i) = -log(xpfb(i)/(p0su(i)+p0(i)-xpfa(i)) - 1) + log(xpfb(i)/(p0(i)-xpfa(i)) - 1)
     ELSE
        a(i) = 0.
     ENDIF
  END DO
  CLOSE(1)

  ifail=0
  OPEN (39,file='xf.b',status='old',form='unformatted',iostat=ifail)
  if (ifail.gt.0) then
     print*,'cost: parameter values set by initmod'
  else
     read(39) x
     print*,'cost: parameter values overwritten by those from x.b'
  endif
  CLOSE(39)

  ALLOCATE( hessinverse(n, n),  aux(n, n), hess_denoised(n, n), hessinverse_denoised(n, n), evec(n, n), e(n, n), einverse(n, n), &
       & evals(n), hess(n,n))

! ... read hessian
  OPEN(unit=2,file='output/hess.bin',form='unformatted')
  READ (2) hess
  CLOSE(2)

! to check neg ev for parameter 31
!  hess(31,:)=0.
!  hess(:,31)=0.
!  hess(31,31)=1.

  PRINT*,'symmetry check ',MINVAL(hess -TRANSPOSE(hess)),MAXVAL(hess -TRANSPOSE(hess))

  ! calculate eigen spectrum 
  PRINT*,'eigenvalues of the Hessian in x-space'

  aux = hess
  CALL dsyev('v', 'u', n, aux, n, evals, work, 100000, info)
  evec = aux

  ! least resolved direction
  OPEN(unit=2,file='pert-smallest-ev.b',form='unformatted')
  write (2) x + evec(:,1)
  CLOSE(2)
  WRITE(6,'(a6,2a16)') 'Least Resolved Direction'
  WRITE(6,'(a6,2a16)') '#','Component'
  DO i=1,n
     WRITE(6,'(i6,2f16.8)') i, evec(i,1)
  END DO

  ! write out e-vals
  WRITE(6,'(a6,2a16)') '#','Eval','Unc red (%)'
  DO i=1,n
     WRITE(6,'(i6,2f16.8)') i,evals(i),100.*(1.-1./sqrt(evals(i)))
  END DO

  einverse=0.
  DO i = 1, n
     if (abs(evals( i)).lt.eps) then 
        print *, 'eigenvalue # ',i,' too close to 0 = ', evals(i)
        print *, 'cannot invert matrix '
        goto 100
     endif
     einverse( i, i) = 1/evals( i)
  END DO
  hessinverse =  MATMUL( evec, MATMUL( einverse, TRANSPOSE( evec)))
  aux = MATMUL( hess, hessinverse) -RESHAPE ( (/ (1.0, (0.0, k= 1,n), j= 1,n-1), 1.0 /), (/ n, n /) )
  PRINT*,'inversion check',MINVAL(aux),MAXVAL(aux)
100 continue

  PRINT*,'sqrt of diagonal of inverse Hessian (parameter uncertainties)'
  DO i=1,n
     WRITE(6,*) i, SQRT(ABS(hessinverse( i,i)))
  END DO

  ! make hessian positive definite, assume that largest neg eigenvalue represent noise level
! cutoff = SUM( MINLOC( evals, (evals > ABS(evals(1)))))    ! noise floor set to abs(lambda_1)
 cutoff = SUM( MINLOC( evals, (evals > 1.+(1.-evals(1))))) ! noise floor set to 1. 1- lambda_1
! cutoff = SUM( MINLOC( evals, (evals > 1.)))               ! noise floor set to 1. 
! if (minval(evals) .gt. 0.) cutoff=1                       ! no noise
  PRINT*,'noise floor at eval ',cutoff
  ! now reconstruct hessian as observed part (> cutoff) + prior part *< cutoff)
  e = 0.
  ! insert eigen-values into diagonal matrix
  DO i = cutoff, n
     e( i, i) = evals( i)
  END DO
  ! and now add prior part, evals here are always 1
  DO i = 1, cutoff -1
     e( i, i) = 1.
  END DO
  hess_denoised =  MATMUL( evec, MATMUL( e, TRANSPOSE( evec)))
 
  OPEN(unit=1, file='output/denoised_hess_x.bin', form='unformatted')
  REWIND 1
  WRITE (1) hess_denoised
  CLOSE(1)

! overwrite eigenvalues by their reciprocals to build inverse
  e = 0.
  ! insert eigen-values into diagonal matrix
  DO i = cutoff, n
     e( i, i) = 1/evals( i)
  END DO
  ! and now add prior part, evals here are always 1
  DO i = 1, cutoff -1
     e( i, i) = 1.
  END DO
  hessinverse_denoised =  MATMUL( evec, MATMUL( e, TRANSPOSE( evec)))

  OPEN(unit=1, file='output/denoised_inverse_hess_x.bin', form='unformatted')
  REWIND 1
  WRITE (1) hessinverse_denoised
  CLOSE(1)

  DO i=1,n
     xs(i) = SQRT(ABS(hessinverse_denoised( i,i)))
  END DO

! produce prior value
  x0 = p2x(p0,xpf,xpfa,xpfb,p0,p0su,a,n)
  
  p = x2p(x, x0, xpf, xpfa, xpfb, p0, p0su, a, n)
  pps = x2p(x+xs, x0, xpf, xpfa, xpfb, p0, p0su, a, n)
  pms = x2p(x-xs, x0, xpf, xpfa, xpfb, p0, p0su, a, n)

  PRINT*,'sqrt diagonal elements of inverse denoised Hessian (parameter unc)'
  OPEN(unit=2, file='params_post.txt', form='formatted')
  REWIND 2
  OPEN(unit=3, file='x_post.txt', form='formatted')
  REWIND 3
  WRITE(2,'(a,2x,3a18)') '#','Parameter ','upper 67%ile','lower 67%ile'
  WRITE(3,'(a,5x,3a17)') '#','Parameter (x)','increment (x)','Unc red (%)'
!  WRITE(3,'(a6,3a16)') '#','Parameter (x)','increment (x)','Unc red (%)'
  WRITE(6,'(a6,3a17)') '#','Parameter (x)','increment (x)','Unc red (%)'
  DO i=1,n
     WRITE(6,'(i6,3f17.8)') i,x(i),x(i)-x0(i),100.*(1.-xs(i))
     WRITE(3,'(i6,3f17.8)') i,x(i),x(i)-x0(i),100.*(1.-xs(i))
     write(2,'(i3,3(e18.10))') i,p(i),pms(i),pps(i)
  END DO
  CLOSE(2)
  CLOSE(3)

  STOP
  
CONTAINS 

  ! syminv calls inverting routines from LAPACK ??
  SUBROUTINE syminv(a, info)
    IMPLICIT NONE
    REAL :: a(:,:)
    INTEGER :: info ! status return
    ! local variables
    INTEGER :: n ! dimension
    INTEGER :: i,j ! index variables
    REAL, DIMENSION(:), ALLOCATABLE :: work    
    INTEGER, DIMENSION(:), ALLOCATABLE :: iwork,ipiv
    
    REAL anorm ! 1-norm of A from slansy
    REAL rcond ! reciprocal condition number
    REAL slansy ! anorm function

    n = SIZE(a, 1)
    ALLOCATE(work(n*n),ipiv(n*n),iwork(n*n))
    
    ! routine for inverting a real symmetric positive definite matrix
    ! e.g. a covariance matrix (I hope)
    
    CALL dsytrf ('u', n, A, n, IPIV, WORK, n*n, info)
    CALL dsytri ('u', n, a, n, IPIV, WORK, info)
    DO j =2,n
       DO i=1,j
          a(j,i) = a(i,j)
       END DO
    END DO
    RETURN
    
  END SUBROUTINE syminv
END PROGRAM inv_hess
