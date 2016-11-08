program dline
  implicit none
  integer :: n

! get the number of control variables
  call numbmod( n )

! action
  call dlinerun( n )

end program dline

subroutine dlinerun( n )
  implicit none

  integer :: n
  integer, parameter  :: uopt = 7        ! unit for options
  integer, parameter  :: usect = 22      ! unit for section
  integer, parameter  :: upar1 = 23      ! unit for parameter values 
  integer, parameter  :: upar2 = 24      ! unit for parameter values 
  integer, parameter  :: udx = 25        ! unit for search direction
  real    :: fc, fc0
  real    :: x0(n), xend(n), x(n), t, dx(n), ndx
  real    :: x_tl(n)
  real    :: fc_tl, fc0_tl

  integer :: nt = 10                     ! number of steps in search direction
  integer :: its=-2,ite=10               ! start and end point 
  integer :: ixs=1,ixe                   ! use search direction projected onto range ixs to ixe
  real :: a = 1.                         ! scaling factor default: 1.

  integer :: i, it, ix
  namelist/nlist/a,nt,its,ite,ixs,ixe

! check for too large n
  if (n.gt.99) print *, 'Line: n = ',n,' too large, change print format'

! default values for options
  ixe = n

! read namelist for options
      open( uopt, file   = 'line.par', access = 'sequential', form = 'formatted'  )
      read( uopt, nlist, end=9000, err=8000 )
      print *, ' Line : options have been read'
 8000 continue
 9000 continue
      close( uopt )

! initialisize the model
  call initmod( n, x0 )
  
  print*, 'Read previous point from line search '
  OPEN (39,file='lx.b',status='unknown',form='unformatted')
  READ(39) xend
  CLOSE(39)

  print*, 'Read final point '
  OPEN (39,file='xf.b',status='unknown',form='unformatted')
  READ(39) x0
  CLOSE(39)

  ixe = min(n,ixe)
  ixs = max(1,ixs)

! compute norm of delta x
  ndx =0.
  dx = 0.
  do ix=ixs,ixe
     ndx=ndx+(xend(ix)-x0(ix))**2
     dx(ix) = xend(ix)-x0(ix)
  enddo
  ndx = sqrt(ndx)

!  open units
  open (usect,file='dline.set',form='formatted')
  open (upar1,file='parset.dat',form='formatted')
  open (upar2,file='parlist.dat',form='formatted')
  open (udx,file='dx.dat',form='formatted')

!  dx
  do ix=1,n
     write (udx,*) ix, dx(ix)
  enddo

!  xmgrace header
  open (usect,file='dline.set',form='formatted')
  write (usect,'(A)') '@WITH G00'
  write (usect,'(A)') '@G00 ON'
  write (usect,'(A)') '@TYPE xy'
  write (usect,'(A,E16.8)') '@WORLD XMAX', ndx
  write (usect,'(A)') '@s0 symbol 1'
  write (usect,'(A)') '@s0 line type 0'
  write (usect,'(A)') '@xaxis  label "norm of difference of control vector to base point"'
  write (usect,'(A)') '@TITLE "directional derivative of cost function"'

  write (upar1,'(A,A19,99I20)') '#', 't', (i,i=1,n)

! initialise tangent linear
  x_tl = dx/ndx
  fc0_tl = 0.

! this is the loop that does the tests
do it=its,ite
     t=a*float(it)/nt
     x=x0+t*dx
     x_tl = dx/ndx
     fc = 0.
     fc_tl = 0.
     call model_tl( n, x, x_tl, fc, fc_tl )
     write (usect,*) t*ndx,fc_tl-fc0_tl
     write (upar1,'(100e20.6)') t*ndx,x(:)
     write (upar2,*) '# parameter vector for it = ',it
     do i=1,n
        write (upar2,*) i, x(i)
     enddo
  enddo

!  close units
  close (usect)
  close (upar1)
  close (upar2)
  close (udx)
end subroutine dlinerun
