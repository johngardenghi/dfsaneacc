program dfsaneaccma

  use bmdfsaneacc, only: dfsaneacc
  use iso_c_binding, only: c_ptr, c_loc, c_f_pointer

  implicit none

  type :: pdata_type
     integer :: counters(1) = 0
     character(len=50) :: pname
     real(kind=8), allocatable :: crhs(:)
  end type pdata_type

  ! LOCAL SCALARS
  character(len=50) :: pname
  integer :: allocstat,fcnt,ndiis,istop,iter,maxit,n,nhlim
  real(kind=8) :: epsf
  real :: finish,start
  type(pdata_type), target :: pdata

  ! LOCAL ARRAYS
  real(kind=8), allocatable :: r(:),x(:)

  call iniprob(n,x,c_loc(pdata))

  allocate(r(n), stat=allocstat)
  if ( allocstat .ne. 0 ) then
     write(*,*) 'Error in allocate.'
     stop
  end if

  epsf = 1.0d-06 * sqrt( real( n ) )

  maxit = huge( 1 ) ! In this case, A CPU time limit may be imposed
  ndiis = 1         ! accelerate every 'ndiis' iterations;
  nhlim = 6         ! points to enter in history

  call evalr(n,x,r,c_loc(pdata))

  open(10,file='tabline.txt')
  write(10,9000) pdata%pname,n,9,norm2(r(1:n)),0,0,0.0
  close(10)

  call cpu_time(start)
  call dfsaneacc(evalr,n,x,r,epsf,maxit,ndiis,nhlim,iter,fcnt,istop,c_loc(pdata))
  call cpu_time(finish)

  write(*,*) 'Satisfied stopping criterion = ',istop
  write(*,*) 'Residual Euclidian norm = ',norm2( r(1:n) )
  write(*,*) 'Number of iterations = ',iter
  write(*,*) 'Number of functional evaluations = ',fcnt
  write(*,*) 'CPU time in seconds = ',finish - start

  open(10,file='tabline.txt')
  write(10,9000) pdata%pname,n,istop,norm2(r(1:n)),iter,fcnt,finish-start
  close(10)

  deallocate(x,r,stat=allocstat)
  if ( allocstat .ne. 0 ) then
     write(*,*) 'Deallocation error.'
     stop
  end if

  call endprob(c_loc(pdata))

9000 format(A10,1X,I6,1X,I2,1X,1P,D7.1,2(1X,I10),1X,0P,F12.6)

contains

  ! ******************************************************************
  ! ******************************************************************

  subroutine iniprob(n,x,pdataptr)

    ! SCALAR ARGUMENTS
    integer, intent(out) :: n
    type(c_ptr), optional, intent(in) :: pdataptr

    ! ARRAY ARGUMENTS
    real(kind=8), allocatable :: x(:)

    ! LOCAL SCALARS
    integer :: allocstat,m,status,e_order,l_order,v_order

    n = 3
    allocate( x(n), stat=allocstat )
    if ( allocstat .ne. 0 ) then
       write(*,*) 'Error in allocate statement!'
       stop
    end if

    x(1:n) = 1.0d0 / ( real(n)**2 )

  end subroutine iniprob

  ! ******************************************************************
  ! ******************************************************************

  subroutine endprob(pdataptr)

    ! SCALAR ARGUMENTS
    type(c_ptr), optional, intent(in) :: pdataptr

    ! LOCAL SCALARS
    integer :: allocstat,status
    type(pdata_type), pointer :: pdata

  end subroutine endprob

  ! ******************************************************************
  ! ******************************************************************

  subroutine evalr(n,x,r,pdataptr)

    ! SCALAR ARGUMENTS
    integer, intent(in) :: n
    type(c_ptr), optional, intent(in) :: pdataptr

    ! ARRAY ARGUMENTS
    real(kind=8), intent(in) :: x(n)
    real(kind=8), intent(out) :: r(n)

    ! LOCAL SCALARS
    integer :: i

    r(1) = exp( x(1) ) - 1.0d0
    do i = 2, n
       r(i) = real(i)/1.0d+1 * ( exp( x(i) ) + x(i-1) - 1.0d0 )
    end do

  end subroutine evalr

end program dfsaneaccma
