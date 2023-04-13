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
    type(pdata_type), pointer :: pdata

    ! LOCAL ARRAYS
    logical, allocatable :: equatn(:),linear(:)
    real(kind=8), allocatable :: cl(:),l(:),u(:),lambda(:)

    call c_f_pointer(pdataptr,pdata)

    open(10,file='OUTSDIF.d',form='formatted',status='old')

    pdata%pname(11:50) = ' '

    rewind 10

    call cutest_pname(status,10,pdata%pname)
    if ( status .ne. 0 ) then
       write(*,*) 'INIPROB ERROR: cutest_pname - status not equal to zero. status = ',status
       stop
    end if

    call cutest_cdimen(status,10,n,m)
    if ( status .ne. 0 ) then
       write(*,*) 'INIPROB ERROR: cutest_cdimen - status not equal to zero. status = ',status
       stop
    end if

    if ( n .ne. m ) then
       write(*,*) 'INIPROB ERROR: This program tackles nonlinear systems with n=m only!',pdata%pname
       stop
    end if

    allocate(l(n),u(n),lambda(m),cl(m),equatn(m),linear(m),x(n),pdata%crhs(m),stat=allocstat)
    if ( allocstat .ne. 0 ) then
       write(*,*) 'INIPROB ERROR: Allocation error.'
       stop
    end if

    e_order = 1
    l_order = 0
    v_order = 0

    call CUTEST_csetup(status,10,6,20,n,m,x,l,u,lambda,cl,pdata%crhs,equatn,linear,e_order,l_order,v_order)
    if ( status .ne. 0 ) then
       write(*,*) 'INIPROB ERROR: cutest_csetup - status not equal to zero. status = ',status
       stop
    end if

    if ( any( .not. equatn(1:m) ) ) then
       write(*,*) 'INIPROB ERROR: The problem has inequalities!'
       stop
    end if

    if ( any( cl(1:m) .ne. pdata%crhs(m) ) ) then
       write(*,*) 'INIPROB ERROR: The problem has an equality with cl not equal to cu!'
       stop
    end if

    if ( any( l(1:n) .gt. - 1.0d+20 ) .or. any( u(1:n) .lt. 1.0d+20 ) ) then
       write(*,*) 'INIPROB ERROR: The problem has bounds!'
       stop
    end if

    close(10)

    deallocate(l,u,lambda,cl,equatn,linear,stat=allocstat)
    if ( allocstat .ne. 0 ) then
       write(*,*) 'INIPROB ERROR: Deallocation error.'
       stop
    end if

  end subroutine iniprob

  ! ******************************************************************
  ! ******************************************************************

  subroutine endprob(pdataptr)

    ! SCALAR ARGUMENTS
    type(c_ptr), optional, intent(in) :: pdataptr

    ! LOCAL SCALARS
    integer :: allocstat,status
    type(pdata_type), pointer :: pdata

    call c_f_pointer(pdataptr,pdata)

    deallocate(pdata%crhs,stat=allocstat)
    if ( allocstat .ne. 0 ) then
       write(*,*) 'Error in allocate statement!'
       stop
    end if

    call cutest_cterminate(status)

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
    integer :: jfun(1),jvar(1),jnnz,status
    real(kind=8) :: jval(1)
    type(pdata_type), pointer :: pdata

    call c_f_pointer(pdataptr,pdata)

    pdata%counters(1) = pdata%counters(1) + 1

    call cutest_ccfsg(status,n,n,x,r,jnnz,1,jval,jvar,jfun,.false.)
    if ( status .ne. 0 ) then
       write(*,*) 'INIPROB ERROR: There was a nonnull flag in a call to CUTEst routine cutest_ccfsg'
       stop
    end if

    r(1:n) = r(1:n) - pdata%crhs(1:n)

  end subroutine evalr

end program dfsaneaccma
