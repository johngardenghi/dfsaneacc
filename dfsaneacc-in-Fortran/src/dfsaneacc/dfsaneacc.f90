module bmdfsaneacc

  use bmqr, only: qr,extractQ,extractR,qrupdate,qrsolve
  use iso_c_binding, only: c_ptr

  implicit none

  interface
     subroutine rsub(n,x,res,pdataptr)
       use iso_c_binding, only: c_ptr
       integer, intent(in) :: n
       real(kind=8), intent(in) :: x(n)
       real(kind=8), intent(out) :: res(n)
       type(c_ptr), optional, intent(in) :: pdataptr
     end subroutine rsub
  end interface
  
  public :: dfsaneacc

contains

  ! ******************************************************************
  ! ******************************************************************

  subroutine dfsaneacc(evalr,n,x,res,epsf,maxit,ndiis,nhlim,iter,fcnt, &
       istop,pdataptr)

    ! PROCEDURE
    procedure(rsub) :: evalr

    ! SCALAR ARGUMENTS
    integer, intent(in) :: ndiis,maxit,n,nhlim
    integer, intent(out) :: fcnt,istop,iter
    real(kind=8), intent(in) :: epsf
    type(c_ptr), optional, intent(in) :: pdataptr

    ! ARRAY ARGUMENTS
    real(kind=8), intent(inout) :: x(n)
    real(kind=8), intent(out) :: res(n)
    
    ! PARAMETERS
    logical, parameter :: diis = .true., tabline = .false.
    integer, parameter :: M = 10, iprint = -1, maxitnoprogr = huge(1)
    real(kind=8), parameter :: gamma = 1.0d-04, sigma1 = 0.1d0, sigma2 = 0.5d0, &
         macheps = epsilon( 1.0d0 ), lamsgmin = sqrt( macheps ), lamsgmax = 1.0d0 / lamsgmin, &
         hsmall = 1.0d-01, hlarge = 0.1d0, hinit = 0.01d0
    
    ! LOCAL SCALARS
    logical :: extra
    integer :: backcnt,itnoprogr,jtbr,maxrank,ncol,nh,rank,nexti
    real(kind=8) :: eta,f,facc,fbest,fmax,ftrial,ftrial1,ftrial2,lambda1,lambda2,lamsg,ltmp
    
    ! LOCAL ARRAYS
    integer :: perm(nhlim-1)
    real(kind=8) :: col(n),lastfv(0:M-1),s(n),resacc(n),restmp(n),restrial(n),xacc(n), &
         xtmp(n),xtrial(n),xhist(n,nhlim),rhist(n,nhlim),q(n,nhlim-1),r(nhlim-1,nhlim-1), &
         gam(nhlim-1)

    iter = 0
    fcnt = 0
    itnoprogr = 0
    
    lastfv(0:M-1) = - huge( 1.0d0 )

    call evalr(n,x,res,pdataptr)
    fcnt = fcnt + 1

    f = dot_product( res(1:n), res(1:n) )
    fbest = f

    nh = 1
    xhist(1:n,1) = x(1:n)
    rhist(1:n,1) = res(1:n)

    maxrank = 0

    nexti = 1
    
    eta = min( 0.5d0 * f, sqrt( f ) )
    
    lamsg = 1.0d0

100 continue

    if ( iprint .ge. 0 ) then
       write(*,*)
     ! write(*,*) 'iter = ',iter,' f = ',f,' x = ',x(1:n),' res = ',res(1:n)
       write(*,*) 'iter = ',iter,' f = ',f
    end if

    if ( tabline ) then
       open(10,file='tabline-interrupted.txt')
       write(10,9000) n,9,sqrt(f),iter,fcnt
       close(10)
    end if
    
    if ( f .le. epsf ** 2 ) then
       if ( iprint .ge. 0 ) then
          write(*,*) 'success!'
       end if
       istop = 0
       return
    end if

    if ( iter .ge. maxit ) then
       if ( iprint .ge. 0 ) then
          write(*,*) 'maximum number of iterations reached!'
       end if
       istop = 1
       return
    end if
    
    if ( itnoprogr .ge. maxitnoprogr ) then
       if ( iprint .ge. 0 ) then
          write(*,*) 'no progress during maxitnoprogr consecutive iterations!'
       end if
       istop = 2
       return
    end if
    
    lastfv(mod(iter,M)) = f
    fmax = maxval( lastfv(0:M-1) )

    iter = iter + 1

    backcnt = 0
    
200 continue

    if ( backcnt .eq. 0 ) then
       lambda1 = lamsg
       
    else ! if ( backcnt .ge. 1 ) then
       ! minimizer of the quadratic q that satisfies q(0) = f,
       ! q(lambda1) = ftrial1, and q'(0) = f'(0). Assuming Jk =
       ! identity, f'(0) = 2 dk^t Jk Fk is approximated by - 2 f.
       ltmp = lambda1 ** 2 * f / ( ftrial1 + ( 2.0d0 * lambda1 - 1.0d0 ) * f )

       if ( isnan( ltmp ) ) then
          lambda1 = sigma1 * lambda1
       else
          lambda1 = max( sigma1 * lambda1, min( ltmp, sigma2 * lambda1 ) )
       end if
    end if
    
    xtrial(1:n) = x(1:n) - lambda1 * res(1:n)

    if ( any( xtrial(1:n) .ne. x(1:n) ) ) then
       call evalr(n,xtrial,restrial,pdataptr)
       fcnt = fcnt + 1

       ftrial1 = dot_product( restrial(1:n), restrial(1:n) )

       if ( iprint .ge. 1 ) then
          ! write(*,*) '    lambda = ',lambda1,' ftrial+ = ',ftrial1,' xtrial = ',xtrial(1:n),' restrial = ',restrial(1:n)
          write(*,*) '    lambda = ',lambda1,' ftrial+ = ',ftrial1
       end if
    
       if ( ftrial1 .le. fmax + eta - 2.0d0 * gamma * lambda1 ** 2 * f ) then
          ftrial = ftrial1
          go to 300
       end if
    end if

    if ( backcnt .eq. 0 ) then
       lambda2 = lamsg
       
    else ! if ( backcnt .ge. 1 ) then
       ! minimizer of the quadratic q that satisfies q(0) = f,
       ! q(lambda2) = ftrial2, and q'(0) = f'(0). Assuming Jk =
       ! identity, f'(0) = 2 dk^t Jk Fk is approximated by - 2 f.
       ltmp = lambda2 ** 2 * f / ( ftrial2 + ( 2.0d0 * lambda2 - 1.0d0 ) * f )

       if ( isnan( ltmp ) ) then
          lambda2 = sigma1 * lambda2
       else
          lambda2 = max( sigma1 * lambda2, min( ltmp, sigma2 * lambda2 ) )
       end if

       ! lambda2 = max( sigma1 * lambda2, min( ltmp, sigma2 * lambda2 ) )
    end if

    xtrial(1:n) = x(1:n) + lambda2 * res(1:n)

    if ( all( xtrial(1:n) .eq. x(1:n) ) ) then
       if ( iprint .ge. 0 ) then
          write(*,*) 'too small step in the line search'
       end if
       istop = 3
       return
    end if

    call evalr(n,xtrial,restrial,pdataptr)
    fcnt = fcnt + 1

    ftrial2 = dot_product( restrial(1:n), restrial(1:n) )
       
    if ( iprint .ge. 1 ) then
     ! write(*,*) '    lambda = ',lambda2,' ftrial- = ',ftrial2,' xtrial = ',xtrial(1:n),' restrial = ',restrial(1:n)
       write(*,*) '    lambda = ',lambda2,' ftrial- = ',ftrial2
    end if
    
    if ( ftrial2 .le. fmax + eta - 2.0d0 * gamma * lambda2 ** 2 * f ) then
       ftrial = ftrial2       
       go to 300
    end if

    backcnt = backcnt + 1
    
    go to 200

300 continue

    if ( .not. ( diis .and. mod( iter, ndiis ) .eq. 0 ) ) then
       go to 500
    end if
    
    if ( iprint .ge. 2 ) then
       write(*,*) 'Including xtrial in history.'
    end if
       
    if ( nh + 1 .le. min( nhlim, n + 1 ) ) then
       if ( iprint .ge. 2 ) then
          write(*,*) 'There is space to include xtrial without removing any other point ', &
               '(nh = ',nh,' min( nhlim, n + 1 ) = ',min( nhlim, n + 1 ),').'
       end if
       
       nh = nh + 1
       jtbr = 0
       
    else
       if ( iprint .ge. 2 ) then
          write(*,*) 'History is full. Oldest point is being removed ', &
               '(nh = ',nh,' min( nhlim, n + 1 ) = ',min( nhlim, n + 1 ),').'
       end if
       
       xhist(1:n,1:nh-1) = xhist(1:n,2:nh)
       rhist(1:n,1:nh-1) = rhist(1:n,2:nh)
       jtbr = 1
    end if

    xhist(1:n,nh) = xtrial(1:n)
    rhist(1:n,nh) = restrial(1:n)
       
    col(1:n) = rhist(1:n,nh) - rhist(1:n,nh-1)

    extra = .false.
       
    if ( nh .eq. 2 ) then
       ncol = 1
       call qr(n,ncol,n,col,gam,perm,rank)
       call extractQ(n,ncol,n,col,gam,rank,n,Q)
       call extractR(ncol,n,col,rank,nhlim-1,R)
       maxrank = max( maxrank, rank )
    else
       call qrupdate(n,ncol,rank,nhlim-1,jtbr,n,Q,nhlim-1,R,perm,col)
       maxrank = max( maxrank, rank )

       if ( iprint .ge. 2 ) then
          write(*,*) 'Number of columns = ',ncol,' rank = ',rank,' (largest achieved rank = ',maxrank,')'
       end if

       if ( ncol .ne. nh - 1 ) then
          write(*,*) 'ERROR in subroutine dfsane: ncol is not equal to nh - 1! ncol = ',ncol,' nh = ',nh
          stop
       end if
          
       if ( rank .lt. maxrank ) then
          if ( iprint .ge. 2 ) then
             write(*,*) 'Current rank is smaller that largest achieved rank. Thus, an extra point will be considerd.'
          end if
             
          extra = .true.
             
          xtmp(1:n) = x(1:n)
          xtmp(nexti) = x(nexti) + hsmall * max( 1.0d0, abs( x(nexti) ) )
             
          nexti = mod( nexti, n ) + 1
             
          call evalr(n,xtmp,restmp,pdataptr)
          fcnt = fcnt + 1
             
          if ( iprint .ge. 2 ) then
             write(*,*) 'Including extra point in history.'
          end if
             
          if ( nh + 1 .le. min( nhlim, n + 1 ) ) then
             if ( iprint .ge. 2 ) then
                write(*,*) 'There is space to include the extra point without removing any other point ', &
                     '(nh = ',nh,' min( nhlim, n + 1 ) = ',min( nhlim, n + 1 ),').'
             end if
             
             nh = nh + 1
             jtbr = 0
          else
             if ( iprint .ge. 2 ) then
                write(*,*) 'History is full. Oldest point is being removed ', &
                     '(nh = ',nh,' min( nhlim, n + 1 ) = ',min( nhlim, n + 1 ),').'
             end if
                
             xhist(1:n,1:nh-1) = xhist(1:n,2:nh)
             rhist(1:n,1:nh-1) = rhist(1:n,2:nh)
             jtbr = 1
          end if

          xhist(1:n,nh) = xtmp(1:n)
          rhist(1:n,nh) = restmp(1:n)

          col(1:n) = rhist(1:n,nh) - rhist(1:n,nh-1)

          call qrupdate(n,ncol,rank,nhlim-1,jtbr,n,Q,nhlim-1,R,perm,col)
          maxrank = max( maxrank, rank )

          if ( iprint .ge. 2 ) then
             write(*,*) 'Number of columns = ',ncol,' rank = ',rank,' (largest achieved rank = ',maxrank,')'
          end if
          
          if ( ncol .ne. nh - 1 ) then
             write(*,*) 'ERROR in subroutine dfsane: ncol is not equal to nh - 1! ncol = ',ncol,' nh = ',nh
             stop
          end if
          
       end if
    end if

    if ( rank .gt. 0 ) then
       call qrsolve(n,ncol,rank,n,Q,nhlim-1,R,perm,restrial,xacc)
       
       xacc(1:n) = xtrial(1:n) - matmul( xhist(1:n,2:nh) - xhist(1:n,1:nh-1), xacc(1:nh-1) )

       if ( extra ) then
          if ( iprint .ge. 2 ) then
             write(*,*) 'The extra point is being removed.'
          end if
          
          nh = nh - 1
          
          jtbr = ncol
          
          call qrupdate(n,ncol,rank,nhlim-1,jtbr,n,Q,nhlim-1,R,perm)
          maxrank = max( maxrank, rank )
          
          if ( iprint .ge. 2 ) then
             write(*,*) 'Number of columns = ',ncol,' rank = ',rank,' (largest achieved rank = ',maxrank,')'
          end if
          
          if ( ncol .ne. nh - 1 ) then
             write(*,*) 'ERROR in subroutine dfsane: ncol is not equal to nh - 1! ncol = ',ncol,' nh = ',nh
             stop
          end if
       end if

       if ( norm2( xacc(1:n) ) .le. 10.0d0 * max( 1.0d0, norm2( x(1:n) ) ) ) then
          call evalr(n,xacc,resacc,pdataptr)
          fcnt = fcnt + 1

          facc = dot_product( resacc(1:n), resacc(1:n) )
          
          if ( iprint .ge. 1 ) then
             ! write(*,*) '    facc = ',facc,' xacc = ',xacc(1:n),' resacc = ',resacc(1:n)
             write(*,*) '    facc = ',facc
          end if
          
          if ( facc .lt. ftrial .and. norm2( xacc(1:n) - x(1:n) ) .gt. macheps * max( 1.0d0, norm2( x(1:n) ) ) ) then
             if ( iprint .ge. 2 ) then
                write(*,*) 'Accelerated point is being accepted; and history will be updated substituting xtrial by xacc.'
             end if
             
             xhist(1:n,nh) = xacc(1:n)
             rhist(1:n,nh) = resacc(1:n)
             
             jtbr = ncol
             
             col(1:n) = rhist(1:n,nh) - rhist(1:n,nh-1)
             
             call qrupdate(n,ncol,rank,nhlim-1,jtbr,n,Q,nhlim-1,R,perm,col)
             maxrank = max( maxrank, rank )
             
             if ( iprint .ge. 2 ) then
                write(*,*) 'Number of columns = ',ncol,' rank = ',rank,' (largest achieved rank = ',maxrank,')'
             end if
             
             if ( ncol .ne. nh - 1 ) then
                write(*,*) 'ERROR in subroutine dfsane: ncol is not equal to nh - 1! ncol = ',ncol,' nh = ',nh
                stop
             end if
             
             xtrial(1:n) = xacc(1:n)
             restrial(1:n) = resacc(1:n)
             ftrial = facc
          end if
       end if
    end if
    
    if ( rank .eq. 0 ) then
       if ( iprint .ge. 0 ) then
          write(*,*) 'Acceleration matrix has null rank; thus, it will be rebuild from scratch.'
          write(*,*) 'It will be populated with neighbours of the current point taking a small step into coordinate directions.'
       end if
       
       nh = 0
          
400    continue
       if ( nh + 1 .lt. min( nhlim, n + 1 ) ) then
          xtmp(1:n) = xtrial(1:n)
          xtmp(nexti) = xtrial(nexti) + hlarge * max( 1.0d0, abs( xtrial(nexti) ) )
             
          nexti = mod( nexti, n ) + 1
          
          call evalr(n,xtmp,restmp,pdataptr)
          fcnt = fcnt + 1
             
          nh = nh + 1
             
          xhist(1:n,nh) = xtmp(1:n)
          rhist(1:n,nh) = restmp(1:n)

          if ( nh .ge. 2 ) then
             col(1:n) = rhist(1:n,nh) - rhist(1:n,nh-1)

             if ( nh .eq. 2 ) then
                ncol = 1
                call qr(n,ncol,n,col,gam,perm,rank)
                call extractQ(n,ncol,n,col,gam,rank,n,Q)
                call extractR(ncol,n,col,rank,nhlim-1,R)
                maxrank = max( maxrank, rank )
                   
             else
                jtbr = 0
                call qrupdate(n,ncol,rank,nhlim-1,jtbr,n,Q,nhlim-1,R,perm,col)
                maxrank = max( maxrank, rank )
          
                if ( iprint .ge. 2 ) then
                   write(*,*) 'Number of columns = ',ncol,' rank = ',rank,' (largest achieved rank = ',maxrank,')'
                end if
          
                if ( ncol .ne. nh - 1 ) then
                   write(*,*) 'ERROR in subroutine dfsane: ncol is not equal to nh - 1! ncol = ',ncol,' nh = ',nh
                   stop
                end if
             end if
          end if
          
          go to 400
       end if

       nh = nh + 1

       xhist(1:n,nh) = xtrial(1:n)
       rhist(1:n,nh) = restrial(1:n)

       jtbr = 0

       col(1:n) = rhist(1:n,nh) - rhist(1:n,nh-1)

       call qrupdate(n,ncol,rank,nhlim-1,jtbr,n,Q,nhlim-1,R,perm,col)
       maxrank = max( maxrank, rank )

       call qrsolve(n,ncol,rank,n,Q,nhlim-1,R,perm,restrial,xacc)
       
       xacc(1:n) = xtrial(1:n) - matmul( xhist(1:n,2:nh) - xhist(1:n,1:nh-1), xacc(1:nh-1) )

       call evalr(n,xacc,resacc,pdataptr)
       fcnt = fcnt + 1
       
       facc = dot_product( resacc(1:n), resacc(1:n) )

       if ( iprint .ge. 1 ) then
        ! write(*,*) '    facc = ',facc,' xacc = ',xacc(1:n),' resacc = ',resacc(1:n)
          write(*,*) '    facc = ',facc
       end if

       if ( facc .lt. ftrial .and. norm2( xacc(1:n) - x(1:n) ) .gt. macheps * max( 1.0d0, norm2( x(1:n) ) ) ) then
          xhist(1:n,nh) = xacc(1:n)
          rhist(1:n,nh) = resacc(1:n)
       
          jtbr = ncol

          col(1:n) = rhist(1:n,nh) - rhist(1:n,nh-1)

          call qrupdate(n,ncol,rank,nhlim-1,jtbr,n,Q,nhlim-1,R,perm,col)
          maxrank = max( maxrank, rank )
                    
          xtrial(1:n) = xacc(1:n)
          restrial(1:n) = resacc(1:n)
          ftrial = facc
       end if
    end if

500 continue
    
    ! Formula of lamsg below is, in theory, equal to: 
    ! lamsg = lambda * f / ( dot_product( res(1:n), restrial(1:n) ) - f )
    ! However, in practive, it presented better results.

    s(1:n) = xtrial(1:n) - x(1:n)
    ! lamsg = abs( dot_product( s(1:n), s(1:n) ) / dot_product( s(1:n), restrial(1:n) - res(1:n) ) )
    lamsg = dot_product( s(1:n), s(1:n) ) / dot_product( s(1:n), restrial(1:n) - res(1:n) )
    if ( .not. ( lamsgmin .le. abs( lamsg ) .and. abs( lamsg ) .le. 1.0 ) ) then
       lamsg = max( lamsgmin, min( norm2( xtrial(1:n) ) / sqrt( ftrial ), lamsgmax ) )
       ! lamsg = max( lamsgmin, min( max( 0.01d0 * norm2( xtrial(1:n) ), norm2( s(1:n) ) ) / sqrt( ftrial ), lamsgmax ) )
       ! lamsg = max( lamsgmin, min( max( 0.1d0 * norm2( xtrial(1:n) ), norm2( s(1:n) ) ) / sqrt( ftrial ), lamsgmax ) )
       ! lamsg = max( lamsgmin, min( max( norm2( xtrial(1:n) ), norm2( s(1:n) ) ) / sqrt( ftrial ), lamsgmax ) )
    end if

    ! s(1:n) = xtrial(1:n) - x(1:n)
    ! lamsg = abs( dot_product( s(1:n), s(1:n) ) / dot_product( s(1:n), restrial(1:n) - res(1:n) ) )
    ! if ( .not. ( lamsgmin .le. lamsg .and. lamsg .le. 1.0d0 ) ) then
    !    lamsg = max( lamsgmin, min( norm2( xtrial(1:n) ) / sqrt( ftrial ), lamsgmax ) )
    ! end if

    ! s(1:n) = xtrial(1:n) - x(1:n)
    ! lamsg = max( sqrt( macheps ), sqrt( macheps ) * hinit * norm2( xtrial(1:n) ), hinit * norm2( s(1:n) ) )
    ! lamsg = max( macheps, min(  lamsg / sqrt( ftrial ), 1.0d0 / macheps ) )

    ! print *, "ftrial = ", ftrial, ", f = ", f, ", ftrial - f = ", ftrial - f, &
    !      "frial - f .lt epsf? ", ftrial - f .lt. epsf
    
    if ( .not. ( sqrt(ftrial) .lt. sqrt(f) - 1.0d-8 * sqrt(f) ) ) then
       itnoprogr = itnoprogr + 1
    else
       itnoprogr = 0
    end if
    
    eta = 0.5d0 * eta
    
    x(1:n) = xtrial(1:n)
    res(1:n) = restrial(1:n)
    f = ftrial

    go to 100

9000 format(I6,1X,I2,1X,1P,D7.1,2(1X,I10))

  end subroutine dfsaneacc

end module bmdfsaneacc
