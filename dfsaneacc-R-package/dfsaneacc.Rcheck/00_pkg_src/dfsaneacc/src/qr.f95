module bmqrf
  use, intrinsic :: iso_c_binding

  implicit none

contains

  subroutine bmqr(n,m,mlim,lda,A,gam,p,rank) bind(C, name="bmqr_")

    ! SCALAR ARGUMENTS
    integer(c_int), intent(in) :: lda,m,mlim,n
    integer(c_int), intent(out) :: rank

    ! ARRAY ARGUMENTS
    integer(c_int), intent(out) :: p(mlim)
    real(c_double), intent(inout) :: A(lda,m)
    real(c_double), intent(out) :: gam(m)

    ! PARAMETERS
    real(c_double) :: tol = 100.0d0 * epsilon( 1.0d0 )

    ! LOCAL SCALARS
    integer(c_int) :: imax,itmp,j,k,kmax
    real(c_double) :: tau,tmp

    ! LOCAL ARRAYS
    real(c_double) :: c(m),cnrm(m+1),cnrmfull(m+1),v(n)

    p(1:m) = (/ (j,j=1,m) /)

    do k = 1,m
       c(k) = maxval( abs( A(1:n,k) ) )
       if ( c(k) .gt. tol ) then
          A(1:n,k) = A(1:n,k) / c(k)
       end if
       cnrm(k) = dot_product( A(1:n,k), A(1:n,k) )
       cnrmfull(k) = cnrm(k)
    end do
    cnrm(m+1) = 0.0d0
    cnrmfull(m+1) = 0.0d0

    k = 0

    kmax = maxloc( cnrm(1:m+1), dim=1 )

    do while ( k .lt. m .and. cnrm(kmax) .ge. tol * max( cnrmfull(kmax), 1.0d0 ) )

       k = k + 1

       if ( kmax .ne. k ) then
          itmp = p(k)
          p(k) = p(kmax)
          p(kmax) = itmp

          v(1:n) = A(1:n,k)
          A(1:n,k) = A(1:n,kmax)
          A(1:n,kmax) = v(1:n)

          tmp = c(k)
          c(k) = c(kmax)
          c(kmax) = tmp

          tmp = cnrm(k)
          cnrm(k) = cnrm(kmax)
          cnrm(kmax) = tmp
       end if

       tau = sqrt( cnrm(k) )

       if ( A(k,k) .lt. 0.0d0 ) then
          tau = - tau
       end if

       A(k,k) = A(k,k) + tau

       A(k+1:n,k) = A(k+1:n,k) / A(k,k)

       gam(k) = A(k,k) / tau

       A(k,k) = - tau

       do j = k + 1,m
          tmp = gam(k) * dot_product( (/ 1.0d0, A(k+1:n,k) /), A(k:n,j) )
          A(k:n,j) = A(k:n,j) - tmp * (/ 1.0d0, A(k+1:n,k) /)

          cnrm(j) = cnrm(j) - A(k,j) ** 2
       end do

       kmax = k + maxloc( cnrm(k+1:m+1), dim=1 )
    end do

    rank = k

    do k = 1,m
       imax = min( k, rank )
       if ( c(k) .gt. tol ) then
          A(1:imax,k) = A(1:imax,k) * c(k)
       end if
    end do

    do k = 1,m
       if ( .not. any( p(1:m) .eq. k ) ) then
          ! write(*,*) 'ERROR in subroutine qr: (A) Permutation array p is wrong. There is no k = ',k
          ! write(*,*) 'm = ',m,' p = ',p(1:m)
          call rexit('ERROR in subroutine qr: (A) Permutation array p is wrong.')
       end if
    end do

  end subroutine bmqr

  ! ******************************************************************
  ! ******************************************************************

  subroutine extractQ(n,m,mlim,lda,A,gam,rank,ldq,Q) bind(C, name="extractQ_")

    ! SCALAR ARGUMENTS
    integer(c_int), intent(in) :: lda,ldq,m,mlim,n,rank

    ! ARRAY ARGUMENTS
    real(c_double), intent(in) :: A(lda,m),gam(m)
    real(c_double), intent(out) :: Q(ldq,mlim)

    ! First rank columns of Q

    ! LOCAL SCALARS
    integer(c_int) :: i,j,k
    real(c_double) :: tmp

    do j = 1,rank
       Q(1:n,j) = (/ (0.0d0,i=1,j-1), 1.0d0, (0.0d0,i=j+1,n) /)
       do k = rank,1,-1
          tmp = gam(k) * dot_product( (/ 1.0d0, A(k+1:n,k) /), Q(k:n,j) )
          Q(k:n,j) = Q(k:n,j) - tmp * (/ 1.0d0, A(k+1:n,k) /)
       end do
    end do

  end subroutine extractQ

  ! ******************************************************************
  ! ******************************************************************

  subroutine extractR(m,mlim,lda,A,rank,ldr,R) bind(C, name="extractR_")

    ! SCALAR ARGUMENTS
    integer(c_int), intent(in) :: lda,ldr,m,mlim,rank

    ! ARRAY ARGUMENTS
    real(c_double), intent(in) :: A(lda,m)
    real(c_double), intent(out) :: R(ldr,mlim)

    ! LOCAL SCALARS
    integer(c_int) :: j

    R(1:rank,1:m) = A(1:rank,1:m)

    do j = 1,rank
       R(j+1:rank,j) = 0.0d0
    end do

  end subroutine extractr

  ! ******************************************************************
  ! ******************************************************************

  subroutine bmqrupdate(n,m,rank,mlim,jtbr,newcolp,ldq,Q,ldr,R,p,newcol) bind(C, name="bmqrupdate_")

    ! SCALAR ARGUMENTS
    integer(c_int), intent(in) :: jtbr,ldq,ldr,mlim,n,newcolp
    integer(c_int), intent(inout) :: m,rank

    ! ARRAY ARGUMENTS
    integer(c_int), intent(inout) :: p(mlim)
    real(c_double), intent(in) :: newcol(n)
    real(c_double), intent(inout) :: Q(ldq,mlim),R(ldr,mlim)

    ! THE DESCRIPTION BELOW NEEDS TO BE UPDATED SINCE IT DOES NOT
    ! CORRESPOND TO THE CURRENT VERSION OF THE ROUTINE

    ! This routine starts with the decomposition A P = Q R.
    ! A is n times m with rank rank.
    ! P is an m times m permutation matrix.
    ! Q is an n times rank matrix with orthonormal columns.
    ! R is a rank times m matrix of the form ( R11 R12 ) with
    ! R11 rank times rank, upper-triangular, and non-singular.
    !
    ! The routine receives a new column c in R^n that substitues the
    ! first column a1 of A and updates the given QR decomposition. As a
    ! result, the rank may remain the same or to be increased or reduced
    ! by one.
    !
    ! Let A P = ( a_i1 a_i2 ... a_im ); and let k such that a_ik is the
    ! first column of A. Let P1 be the permutation matrix such that A P
    ! P1 = ( a_i1 ... a_ik-1 a_ik+1 ... a_im a_ik ); i.e. P1 moves the
    ! first column of A to the last place in A P P1 and makes a shift to
    ! close the gap.
    !
    ! Then A ( P P1 ) = Q ( R P1 ).
    !
    ! Let S in m times (m-1) be the identity matrix with the last column
    ! removed. Then, we remove the last column of ( A P P1 ) and of ( R
    ! P1 ) by applying S from the right:
    !
    ! A ( P P1 ) S = Q ( R P1 ) S.
    !
    ! Now, there are two cases to be considered: (a) the first rank
    ! columns of R P1 S continued being equal to R11, and (b) applying
    ! P1 S removed a column from R11 and, by making the shift to close
    ! the gap, let the first column of R12 as the last colum of R11.
    !
    ! Let Q1^t be the matrix that (using Givens rotations) nullifies the
    ! subdiagonal elements of R11 P1 S, in case there is any. (In case
    ! (a), we have Q1 = Identity.)
    !
    ! Then Q1^t Q^t A P P1 S = Q1^t R P1 S or, equivalently,
    ! A ( P P1 ) S = ( Q Q1 ) ( Q1^t R P1 S ).
    !
    ! In case (b), we must find 'the column of R12' (currently allocated
    ! from column rank to column m-1 in Q1^t R P1 S) with largest
    ! rank-th element in absolute value. Let j be this column. Again,
    ! there are two cases to be considered. In case (i), the absolute
    ! value of the rank-th element is considered negligible, no
    ! permutation is done and the rank is reduced by one. (In this case,
    ! the last column of ( Q Q1 ) can be deleted since the last row of R
    ! is null.) In case (ii), columns rank and j must be permuted and
    ! the rank is preserved.
    !
    ! If we are in case (i), let T in rank times ( rank - 1 ) be the
    ! Identity matrix with the last column removed. If we are in case
    ! (ii), let T be rank times rank Identity matrix. Then
    !
    ! A ( P P1 S ) = ( Q Q1 T ) ( T^t Q1^t R P1 S ).
    !
    ! If we are in case (i), P2 is an Identity matrix. If we are in case
    ! (ii), it is the permutation matrix that permutes columns j and
    ! rank in ( T^t Q1^t R P1 S ). Then
    !
    ! A ( P P1 S P2 ) = ( Q Q1 T ) ( T^t Q1^t R P1 S P2 ).
    !
    ! We now compute
    !
    ! q = 0 if ||qtilde|| = 0 and q = qtilde / ||qtilde||, otherwise,
    !
    ! where
    !
    ! qtilde = c - Qtilde Qtilde^t c
    !
    ! and
    !
    ! Qtilde = Q Q1 T.
    !
    ! And there are two case, depending on whether ||qtilde|| = 0 or
    ! ||qtilde|| > 0.
    !
    ! If ||qtilde||=0, we add c as a new last column in A ( P P1 S P2 )
    ! and ( Qtilde^t c ) as a new last column in ( T^t Q1^t R P1 S P2 ).
    !
    ! If ||qtilde||>0, we add c as a new last column in A ( P P1 S P2 ),
    ! we add q as a new last column in ( Q Q1 T ), and we add ( Qtilde^t
    ! c ) concatenated with ||qtilde|| as a new last column in ( T^t
    ! Q1^t R P1 S P2 ). Finally we apply a permutation matrix P3 that
    ! permutes this last new column with column rank+1 and increase the
    ! rank by one.

    ! PARAMETERS
    real(c_double) :: tolsqr = 10.0d0 * sqrt( epsilon( 1.0d0) )

    ! LOCAL SCALARS
    integer(c_int) :: i,itmp,j,jmax,jtbrpos,k
    real(c_double) :: cc,ss,tau,tmp

    ! LOCAL ARRAYS
    real(c_double) :: v(mlim)

    if ( n .lt. m .or. rank .gt. m .or. m .gt. mlim .or. ldq .lt. n ) then ! .or. ldr .lt. mlim ) then
       ! write(*,*) 'ERROR in subroutine qrupdate: Check n, m, rank, mlim, ldq, and ldr!', &
       !     n, m, rank, mlim, ldq, ldr, n .lt. m, rank .gt. m, m .gt. mlim, ldq .lt. n, ldr .lt. mlim
       call rexit('ERROR in subroutine qrupdate: Check n, m, rank, mlim, ldq, and ldr!')
    end if

    if ( .not. ( 0 .le. jtbr .and. jtbr .le. m ) ) then
       ! write(*,*) 'ERROR in subroutine qrupdate: jtbr must be between 0 and m! m = ',m,' jtbr = ',jtbr
       call rexit('ERROR in subroutine qrupdate: jtbr must be between 0 and m!')
    end if

    if ( jtbr .eq. 0 .and. m .eq. mlim ) then
       ! write(*,*) 'ERROR in subroutine qrupdate: There is no space for adding a new column without removing another column!'
       call rexit('ERROR in subroutine qrupdate: There is no space for adding a new column without removing another column!')
    end if

    if ( any( Q(1:n,1:rank) /= Q(1:n,1:rank) ) .or. any( R(1:rank,1:m) /= R(1:rank,1:m) ) ) then
       ! write(*,*) 'ERROR in subroutine qrupdate: NaN in Q or R!'
       call rexit('ERROR in subroutine qrupdate: NaN in Q or R!')
    end if

    if ( newcolp .ne. 0 ) then
       if ( any( newcol(1:n) /= newcol(1:n) ) ) then
          ! write(*,*) 'ERROR in subroutine qrupdate: NaN in newcol!'
          call rexit('ERROR in subroutine qrupdate: NaN in newcol!')
       end if
    end if

    do k = 1,m
       if ( .not. any( p(1:m) .eq. k ) ) then
          ! write(*,*) 'ERROR in subroutine qrupdate: (B) Permutation array p is wrong. There is no k = ',k
          ! write(*,*) 'm = ',m,' p = ',p(1:mlim)
          call rexit('ERROR in subroutine qrupdate: (B) Permutation array p is wrong.')
       end if
    end do

    if ( jtbr .ne. 0 ) then
       ! Compute P P1 (saved in perm)

       jtbrpos = 1; do while ( p(jtbrpos) .ne. jtbr ); jtbrpos = jtbrpos + 1; end do

       p(jtbrpos:m-1) = p(jtbrpos+1:m)

       where ( p(1:m-1) .ge. jtbr )
          p(1:m-1) = p(1:m-1) - 1
       end where

       ! Compute R P1 (saved in R). Last column is ignored because it will
       ! be deleted.

       R(1:rank,jtbrpos:m-1) = R(1:rank,jtbrpos+1:m)

       ! If the removed column were part of R11, compute Q1 (with Given
       ! rotations) to nullify the subdiagonal elements in ( R P1 ); and
       ! apply it from left to ( R P1 ) and from right to Q (ignoring the
       ! last column of ( R P1 )).

       if ( jtbrpos .le. rank ) then

          ! Givens rotations to nullify the subdiagonal element of ( R P1 )
          ! that are located from column jtbrpos to column rank - 1.
          !
          ! I f x not equal zero then
          !
          ! |  cc ss | | x1 | = | y |,
          ! | -ss cc | | x2 |   | 0 |
          !
          ! where
          !
          ! y = ||x||
          ! cc = x1 / ||x||
          ! ss = x2 / ||x||

          ! Compute ( Q1^t R P1 ) and ( Q Q1 )

          do k = jtbrpos,rank - 1
             tmp = norm2( R(k:k+1,k) )

             if ( tmp .ne. 0.0d0 ) then
                cc = R(k,k)   / tmp
                ss = R(k+1,k) / tmp

                R(k,k) = tmp
                R(k+1,k) = 0.0d0

                do j = k + 1,m - 1
                   tmp = R(k,j)
                   R(k,j)   =   cc * tmp + ss * R(k+1,j)
                   R(k+1,j) = - ss * tmp + cc * R(k+1,j)
                end do

                do i = 1,n
                   tmp = Q(i,k)
                   Q(i,k)   =   cc * tmp + ss * Q(i,k+1)
                   Q(i,k+1) = - ss * tmp + cc * Q(i,k+1)
                end do
             end if
          end do

          ! Compute ( P P1 P2 ) (saved in perm)

          ! Check whether the rank of ( Q1^t R P1 ) is rank or rank -
          ! 1. (If the rank were m, there is nothing to be checked since by
          ! removing a column for sure we have that the rank was reduced by
          ! one.)

          if ( rank .eq. m ) then
             rank = rank - 1

          else
             jmax = rank - 1 + maxloc( abs( R(rank,rank:m-1) ), dim=1 )

             if ( abs( R(rank,jmax) ) .ge. tolsqr * max( norm2( R(1:rank,jmax) ), 1.0d0 )  ) then
                itmp = p(rank)
                p(rank) = p(jmax)
                p(jmax) = itmp

                v(1:rank) = R(1:rank,rank)
                R(1:rank,rank) = R(1:rank,jmax)
                R(1:rank,jmax) = v(1:rank)

             else
                rank = rank - 1
             end if
          end if
       end if

       m = m - 1

       if ( any( Q(1:n,1:rank) /= Q(1:n,1:rank) ) .or. any( R(1:rank,1:m) /= R(1:rank,1:m) ) ) then
          ! write(*,*) 'ERROR in subroutine qrupdate: NaN was created! (place 1)'
          call rexit('ERROR in subroutine qrupdate: NaN was created! (place 1)')
       end if
    end if

    if ( newcolp .eq. 1 ) then
       ! Compute ( Q Q1 )^t c (saved in colum m+1 of ( Q1^t R P1 )) and
       ! q = c - (Q Q1) ( Q Q1 )^t c

       R(1:rank,m+1) = matmul( transpose( Q(1:n,1:rank) ), newcol(1:n) )
       p(m+1) = m + 1

       Q(1:n,rank+1) = newcol(1:n) - matmul( Q(1:n,1:rank), R(1:rank,m+1) )

       tau = norm2( Q(1:n,rank+1) )

       if ( tau .ge. tolsqr * max( norm2( (/ R(1:rank,m+1), tau /) ), maxval( abs( (/ (R(k,k),k=1,rank) /) ) ), 1.0d0 ) ) then
          Q(1:n,rank+1) = Q(1:n,rank+1) / tau

          v(1:rank) = R(1:rank,rank+1)
          R(1:rank,rank+1) = R(1:rank,m+1)
          R(1:rank,m+1) = v(1:rank)

          R(rank+1,1:m+1) = 0.0d0
          R(rank+1,rank+1) = tau

          itmp = p(m+1)
          p(m+1) = p(rank+1)
          p(rank+1) = itmp

          rank = rank + 1
       end if

       m = m + 1
    end if

    do k = 1,m
       if ( .not. any( p(1:m) .eq. k ) ) then
          ! write(*,*) 'ERROR in subroutine qrupdate: (D) Permutation array p is wrong. There is no k = ',k
          ! write(*,*) 'm = ',m,' p = ',p(1:m)
          call rexit('ERROR in subroutine qrupdate: (D) Permutation array p is wrong.')
       end if
    end do

    if ( rank .gt. m ) then
       ! write(*,*) 'ERROR in subroutine qrupdate: There is something wrong because we ended with rank > m!'
       ! write(*,*) 'rank = ',rank,' m = ',m
       call rexit('ERROR in subroutine qrupdate: There is something wrong because we ended with rank > m!')
    end if

    if ( any( Q(1:n,1:rank) /= Q(1:n,1:rank) ) .or. any( R(1:rank,1:m) /= R(1:rank,1:m) ) ) then
       ! write(*,*) 'ERROR in subroutine qrupdate: NaN was created! (place 2)'
       call rexit('ERROR in subroutine qrupdate: NaN was created! (place 2)')
    end if

  end subroutine bmqrupdate

  ! ******************************************************************
  ! ******************************************************************

  subroutine bmqrsolve(n,m,rank,ldq,Q,ldr,R,p,b,x) bind(C, name="bmqrsolve_")

    ! SCALAR ARGUMENTS
    integer(c_int), intent(in) :: ldq,ldr,m,n,rank

    ! ARRAY ARGUMENTS
    integer(c_int), intent(in) :: p(m)
    real(c_double), intent(in) :: Q(ldq,rank),R(ldr,m),b(n)
    real(c_double), intent(out) :: x(m)

    ! LOCAL SCALARS
    integer(c_int) :: j,k
    real(c_double) :: tau,tmp

    ! LOCAL ARRAYS
    real(c_double) :: kap(rank),Rt(m,rank)

    if ( any( b(1:n) /= b(1:n) ) .or. any( Q(1:n,1:rank) /= Q(1:n,1:rank) ) .or. any( R(1:rank,1:m) /= R(1:rank,1:m) ) ) then
       ! write(*,*) 'ERROR in subroutine qrsolve: NaN in Q, R, or b!'
       call rexit('ERROR in subroutine qrsolve: NaN in Q, R, or b!')
    end if

    if ( any( (/  ( R(k,k), k=1,rank ) /) .eq. 0.0d0 ) ) then
       ! write(*,*) 'ERROR in subroutine qrsolve: R with null diagonal!'
       call rexit('ERROR in subroutine qrsolve: R with null diagonal!')
    end if

    ! Q^t b
    x(1:rank) = matmul( transpose( Q(1:n,1:rank) ), b(1:n) )

    if ( rank .eq. m ) then

       ! Solve R x = Q^t b by back substitution
       do k = rank,1,-1
          if ( R(k,k) .eq. 0.0d0 ) then
             ! write(*,*) 'ERROR in subroutine qrsolve: Division by zero (place 1)'
             call rexit('ERROR in subroutine qrsolve: Division by zero (place 1)')
          end if

          x(k) = x(k) / R(k,k)
          x(1:k-1) = x(1:k-1) - x(k) * R(1:k-1,k)
       end do

    else
       ! R^t = U T (R^t has full column rank, UT is a QR decomposition of R^t)

       Rt(1:m,1:rank) = transpose( R(1:rank,1:m) )

       do k = 1,rank
          tau = norm2( Rt(k:m,k) )

          if ( Rt(k,k) .lt. 0.0d0 ) then
             tau = - tau
          end if

          Rt(k,k) = Rt(k,k) + tau

          if ( Rt(k,k) .eq. 0.0d0 ) then
             ! write(*,*) 'ERROR in subroutine qrsolve: Division by zero (place 2)'
             call rexit('ERROR in subroutine qrsolve: Division by zero (place 2)')
          end if

          Rt(k+1:m,k) = Rt(k+1:m,k) / Rt(k,k)

          if ( tau .eq. 0.0d0 ) then
             ! write(*,*) 'ERROR in subroutine qrsolve: Division by zero (place 3)',k,rank
             call rexit('ERROR in subroutine qrsolve: Division by zero (place 3)')
          end if

          kap(k) = Rt(k,k) / tau

          Rt(k,k) = - tau

          do j = k + 1,rank
             tmp = kap(k) * dot_product( (/ 1.0d0, Rt(k+1:m,k) /), Rt(k:m,j) )
             Rt(k:m,j) = Rt(k:m,j) - tmp * (/ 1.0d0, Rt(k+1:m,k) /)
          end do
       end do

       ! R x = (U T)^t x = Q^t b

       ! x = T^-t Q^t b
       do k = 1,rank
          if ( Rt(k,k) .eq. 0.0d0 ) then
             ! write(*,*) 'ERROR in subroutine qrsolve: Division by zero (place 4)'
             call rexit('ERROR in subroutine qrsolve: Division by zero (place 4)')
          end if

          x(k) = ( x(k) - dot_product( Rt(1:k-1,k), x(1:k-1) ) ) / Rt(k,k)
       end do

       x(rank+1:m) = 0.0d0

       ! x = U T^-t Q^t b
       do k = rank,1,-1
          tmp = kap(k) * dot_product( (/ 1.0d0, Rt(k+1:m,k) /), x(k:m) )
          x(k:m) = x(k:m) - tmp * (/ 1.0d0, Rt(k+1:m,k) /)
       end do
    end if

    ! Permutation
    x(p(1:m)) = x(1:m)

  end subroutine bmqrsolve

end module bmqrf
