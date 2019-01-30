MODULE mmisc_opr

  ! This module contains subroutines and functions of operators

  USE parameters
  USE dts
  USE grid
  IMPLICIT NONE

CONTAINS

! *****************************************************************************************
  SUBROUTINE vort2flux( dtype_data, nlev )

    !***************************************************************!
    !*    Multiscale method to solve C^T C s = omega               *!
    !*    and return the velocity, C s.                            *!
    !*    Results are returned in vel on each of the first nlev    *!
    !*    grids.                                                   *!
    !*    Warning: the vorticity field on all but the finest mesh  *!
    !*     is modified by the routine in the following way:        *!
    !*     the value in the center of the domain is interpolated   *!
    !*     from the next finer mesh (the value near the edge is    *!
    !*     not changed.                                            *!
    !***************************************************************!
    USE myfft
    USE mmisc_bc
    IMPLICIT NONE
    INTEGER, INTENT(in) :: nlev
    TYPE(data_t), INTENT(inout) :: dtype_data
    TYPE(circulation_t) :: vort
    REAL(KIND(0.D0)) :: sbc(2*(m+1)+2*(n+1))  ! streamfun bc
    INTEGER :: k

    ! find coarsifieid omegas
    DO k = 2, mgridlev
        dtype_data%levels(k)%omega%values = &
                  coarsify( crhs = dtype_data%levels(k)%omega%values, &
                             rhs = dtype_data%levels(k-1)%omega%values )
    END DO

    CALL dts_alloc_circulation_t( dtype = vort )

    ! invert laplacian on largest grid and then telescope in
    DO k = mgridlev, 1, -1

      ! don't mess up with omega. Make a copy and use it as the right-hand side
      vort%values = dtype_data%levels(k)%omega%values

       IF (k.eq.mgridlev) THEN
         ! zero bc's on stfn assumed for largest grid
         sbc = 0.d0
       ELSE
         ! get bc for stream function
         CALL mmisc_bc_getbc_spline( r = dtype_data%levels(k+1)%stfn%values, &
                                     rbc = sbc, fac = 1.d0)       ! get bc's from pervious s
       END IF

       ! apply bc to right-hand-side
       CALL mmisc_bc_applybc( r = vort%values, rbc = sbc, fac = 1.d0)  ! apply them

       ! solve poisson eq
       dtype_data%levels(k)%stfn%values = ctci( vort%values )

       ! compute vel if k <= nlev
       IF (k.le.nlev) THEN
         dtype_data%levels(k)%q%values = &
                        curl( dtype_data%levels(k)%stfn%values, sbc ) !dtype_data%levels(k)%stfn%values, sbc )
       END IF

    END DO

    CALL dts_dealloc_circulation_t( dtype = vort )

  END SUBROUTINE vort2flux

  !********************************************************************!
    FUNCTION a_times(x)
      !***************************************************************!
      !*   Cholesky factorization of EC(C^t C)^-1 C^t E^t            *!
      !*      performed once only for stationary bodies              *!
      !***************************************************************!
      USE myfft
      IMPLICIT NONE
      INTEGER :: i
      REAL(KIND(0.D0)) :: x(Nf), a_times(Nf)
      TYPE(data_t) :: dtype_temp


      ! zero all levels of temporary vorticity
      CALL dts_alloc_data_t( dtype = dtype_temp, &
                             nb = nb, alloc_ibd =.FALSE.) ! nct is a dummy const. here
      DO i = 1, mgridlev
        dtype_temp%levels(i)%omega%values = 0.D0
      END DO
      ! regularize immersed body forces for level 1 vorticity
      dtype_temp%levels(1)%omega%values = ainv( curlT( reg(x) ) )
      ! coarsify vorticity and find curl(laplace^-1)
      CALL vort2flux( dtype_data = dtype_temp, nlev = 1 )
      ! regularize
      a_times = regT( dtype_temp%levels(1)%q%values )

      CALL dts_dealloc_data_t( dtype = dtype_temp, &
                               alloc_ibd =.FALSE. )

    END FUNCTION a_times

  ! *********************************************************************!
  FUNCTION coarsify( crhs, rhs ) RESULT( arhs )
     !***************************************************************!
     !*   given vorticity on a smaller, fine mesh, (rhs) interp.    *!
     !*   values to the center region of a larger, coarser mesh     *!
     !*   (crhs).  The values outside the center region are         *!
     !*   not unmodified. Result is placed in arhs                  *!
     !***************************************************************!
     IMPLICIT NONE

     REAL(KIND(0.D0)) :: crhs(:,:), rhs(:,:)
     REAL(KIND(0.D0)) :: arhs(SIZE(rhs,1),SIZE(rhs,2))
     INTEGER          :: i, j, indi, indj

     arhs = crhs
     DO j=-n/4+1,n/4-1
        indj = n/2+2*j
        DO i=-m/4+1,m/4-1
           indi = m/2+2*i
           arhs(m/2+i,n/2+j) = rhs(indi  ,indj)   + &
                       0.5d0*( rhs(indi+1,indj)   + rhs(indi  ,indj+1)   + &
                               rhs(indi-1,indj)   + rhs(indi  ,indj-1) ) + &
                      0.25d0*( rhs(indi+1,indj+1) + rhs(indi+1,indj-1)   + &
                               rhs(indi-1,indj-1) + rhs(indi-1,indj+1) )
         END DO
      END DO

 END FUNCTION coarsify

 ! *****************************************************************************************

 FUNCTION nonlinear_qcrossomega( q_conv, qbc, omega, wbc ) RESULT(qcgamma)
   !***************************************************************!
   !*   nonlinear terms in rotational form  : q_conv x omega      *!
   !***************************************************************!
   IMPLICIT NONE
   REAL(KIND(0.D0)) :: q_conv(Nq), omega(2:m,2:n)
   REAL(KIND(0.D0)) :: wbc(2*(m+1)+2*(n+1)), qbc(2*(m+1)+2*(n+1))
   REAL(KIND(0.D0)) :: qcgamma(Nq)
   INTEGER :: i,j

   qcgamma = 0.D0

    ! get ( q x omega )_x
   DO i=2,m
      qcgamma(u(i,1)) = 0.25D0*(  &
                          (q_conv(v(i,2))+q_conv(v(i-1,2))) * omega(i,2) &
                        + (q_conv(v(i,1))+q_conv(v(i-1,1))) * wbc(bottom+i) )
      DO j=2,n-1
        qcgamma(u(i,j)) = 0.25D0*(  &
                            (q_conv(v(i,j+1))+q_conv(v(i-1,j+1))) * omega(i,j+1) &
                          + (q_conv(v(i,j  ))+q_conv(v(i-1,j  ))) * omega(i,j  ) )
      END DO
      qcgamma(u(i,n)) = 0.25D0*(  &
                          (q_conv(v(i,n+1))+q_conv(v(i-1,n+1))) * wbc(top+i) &
                        + (q_conv(v(i,n))+q_conv(v(i-1,n))) * omega(i,n) )
   END DO
   ! for i = 1
   qcgamma(u(1,1)) = 0.25D0*(  &
                       (q_conv(v(1,2))+qbc(left+2)) * wbc(left+2) &
                     + (q_conv(v(1,1))+qbc(left+1)) * wbc(left+1) )
   DO j=2,n-1
     qcgamma(u(1,j)) = 0.25D0*(  &
                         (q_conv(v(1,j+1))+qbc(left+j+1)) * wbc(left+j+1) &
                       + (q_conv(v(1,j  ))+qbc(left+j  )) * wbc(left+j  ) )
   END DO
   qcgamma(u(1,n)) = 0.25D0*(  &
                       (q_conv(v(1,n+1))+qbc(left+n+1)) * wbc(left+n+1) &
                     + (q_conv(v(1,n  ))+qbc(left+n  )) * wbc(left+n  ) )
   ! for i = m+1
   qcgamma(u(m+1,1)) = 0.25D0*(  &
                         (qbc(right+2)+q_conv(v(m,2))) * wbc(right+2) &
                       + (qbc(right+1)+q_conv(v(m,1))) * wbc(right+1) )
   DO j=2,n-1
     qcgamma(u(m+1,j)) = 0.25D0*(  &
                           (qbc(right+j+1)+q_conv(v(m,j+1))) * wbc(right+j+1) &
                         + (qbc(right+j  )+q_conv(v(m,j  ))) * wbc(right+j  ) )
   END DO
   qcgamma(u(m+1,n)) = 0.25D0*(  &
                         (qbc(right+n+1)+q_conv(v(m,n+1))) * wbc(right+n+1) &
                       + (qbc(right+n  )+q_conv(v(m,n  ))) * wbc(right+n  ) )
    !note : i=1 or i=m+1 are not needed by rot but by div

    ! get ( q x omega )_y
    DO j=2,n
       qcgamma(v(1,j))    = -0.25D0*( &
                           (q_conv(u(2,j))+q_conv(u(2,j-1))) * omega(2,j) &
                         + (q_conv(u(1,j))+q_conv(u(1,j-1))) * wbc(left+j) )
       DO i=2,m-1
         qcgamma(v(i,j)) = -0.25D0*( &
                             (q_conv(u(i+1,j))+q_conv(u(i+1,j-1))) * omega(i+1,j) &
                           + (q_conv(u(i  ,j))+q_conv(u(i  ,j-1))) * omega(i  ,j) )
       END DO
       qcgamma(v(m,j))    = -0.25D0*(  &
                           (q_conv(u(m+1,j))+q_conv(u(m+1,j-1))) * wbc(right+j) &
                         + (q_conv(u(m  ,j))+q_conv(u(m  ,j-1))) * omega(m,j) )
    END DO
    ! j = 1
    qcgamma(v(1,1))    = -0.25D0*( &
                           (qbc(bottom+2)+q_conv(u(2,1))) * wbc(bottom+2) &
                         + (qbc(bottom+1)+q_conv(u(1,1))) * wbc(bottom+1) )
    DO i=2,m-1
      qcgamma(v(i,1)) = -0.25D0*( &
                          (qbc(bottom+i+1)+q_conv(u(i+1,1))) * wbc(bottom+i+1) &
                        + (qbc(bottom+i  )+q_conv(u(i  ,1))) * wbc(bottom+i  ) )
    END DO
    qcgamma(v(m,1))    = -0.25D0*(  &
                           (qbc(bottom+m+1)+q_conv(u(m+1,1)))*wbc(bottom+m+1) &
                         + (qbc(bottom+m  )+q_conv(u(m  ,1)))*wbc(bottom+m  ) )
    ! j = n+1
    qcgamma(v(1,n+1))    = -0.25D0*( &
                             (q_conv(u(2,n))+qbc(top+2)) * wbc(top+2) &
                           + (q_conv(u(1,n))+qbc(top+1)) * wbc(top+1) )
    DO i=2,m-1
      qcgamma(v(i,n+1)) = -0.25D0*( &
                            (q_conv(u(i+1,n))+qbc(top+i+1)) * wbc(top+i+1) &
                          + (q_conv(u(i  ,n))+qbc(top+i  )) * wbc(top+i  ) )
    END DO
    qcgamma(v(m,n+1))    = -0.25D0*(  &
                             (q_conv(u(m+1,n))+qbc(top+m+1))*wbc(top+m+1) &
                           + (q_conv(u(m  ,n))+qbc(top+m  ))*wbc(top+m  ) )
    !note : j=1 or j=n+1 are not needed by rot but by div

  END FUNCTION nonlinear_qcrossomega

! *****************************************************************************************
  FUNCTION curl( x, xbc )
    !***************************************************************!
    !*   returns curl(x) given x and the values of s'fun on bdy    *!
    !***************************************************************!
    IMPLICIT NONE

    REAL(KIND(0.D0)) :: x(2:m,2:n)
    REAL(KIND(0.D0)) :: xbc(2*(m+1)+2*(n+1))
    REAL(KIND(0.D0)) :: curl(Nq)
    INTEGER :: i,j

    DO j=2,n-1
       DO i=2,m
          curl(u(i,j)) = x(i,j+1) - x(i,j)
       ENDDO
    ENDDO

    DO i=2,m
       j=1
       curl(u(i,j)) = x(i,j+1) - xbc( bottom + i )
       j=n
       curl(u(i,j)) = xbc( top + i ) - x(i,j)
    ENDDO

    DO j=1,n
       i = 1
       curl(u(i,j)) = xbc( left + j + 1 ) -xbc( left + j )
       i = m+1
       curl(u(i,j)) = xbc( right + j + 1 ) -xbc( right + j )
    ENDDO

    DO j=2,n
       i = 1
       curl(v(i,j)) = - x(i+1,j) + xbc( left + j )
       DO i=2,m-1
          curl(v(i,j)) = x(i,j) - x(i+1,j)
       ENDDO
       i = m
       curl(v(i,j)) = - xbc( right + j) +  x(i,j)
    ENDDO

    DO i=1,m
       j = 1
       curl(v(i,j)) = xbc( bottom + i ) -xbc( bottom + i + 1 )
       j = n+1
       curl(v(i,j)) = xbc( top + i ) -xbc( top + i + 1 )
    ENDDO

  END FUNCTION curl

! *****************************************************************************************
 FUNCTION curlT( x )
   !***************************************************************!
   !*   Transpose of curl                                         *!
   !***************************************************************!
   IMPLICIT NONE

   REAL(KIND(0.D0)) :: x(Nq)
   REAL(KIND(0.D0)) :: curlT(2:m,2:n)
   INTEGER :: i,j

   DO j=2,n
      DO i=2,m
         curlT(i,j) = x(v(i,j)) - x(v(i-1,j)) - x(u(i,j)) + x(u(i,j-1))
      ENDDO
   ENDDO

 END FUNCTION curlT

 !*****************************************************************!
  FUNCTION reg( h0 ) RESULT( h )

    !***************************************************************!
    !*   regularization of body force                              *!
    !***************************************************************!
    IMPLICIT NONE
    REAL(KIND(0.D0)), DIMENSION(Nf), INTENT(IN) :: h0
    REAL(KIND(0.D0)), DIMENSION(Nq)             :: h
    INTEGER :: i, j, k, l, p, next

    ! initalize regularization field
    h = 0.D0

    DO k = 1, nb
       i = indexx(k)
       j = indexx(k+nb)
       next = 0
       DO l = -support, support
          DO p = -support, support
             next = next+1
             h(u(i+p,j+l)) = h(u(i+p,j+l)) &
                  + weight(next,k)*h0(k)
             h(v(i+p,j+l)) = h(v(i+p,j+l)) &
                  + weight(next,k+nb)*h0(k+nb)
          END DO
       END DO
    END DO

  END FUNCTION reg

! *****************************************************************************************
  FUNCTION regT( h ) RESULT( h0 )

    !***************************************************************!
    !*  interpolation to body point  (Transpose of reg)            *!
    !***************************************************************!
    IMPLICIT NONE
    REAL(KIND(0.D0)), DIMENSION(Nq), INTENT(IN) :: h
    REAL(KIND(0.D0)), DIMENSION(Nf)             :: h0
    INTEGER :: i, j, k, l, p, next

    ! initalize interpolation field
    h0 = 0.D0

    DO k = 1, nb
       i = indexx(k)
       j = indexx(k+nb)
       next = 0
       DO l = -support, support
          DO p = -support, support
             next = next+1
             h0(k) = h0(k) + weight(next,k)*h(u(i+p,j+l))
             h0(k+nb) = h0(k+nb) + weight(next,k+nb)*h(v(i+p,j+l))
          END DO
       END DO
    END DO

  END FUNCTION regT

 ! ***********************************************************************!
! ********************************************************************!
! subroutines and functions below are pressure calculation related    !
! ********************************************************************!
 ! **********************************************************************!

 FUNCTION coarsify_pressure( crhs, rhs ) RESULT( arhs )
    IMPLICIT NONE
    !***************************************************************!
    !*   given mass source on a smaller, fine mesh, (rhs) interp.  *!
    !*   values to the center region of a larger, coarser mesh     *!
    !*   (crhs).  The values outside the center region are         *!
    !*   not unmodified. Result is placed in arhs                  *!
    !***************************************************************!
    REAL(KIND(0.D0)) :: crhs(:,:), rhs(:,:)
    REAL(KIND(0.D0)) :: arhs(SIZE(rhs,1),SIZE(rhs,2))
    INTEGER :: i,j,indi,indj, ischeme
    REAL(KIND(0.D0)), ALLOCATABLE :: temp(:,:)

    ischeme = 2
    IF (ischeme.eq.1) THEN
      ! *************** scheme 1: naive averging (old) ***********************
      ! this coarsifying scheme is too naive, diverence wll not be perserve
      arhs = crhs
      DO j=-n/4+1, n/4
        indj = n/2+2*j
        DO i=-m/4+1, m/4
          indi = m/2+2*i
          arhs(m/2+i,n/2+j) = ( rhs(indi-1, indj-1) + &  ! Left Bottom
                                rhs(indi-1, indj  ) + &  ! Left Top
                                rhs(indi,   indj-1) + &  ! Right Bottom
                                rhs(indi,   indj  )   )  ! Right Top
        END DO
      END DO
    ! *************** scheme 1: naive averging (old) ***********************
    ELSEIF(ischeme.eq.2) THEN
      ! ************** scheme 2: Divergence perserving (new) ***************
      ! The new coarsifying scheme will perserve diverence
      ALLOCATE( temp(1:m/2-1,1:n/2-1) )
      DO j = 1, n/2-1
        indj = 2*j
        DO i = 1, m/2-1
          indi = 2*i
          temp(i,j) =          rhs(indi  ,indj)   + &
                      0.5d0*(  rhs(indi+1,indj)   + rhs(indi  ,indj+1)   + &
                               rhs(indi-1,indj)   + rhs(indi  ,indj-1) ) + &
                      0.25d0*( rhs(indi+1,indj+1) + rhs(indi+1,indj-1)   + &
                               rhs(indi-1,indj-1) + rhs(indi-1,indj+1) )
        END DO
      END DO

      arhs = crhs
      DO j=-n/4+2,n/4-1
        indj = n/4+j-1
        DO i=-m/4+2,m/4-1
          indi = m/4+i-1
          arhs(m/2+i,n/2+j) = (      temp(indi,indj)   + 3.D0*temp(indi+1,indj) + &
                                3.D0*temp(indi,indj+1) + 9.D0*temp(indi+1,indj+1) )/16.D0
        END DO
      END DO

      DEALLOCATE( temp )
      ! ************** scheme 2: Divergence perserving (new) ***************
    END IF

END FUNCTION coarsify_pressure

 ! **********************************************************************!

FUNCTION divergence( x )
  !***************************************************************!
  !*   divergence(x), x at edge & divergence at center           *!
  !***************************************************************!
   IMPLICIT NONE
   REAL(KIND(0.D0)) :: x(Nq)
   REAL(KIND(0.D0)) :: divergence(1:m,1:n)
   INTEGER          :: i,j

   DO j = 1, n
      DO i = 1, m
         divergence(i,j) = x(u(i+1,j)) - x(u(i,j)) + x(v(i,j+1)) - x(v(i,j))
      END DO
   END DO

END FUNCTION divergence

 ! **********************************************************************!
! ********************************************************************!
! subroutines and functions above are pressure calculation related    !
! ********************************************************************!
 ! **********************************************************************!

END MODULE mmisc_opr
