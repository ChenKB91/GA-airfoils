MODULE myfft

   IMPLICIT NONE

   INCLUDE 'fftw3.f'           ! needed for defining the plan

   REAL(KIND(0.D0)), ALLOCATABLE :: viscfac(:),vfac(:),con1(:),con2(:)
   REAL(KIND(0.D0)), ALLOCATABLE :: lam(:,:), in(:,:), laminv(:,:), in_ddti(:,:), laminv_ddti(:,:)
   REAL(KIND(0.D0)), ALLOCATABLE :: lam1(:,:,:), lam1i(:,:,:), lam1inv(:,:,:)

!   REAL(KIND(0.D0)), DIMENSION(:,:,:), ALLOCATABLE :: intfac1,intfac2,intfac3
   INTEGER*8 :: forward, inverse, forward_ddti
   INTEGER :: mm, nn
   REAL*8  :: normalize

CONTAINS

! *****************************************************************************************
  SUBROUTINE myfft_setup_fft

    USE grid
    USE parameters
    IMPLICIT NONE
    INTEGER  :: i,j,k
    REAL*8   :: del2, del22, normalize_ddti
    REAL*8   :: lam_ddti(m,n)

    mm = m
    nn = n
    ALLOCATE( viscfac(mgridlev), con1(mgridlev), con2(mgridlev) )
    ALLOCATE( vfac(mgridlev) )
    ALLOCATE( in(mm-1,nn-1), laminv(mm-1,nn-1) )
    ALLOCATE( lam(mm-1,nn-1) )
    ALLOCATE( lam1(mm-1,nn-1,mgridlev), &
              lam1i(mm-1,nn-1,mgridlev), &
              lam1inv(mm-1,nn-1,mgridlev) )
!    ALLOCATE( intfac1(mm-1,nn-1,mgridlev), intfac2(mm-1,nn-1,mgridlev), intfac3(mm-1,nn-1,mgridlev)  )
    ALLOCATE( in_ddti(mm,nn), laminv_ddti(mm,nn)  )

    CALL dfftw_plan_r2r_2d(forward, mm-1,nn-1,in,in, FFTW_RODFT00, FFTW_RODFT00, FFTW_EXHAUSTIVE)
    CALL dfftw_plan_r2r_2d(forward_ddti, mm,nn,in_ddti,in_ddti, FFTW_RODFT00, FFTW_RODFT00, FFTW_EXHAUSTIVE)

    normalize_ddti = 4.d0*REAL( (mm+1)*(nn+1) )
    normalize      = 4.d0*REAL(     mm*nn     )

    ! eigenvalues for inverse of C^T C
    del2 = delta*delta
    DO k=1,mgridlev
       del22      =  del2*4.d0**(REAL(k)-1)
       viscfac(k) =  dt/Re/del22
       vfac(k)    =  0.5d0*dt/Re/del22/normalize
       con1(k)    =  1.5d0*dt/del22/normalize
       con2(k)    = -0.5d0*dt/del22/normalize
    ENDDO

    DO j=1,nn-1
       DO i=1,mm-1

          lam(i,j) = -2.d0*( COS( pi*REAL(i)/REAL(mm) ) + &
                             COS( pi*REAL(j)/REAL(nn) ) - 2.d0 )
          laminv(i,j) = 1.d0/lam(i,j)/normalize
          ! integrating factors for viscous terms
          DO k=1,mgridlev
             del22 = del2* 4.d0**(REAL(k)-1.D0)
             lam1(i,j,k)    =      (1.d0 - 0.5d0*dt*lam(i,j)/del22/Re)/normalize ! original
             lam1i(i,j,k)   = 1.d0/(1.d0 + 0.5d0*dt*lam(i,j)/del22/Re) ! original

! it appears that this HUGE array is NEVER used anywhere, so why compute/store it????
!             lam1inv(i,j,k) = 1.d0/(1.d0 - 0.5d0*dt*lam(i,j)/del22/Re)

!computed, but never USED !!????
!             intfac2(i,j,k) =  1.5*dt* EXP(      -dt*lam(i,j) / del22 / Re )/del22/normalize
!             intfac3(i,j,k) = -0.5*dt* EXP( -2.d0*dt*lam(i,j) / del22 / Re )/del22/normalize
!             intfac1(i,j,k) =          EXP(      -dt*lam(i,j) / del22 / Re )      /normalize
          END DO
      END DO
    END DO
    DO j=1,nn
       DO i=1,mm
          lam_ddti(i,j)    = 2.d0*( COS( pi*REAL(i)/REAL(mm+1) ) + &
                                    COS( pi*REAL(j)/REAL(nn+1) ) - 2.d0 )
          laminv_ddti(i,j) = 1.d0/lam_ddti(i,j)/normalize_ddti
       END DO
    END DO

  END SUBROUTINE myfft_setup_fft
! *****************************************************************************************
  SUBROUTINE myfft_destroy_fft
    IMPLICIT NONE
    CALL dfftw_destroy_plan(forward)
    DEALLOCATE( lam,in, laminv, viscfac,vfac,con1,con2 )
!    deallocate( intfac1,intfac2,intfac3 )
    DEALLOCATE( lam1, lam1i, lam1inv )

    CALL dfftw_destroy_plan(forward_ddti)
    DEALLOCATE( in_ddti, laminv_ddti )

  END SUBROUTINE myfft_destroy_fft
! *****************************************************************************************
  FUNCTION dst( psi)
    IMPLICIT NONE
    REAL(KIND(0.D0)) :: psi(:,:)
    REAL(KIND(0.D0)) :: dst(2:mm,2:nn)

    ! discrete sine transform
    ! careful...two calls of dst need to be divided by "normalize"
    ! to return original vector

    in  = psi
    CALL dfftw_execute_r2r(forward,in,in)
    dst = in

  END FUNCTION dst
! *****************************************************************************************
  FUNCTION ctci( omega )
    IMPLICIT NONE
    REAL(KIND(0.D0)) :: omega(:,:)
    REAL(KIND(0.D0)) :: ctci(2:mm,2:nn)

    in =  omega
    CALL dfftw_execute_r2r(forward,in,in)
    in = laminv * in
    CALL dfftw_execute_r2r(forward,in,in)
    ctci =  in

  END FUNCTION ctci
! *****************************************************************************************
  FUNCTION ainv( omega )
    IMPLICIT NONE
    REAL(KIND(0.D0)) :: omega(:,:)
    REAL(KIND(0.D0)) :: ainv(2:mm,2:nn)

    in = omega
    CALL dfftw_execute_r2r(forward,in,in)
    in = lam1i(:,:,1) * in / normalize
    CALL dfftw_execute_r2r(forward,in,in)
    ainv = in

  END FUNCTION ainv
! *****************************************************************************************
  FUNCTION ddti( phi )
    IMPLICIT NONE
    REAL(KIND(0.D0)) :: phi(:,:)
    REAL(KIND(0.D0)) :: ddti(1:mm,1:nn)

    in_ddti = phi
    CALL dfftw_execute_r2r(forward_ddti,in_ddti,in_ddti)
    in_ddti = laminv_ddti * in_ddti
    CALL dfftw_execute_r2r(forward_ddti,in_ddti,in_ddti)
    ddti = in_ddti

  END FUNCTION ddti
! *****************************************************************************************

END MODULE myfft
