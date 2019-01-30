MODULE ibpm

  USE parameters
  USE dts
  USE grid
  USE mmisc
  USE controls
  USE variables  ! use of ibd_init in SUBROUTINE ibpm_pcm_solve_force_unsteady
                 ! use of chol_data in SUBROUTINE ibpm_pcm_solve_force_steady

  IMPLICIT NONE

CONTAINS

! *****************************************************************************************
! Steps the forward nonlinear simulations forward by one time step

    SUBROUTINE ibpm_solver( dtype_data )
      USE myfft
      USE user
      IMPLICIT NONE
      TYPE(data_t), INTENT(inout) :: dtype_data

      ! immersed boundary projection method sovler
      !WRITE(*,*) 'Solving intermediate vorticity at it =', dtype_data_1%it ! for debugging
      ! Step 1: solve intermediate curl(momentum) eq. on each grid
      CALL ibpm_compute_intermed_vort( dtype_data = dtype_data )

      IF (alloc_ibd) THEN
        ! we now have vorticity evolution on all domains.
        ! For smallest domain, we add the vorticity created by the body
        ! Step 2: solve the lagrangian force that satisfies no slip bc.
        !WRITE(*,*) 'Solving immersed boundary forces at it =', dtype_data_1%it ,', icpu =', icpu   ! for debugging
        CALL ibpm_solve_force( dtype_data = dtype_data )

        ! Step 3: correct the vorticity field with the lagrangian force on level 1
        !WRITE(*,*) 'Correcting vorticity field with immersed boundary forces ....'   ! for debugging
        dtype_data%levels(1)%omega%values = &
                       dtype_data%levels(1)%omega%values &
                       - ainv( curlT( reg( dtype_data%ibd%fb ) ) )

      END IF

      !WRITE(*,*) 'Computing the correct velocity field ....'   ! for debugging
      ! coarsify final omega and correct velocities on all grids
      CALL vort2flux( dtype_data = dtype_data, nlev = mgridlev )

      ! for steady motion (n_constraint=0), q0 is update here
      ! If no constraints, set q0 based on user-defined motion
      CALL user_center_of_mass_motion_grid( it = dtype_data%it, &
                                            dtype_data = dtype_data )
      ! assign q0 based on cm motion
      CALL mmisc_motion_cm2q0( dtype_data = dtype_data )
      
    END SUBROUTINE ibpm_solver

! *****************************************************************************************
  SUBROUTINE ibpm_compute_intermed_vort(dtype_data)
    ! Step 1: solve intermediate curl(momentum) eq. on each grid
    USE myfft
    IMPLICIT NONE
    TYPE(data_t), INTENT(inout) :: dtype_data
    INTEGER :: i, k
    REAL(KIND(0.D0)) :: rhs(2:m,2:n), rhsbc(2:m,2:n)
    REAL(KIND(0.D0)) :: rhs_old(2:m,2:n,mgridlev)
    REAL(KIND(0.D0)) :: omega_bc(2*(m+1)+2*(n+1))
    REAL(KIND(0.D0)) :: lastbc(2*(m+1)+2*(n+1),mgridlev), qlastbc(2*(m+1)+2*(n+1))
    REAL(KIND(0.D0)) :: nl_temp(Nq)

    ! save current vorticity bc as lastbc for the next evolution
    DO k = 1, mgridlev-1
       CALL mmisc_bc_getbc_spline( &
          r   = dtype_data%levels(k+1)%omega%values, &
          rbc = lastbc(:,k), &
          fac = 0.25d0 )
    END DO
    lastbc(:,mgridlev) = 0.d0

    !WRITE(*,*) 'get lastbc at it =', dtype_data_1%it ,', icpu =', icpu ! for debugging
    ! telescoping in from the coarsest domain
    DO k = mgridlev,1,-1
      ! *******************  telescoping in (begins) *********************
       ! get bc for intermediate vorticity during telescoping in
       IF (k.eq.mgridlev) THEN
          omega_bc = 0.d0
          qlastbc = 0.d0
       ELSE
          CALL mmisc_bc_getbc_spline( &
                      r   = dtype_data%levels(k+1)%omega%values, &
                      rbc = omega_bc, &
                      fac = 0.25d0 )
          CALL mmisc_bc_getbc_flux(&
                      r   = dtype_data%levels(k+1)%q%values + &
                            dtype_data%levels(k+1)%q0%values, &
                      rbc = qlastbc, &
                      fac = 0.5d0 )
       END IF

       ! Computing the nonliner term (nonlinear term mult. by dt/2)
         ! OLD: nl_temp = nonlinear( omega(:,:,k), q(:,k), q0(:,k), lastbc(:,k) )
         nl_temp = nonlinear_qcrossomega( &
                                q_conv = dtype_data%levels(k)%q%values + &
                                         dtype_data%levels(k)%q0%values, &
                                qbc    = qlastbc, &
                                omega  = dtype_data%levels(k)%omega%values, &
                                wbc    = lastbc(:,k) )

       ! add user-defined RHS forcing term to momentum eq.
       IF (k.eq.1) THEN
         ! specified your rhs_forcing
           DO i = 0, n_control ! i = 0 is for user-defined body force
             nl_temp = nl_temp + rhs_forcing( dtype_data = dtype_data, &
                                              i_ctrl = i )
           END DO
       END IF
       !rhs = rot( nonlinear( q_conv, omega, lastbc(:,k) ) )
       rhs = curlT( nl_temp )

       ! If this is the very first time step, we need to use explicit Euler
       IF ( dtype_data%it .eq. 1 ) THEN
          rhs_old(:,:,k) = rhs(:,:)
       ELSE
          rhs_old(:,:,k) = curlT( dtype_data%levels(k)%nl_bf_old%values )
       END IF

       ! update for next time step
       dtype_data%levels(k)%nl_bf_old%values = nl_temp

       ! apply bc to right-hand-side, we apply twice because of CN scheme
       rhsbc = 0.d0
       CALL mmisc_bc_applybc( r = rhsbc, rbc = lastbc(:,k), fac = vfac(k) )
       CALL mmisc_bc_applybc( r = rhsbc, rbc = omega_bc,    fac = vfac(k) )

       ! compute intermediate vorticity on the kth level
       dtype_data%levels(k)%omega%values = &
           dst( lam1i(:,:,k) * &
                ( dst( con1(k)*rhs(:,:)       + &
                       con2(k)*rhs_old(:,:,k) + &
                               rhsbc(:,:)          ) + &
                  lam1(:,:,k)*dst(dtype_data%levels(k)%omega%values)  )  )

    ! *******************  telescoping in (ends) *********************
    END DO

  END SUBROUTINE ibpm_compute_intermed_vort

! *****************************************************************************************
  SUBROUTINE ibpm_solve_force(dtype_data)
    ! Step 2: solve the lagrangian force that satisfies no slip bc.
    IMPLICIT NONE
    TYPE(data_t), INTENT(inout) :: dtype_data
    TYPE(ib_data_t) :: ibd_temp
    REAL(KIND(0.D0)), DIMENSION(Nf) :: rhsf, motion

    ! get intermediate velocity on the first grid level to compute the right-hand-side
    CALL vort2flux( dtype_data = dtype_data, nlev = 1 )

    ! apply actuators
    CALL dts_alloc_ib_data_t( dtype = ibd_temp, nf = nf )
    ibd_temp = dtype_data%ibd
    CALL ibpm_actuators( it = dtype_data%it, &
                         dtype_data_old = dtype_data, &
                         dtype_ibd = ibd_temp )

    dtype_data%ibd = ibd_temp
    CALL dts_dealloc_ib_data_t( dtype = ibd_temp )

    ! compute force
    dtype_data%ibd%slip_vel = delta * dtype_data%ibd%vb
    !WRITE(*,*) 'slip velocity =',  dtype_data%ibd%slip_vel/delta
    motion = dtype_data%ibd%slip_vel - regT( dtype_data%levels(1)%q0%values )

    rhsf = regT( dtype_data%levels(1)%q%values ) - motion

    dtype_data%ibd%fb = cholsl(rhsf,chol_data%cholmat,chol_data%cholvec)

  END SUBROUTINE ibpm_solve_force

  ! *****************************************************************************************
  ! Computes the velocity of the actuators (which may potentially be slaved to a moving body)
    SUBROUTINE ibpm_actuators(it, dtype_data_old, dtype_ibd)
      USE user
      IMPLICIT NONE

      TYPE(ib_data_t), INTENT(inout) :: dtype_ibd
      TYPE(data_t), INTENT(in) :: dtype_data_old

      INTEGER :: iact,istrt,iend
      INTEGER :: i, i_ctrl,islv
      INTEGER, INTENT(in) :: it
      REAL(KIND(0.D0)), DIMENSION(:), ALLOCATABLE :: vx,vy
      REAL(KIND(0.D0)) :: xmove, ymove, theta_g2l

      DO iact=1,n_actuator

         istrt = act(iact)%ipos
         iend = act(iact)%ipos+act(iact)%npts - 1

         ALLOCATE( vx(iend-istrt+1), vy(iend-istrt+1) )
         theta_g2l = dtype_data_old%cm_motion%theta_g2l
         dtype_ibd%vb(f(istrt:iend,1)) = 0.D0
         dtype_ibd%vb(f(istrt:iend,2)) = 0.D0
         DO i_ctrl = 1,n_control
           CALL user_actuator_vel( actno = iact, &
                                   it = it, &
                                   theta_g2l = theta_g2l, &
                                   x = dtype_ibd%xb(f(istrt:iend,1)), &
                                   y = dtype_ibd%xb(f(istrt:iend,2)), &
                                   vx = vx, &
                                   vy = vy, &
                                   i_ctrl = i_ctrl )
           dtype_ibd%vb(f(istrt:iend,1)) = dtype_ibd%vb(f(istrt:iend,1)) + vx
           dtype_ibd%vb(f(istrt:iend,2)) = dtype_ibd%vb(f(istrt:iend,2)) + vy
         END DO
         DEALLOCATE( vx,vy )
      END DO

    END SUBROUTINE ibpm_actuators

    ! *****************************************************************************************
    !DO Not Change fucntion rhs_forcing
      FUNCTION rhs_forcing( dtype_data, i_ctrl ) RESULT(dq)
        USE user
       ! Computes the added force term to the RHS of the Navier-Stokes equation, according to
       ! user specified controls
        IMPLICIT NONE
        INTEGER, INTENT(in) :: i_ctrl
        TYPE(data_t), INTENT(in) :: dtype_data
        INTEGER :: it, i, j, ilev
        REAL(KIND(0.D0)) :: q(Nq), dq(Nq)
        REAL(KIND(0.D0)), ALLOCATABLE :: dqx(:,:), dqy(:,:)
        REAL(KIND(0.D0)) :: xx,yy, uu,vv, theta_g2l, del, bf(2)

        dq = 0.d0

        ilev = 1
        q = dtype_data%levels(ilev)%q%values &
            + dtype_data%levels(ilev)%q0%values

        it = dtype_data%it
        theta_g2l = dtype_data%cm_motion%theta_g2l

        ALLOCATE( dqx(1:m,1:n), dqy(1:m,1:n) )
        dqx = 0.d0
        dqy = 0.d0
        ! initial vorticity field
        del = delta*2.d0**(REAL(ilev)-1.D0)  ! cell face length on each grid
        bf = 0.D0
        DO j = 1, n
          DO i = 1, m
            ! find coord where x-comp velocity lives
             xx =  (REAL(i)-0.5D0-REAL(m)/2.D0) * del + &
                                        REAL(m)/2.D0*delta - offsetx
             yy =  (REAL(j)-0.5D0-REAL(n)/2.D0) * del + &
                                        REAL(n)/2.D0*delta - offsety
             uu = 0.5d0 * ( q(u(i,j)) + q(u(i+1,j)) ) / del
             vv = 0.5d0 * ( q(v(i,j)) + q(v(i,j+1)) ) / del

             bf = bodyforce(it,theta_g2l,xx,yy,uu,vv,i_ctrl)

             dqx(i,j) = dqx(i,j) + (del**3.d0) * bf(1)
             dqy(i,j) = dqx(i,j) + (del**3.d0) * bf(2)
          END DO
       END DO

       DO j = 2,n-1
         DO i = 2,m
            dq(u(i,j)) = 0.5d0 * ( dqx(i-1,j)+dqx(i,j) )
         END DO
      END DO
      DO i = 2,m-1
         DO j = 2,n
           dq(v(i,j)) = 0.5d0 * ( dqy(i,j-1)+dqy(i,j) )
         END DO
      END DO

      DEALLOCATE( dqx, dqy )

      END FUNCTION rhs_forcing

! *****************************************************************************************

END MODULE ibpm
