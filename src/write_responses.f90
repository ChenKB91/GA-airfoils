MODULE write_responses
  ! This moduel contains subroutines that compute the response you want
  USE parameters
  USE dts
  USE grid

  IMPLICIT NONE
  REAL :: now_lift, now_drag

CONTAINS
! *****************************************************************************************
SUBROUTINE write_responses_write_file( dtype_data )
  !***************************************************************!
  !*   write maximum u*dt/dx and slip error                      *!
  !***************************************************************!
  IMPLICIT NONE
  TYPE(data_t), INTENT(in) :: dtype_data
  TYPE(response_t) :: dtype_resp
  LOGICAL :: error
  INTEGER :: i, j, it

  ! allocate response metrix
  CALL dts_alloc_response_t( dtype = dtype_resp, &
                             nbdy = n_body+n_actuator, &
                             ncr = Ncr )
  ! compute responses (for all fwd, bwd, and opt)
  CALL write_responses_compute_response( dtype_data = dtype_data, &
                                         dtype_resp = dtype_resp, error = error )

  ! record current drag & lift
  now_drag = dtype_resp%force_lab(1,1)
  now_lift = dtype_resp%force_lab(1,2)

  IF (error.and.(dtype_resp%it.eq.1000)) THEN
        WRITE(*,*) '-----------------WRITING GEN_FORCE FILE-----------------'
        OPEN(unit=250,file="output/responds/gen_force.dat",form="formatted",status="unknown",position="append")
        WRITE(250,*) 0,0
        CLOSE(250)
        return
  END IF

  !WRITE(*,*) 'Responses computed ...' ! for debugging
  ! write response only in fwd and bwd, not opt

    ! open/create files
    IF (.not.alloc_ibd) THEN
      ! ************************  if there is no body (begins) **************************
      ! open/create files
      IF (dtype_data%it.eq.1) THEN
        OPEN(unit=100,file="output/responds/cfl.dat",form="formatted",status="replace")
      ELSE
        OPEN(unit=100,file="output/responds/cfl.dat",form="formatted",status="unknown",position="append")
      END IF
      WRITE(100,*) dtype_resp%it, dtype_resp%cfl, dt
      CLOSE(100)
      ! ************************  if there is no body (ends) **************************
    ELSE
      ! ************************  ! if ther are bodies (begins) **************************
      ! open/create files
      IF (dtype_resp%it.eq.1) THEN
        OPEN(unit=100,file="output/responds/cfl_slip.dat",form="formatted",status="replace")
        OPEN(unit=200,file="output/responds/force.dat",form="formatted",status="replace")
        OPEN(unit=300,file="output/responds/motion.dat",form="formatted",status="replace")
        OPEN(unit=400,file="output/responds/moment_power.dat",form="formatted",status="replace")
        OPEN(unit=150,file="output/responds/simple_force.dat",form="formatted",status="replace")
      ELSE
        OPEN(unit=100,file="output/responds/cfl_slip.dat",form="formatted",status="unknown",position="append")
        OPEN(unit=200,file="output/responds/force.dat",form="formatted",status="unknown",position="append")
        OPEN(unit=300,file="output/responds/motion.dat",form="formatted",status="unknown",position="append")
        OPEN(unit=400,file="output/responds/moment_power.dat",form="formatted",status="unknown",position="append")
        OPEN(unit=150,file="output/responds/simple_force.dat",form="formatted",status="unknown",position="append")
      END IF

      IF (dtype_resp%it.eq.1000) THEN
        WRITE(*,*) '-----------------WRITING GEN_FORCE FILE-----------------'
        OPEN(unit=250,file="output/responds/gen_force.dat",form="formatted",status="unknown",position="append")
        WRITE(250,*) ((dtype_resp%force_lab(i,j), j=1,2), i=1,n_body+n_actuator)
        CLOSE(250)
      END IF

      WRITE(100,*) dtype_resp%it, dtype_resp%cfl, dtype_resp%slip

      IF (Ncr.eq.0) THEN
        WRITE(200,*) dtype_resp%it, &
                     ((dtype_resp%force(i,j), j=1,2), i=1,n_body+n_actuator)

        WRITE(150,*) dtype_resp%it, &
                     ((dtype_resp%force_lab(i,j), j=1,2), i=1,n_body+n_actuator)
      ELSE
        WRITE(200,*) dtype_resp%it, &
                     ((dtype_resp%force(i,j), j=1,2), i=1,n_body+n_actuator),&
                     ((dtype_resp%force_lab(i,j), j=1,2), i=1,n_body+n_actuator)
                     
        WRITE(150,*) dtype_resp%it, &
                     ((dtype_resp%force_lab(i,j), j=1,2), i=1,n_body+n_actuator)

      END IF

      IF (Ncr.eq.0) THEN
        WRITE(300,*) dtype_data%it, &
                     dtype_data%cm_motion%trans_cm%xcm, &
                     dtype_data%cm_motion%trans_cm%vcm
      ELSE
        WRITE(300,*) dtype_data%it, &
                     dtype_data%cm_motion%theta_g2l, &
                     dtype_data%cm_motion%omega_g2l, &
                     dtype_data%cm_motion%trans_cm%xcm, &
                     dtype_data%cm_motion%trans_cm%vcm, &
                     (dtype_data%cm_motion%rot_cm(j)%rox, j=1,ncr), &
                     (dtype_data%cm_motion%rot_cm(j)%roy, j=1,ncr), &
                     (dtype_data%cm_motion%rot_cm(j)%theta, j=1,ncr), &
                     (dtype_data%cm_motion%rot_cm(j)%omega, j=1,ncr)
      END IF

      IF (Ncr.eq.0) THEN
        WRITE(400,*) dtype_resp%it, &
                     ((dtype_resp%moment(i,j), i=1,n_body+n_actuator), j=1,ncr), &
                     (dtype_resp%power_trans(i), i=1,n_body+n_actuator), &
                     ((dtype_resp%power_rot(i,j), i=1,n_body+n_actuator), j=1,ncr)
      ELSE
        WRITE(400,*) dtype_resp%it, (dtype_resp%power_trans(i), i=1,n_body+n_actuator)
      END IF

      CLOSE(100)
      CLOSE(200)
      CLOSE(300)
      CLOSE(400)
      CLOSE(150)
      ! ************************  ! if ther are bodies (begins) **************************
    END IF


  ! deallocate response metrix
  CALL dts_dealloc_response_t( dtype = dtype_resp, &
                               nbdy = n_body+n_actuator, &
                               ncr = Ncr )

END SUBROUTINE write_responses_write_file

! *****************************************************************************************
  SUBROUTINE write_responses_compute_response( dtype_data,dtype_resp, error)
    !*****************************************************************!
    !*   write immersed body forces to file                          *!
    !*****************************************************************!
    USE mmisc
    IMPLICIT NONE
    TYPE(data_t), INTENT(in) :: dtype_data
    TYPE(response_t), INTENT(inout) :: dtype_resp
    LOGICAL, INTENT(out) :: error
    INTEGER          :: i, j, k
    REAL(KIND(0.D0)), ALLOCATABLE :: xb(:), vb(:), fb(:)
    REAL(KIND(0.D0)), ALLOCATABLE :: rox(:), roy(:), rot_ang(:), rot_rate(:)
    REAL(KIND(0.D0)), ALLOCATABLE :: forcex(:), forcey(:)
    REAL(KIND(0.D0)), ALLOCATABLE :: forcex_lab(:), forcey_lab(:)
    REAL(KIND(0.D0)), ALLOCATABLE :: moment(:,:)
    REAL(KIND(0.d0)), ALLOCATABLE :: moment_total(:)
    REAL(KIND(0.D0)) :: ftemp(2)
    REAL(KIND(0.D0)) :: qb(Nq)
    REAL(KIND(0.D0)) :: vavg(1:m,1:n)
    REAL(KIND(0.D0)) :: theta_g2l, theta_rot

    dtype_resp%it = dtype_data%it
    theta_g2l = dtype_data%cm_motion%theta_g2l

    !WRITE(*,*) 'computing CFL ...' ! for debugging
    ! qb = q + q0
    qb = dtype_data%levels(1)%q%values + dtype_data%levels(1)%q0%values
    ! computing cfl = u*dt/dx
    DO j=1,n
       DO i=1,m
          vavg(i,j) = 0.5D0/delta*( ( qb(u(i,j))+qb(u(i+1,j)) )**2 + &
                      ( qb(v(i,j))+qb(v(i,j+1)) )**2 )**0.5
       END DO
    END DO
    dtype_resp%cfl = MAXVAL(vavg)*dt/delta
    error = .FALSE.
    IF(abs(dtype_resp%cfl).gt.10**4)THEN
      CFL_FLAG=.TRUE.
      error = .TRUE.
      return
      ! STOP 'ERROR: CFL is larger than 10^4'
    ENDIF

    IF ((n_body+n_actuator).ne.0) THEN
      ALLOCATE( xb(Nf), vb(Nf), fb(Nf) )
      xb = 0.D0
      vb = 0.D0
      fb = 0.D0
      IF (Ncr.ne.0) THEN
        ALLOCATE( rox(Ncr), roy(Ncr), rot_ang(Ncr), rot_rate(Ncr) )
      END IF

      xb = dtype_data%ibd%xb
      vb = dtype_data%ibd%vb
      fb = dtype_data%ibd%fb/dt
      DO j = 1, Ncr
        rox(j) = dtype_data%cm_motion%rot_cm(j)%rox
        roy(j) = dtype_data%cm_motion%rot_cm(j)%roy
        rot_ang(j) = dtype_data%cm_motion%rot_cm(j)%theta
        rot_rate(j) = dtype_data%cm_motion%rot_cm(j)%omega
      END DO

      !WRITE(*,*) 'computing slip error ...' ! for debugging
      ! computing slip
      dtype_resp%slip = MAXVAL(ABS( regT( qb ) - dtype_data%ibd%slip_vel ))

      DO i = 1, n_body
         dtype_resp%force(i,1) = 2.d0*delta*SUM( getbody(fb, 1, i, bdy(i)%npts) )
         dtype_resp%force(i,2) = 2.d0*delta*SUM( getbody(fb, 2, i, bdy(i)%npts) )
         dtype_resp%force_lab(i,:) = dtype_resp%force(i,:)
         dtype_resp%power_trans(i) = 2.d0*delta*( &
                               dot_product( getbody(fb, 1, i, bdy(i)%npts), &
                                            getbody(vb, 1, i, bdy(i)%npts) ) &
                             + dot_product( getbody(fb, 2, i, bdy(i)%npts), &
                                            getbody(vb, 2, i, bdy(i)%npts) ) ) &
                             + dot_product( dtype_resp%force(i,:) , &
                                            dtype_data%cm_motion%trans_cm%vcm )
         DO j = 1, Ncr
           IF (j.eq.1) THEN
              ftemp(1) = dtype_resp%force(i,1)
              ftemp(2) = dtype_resp%force(i,2)
           END IF
           theta_rot = theta_g2l
           !theta_rot = atan( sin(dtype_data%cm_motion%rot_cm(2)%theta), &
          !                   tsr+cos(dtype_data%cm_motion%rot_cm(2)%theta) ) &
          !             + dtype_data%cm_motion%rot_cm(1)%theta
           dtype_resp%force_lab(i,1) = &
                          ftemp(1)*cos(theta_rot) - ftemp(2)*sin(theta_rot)
           dtype_resp%force_lab(i,2) = &
                          ftemp(1)*sin(theta_rot) + ftemp(2)*cos(theta_rot)
           ftemp(1) = dtype_resp%force_lab(i,1)
           ftemp(2) = dtype_resp%force_lab(i,2)

           dtype_resp%moment(i,j) = 2.d0*delta*( &
                        dot_product( getbody(fb, 2, i, bdy(i)%npts), &
                                     getbody(xb, 1, i, bdy(i)%npts)-rox(j) ) &
                       -dot_product( getbody(fb, 1, i, bdy(i)%npts), &
                                     getbody(xb, 2, i, bdy(i)%npts)-roy(j) ) )
           dtype_resp%power_rot(i,j) = dtype_resp%moment(i,j) * rot_rate(j)
         END DO
      END DO

      DO i=1,n_actuator
        dtype_resp%force(n_body+i,1) = 2.d0*delta*SUM( getact(fb, 1, i, act(i)%npts) )
        dtype_resp%force(n_body+i,2) = 2.d0*delta*SUM( getact(fb, 2, i, act(i)%npts) )
        dtype_resp%force_lab(n_body+i,:) = dtype_resp%force(n_body+i,:)
        dtype_resp%power_trans(n_body+i) = 2.d0*delta*( &
                              dot_product( getact(fb, 1, i, act(i)%npts), &
                                           getact(vb, 1, i, act(i)%npts) ) &
                            + dot_product( getact(fb, 2, i, act(i)%npts), &
                                           getact(vb, 2, i, act(i)%npts) ) )
        DO j = 1, Ncr
          IF (j.eq.1) THEN
             ftemp(1) = dtype_resp%force(n_body+i,1)
             ftemp(2) = dtype_resp%force(n_body+i,2)
          END IF
          theta_rot = theta_g2l
          !theta_rot = atan( sin(dtype_data%cm_motion%rot_cm(2)%theta), &
          !                  tsr+cos(dtype_data%cm_motion%rot_cm(2)%theta) ) &
          !            + dtype_data%cm_motion%rot_cm(1)%theta
          dtype_resp%force_lab(n_body+i,1) = &
                          ftemp(1)*cos(theta_rot) - ftemp(2)*sin(theta_rot)
          dtype_resp%force_lab(n_body+i,2) = &
                          ftemp(1)*sin(theta_rot) + ftemp(2)*cos(theta_rot)
          ftemp(1) = dtype_resp%force_lab(n_body+i,1)
          ftemp(2) = dtype_resp%force_lab(n_body+i,2)

          dtype_resp%moment(n_body+i,j) = 2.d0*delta*( &
                       dot_product( getact(fb, 2, i, act(i)%npts), &
                                    getact(xb, 1, i, act(i)%npts)-rox(j) ) &
                      -dot_product( getact(fb, 1, i, act(i)%npts), &
                                    getact(xb, 2, i, act(i)%npts)-roy(j) ) )
          dtype_resp%power_rot(n_body+i,j) = &
                               dtype_resp%moment(n_body+i,j) * rot_rate(j)
        END DO
      END DO

      DEALLOCATE( xb, vb, fb )
      IF (Ncr.ne.0) THEN
        DEALLOCATE( rox, roy, rot_ang, rot_rate )
      END IF
    END IF

  END SUBROUTINE write_responses_compute_response

! *****************************************************************************************

END MODULE write_responses
