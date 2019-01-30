MODULE user

  USE parameters
  USE dts
  USE grid
  USE mmisc
  USE controls

  IMPLICIT NONE

  ! module for user-defined routines to move bodies relative to the grid, add actuators, or move bodies with the grid.
  !
  ! see specific examples in the test cases

CONTAINS
!*****************************************************************!
 !*****************************************************************!
 !                  Assign user-defined motion                     !
 !*****************************************************************!
  SUBROUTINE user_center_of_mass_motion_grid( it, dtype_data )
    ! only be used in fwd sim
    IMPLICIT NONE
    INTEGER, INTENT(in) :: it
    TYPE(data_t), INTENT(inout) :: dtype_data
    INTEGER  :: i
    REAL(KIND(0.D0)) :: u_grid(2), u_lab(2), x_grid(2), x_lab(2)
    REAL(KIND(0.D0)) :: ang, rate, rox, roy, theta_g2l, omega_g2l

    ! motion_grid specifies the (possibly time-varying)
    ! u component, v component of the velocity, angular velocity,
    ! and the coordinates of the center of mass

    ! ***********  rotational motion  ***********
    ! Here assign the rotation angle and rotating rate and update theta_g2l, omega_g2l
    !do not change
    theta_g2l = 0.D0; omega_g2l = 0.D0
    DO i = 1, Ncr
      !do not change
      rate = 0.D0
      ang = 0.D0
      rox = 0.D0
      roy = 0.D0
      ! *********** define the user-defined rotational motion (starts) ! ***********
      IF (i.eq.1) THEN
          rate = 0.D0!1.D0/vawt_rad
          ang = - aoa !-15.D0*pi/180.D0 !rate * it * dt    ! integration of rate
          rox = 0.D0
          roy = 0.D0
      END IF
      ! *********** define the user-defined rotational motion (ends) ! ***********
      ! assign values
      dtype_data%cm_motion%rot_cm(i)%theta = ang
      dtype_data%cm_motion%rot_cm(i)%omega = rate
      dtype_data%cm_motion%rot_cm(i)%rox = rox
      dtype_data%cm_motion%rot_cm(i)%roy = roy
      theta_g2l = theta_g2l + dtype_data%cm_motion%rot_cm(i)%theta
      omega_g2l = omega_g2l + dtype_data%cm_motion%rot_cm(i)%omega
    END DO
    
    dtype_data%cm_motion%theta_g2l = theta_g2l
    dtype_data%cm_motion%omega_g2l = omega_g2l
    
! ***********  translational motion  ***********
    ! user-define labframe translational velocity and motion
    ! **** define the user-defined labframe translational motion (starts) ! ****
    u_lab(1) = -1.D0
    u_lab(2) = 0.d0 !0.3D0*sin(2.D0*pi*0.2D0*dt*it)
    ! this is use for initial condition of xcm
    x_lab(1) = 0.D0 !-1.D0 * it * dt
    x_lab(2) = 0.D0 !0.3D0/(2.D0*pi*0.2D0)*(-cos(2.D0*pi*0.2D0*dt*it))
    ! **** define the user-defined labframe translational motion (ends) ! ****

    !WRITE(*,*) 'theta_g2l = ' ,theta_g2l*180.D0/pi, ' deg'    ! for debugging
    x_grid = express_vector_lab2grid( x_lab, theta_g2l ) ! coordinate in the body-fixed frame
    u_grid = express_vector_lab2grid( u_lab, theta_g2l ) ! velocity in the body-fixed frame
    ! assign value
    dtype_data%cm_motion%trans_cm%xcm = x_grid
    dtype_data%cm_motion%trans_cm%vcm = u_grid

  END SUBROUTINE user_center_of_mass_motion_grid

 ! ****************************************************************!
 !*****************************************************************!
 !             Assign user-defined actuator velocity               !
 !*****************************************************************!
 SUBROUTINE user_actuator_vel( actno, it, theta_g2l, x, y, vx, vy, i_ctrl )
   IMPLICIT NONE
   INTEGER, INTENT(in) :: actno, it, i_ctrl
   REAL(KIND(0.D0)), INTENT(in) :: theta_g2l
   REAL(KIND(0.D0)), DIMENSION(:) :: x,y,vx,vy
   INTEGER :: i
   REAL(KIND(0.D0)) :: u_grid(2), u_lab(2), u_gridp(2)
   LOGICAL :: assignvalue

   ! DO NOT MODIFY x,y

   vx = 0.D0 ! do not remove; overwrite values below inside if block
   vy = 0.D0 ! do not remove; overwrite values below inside if block

   ! for each actuator, specify the imposed velocity as a function of time

   ! inputs
   !  actno : the actuator number (the number x in the filename actuator.00x.inp
   !  it : the time step (multiply by dt to get physical time)
   !  x,y  : the positions of each point defining the actuator
   !
   ! outputs
   !  vx, vy : the velocity components

   ! contribution from controls
   u_grid = 0.D0 ; u_lab = 0.D0
   assignvalue =.FALSE.
   CALL constrols_assign_actuator_vel( actno = actno, &
                                       it = it, &
                                       theta_g2l = theta_g2l, &
                                       i_ctrl = i_ctrl, &
                                       u_grid = u_grid, &
                                       assignvalue = assignvalue )

   ! if no control on actuator velocity or not use on computing control gradient,
   ! use user-defined velocity
   IF ( .not.assignvalue ) THEN
     ! *********** define the user-defined actuator velocity (starts) ! ***********
     IF (actno.eq.1) THEN ! for actuator 1
        u_grid(1) = 0.D0 ! Steady actuation with speed 0 in x dirn.
        u_grid(2) = 0.D0
     END IF
     ! *********** define the user-defined actuator velocity (ends) ! ***********
   END IF

   ! do not remove
   vx = u_grid(1)
   vy = u_grid(2)

 END SUBROUTINE user_actuator_vel

 !*****************************************************************!
 !*****************************************************************!
 !             Assign user-defined bodyforce                       !
 !*****************************************************************!
FUNCTION bodyforce(it, theta_g2l, xin, yin, uin, vin, i_ctrl)
   IMPLICIT NONE
   INTEGER, INTENT(in) :: it, i_ctrl
   REAL(KIND(0.D0)), INTENT(in) :: theta_g2l, xin, yin, uin, vin
   REAL(KIND(0.D0)):: bodyforce(2)
   INTEGER :: i
   REAL(KIND(0.D0)) :: A_grid(2), xbf_grid(2), width

   ! specify a body force (per unit mass) on RHS of u-momentum equation,
   ! as a function of location (x,y) and possibly velocity (u,v) and possibly time

   ! the specified function should decay to zero near the edges of the (mgrid=1)
   ! inner domain
   !
   ! units are in terms of velocity (result should be consistent
   ! with U^2/L units

   bodyforce = 0.d0 ! do not remove; overwrite values below inside if block

   A_grid = 0.D0 ; xbf_grid = 0.D0; width = 1.D0
   IF (i_ctrl.gt.0) THEN
     IF ( (ctrl(i_ctrl)%type.ge.1).and.(ctrl(i_ctrl)%type.le.4) ) THEN
       ! contribution from controls
       CALL constrols_assign_bodyforces( it = it, &
                                         theta_g2l = theta_g2l, &
                                         i_ctrl = i_ctrl, &
                                         A_grid = A_grid, &
                                         xbf_grid = xbf_grid, &
                                         width = width )
     END IF
   ! *********** define the user-defined bodyforces (starts) ! ***********
   ELSEIF (i_ctrl.eq.0) THEN
       A_grid(1) = 0.D0 ! Steady blowing with a magnitude 0 in x dirn.
       A_grid(2) = 0.D0 ! Steady blowing with a magnitude 0 in y dirn.
       xbf_grid(1) = 0.D0
       xbf_grid(2) = 0.D0
       width = 1.D0
   ! *********** define the user-defined bodyforces (ends) ! ***********
   END IF

   bodyforce = bodyforce + A_grid * &
                   gaussianforce( xin, yin, xbf_grid(1), xbf_grid(2), width )

 END FUNCTION bodyforce

 !*****************************************************************!
 !*****************************************************************!
 !        Assign user-defined initial vorticity field              !
 !*****************************************************************!
  FUNCTION initial_vorticity( xx, yy ) RESULT(vort)

    !***************************************************************!
    !*  Initialize vorticty and velocity for the first time step   *!
    !***************************************************************!

   ! specify a inital vorticity and velocity flow field at it=0
   ! vorticity and velocity should decay to zero at the edge of the largest domain
   ! om (omega) is the cell circulation, and q (vel) velocity flux
   ! commented line shows the case for the Taylor vortex

   IMPLICIT NONE
   REAL(KIND(0.D0)), INTENT(in) :: xx, yy
   REAL(KIND(0.D0)) :: vort

   vort = 0.D0 ! do not remove

  ! ******* define the user-defined initial vorticity field   (starts) ! *********
  ! initialize taylor vortex at (xc,yc)
  ! vort = taylor_vortex(x = xx,y = yy,r = r,xc = xc,yc = yc, t = t_initial)
  ! ******* define the user-defined initial vorticity field   (ends) ! *********

  END FUNCTION initial_vorticity

!*****************************************************************!

END MODULE user
