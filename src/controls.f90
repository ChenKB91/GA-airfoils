MODULE controls

  ! Control types :
  !   1 : x-component body force in the lab frame
  !   2 : y-component body force in the lab frame
  !   3 : x-component body force in the body-fixed frame
  !   4 : y-component body force in the body-fixed frame

  !   5 : x-component actuator velocity in the lab frame
  !   6 : y-component actuator velocity in the lab frame
  !   7 : x-component actuator velocity in the body-fixed frame
  !   8 : y-component actuator velocity in the body-fixed frame

  ! controls input example :
  ! for each control_input.xxx.inp (xxx is the number) in input directory :
  ! ctrl(xxx)%type
  ! ctrl(xxx) parameters  (if there is one)
  ! ctrl%ncsteps
  ! ( ihstart  ,  ctrl%values(ihstart)
  !  ihstart+1 , ctrl%values(ihstart+1)
  !            .
  !            .
  !            .
  !   ihstop   ,  ctrl%values(ihstop) ) if ctrl%ncsteps =/= 0
  !  Note : If ctrl%ncsteps = 0, the code will automatically assume ctrl%values(:) = 0
  !
  ! Note : parameters for control(xxx) has different chioces depend on its type
  ! for control(xxx)%type = 1 ~ 4:
  !  ctrl(xxx) parameters = ctrl%xbf, ctrl%ybf, ctrl%width
  ! for control(xxx)%type = 5 ~ 8 :
  !  ctrl(xxx) parameters = ctrl%ibody


  USE parameters
  USE dts

  IMPLICIT NONE

  !  for bodies
  INTEGER, PARAMETER :: maxcontrols = 999 ! a large number
  TYPE(control_t) :: ctrl(maxcontrols)

CONTAINS
! *****************************************************************************************

  SUBROUTINE controls_setup_control
    IMPLICIT NONE
    LOGICAL :: readinput, existctrl
    INTEGER :: i
    CHARACTER(3) :: file_num

    ! look for controls in input directory
    readinput = .TRUE.
    LABFRAME = .FALSE.
    n_control = 0
    DO WHILE (readinput)
       WRITE(file_num,"(I3.3)")  n_control+1
       INQUIRE(file="input/control_input."//file_num//".inp",exist=readinput)
       IF (readinput) THEN
          n_control=n_control+1
          ctrl(n_control)%id = n_control
          ! in each control input files
          OPEN(unit=80,file="input/control_input."//file_num//".inp",form='formatted',status='old')
          ! read control type
          READ(80,*) ctrl(n_control)%type
          ctrl(n_control)%icr = 1001 ! dummy number; shouldn't be the same as cost/constraint dummy number
          ctrl(n_control)%ibody = 0   ! which actuator you are controling
          ctrl(n_control)%xbf = 0.D0 ! dummy number
          ctrl(n_control)%ybf = 0.D0 ! dummy number
          ctrl(n_control)%width = 0.D0 ! dummy number
          ! read correspoding icr if neccessary
          IF ( (ctrl(n_control)%type.ge.1).and.(ctrl(n_control)%type.le.4) ) THEN
            READ(80,*) ctrl(n_control)%xbf, ctrl(n_control)%ybf, ctrl(n_control)%width
            IF (ctrl(n_control)%width.eq.0.D0) STOP 'ERROR : control to bodyforce should have non-zero width'
          ELSEIF ( (ctrl(n_control)%type.ge.5).and.(ctrl(n_control)%type.le.8) ) THEN
            READ(80,*) ctrl(n_control)%ibody
            IF (ctrl(n_control)%ibody.eq.0) STOP 'ERROR : control should be applied to a actuator'
          END IF
          CLOSE(80)

          ! allocate values
          ALLOCATE( ctrl(n_control)%values(istart:istop) )
          ! initialize values
          ctrl(n_control)%values = 0.d0        ! final control update

            ! show the controls
            SELECT CASE (ctrl(n_control)%type)
            CASE (1)
              WRITE(*,*) '=> Control no.', n_control, &
                         ' is controlling x-comp. bodyforce in the lab frame'
              LABFRAME = .TRUE.
            CASE (2)
              WRITE(*,*) '=> Control no.', n_control, &
                         ' is cntrolling y-comp. bodyforce in the lab frame'
              LABFRAME = .TRUE.
            CASE (3)
              WRITE(*,*) '=> Control no.', n_control, &
                         ' is controlling x-comp. bodyforce in the body-fixed frame'
            CASE (4)
              WRITE(*,*) '=> Control no.', n_control, &
                         ' is controlling y-comp. bodyforce in the body-fixed frame'
            CASE (5)
              WRITE(*,*) '=> Control no.', n_control, &
                         ' is controlling x-comp. actuator velocity in the lab frame'
              LABFRAME = .TRUE.
            CASE (6)
              WRITE(*,*) '=> Control no.', n_control, &
                         ' is controlling y-comp. actuator velocity in the lab frame'
              LABFRAME = .TRUE.
            CASE (7)
              WRITE(*,*) '=> Control no.', n_control, &
                         ' is controlling x-comp. actuator velocity in the body-fixed frame'
            CASE (8)
              WRITE(*,*) '=> Control no.', n_control, &
                         ' is controlling y-comp. actuator velocity in the body-fixed frame'
            CASE DEFAULT
              STOP '=> ERROR: Control type should be within 1 ~ 8.'
            END SELECT
         END IF
    END DO
    WRITE(*,*) '=> Import all controls. There are',n_control,'control variables.'

    ! reading in controls for all run status
    WRITE(*,*) 'importing control ...'   ! for debugging
    CALL controls_read_control

    ! copy input control to control.dat if it was not copied before
    IF ( n_control.gt.0 ) THEN
      INQUIRE(file="output/control.dat", exist = existctrl)
      IF (.not.existctrl) CALL controls_write_control
    END IF

    ! sanity check

  END SUBROUTINE controls_setup_control
! **************************************
  SUBROUTINE controls_destroy_control
    IMPLICIT NONE
    INTEGER :: i
      DO i = 1, n_control
         DEALLOCATE( ctrl(i)%values )
      END DO
  END SUBROUTINE controls_destroy_control

! *****************************************************************************************
!                 Controls, gradients, and conjugate gradients I/O                        !
! *****************************************************************************************

  SUBROUTINE controls_read_control
    IMPLICIT NONE
    INTEGER :: i_ctrl, i, itemp, icr_temp
    CHARACTER(3) :: file_num, file_num_output
    LOGICAL :: readinput

      DO i_ctrl = 1, n_control
        WRITE(file_num,"(I3.3)")  i_ctrl
        OPEN(unit=80,file="input/control_input."//file_num//".inp",form='formatted',status='old')
        ! in each control input file, type should already been read by control setup
        READ(80,*) itemp
        ! read number of steps for controls
        READ(80,*) ctrl(i_ctrl)%ncsteps

        IF ( ctrl(i_ctrl)%ncsteps .eq. nhoriz ) THEN
          ! ****  read control input depends on whether ncsteps == nhoriz (begins) ******
          ! read given control inputs
          DO i = istart, istop
            READ(80,*) itemp, ctrl(i_ctrl)%values(i)
          END DO
          CLOSE(80)
          WRITE(*,*) '=> Read all controls values over the horizon.'
        ! if ncsteps!=nhoriz
        ELSEIF ( ctrl(i_ctrl)%ncsteps.eq.0 ) THEN
          ! if ncsteps = 0 , assume inputs are zeros and rewrite input files
          CLOSE(80)
          ctrl(i_ctrl)%ncsteps = nhoriz
          WRITE(*,*) '=> assume control No.', i_ctrl,' are zeros over the horizon'
        ELSE
          CLOSE(80)
          ! if ncsteps!= nhoriz and ncsteps != 0, then stop
          WRITE(*,*) '=> Control no.',i_ctrl,'has',ctrl(i_ctrl)%ncsteps, &
                     'time steps and Type : ', ctrl(i_ctrl)%type
          WRITE(*,*) '   but horizon is with', nhoriz, 'time steps.'
          STOP 'ERROR: control length should be the same as the length of the horizon'
        END IF
        ! ****  read control input depends on whether ncsteps == nhoriz (ends) ******
      END DO

  END SUBROUTINE controls_read_control

  ! *************************************
    SUBROUTINE controls_write_control
      IMPLICIT NONE
      INTEGER :: i, i_ctrl
      CHARACTER(3) :: file_num

      OPEN(unit=80,file="output/control.dat",form='formatted',status='replace')
      DO i = istart, istop
        WRITE(80,*) i, ( ctrl(i_ctrl)%values(i), i_ctrl = 1, n_control )
      END DO
      CLOSE(80)

    END SUBROUTINE controls_write_control

! ******************************************************************************
!  Assign controls
! *******************************************************************************
  ! type 1~4 : bodyforce
  SUBROUTINE constrols_assign_bodyforces( it, theta_g2l, i_ctrl, &
                                          A_grid, xbf_grid, width )
    USE mmisc
    IMPLICIT NONE
    INTEGER, INTENT(in) :: it, i_ctrl
    REAL(KIND(0.D0)), INTENT(in) :: theta_g2l
    REAL(KIND(0.D0)), INTENT(inout) :: A_grid(2), xbf_grid(2), width
    INTEGER :: i
    REAL(KIND(0.D0)) :: A_lab(2), A_gridp(2), xbf_lab(2)

    A_grid = 0.D0 ; A_lab = 0.D0; xbf_grid = 0.D0; xbf_lab = 0.D0; width = 1.D0
    ! contribution from controls
    IF (ctrl(i_ctrl)%type.eq.1) THEN  ! ctrl type 1 : x-component body force in the lab frame
      A_lab(1) = 1.D0 ; A_lab(2) = 0.D0   ! \hat{gL}
      A_lab = A_lab * ctrl(i_ctrl)%values(it)  ! g_L
      ! note : in fwd sim, values_temp = values ; in opt, values_temp = values update by a certain distance
      A_grid = express_vector_lab2grid( A_lab, theta_g2l ) ! g_B = B g_L
      xbf_lab(1) = ctrl(i_ctrl)%xbf ; xbf_lab(2) = ctrl(i_ctrl)%ybf
      width = ctrl(i_ctrl)%width
      xbf_grid = express_vector_lab2grid( xbf_lab, theta_g2l ) ! coord in the body-fixed frame
    ELSEIF (ctrl(i_ctrl)%type.eq.2) THEN  ! ctrl type 2 : y-component body force in the lab frame
      A_lab(1) = 0.D0 ; A_lab(2) = 1.D0  ! \hat{gL}
      A_lab = A_lab * ctrl(i_ctrl)%values(it)  ! g_L
      A_grid = express_vector_lab2grid( A_lab, theta_g2l ) ! g_B = B g_L
      xbf_lab(1) = ctrl(i_ctrl)%xbf ; xbf_lab(2) = ctrl(i_ctrl)%ybf
      width = ctrl(i_ctrl)%width
      xbf_grid = express_vector_lab2grid( xbf_lab, theta_g2l ) ! coord in the body-fixed frame
    ELSEIF (ctrl(i)%type.eq.3) THEN  ! ctrl type 3 : x-component body force in the body-fixed frame
      A_grid(1) = 1.D0 ; A_grid(2) = 0.D0  ! \hat{g_B}
      A_grid = A_grid * ctrl(i_ctrl)%values(it) ! g_B
      xbf_grid(1) = ctrl(i_ctrl)%xbf ; xbf_grid(2) = ctrl(i_ctrl)%ybf
      width = ctrl(i_ctrl)%width
    ELSEIF (ctrl(i_ctrl)%type.eq.4) THEN  ! ctrl type 4 : y-component body force in the body-fixed frame
      A_grid(1) = 0.D0 ; A_grid(2) = 1.D0   ! \hat{g_B}
      A_grid = A_grid * ctrl(i_ctrl)%values(it)  ! g_B
      xbf_grid(1) = ctrl(i_ctrl)%xbf ; xbf_grid(2) = ctrl(i_ctrl)%ybf
      width = ctrl(i_ctrl)%width
    END IF

  END SUBROUTINE constrols_assign_bodyforces

! *****************************************************************************************
! type 5~8 : actuator svelocity
  SUBROUTINE constrols_assign_actuator_vel( actno, it, theta_g2l, &
                                            i_ctrl, u_grid, assignvalue )
    USE mmisc
    IMPLICIT NONE
    INTEGER, INTENT(in) :: actno, it, i_ctrl
    REAL(KIND(0.D0)), INTENT(in) :: theta_g2l
    REAL(KIND(0.D0)), INTENT(inout) :: u_grid(2)
    LOGICAL, INTENT(inout) :: assignvalue
    INTEGER :: i
    REAL(KIND(0.D0)) :: u_lab(2), u_gridp(2)

    u_grid = 0.D0 ; u_lab = 0.D0
    IF (ctrl(i_ctrl)%ibody.eq.actno) THEN
      IF (ctrl(i_ctrl)%type.eq.5) THEN  ! ctrl type 5 : x-component actuator velocity in the lab frame
        assignvalue = .TRUE.
        u_lab(1) = 1.D0 ; u_lab(2) = 0.D0  ! hat{qaL}
        u_lab = u_lab * ctrl(i_ctrl)%values(it) ! qaL
        ! Note: in fwd sim, values_temp = values ; in opt, values_temp = values update by a certain distance
        u_grid = express_vector_lab2grid( u_lab, theta_g2l ) ! qaB ! velocity in the body-fixed frame
      ELSEIF (ctrl(i_ctrl)%type.eq.6) THEN  ! ctrl type 6 : y-component actuator velocity in the lab frame
        assignvalue = .TRUE.
        u_lab(1) = 0.D0 ; u_lab(2) = 1.D0  ! hat{qaL}
        u_lab = u_lab * ctrl(i_ctrl)%values(it)  ! qaL
        u_grid = express_vector_lab2grid( u_lab, theta_g2l ) ! qaB
      ELSEIF (ctrl(i_ctrl)%type.eq.7) THEN  ! ctrl type 7 : x-component actuator velocity in the body-fixed frame
        assignvalue = .TRUE.
        u_grid(1) = 1.D0; u_grid(2) = 0.D0   ! \hat{qaB}
        u_grid = u_grid * ctrl(i_ctrl)%values(it)   ! qaB
      ELSEIF (ctrl(i_ctrl)%type.eq.8) THEN  ! ctrl type 8 : y-component actuator velocity in the body-fixed frame
        assignvalue = .TRUE.
        u_grid(1) = 0.D0; u_grid(2) = 1.D0  ! \hat{qaB}
        u_grid = u_grid * ctrl(i_ctrl)%values(it)  ! qaB
      END IF
    END IF

  END SUBROUTINE constrols_assign_actuator_vel

! *****************************************************************************************

END MODULE controls
