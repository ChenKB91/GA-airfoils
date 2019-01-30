MODULE mmisc_motion

  ! This module contains subroutines and functions relate to motion
  USE parameters
  USE dts
  USE grid
  IMPLICIT NONE

CONTAINS

  !*****************************************************************!

  FUNCTION motion_translation( ux_in, uy_in, ilev ) RESULT(qt)
    IMPLICIT NONE
    INTEGER, INTENT(in) :: ilev
    REAL(KIND(0.D0)), INTENT(in) :: ux_in, uy_in
    REAL(KIND(0.D0)) :: fac
    REAL(KIND(0.D0)) :: qt(Nq)

    qt = 0.d0
    fac = delta*2.d0**(REAL(ilev)-1.D0)  ! cell face length on each grid
    ! qt = - U_b
    qt(1:(m+1)*n)    = - ux_in
    qt((m+1)*n+1:Nq) = - uy_in
    qt = fac * qt

  END FUNCTION motion_translation

  !*****************************************************************!

  FUNCTION motion_rotation( rate_in, rox_in, roy_in, ilev ) RESULT(qr)
    IMPLICIT NONE
    INTEGER, INTENT(in) :: ilev
    REAL(KIND(0.D0)), INTENT(in) :: rate_in, rox_in, roy_in
    INTEGER :: i, j
    REAL(KIND(0.D0)) :: xx, yy, fac
    REAL(KIND(0.D0)) :: qr(Nq)

    qr = 0.d0
    fac = delta*2.d0**(REAL(ilev)-1.D0)  ! cell face length on each grid
    ! q0r = - omega x (r-ro)
    DO i=1,m+1
       DO j=1,n
         ! find coord where x-comp velocity lives
          xx =  (REAL(i)-1.D0-REAL(m)/2.D0) * fac + REAL(m)/2.D0*delta - offsetx
          yy =  (REAL(j)-0.5D0-REAL(n)/2.D0) * fac + REAL(n)/2.D0*delta -offsety
          qr(u(i,j)) = + rate_in *( yy - roy_in )
       END DO
    END DO
    DO i=1,m
       DO j=1,n+1
         ! find coord where y-comp velocity lives
          xx =  (REAL(i)-0.5d0-REAL(m)/2.D0)*delta*(2.D0**(REAL(ilev)-1.D0)) &
               + REAL(m)/2.D0*delta - offsetx
          yy =  (REAL(j)-1.D0-REAL(n)/2.D0)*delta*(2.D0**(REAL(ilev)-1.D0)) &
               + REAL(n)/2.D0*delta -offsety
          qr(v(i,j))  = - rate_in *( xx - rox_in )
       END DO
    END DO
    qr = fac * qr

  END FUNCTION motion_rotation

  !*****************************************************************!

  SUBROUTINE mmisc_motion_cm2q0( dtype_data )
    IMPLICIT NONE
    TYPE(data_t), INTENT(inout) :: dtype_data
    INTEGER :: i, k

    DO k = 1, mgridlev
      ! q0p = - U_b
      dtype_data%levels(k)%q0p%values = &
       motion_translation( ux_in = dtype_data%cm_motion%trans_cm%vcm(1), &
                           uy_in = dtype_data%cm_motion%trans_cm%vcm(2), &
                           ilev = k )
     ! add q0p to q0
      dtype_data%levels(k)%q0%values = dtype_data%levels(k)%q0p%values
      ! q0r = - U_omeag
      DO i = 1, Ncr
        dtype_data%levels(k)%q0r(i)%values = &
        motion_rotation( rate_in = dtype_data%cm_motion%rot_cm(i)%omega, &
                         rox_in =  dtype_data%cm_motion%rot_cm(i)%rox, &
                         roy_in =  dtype_data%cm_motion%rot_cm(i)%roy, &
                         ilev = k )
        ! add q0r to q0
        dtype_data%levels(k)%q0%values = &
                dtype_data%levels(k)%q0%values + &
                dtype_data%levels(k)%q0r(i)%values
      END DO
    END DO

  END SUBROUTINE mmisc_motion_cm2q0

  !*****************************************************************!

  FUNCTION rotate_vector( v_in, theta ) RESULT(v_out)
    IMPLICIT NONE
    REAL(KIND(0.D0)), INTENT(in) :: v_in(2), theta
    REAL(KIND(0.D0)) :: v_out(2)

    v_out(1) = v_in(1) * cos(theta) - v_in(2) * sin(theta)
    v_out(2) = v_in(1) * sin(theta) + v_in(2) * cos(theta)

  END FUNCTION rotate_vector

  !*****************************************************************!

  FUNCTION express_vector_grid2lab( v_grid_in, theta_g2l ) RESULT(v_lab_out)
    IMPLICIT NONE
    REAL(KIND(0.D0)), INTENT(in) :: v_grid_in(2), theta_g2l
    REAL(KIND(0.D0)) :: v_lab_out(2)

    v_lab_out = rotate_vector( v_in = v_grid_in, theta = theta_g2l )

  END FUNCTION express_vector_grid2lab

  !*****************************************************************!

  FUNCTION express_vector_lab2grid( v_lab_in, theta_g2l ) RESULT(v_grid_out)
    IMPLICIT NONE
    REAL(KIND(0.D0)), INTENT(in) :: v_lab_in(2), theta_g2l
    REAL(KIND(0.D0)) :: v_grid_out(2)

    v_grid_out = rotate_vector( v_in = v_lab_in, theta = - theta_g2l )

  END FUNCTION express_vector_lab2grid

!*****************************************************************!

END MODULE mmisc_motion
