MODULE dts

  ! This module defines common derived types
  USE parameters
  IMPLICIT NONE
!=============================================================================!
  ! derived types for general flow field data
  TYPE volumetric_t
    REAL(KIND(0.D0)), ALLOCATABLE :: values(:,:) ! pn(1:m,1:n)
  END TYPE volumetric_t

  TYPE flux_t
    REAL(KIND(0.D0)), ALLOCATABLE :: values(:) ! q(nq)
  END TYPE flux_t

  TYPE circulation_t
    REAL(KIND(0.D0)), ALLOCATABLE :: values(:,:)  ! wn(2:m,2:n)
  END TYPE circulation_t

! values -> : = vector
 
  TYPE data_level_t
    INTEGER :: level_id
    TYPE(flux_t) :: q, q0, q0p, nl_bf_old
    TYPE(flux_t), ALLOCATABLE :: q0r(:)
    TYPE(circulation_t) :: omega, stfn
    TYPE(volumetric_t) :: pressure, pressure_dynamic
  END TYPE data_level_t
!=============================================================================!
  ! derived types for general immersed boundary data
  TYPE ib_code_t
     INTEGER, ALLOCATABLE :: code(:) ! code(nf)
  END TYPE ib_code_t

  TYPE ib_data_t
     REAL(KIND(0.D0)), ALLOCATABLE :: xb(:), vb(:), fb(:), slip_vel(:)  ! xb/vb/ab/fb(nf)
  END TYPE ib_data_t

!=============================================================================!
  ! special types for Chen's vawt simulation
  TYPE translation_t
    ! position/velocity/acceleration of center of mass
    REAL(KIND(0.D0)) :: xcm(2), vcm(2)
  END TYPE translation_t

  TYPE rotation_t
    REAL(KIND(0.D0)) :: rox, roy  ! center of rotation
    REAL(KIND(0.D0)) :: theta, omega ! rotational angle/velocity/acceleration
  END TYPE rotation_t

  TYPE motion_t
    TYPE(translation_t) :: trans_cm
    TYPE(rotation_t), ALLOCATABLE :: rot_cm(:)
    REAL(KIND(0.D0)) :: theta_g2l   ! The angle to rotate vectors in the grid-fixed frame to the lab frame
    REAL(KIND(0.D0)) :: omega_g2l   ! d(theta_g2l)/dt
  END TYPE motion_t
!=============================================================================!
  ! General data
  TYPE data_t
    INTEGER :: it
    TYPE(data_level_t), ALLOCATABLE :: levels(:)  ! levels(mgridlev)
    TYPE(ib_data_t) :: ibd
    TYPE(ib_code_t) :: codeb
    TYPE(motion_t)  :: cm_motion
  END TYPE data_t
!=============================================================================!
  !  for Cholesky (stationary body)
  TYPE chol_data_t
    LOGICAL :: recompute
    REAL(KIND(0.D0)), ALLOCATABLE :: cholmat(:,:)  ! cholmat(Nf,Nf)
    REAL(KIND(0.D0)), ALLOCATABLE :: cholvec(:)    ! cholvec(Nf)
  END TYPE chol_data_t
!=============================================================================!
  ! a special type for bodies
  TYPE body_t
     INTEGER :: npts   ! number of points in this body
     INTEGER :: ipos   ! position in overall array
     REAL(KIND(0.D0)), ALLOCATABLE :: x(:), y(:)  ! x/y(npts)
  END TYPE body_t

  TYPE actuator_t
     INTEGER :: npts   ! number of points in this actuator
     INTEGER :: ipos ! position in overall array
     ! s is the strength of actuator for multipoint actuators
     REAL(KIND(0.D0)), ALLOCATABLE :: x(:), y(:), s(:) ! x/y/s(npts)
  END TYPE actuator_t
!=============================================================================!
! type for outputing response
  TYPE response_t
    INTEGER :: it
    REAL(KIND(0.D0)) :: cfl, slip
    REAL(KIND(0.D0)), ALLOCATABLE :: force(:,:)      ! force(nbdy,2)
    REAL(KIND(0.D0)), ALLOCATABLE :: force_lab(:,:)  ! force(nbdy,2)
    REAL(KIND(0.D0)), ALLOCATABLE :: moment(:,:)     ! moment(nbdy,ncr)
    REAL(KIND(0.D0)), ALLOCATABLE :: power_trans(:)  ! power_trans(nbdy)
    REAL(KIND(0.D0)), ALLOCATABLE :: power_rot(:,:)  ! power_rot(nbdy,ncr)
  END TYPE response_t
!=============================================================================!
  ! derived types for control variables and parameters
    TYPE control_t
        INTEGER :: id
        INTEGER :: ncsteps
        INTEGER :: type
        INTEGER :: icr
        INTEGER :: ibody
        REAL(KIND(0.D0)) :: xbf, ybf, width
        REAL(KIND(0.D0)), ALLOCATABLE :: values(:) ! control(nh)
    END TYPE control_t

!=============================================================================!

CONTAINS ! all global before this

!=============================================================================!
SUBROUTINE dts_alloc_volumetric_t(dtype)
  IMPLICIT NONE
  TYPE(volumetric_t), INTENT(inout) :: dtype
  ALLOCATE(dtype%values(1:m,1:n))
  dtype%values = 0.D0
END SUBROUTINE dts_alloc_volumetric_t
!=======================================!
SUBROUTINE dts_dealloc_volumetric_t(dtype)
  IMPLICIT NONE
  TYPE(volumetric_t), INTENT(inout) :: dtype
  DEALLOCATE(dtype%values)
END SUBROUTINE dts_dealloc_volumetric_t
!=============================================================================!

SUBROUTINE dts_alloc_flux_t(dtype)
  IMPLICIT NONE
  TYPE(flux_t), INTENT(inout) :: dtype
  INTEGER :: Nq
  Nq = (m+1)*n + m*(n+1)
  ALLOCATE(dtype%values(Nq))
  dtype%values = 0.D0
END SUBROUTINE dts_alloc_flux_t
!=======================================!
SUBROUTINE dts_dealloc_flux_t(dtype)
  IMPLICIT NONE
  TYPE(flux_t), INTENT(inout):: dtype
  DEALLOCATE(dtype%values)
END SUBROUTINE dts_dealloc_flux_t
!=============================================================================!

SUBROUTINE dts_alloc_circulation_t(dtype)
  IMPLICIT NONE
  TYPE(circulation_t), INTENT(inout) :: dtype
  ALLOCATE(dtype%values(2:m,2:n))
  dtype%values = 0.D0
END SUBROUTINE dts_alloc_circulation_t
!=======================================!
SUBROUTINE dts_dealloc_circulation_t(dtype)
  IMPLICIT NONE
  TYPE(circulation_t), INTENT(inout) :: dtype
  DEALLOCATE(dtype%values)
END SUBROUTINE dts_dealloc_circulation_t
!=============================================================================!

SUBROUTINE dts_alloc_data_level_t(dtype)
  IMPLICIT NONE
  TYPE(data_level_t), INTENT(inout) :: dtype
  INTEGER :: i
  CALL dts_alloc_flux_t( dtype = dtype%q         )
  CALL dts_alloc_flux_t( dtype = dtype%q0        )
  CALL dts_alloc_flux_t( dtype = dtype%q0p       )
  CALL dts_alloc_flux_t( dtype = dtype%nl_bf_old )
  ALLOCATE( dtype%q0r(Ncr) )
  DO i = 1, Ncr
    CALL dts_alloc_flux_t( dtype = dtype%q0r(i) )
  END DO
  CALL dts_alloc_circulation_t( dtype = dtype%omega )
  CALL dts_alloc_circulation_t( dtype = dtype%stfn  )
  IF (compute_pressure) THEN
    CALL dts_alloc_volumetric_t( dtype = dtype%pressure )
    CALL dts_alloc_volumetric_t( dtype = dtype%pressure_dynamic )
  END IF
END SUBROUTINE dts_alloc_data_level_t
!=======================================!
SUBROUTINE dts_dealloc_data_level_t(dtype)
  IMPLICIT NONE
  TYPE(data_level_t), INTENT(inout) :: dtype
  INTEGER :: i
  CALL dts_dealloc_flux_t( dtype = dtype%q         )
  CALL dts_dealloc_flux_t( dtype = dtype%q0        )
  CALL dts_dealloc_flux_t( dtype = dtype%q0p       )
  CALL dts_dealloc_flux_t( dtype = dtype%nl_bf_old )
  DO i = 1, Ncr
    CALL dts_dealloc_flux_t( dtype = dtype%q0r(i) )
  END DO
  DEALLOCATE(dtype%q0r)
  CALL dts_dealloc_circulation_t( dtype = dtype%omega )
  CALL dts_dealloc_circulation_t( dtype = dtype%stfn  )
  IF (compute_pressure) THEN
    CALL dts_dealloc_volumetric_t( dtype = dtype%pressure )
    CALL dts_dealloc_volumetric_t( dtype = dtype%pressure_dynamic )
  END IF
END SUBROUTINE dts_dealloc_data_level_t
!=============================================================================!

SUBROUTINE dts_alloc_ib_code_t(dtype,nb)
  IMPLICIT NONE
  TYPE(ib_code_t), INTENT(inout) :: dtype
  INTEGER, INTENT(in) :: nb
  ALLOCATE(dtype%code(nb))
  dtype%code = 0
END SUBROUTINE dts_alloc_ib_code_t
!=======================================!
SUBROUTINE dts_dealloc_ib_code_t(dtype)
  IMPLICIT NONE
  TYPE(ib_code_t) :: dtype
  DEALLOCATE(dtype%code)
END SUBROUTINE dts_dealloc_ib_code_t
!=============================================================================!

SUBROUTINE dts_alloc_ib_data_t(dtype,nf)
  IMPLICIT NONE
  TYPE(ib_data_t), INTENT(inout) :: dtype
  INTEGER, INTENT(in) :: nf
  ALLOCATE(dtype%xb(nf), dtype%vb(nf), dtype%fb(nf), dtype%slip_vel(nf) )
  dtype%xb = 0.D0
  dtype%vb = 0.D0
  dtype%fb = 0.D0
  dtype%slip_vel = 0.D0
END SUBROUTINE dts_alloc_ib_data_t
!=======================================!
SUBROUTINE dts_dealloc_ib_data_t(dtype)
  IMPLICIT NONE
  TYPE(ib_data_t), INTENT(inout) :: dtype
  DEALLOCATE(dtype%xb, dtype%vb, dtype%fb, dtype%slip_vel)
END SUBROUTINE dts_dealloc_ib_data_t

!=============================================================================!

SUBROUTINE dts_alloc_motion_t(dtype)
  IMPLICIT NONE
  TYPE(motion_t), INTENT(inout) :: dtype
  INTEGER :: i
  dtype%trans_cm%xcm = 0.d0
  dtype%trans_cm%vcm = 0.d0
  IF (Ncr.ne.0) THEN
    ALLOCATE(dtype%rot_cm(Ncr))
    DO i = 1, Ncr
      dtype%rot_cm(i)%rox = 0.d0
      dtype%rot_cm(i)%roy = 0.d0
      dtype%rot_cm(i)%theta = 0.d0
      dtype%rot_cm(i)%omega = 0.d0
    END DO
  END IF
  dtype%theta_g2l = 0.D0
  dtype%omega_g2l = 0.D0
END SUBROUTINE dts_alloc_motion_t
!=======================================!
SUBROUTINE dts_dealloc_motion_t(dtype)
  IMPLICIT NONE
  TYPE(motion_t), INTENT(inout) :: dtype
  IF (Ncr.ne.0) THEN
    DEALLOCATE(dtype%rot_cm)
  END IF
END SUBROUTINE dts_dealloc_motion_t
!=============================================================================!

SUBROUTINE dts_alloc_data_t(dtype,nb,alloc_ibd)
  IMPLICIT NONE
  TYPE(data_t), INTENT(inout) :: dtype
  INTEGER, INTENT(in) :: nb
  LOGICAL, INTENT(in) :: alloc_ibd
  INTEGER :: i
  dtype%it = 0
  ALLOCATE(dtype%levels(mgridlev))
  DO i = 1, mgridlev
    CALL dts_alloc_data_level_t(dtype = dtype%levels(i))
    dtype%levels(i)%level_id = i
  END DO
  IF (alloc_ibd) THEN
    CALL dts_alloc_ib_code_t(dtype = dtype%codeb, nb = nb)
    CALL dts_alloc_ib_data_t(dtype = dtype%ibd, nf = 2*nb)
  END IF
  CALL dts_alloc_motion_t(dtype = dtype%cm_motion)
END SUBROUTINE dts_alloc_data_t
!=======================================!
SUBROUTINE dts_dealloc_data_t(dtype,alloc_ibd)
  IMPLICIT NONE
  TYPE(data_t), INTENT(inout) :: dtype
  LOGICAL, INTENT(in) :: alloc_ibd
  INTEGER :: i
  DO i = 1, mgridlev
    CALL dts_dealloc_data_level_t(dtype = dtype%levels(i))
  END DO
  DEALLOCATE(dtype%levels)
  IF (alloc_ibd) THEN
    CALL dts_dealloc_ib_code_t(dtype = dtype%codeb)
    CALL dts_dealloc_ib_data_t(dtype = dtype%ibd)
  END IF
  CALL dts_dealloc_motion_t(dtype = dtype%cm_motion)
END SUBROUTINE dts_dealloc_data_t
!=============================================================================!
SUBROUTINE dts_alloc_chol_data_t(dtype,nf)
  IMPLICIT NONE
  TYPE(chol_data_t), INTENT(inout) :: dtype
  INTEGER, INTENT(in) :: nf
  ALLOCATE(dtype%cholmat(nf,nf), dtype%cholvec(nf))
  dtype%recompute =.FALSE.
  dtype%cholmat = 0.D0
  dtype%cholvec = 0.D0
END SUBROUTINE dts_alloc_chol_data_t
!=======================================!
SUBROUTINE dts_dealloc_chol_data_t(dtype)
  IMPLICIT NONE
  TYPE(chol_data_t), INTENT(inout) :: dtype
  DEALLOCATE(dtype%cholmat, dtype%cholvec)
END SUBROUTINE dts_dealloc_chol_data_t
!=============================================================================!
SUBROUTINE dts_alloc_response_t(dtype,nbdy,ncr)
  IMPLICIT NONE
  TYPE(response_t), INTENT(inout) :: dtype
  INTEGER, INTENT(in) :: nbdy, ncr
  IF (nbdy.ne.0) THEN
    ALLOCATE( dtype%force(nbdy,2), dtype%force_lab(nbdy,2) )
    ALLOCATE( dtype%power_trans(nbdy) )
    dtype%force = 0.D0
    dtype%force_lab = 0.D0
    dtype%power_trans = 0.D0
    IF (ncr.ne.0) THEN
        ALLOCATE( dtype%moment(nbdy,ncr), dtype%power_rot(nbdy,ncr) )
        dtype%moment = 0.D0
        dtype%power_rot = 0.D0
    END IF
  END IF
END SUBROUTINE dts_alloc_response_t
!=======================================!
SUBROUTINE dts_dealloc_response_t(dtype,nbdy,ncr)
  IMPLICIT NONE
  TYPE(response_t), INTENT(inout) :: dtype
  INTEGER, INTENT(in) :: nbdy, ncr
  IF (nbdy.ne.0) THEN
    DEALLOCATE( dtype%force, dtype%force_lab, dtype%power_trans )
    IF (ncr.ne.0) THEN
        DEALLOCATE( dtype%moment, dtype%power_rot )
    END IF
  END IF
END SUBROUTINE dts_dealloc_response_t
!=============================================================================!

END MODULE dts
