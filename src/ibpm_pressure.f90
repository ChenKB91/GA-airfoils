 MODULE ibpm_pressure

  !******************************************************************!
  !*   Module to solve for the total pressure (p+.5*u^2), from a    *!
  !*   given velocity and vorticity field                           *!
  !*   The pressure is defined at the cell centers                  *!
  !*   New features                                                 *!
  !*     - Calculate pressure using Sine FFT                        *!
  !*   New updates 2013.01.17 by Hsieh-Chen Tsai                    *!
  !*     - Calculate pressure in the rotating frame                 *!
  !******************************************************************!

  USE parameters
  USE dts
  USE grid
  IMPLICIT NONE

CONTAINS

!************************************************************************************!

 SUBROUTINE ibpm_pressure_calculate_pressure( dtype_data )
   !***************************************************************!
   !*  Multiscale method to solve D D^T pressure = pressure_rhs   *!
   !*  D D^T pressure = D N(q) - D E^T f - d/dt(mdot) + L'mdot    *!
   !***************************************************************!
   USE myfft
   USE mmisc
   IMPLICIT NONE

   TYPE(data_t), INTENT(inout) :: dtype_data
   INTEGER :: k
   REAL(KIND(0.D0)) :: pressure_rhs_temp(1:m,1:n,mgridlev)  !global variable to store pressure rhs
   REAL(KIND(0.D0)) :: pressure_rhs(1:m,1:n)     ! (at centers)
   REAL(KIND(0.D0)) :: pressure_bc(2*m+2*n)   ! pressure at the boundaries

   ! ======================== DD^T (0.5*|u|^2 + P + 0.5*omegab^2*|x|^2) = RHS    ==================== !
   ! = At boundary, P=0 --> BC = 0.5*|q0/del|^2 - 0.5*|cross(omegab,x)|^2 + dot(U,cross(omegab,x))  = 0.5*|q0p/del|^2!

   ! compute the source term of the pressure equation on each level
   !WRITE(*,*) 'Compute pressure rhs'
   CALL ibpm_pressure_calculate_pressure_rhs( dtype_data = dtype_data, &
                                              pressure_rhs = pressure_rhs_temp )

   ! generate dynamic pressure due to local velocity and rotation
   !WRITE(*,*) 'Compute dynamic pressure'
   CALL ibpm_pressure_set_dynamic_pressure( dtype_data = dtype_data )

   ! before here pressure_rhs_temp should be grid size independent
   ! ==== solve for pressure and telescope in from the coarsest domain ===== !
   DO k = mgridlev,1,-1  ! === Telescope in to solve for pressure on SMALLER DOMAIN

      IF (k.eq.mgridlev) THEN
        ! get boundary condition on the coarsest domain
        CALL getbc_pressurebc_coarsest( dtype_data = dtype_data, &
                                        pressure_bc = pressure_bc  )
      ELSE
        ! get pressure_bc from k+1 level
        CALL mmisc_bc_get_bc_phi( dtype_data%levels(k+1)%pressure%values, &
                                  pressure_bc )
      END IF

      ! pressure_rhs at k level
      pressure_rhs = pressure_rhs_temp(:,:,k)

      ! apply pressure_bc
      CALL mmisc_bc_apply_bc_phi( pressure_rhs, pressure_bc )

      ! compute new pressure at k level  ! pressure=(D D^T)-1 rhs(:,:,mgridlev)
      dtype_data%levels(k)%pressure%values = ddti( pressure_rhs )

   END DO

   ! ====================================================================== !
   !pressure_stag = pressure ! pressure_stag = pressure_static + pressure_dynamic [== 0.5 rho |u|^2 ]
   DO k = 1, mgridlev
     dtype_data%levels(k)%pressure%values = &
                            dtype_data%levels(k)%pressure%values &
                          - dtype_data%levels(k)%pressure_dynamic%values
   END DO

END SUBROUTINE ibpm_pressure_calculate_pressure

!************************************************************************************!

SUBROUTINE ibpm_pressure_calculate_pressure_rhs( dtype_data, pressure_rhs )
  USE mmisc
  IMPLICIT NONE
  TYPE(data_t), INTENT(in) :: dtype_data
  REAL(KIND(0.D0)) :: pressure_rhs(1:m,1:n,mgridlev)
  REAL(KIND(0.D0)) :: del
  INTEGER :: k

  ! compute the source term of the pressure equation on each level
  pressure_rhs = 0.0d0
  DO k = 1, mgridlev   ! coarsify onto bigger domain
    IF (k.eq.1) THEN
      pressure_rhs(:,:,1) = - divergence( reg(dtype_data%ibd%fb/dt) ) ! E^T(fb)
    ELSE
      pressure_rhs(:,:,k) = coarsify_pressure( pressure_rhs(:,:,k), &
                                               pressure_rhs(:,:,k-1) )
    END IF
  ENDDO

  DO k = 1, mgridlev   ! Add D(N(q)) to RHS at each level
     del =  delta * (2.d0**REAL(k-1))
     pressure_rhs(:,:,k) = pressure_rhs(:,:,k) &
             + divergence( dtype_data%levels(k)%nl_bf_old%values )/(del**2.D0)
  END DO

END SUBROUTINE ibpm_pressure_calculate_pressure_rhs

!************************************************************************************!

SUBROUTINE ibpm_pressure_set_dynamic_pressure( dtype_data )
  IMPLICIT NONE
  TYPE(data_t), INTENT(inout) :: dtype_data
  INTEGER :: i, j, k
  REAL(KIND(0.D0)) :: del
  REAL(KIND(0.D0)) :: q(Nq), q0p(Nq)
  REAL(KIND(0.D0)) :: uu, vv, u0p, v0p

  ! generate dynamic pressure due to local velocity and rotation
  DO k = 1,mgridlev
     del =  delta * (2.d0**REAL(k-1))

     q = dtype_data%levels(k)%q%values
     q0p = dtype_data%levels(k)%q0p%values

     DO i = 1,m
        DO j = 1,n
          uu = 0.5D0*( q(u(i+1,j))+q(u(i,j)) )/del
          vv = 0.5D0*( q(v(i,j+1))+q(v(i,j)) )/del
          u0p = 0.5D0*( q0p(u(i+1,j))+q0p(u(i,j)) )/del
          v0p = 0.5D0*( q0p(v(i,j+1))+q0p(v(i,j)) )/del

          dtype_data%levels(k)%pressure_dynamic%values(i,j) = &
                0.5D0 * ( (uu+u0p)**2.D0+(vv+v0p)**2.D0 )
          !dtype_data%levels(k)%pressure_dynamic%values(i,j) = &
          !      0.5D0 * ( u0p**2.D0+v0p**2.D0 )

        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE ibpm_pressure_set_dynamic_pressure

!************************************************************************************!
SUBROUTINE getbc_pressurebc_coarsest( dtype_data, pressure_bc )
  IMPLICIT NONE
  TYPE(data_t), INTENT(in) :: dtype_data
  REAL(KIND(0.D0)) :: pressure_bc(2*m+2*n)


  pressure_bc(bottom_phi+1:bottom_phi+m) = &
                    dtype_data%levels(mgridlev)%pressure_dynamic%values(1:m,1)
  pressure_bc(top_phi+1:top_phi+m)       = &
                    dtype_data%levels(mgridlev)%pressure_dynamic%values(1:m,n)
  pressure_bc(left_phi+1:left_phi+n)     = &
                    dtype_data%levels(mgridlev)%pressure_dynamic%values(1,1:n)
  pressure_bc(right_phi+1:right_phi+n)   = &
                    dtype_data%levels(mgridlev)%pressure_dynamic%values(m,1:n)

  ! get boundary condition
  !k           = mgridlev
  !pressure_bc = 0.D0
  !! TOP and BOTTOM
  !DO i=1,m
  !   j = 1
  !   pressure_bc(bottom_phi + i) = &
  !            dtype_data%levels(mgridlev)%pressure_dynamic%values(i,j)
  !   j = n
  !   pressure_bc(top_phi + i)    = &
  !            dtype_data%levels(mgridlev)%pressure_dynamic%values(i,j)
  !END DO
  !! LEFT and RIGHT
  !DO j=1,n
  !   i = 1
  !   pressure_bc(left_phi + j)   = &
  !            dtype_data%levels(mgridlev)%pressure_dynamic%values(i,j)
  !   i = m
  !   pressure_bc(right_phi + j)  = &
  !            dtype_data%levels(mgridlev)%pressure_dynamic%values(i,j)
  !END DO

END SUBROUTINE getbc_pressurebc_coarsest

!************************************************************************************!

SUBROUTINE ibpm_pressure_write_pressure( dtype_data )
  IMPLICIT NONE
  TYPE(data_t), INTENT(in) :: dtype_data
  INTEGER :: k
  CHARACTER(7) :: charit

  WRITE(charit,"(I7.7)") dtype_data%it
  OPEN(unit=100,file="output/snapshots/pressure"//charit//".var",form="unformatted",status="unknown")
  DO k = 1, mgridlev
    WRITE(100) dtype_data%levels(k)%pressure%values
    WRITE(100) dtype_data%levels(k)%pressure_dynamic%values
  END DO
  CLOSE(100)

END SUBROUTINE ibpm_pressure_write_pressure

!************************************************************************************!
SUBROUTINE ibpm_pressure_read_pressure( dtype_data )
  IMPLICIT NONE
  TYPE(data_t), INTENT(inout) :: dtype_data
  INTEGER :: k
  CHARACTER(7) :: charit

  WRITE(charit,"(I7.7)") dtype_data%it
  OPEN(unit=100,file="output/snapshots/pressure"//charit//".var",form="unformatted",status="unknown")
  DO k = 1, mgridlev
    READ(100) dtype_data%levels(k)%pressure%values
    READ(100) dtype_data%levels(k)%pressure_dynamic%values
  END DO
  CLOSE(100)

  END SUBROUTINE ibpm_pressure_read_pressure

!************************************************************************************!

END MODULE ibpm_pressure
