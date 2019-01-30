MODULE parameters

  IMPLICIT NONE
  ! *****************************************************************
  ! general parameters don't need to input in the input file.
  ! they will be determined automatically

  INTEGER :: n_actuator      ! number of actuators (value assigned in grid.f90)
  INTEGER :: n_body          ! number of bodies    (value assigned in grid.f90)
  INTEGER :: n_control       ! number of controls  (value assigned in control.f90)

  ! *****************************************************************
  ! general parameters needed to be input in the input file
  INTEGER :: istart            ! initial time index
  INTEGER :: istop             ! last time index
  INTEGER :: isave             ! save a restart for fwd sim every isave steps
  INTEGER :: nhoriz            ! time steps for optimization horizon
  REAL(KIND(0.0D0)) :: dt      ! time step

  INTEGER :: m                   ! cells in x
  INTEGER :: n                   ! cells in y
  REAL(KIND(0.D0)) :: len        ! length scale for grid
  REAL(KIND(0.D0)) :: offsetx    ! offset for grid in x
  REAL(KIND(0.D0)) :: offsety    ! offset for grid in y
  INTEGER :: mgridlev            ! how many layer of grid

  REAL(KIND(0.0D0)) :: Re        ! Reynolds number
  INTEGER :: ncr                 ! number of centers of rotation

  REAL(KIND(0.D0)) :: pi = 4.d0*atan(1.0d0) ! a useful constant
  INTEGER :: it_advance_show

  ! *****************************************************************
  ! flags don't need to input in the input file.
  ! they will be determined automatically
  LOGICAL :: alloc_ibd         ! whether to allocate ib data (determined in grid.f90)
  LOGICAL :: cfl_flag          ! whether cfl is beyond critical value (determined in write_responses.f90)

  ! *****************************************************************
  ! flags needed to be input in the input file
  LOGICAL :: compute_pressure  ! whether to output pressure
  LOGICAL :: output_sf         ! output stream function ?

  ! flags don't need to input in the input file.
  ! they will be determined automatically
  LOGICAL :: labframe          ! Evaluate everything in the lab frame?? (determined in control.f90)

  !----------------------------------------------------------
  ! NACA0018 parameters
  !----------------------------------------------------------
  REAL(KIND(0.D0)) :: naca_area =  0.123315D0
  REAL(KIND(0.D0)) :: naca_moi_center_chard = 0.007818002863031893D0
  REAL(KIND(0.D0)) :: vawt_rad = 1.5D0
  REAL(KIND(0.D0)) :: naca_moi_vawt_center

  ! angle of attack
  REAL(KIND(0.0D0)) :: aoa  
CONTAINS

! *****************************************************************************************

  SUBROUTINE parameters_input
  IMPLICIT NONE

    LOGICAL :: readinput
    INTEGER :: i

    NAMELIST /read_parameters/  istart, istop, isave, dt,                &
                                m,n,len,offsetx,offsety,mgridlev,        &
                                Re,ncr,output_sf, compute_pressure, aoa
    ! read input
    readinput = .TRUE.
    INQUIRE(file='input/ib.inp',exist=readinput)
    IF (readinput) THEN
      OPEN(unit=3,file='input/ib.inp',form='formatted',status='old')
      READ(unit=3,nml=read_parameters)
      CLOSE(3)
      WRITE(*,*) '=> Read input file'
    ELSE
      STOP 'ERROR: cannot find input file'
    END IF

    nhoriz = istop - istart + 1

    ! ************  Sanity check begins *********************
    IF (istop .lt. istart) THEN
      STOP 'ERROR: istop cannot be less than istart'
    END IF
    ! ************  Sanity check ends *********************

    ! ************ printing information *********************
    IF (Ncr.gt.2) THEN
        WRITE(*,*) '=> Notice : Using more than 2 centers of rotation. Procede with causion'
    END IF

    ! ************ Assigning initial parameters *********************
    ! reset cfl_flag
    CFL_FLAG = .FALSE.

    ! compute moment of inertia w.r.t vawt center
    naca_moi_vawt_center = naca_moi_center_chard + &
                           naca_area*(vawt_rad ** 2.D0)
    !WRITE(*,*) 'Density input = ',  density
    !WRITE(*,*) 'NACA 0018 has an area = ',  naca_area, ', '
    !WRITE(*,*) 'NACA 0018 has a moment of inertia = ',  naca_moi_center_chard, ' w.r.t its center chord '
    !WRITE(*,*) 'With the radius of turbine =', vawt_rad, &
    !           ', the moment of inertia of NACA0018 w.r.t. VAWT rotation center =', naca_moi_vawt_center

    aoa = aoa*pi/180.d0

  END SUBROUTINE parameters_input

  ! *****************************************************************************************

END MODULE parameters
