MODULE variables

  USE parameters
  USE dts
  USE grid
  USE user
  USE ibpm_pressure
  IMPLICIT NONE

  ! fwd and bwd flow field variables
  TYPE(data_t) :: var
  ! variables for Cholesky (stationary body)
  TYPE(chol_data_t) :: chol_data

CONTAINS

! *****************************************************************************************

  SUBROUTINE variables_setup_variables
    USE myfft
    IMPLICIT NONE
    INTEGER :: s, dit
    LOGICAL :: exist_chol, readic

    ! Allocate variables and set constants
    CALL myfft_setup_fft

    WRITE(*,*) '=> Setup multidomain levels and global flow variables'
    CALL dts_alloc_data_t(dtype = var, nb = nb, alloc_ibd = alloc_ibd )

    var%it = istart

    ! we initialize xb to the positions read from the files.  If this is a restart, this will be overwritten later
    WRITE(*,*) '=> Collecting all bodies'
    If ( alloc_ibd ) THEN
      ! when bodies are stationary, use cholesky docomp
      CALL dts_alloc_chol_data_t(dtype = chol_data, nf = nf)
      CALL grid_collect_bodies(xb = var%ibd%xb, code = var%codeb%code)
    END IF


    IF (istart.eq.0) THEN
        ! *********************** istart = 0 (begins) ************************
        CALL variables_initial_condition(dtype_data = var) ! only initial condition has read_type different
        CALL variables_write_variables(it = istart, dtype_data = var)
        ! *********************** istart = 0 (ends) ************************
    ELSE
        ! *********************** istart > 0 (begins) ************************
        ! readin restart file
        WRITE(*,*) '=> Reading flow variables'
        CALL variables_read_variables( it = istart, dtype_data = var)
        WRITE(*,*) '=> Flow variables read'
        ! *********************** istart > 0 (ends) ************************
    END IF


    IF (alloc_ibd) THEN
      WRITE(*,*) '=> Setup regularization of initial geometry'
      CALL grid_setup_reg(xb = var%ibd%xb)
      WRITE(*,*) '=> Regularization setup'

      ! if we start impulsively, and the bodies are stationary,
      ! do Cholesky decomposition
      CALL variables_preprocess
    END IF

  END SUBROUTINE variables_setup_variables
! ***********************************************
SUBROUTINE variables_destroy_variables
 USE myfft
 IMPLICIT NONE

 CALL myfft_destroy_fft
 IF ( alloc_ibd ) THEN
    CALL dts_dealloc_chol_data_t(dtype = chol_data)
 END IF
 CALL dts_dealloc_data_t(dtype = var, alloc_ibd = alloc_ibd )

END SUBROUTINE variables_destroy_variables

! *****************************************************************************************
! Reads or computes the Cholesky matrix for stationary bodies, originating from a
! Cholesky factorization of EC(C^t C)^-1 C^t E^t
SUBROUTINE variables_preprocess
  USE mmisc

  IMPLICIT NONE
  REAL(KIND(0.D0)) :: z(Nf)   ! Temporary immersed body force vector
  INTEGER :: i
  LOGICAL :: exist_chol

  WRITE(*,*) '=> Performing preprocessing ...'
  exist_chol =.FALSE.
  ! INQUIRE(file="output/ib.chd",exist = exist_chol)  ! ignore prev chol
  IF ( exist_chol )THEN
    ! **************** the Cholesky matrix exists (begins) **************
     WRITE(*,*) '=> Cholesky matrix already computed. Checking if we can reuse this matrix'
     CALL variables_read_cholesky
     IF (chol_data%recompute) THEN  ! check if we need to recompute it
       WRITE(*,*) '=> Recomputing body matrix for stationary geometry...please be patient'
       DO i=1,nf   ! Build matrix one column at a time
         z = 0.d0
         z(i) = 1.d0
         chol_data%cholmat(1:Nf,i) = a_times( z )
       END DO
       CALL mmisc_choldc(dtype = chol_data)
       CALL variables_write_cholesky  ! Save results
       chol_data%recompute =.FALSE.
     END IF
     ! **************** the Cholesky matrix exists (ends) **************
  ELSE
    ! ************* the Cholesky matrix does not exist (begins) ***********
    WRITE(*,*) '=> Precomputing body matrix for stationary geometry...please be patient'
    DO i=1,nf   ! Build matrix one column at a time
      z = 0.d0
      z(i) = 1.d0
      chol_data%cholmat(1:Nf,i) = a_times( z )
    END DO
    WRITE(*,*) '=> Performing Cholesky decomposition...please be patient'
    CALL mmisc_choldc(dtype = chol_data)
    WRITE(*,*) '=> Cholesky matrix computed. Saving the Cholesky matrix'
    CALL variables_write_cholesky  ! Save results
    ! ************* the Cholesky matrix does not exist (ends) **************
  END IF
  WRITE(*,*) '=> Done preprocessing ...'

END SUBROUTINE variables_preprocess

! *****************************************************************************************

  SUBROUTINE variables_write_cholesky
    IMPLICIT NONE
    OPEN(unit=100,file="output/ib.chd",form="unformatted",status="unknown")
    WRITE(100) m, n, mgridlev, nb
    WRITE(100) dt, Re, len
    WRITE(100) chol_data%cholmat, chol_data%cholvec
    CLOSE(100)
  END SUBROUTINE variables_write_cholesky
! ********************************************
  SUBROUTINE variables_read_cholesky
    IMPLICIT NONE
    INTEGER :: m_in, n_in, mgridlev_in, nb_in
    REAL(KIND(0.D0)) :: dt_in, Re_in, len_in

    OPEN(unit=100,file="output/ib.chd",form="unformatted",status="unknown")
    READ(100) m_in, n_in, mgridlev_in, nb_in
    READ(100) dt_in, Re_in, len_in

    ! there are 7 parameters that could change the Cholesky matrix
    ! default setting of chol_data%recompute is F
    IF ( (m_in.ne.m).or.(n_in.ne.n).or.(mgridlev_in.ne.mgridlev).or. &
         (nb_in.ne.nb).or.(dt_in.ne.dt).or.(Re_in.ne.Re).or.(len_in.ne.len) ) THEN
      chol_data%recompute =.TRUE.
      WRITE(*,*) '=> Parameters of the Cholesky matrix has changed.'
      CLOSE(100, status='delete')  ! deleting the old file
    ELSE
      WRITE(*,*) '=> Reusing the Cholesky matrix ...'
      READ(100) chol_data%cholmat, chol_data%cholvec
      CLOSE(100)
    END IF

  END SUBROUTINE variables_read_cholesky

! *****************************************************************************************

  SUBROUTINE variables_write_variables(it, dtype_data)
    IMPLICIT NONE
    INTEGER, INTENT(in) :: it
    TYPE(data_t), INTENT(in) :: dtype_data
    CHARACTER(7) :: charit
    INTEGER :: i, k, io_sf

    write(*,*) '=> Writing variables at it=',it
    write(charit,"(I7.7)") it
    OPEN( unit=100, &
          file = "output/snapshots/ib"//charit//".var", &
          form="unformatted", &
          status="unknown")


    IF (output_sf) THEN
      io_sf = 1
    ELSE
      io_sf = 0
    END IF
    ! 1st line : integer parameters
    WRITE(100) dtype_data%it,m,n,mgridlev,nb,ncr,io_sf
    ! 2nd line : real parameters
    WRITE(100) re,dt,len,offsetx,offsety
    DO k = 1, mgridlev
      ! write q, q0p, nl_bf_old on each level
      WRITE(100) dtype_data%levels(k)%q%values, &
                 dtype_data%levels(k)%q0p%values, &
                 dtype_data%levels(k)%nl_bf_old%values
      ! write q0r on each level
      DO i = 1, ncr
        WRITE(100) dtype_data%levels(k)%q0r(i)%values
      END DO
      ! write circulation on each level
      WRITE(100) dtype_data%levels(k)%omega%values
      ! write streamfunction on each level if output_st = T
      IF (output_sf) WRITE(100) dtype_data%levels(k)%stfn%values
    END DO
    ! write ibdata
    IF (alloc_ibd) THEN
      WRITE(100) dtype_data%ibd%xb, dtype_data%ibd%vb, dtype_data%ibd%fb
      WRITE(100) dtype_data%codeb%code
    END IF
    ! write translational cm motion info
    WRITE(100) dtype_data%cm_motion%trans_cm%xcm, &
               dtype_data%cm_motion%trans_cm%vcm, &
               dtype_data%cm_motion%theta_g2l, &
               dtype_data%cm_motion%omega_g2l
    ! write rotational cm motion info
    DO i = 1, ncr
      WRITE(100) dtype_data%cm_motion%rot_cm(i)%rox, &
                 dtype_data%cm_motion%rot_cm(i)%roy, &
                 dtype_data%cm_motion%rot_cm(i)%theta, &
                 dtype_data%cm_motion%rot_cm(i)%omega
    END DO
    CLOSE(100)

  END SUBROUTINE variables_write_variables
! **********************************************
  SUBROUTINE variables_read_variables(it, dtype_data)
    IMPLICIT NONE
    INTEGER, INTENT(in) :: it
    TYPE(data_t), INTENT(inout) :: dtype_data
    TYPE(data_t) :: dtype_data_dummy ! use for dummy state
    TYPE(ib_data_t) :: ibd_old
    TYPE(ib_code_t) :: codeb_old
    CHARACTER(7) :: charit
    INTEGER :: it_temp, m_in, n_in, mgridlev_in, nb_in, ncr_in
    INTEGER :: i, k, io_sf
    LOGICAL :: output_sf_in
    REAL(KIND(0.D0)) :: re_in,dt_in,len_in,offsetx_in,offsety_in
    REAL(KIND(0.D0)) :: temp_flux(Nq), temp

    LOGICAL :: readic

      IF (it.eq.0) THEN
        WRITE(*,*) ' => Inputing initial condition file initial.var'
        OPEN(unit=100,file="input/initial.var",form="unformatted",status="unknown")
      ELSE
       write(*,*) '=> Reading variables at it=',it
       write(charit,"(I7.7)") it
       OPEN( unit=100, &
             file = "output/snapshots/ib"//charit//".var", &
             form="unformatted", &
             status="unknown")
      END IF


    ! 1st line : integer parameters
    READ(100) it_temp, m_in, n_in, mgridlev_in, nb_in, ncr_in, io_sf
    IF ( it.eq.0 ) THEN
      dtype_data%it = 0
    ELSE
      dtype_data%it = it
    END IF
    IF (io_sf.eq.1) THEN
      output_sf_in = .TRUE.
    ELSE
      output_sf_in = .FALSE.
    END IF
    ! Perform a few checks to see restart will work
    IF ((m_in.ne.m).or.(n_in.ne.n).or.(mgridlev.ne.mgridlev_in)) THEN
       WRITE(*,*) 'Initial condition file used different numbers of grid points'
       STOP 'Unable to continue'
    END IF

    ! 2nd line : real parameters
    READ(100) re_in,dt_in,len_in,offsetx_in,offsety_in
    !READ(100) output_sf_in
    ! Perform a few checks to see restart will work
    IF (re_in.ne.re) THEN
       WRITE(*,*) 'Reynolds number has changed from initial condition'
       STOP 'Unable to continue'
    END IF
    IF ((len_in.ne.len).or.(offsetx.ne.offsetx_in).or.(offsety.ne.offsety_in)) THEN
       WRITE(*,*) 'Physical dimensions of grid have changed'
       WRITE(*,*) 'Disregarding old values; proceed with caution!'
    END IF

    DO k = 1, mgridlev
      ! read q, q0p, nl_bf_old on each level
      READ(100) dtype_data%levels(k)%q%values, &
                  dtype_data%levels(k)%q0p%values, &
                  dtype_data%levels(k)%nl_bf_old%values
      dtype_data%levels(k)%q0%values = dtype_data%levels(k)%q0p%values
      ! write q0r on each level
      IF (ncr_in.lt.ncr) THEN
        WRITE(*,*) '=> Warning : number of centers of rotation in snapshot is less than in input'
        WRITE(*,*) '=> Motion that has icr larger than ncr snapshot will not be imported'
        DO i = 1, ncr_in
          READ(100) dtype_data%levels(k)%q0r(i)%values
          dtype_data%levels(k)%q0%values = dtype_data%levels(k)%q0%values + &
                                           dtype_data%levels(k)%q0r(i)%values
        END DO
      ELSEIF (ncr_in.gt.ncr) THEN
        WRITE(*,*) '=> Warning : number of centers of rotation in snapshot is more than in input'
        WRITE(*,*) '=> Motion that has icr larger than ncr input will not be imported'
        DO i = 1, ncr
          READ(100) dtype_data%levels(k)%q0r(i)%values
          dtype_data%levels(k)%q0%values = dtype_data%levels(k)%q0%values + &
                                           dtype_data%levels(k)%q0r(i)%values
        END DO
        DO i = ncr+1, ncr_in
          READ(100) temp_flux
        END DO
      ELSE
        DO i = 1, ncr
          READ(100) dtype_data%levels(k)%q0r(i)%values
          dtype_data%levels(k)%q0%values = dtype_data%levels(k)%q0%values + &
                                           dtype_data%levels(k)%q0r(i)%values
        END DO
      END IF
      ! read circulation on each level
      READ(100) dtype_data%levels(k)%omega%values
      ! read streamfunction on each level if output_st_in = T
      IF (output_sf_in) READ(100) dtype_data%levels(k)%stfn%values
    END DO

    ! read ibdata
    IF (nb_in.ne.0) THEN
      CALL dts_alloc_ib_code_t( dtype = codeb_old, nb = nb_in )
      CALL dts_alloc_ib_data_t( dtype = ibd_old, nf = 2*nb_in )
      READ(100) ibd_old%xb, ibd_old%vb, ibd_old%fb
      READ(100) codeb_old%code
      IF (nb.ne.nb_in) THEN
        WRITE(*,*) 'Geometry has changed from initial condition.'
        WRITE(*,*) 'Disregarding old geometry and forces'
      ELSE
        IF ( ANY(ibd_old%xb .ne. dtype_data%ibd%xb) .or. &
              ANY(codeb_old%code .ne.  dtype_data%codeb%code) ) THEN
          WRITE(*,*) 'Geometry has changed from initial condition.'
          WRITE(*,*) 'Disregarding old geometry and forces'
        ELSE
          dtype_data%ibd%vb = ibd_old%vb
          dtype_data%ibd%fb = ibd_old%fb
        END IF
      END IF
      CALL dts_dealloc_ib_code_t( dtype = codeb_old )
      CALL dts_dealloc_ib_data_t( dtype = ibd_old )
    END IF

    ! read translational cm motion info
    READ(100) dtype_data%cm_motion%trans_cm%xcm, &
              dtype_data%cm_motion%trans_cm%vcm, &
              dtype_data%cm_motion%theta_g2l, &
              dtype_data%cm_motion%omega_g2l
    ! read translational cm motion info
    IF (ncr_in.lt.ncr) THEN
      DO i = 1, ncr_in
        READ(100) dtype_data%cm_motion%rot_cm(i)%rox, &
                  dtype_data%cm_motion%rot_cm(i)%roy, &
                  dtype_data%cm_motion%rot_cm(i)%theta, &
                  dtype_data%cm_motion%rot_cm(i)%omega
      END DO
    ELSEIF (ncr_in.gt.ncr) THEN
      WRITE(*,*) '=> Warning : number of centers of rotation in snapshot is more than in input'
      WRITE(*,*) '=> Motion that has icr larger than ncr input will not be imported'
      DO i = 1, ncr
        READ(100) dtype_data%cm_motion%rot_cm(i)%rox, &
                  dtype_data%cm_motion%rot_cm(i)%roy, &
                  dtype_data%cm_motion%rot_cm(i)%theta, &
                  dtype_data%cm_motion%rot_cm(i)%omega
      END DO
      DO i = ncr+1, ncr_in
        READ(100) temp, temp, temp, temp
      END DO
    ELSE
      DO i = 1, Ncr
        READ(100) dtype_data%cm_motion%rot_cm(i)%rox, &
                  dtype_data%cm_motion%rot_cm(i)%roy, &
                  dtype_data%cm_motion%rot_cm(i)%theta, &
                  dtype_data%cm_motion%rot_cm(i)%omega
      END DO
    END IF

    CLOSE(100)

  END SUBROUTINE variables_read_variables

  !***********************************************************************

  SUBROUTINE variables_initial_condition(dtype_data)
    USE mmisc
    USE user
    IMPLICIT NONE
    TYPE(data_t), INTENT(inout) :: dtype_data
    INTEGER :: i, j, k
    REAL(KIND(0.D0)) :: fac, xx, yy
    LOGICAL :: readic

    ! check whether there is an initial file
    INQUIRE(file="input/initial.var",exist=readic)

    IF (readic) THEN ! read initial file if there is one
      ! ****************** intial file exits (begins) *******************
        WRITE(*,*) '=> Find initial file for fwd sim .... reading'
        CALL variables_read_variables( it = istart, dtype_data = dtype_data )
      WRITE(*,*) '=> Find initial file .... done reading'
      ! ****************** intial file exits (ends) *******************
    ELSE  ! no initial file
      ! ****************** no intial file (begins) *******************
      WRITE(*,*) '=> Cannot find initial file .... initializing the flow'
      ! ********** User-defined initial vorticity field (begins) ***************
      ! initial vorticity field
      DO k = 1, mgridlev
        fac = delta*2.d0**(REAL(k)-1.D0)  ! cell face length on each grid
        DO j = 2, n
          DO i = 2, m
            ! find coord where x-comp velocity lives
             xx =  (REAL(i)-1.D0-REAL(m)/2.D0) * fac + &
                                        REAL(m)/2.D0*delta - offsetx
             yy =  (REAL(j)-1.D0-REAL(n)/2.D0) * fac + &
                                        REAL(n)/2.D0*delta - offsety
             dtype_data%levels(k)%omega%values(i,j) = ( fac**2.D0 ) * &
                                                initial_vorticity( xx, yy )
          END DO
       END DO
     END DO

    CALL vort2flux( dtype_data = dtype_data, nlev = mgridlev )

    ! ********** User-defined initial vorticity field (ends) ***************
      WRITE(*,*) '=> Cannot find initial file .... flow initialized'
      ! ****************** no intial file (ends) *******************
    END IF

  END SUBROUTINE variables_initial_condition

! *****************************************************************************************

END MODULE variables
