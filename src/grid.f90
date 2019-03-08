MODULE grid

  ! bodies/actuators input example :
  ! for each body.xxx.inp (xxx is the number) in input directory :
  ! bdy%npts
  !     bdy%x(1)    ,    bdy%y(1)
  !     bdy%x(2)    ,    bdy%y(2)
  !                 .
  !                 .
  !                 .
  ! bdy%x(bdy%npts) , bdy%y(bdy%npts)
  !
  ! for each actuator.xxx.inp (xxx is the number) in input directory :
  ! act%npts
  !     act%x(1)    ,    act%y(1)
  !     act%x(2)    ,    act%y(2)
  !                 .
  !                 .
  !                 .
  ! act%x(act%npts) , act%y(act%npts)

  USE parameters
  USE dts

  IMPLICIT NONE

  ! Variables for immersed boundary geometry and its motion (if any)
  INTEGER :: support = 3 ! support for smearing delta functions
  REAL(KIND(0.D0)) :: delta ! near field grid spacing

  ! index of each ib points relative to grid
  INTEGER, ALLOCATABLE          :: indexx(:)
  ! arrays for smearing coefficients
  REAL(KIND(0.D0)), ALLOCATABLE :: weight(:,:)

  ! Numbers of faces, ib, and ibfs.
  INTEGER :: nq, nb, nf

  ! Integer pointer for field variables
  INTEGER, ALLOCATABLE :: f(:,:) ! f(i,k) gives kth-comp. of immersed body force (or position or vel) at ith point on body
  INTEGER, ALLOCATABLE :: u(:,:) ! u(i,j) gives point in the q array where u at face of cell i,j lives
  INTEGER, ALLOCATABLE :: v(:,:) ! v(i,j) gives point in the q array where v at face of cell i,j lives

  ! for bcs
  INTEGER :: top,bottom,left,right
  INTEGER :: top_phi,bottom_phi,left_phi,right_phi

  !  for bodies
  INTEGER, PARAMETER :: maxbodies = 999 ! a large number
  TYPE(body_t) :: bdy(maxbodies)
  TYPE(actuator_t) :: act(maxbodies)

CONTAINS

! *****************************************************************************************
  SUBROUTINE grid_setup_eulerian_grid(n_foil)

    IMPLICIT NONE
    INTEGER :: i,j,ns,next
    INTEGER, INTENT(in) :: n_foil

    ! firt order of business is to setup the grid
    delta = len/REAL(m)

    ! a few things to setup for multigrid laplace solution
    ! check that grid is divisible by 4
    IF ((MOD(m,4).ne.0).or.(MOD(n,4).ne.0)) THEN
       STOP 'grid must be divisible by 4'
    END IF
    ! for multigrid boundary conditions
    left = 0
    right = n+1
    bottom = 2*(n+1)
    top = 2*(n+1) + m+1

    left_phi = 0
    right_phi = n
    bottom_phi = 2*n
    top_phi = 2*n+m

    ! indexing for streamfunction
    ns = (m-1)*(n-1)
    ! indexing for velocities (flux components)
    nq = (m+1)*n + (n+1)*m

    ALLOCATE( u(1:m+1,1:n), v(1:m,1:n+1) )

    next = 0
    DO j=1,n
       DO i=1,m+1
          next = next+1
          u(i,j) = next
       END DO
    END DO
    DO j=1,n+1
       DO i=1,m
          next = next+1
          v(i,j) = next
       END DO
    END DO
    IF (next.ne.nq) STOP "ERROR: error in setup_parms - u"

    ! now set up a body
    CALL grid_setup_eulerian_grid_setup_geometry(n_foil = n_foil)

  END SUBROUTINE grid_setup_eulerian_grid

! *****************************************************************************************
  SUBROUTINE grid_setup_eulerian_grid_setup_geometry(n_foil)

    IMPLICIT NONE
    LOGICAL :: readinput
    INTEGER :: i, j, next
    INTEGER, INTENT(in) :: n_foil
    CHARACTER(3) :: file_num

    ! look for bodies in input directory
    readinput = .TRUE.
    n_body = 1
    ! DO WHILE (readinput)
       WRITE(file_num,"(I3.3)") n_foil
       INQUIRE(file="input/ib"//file_num//".inp",exist=readinput)
       IF (readinput) THEN
          ! n_body=n_body+1
          OPEN(unit=8,file="input/ib"//file_num//".inp",form='formatted',status='old')
          READ(8,*) bdy(n_body)%npts
          !READ(8,*) bdy(n_body)%moving
          ALLOCATE( bdy(n_body)%x(bdy(n_body)%npts), bdy(n_body)%y(bdy(n_body)%npts) )
          DO i = 1, bdy(n_body)%npts
             READ(8,*) bdy(n_body)%x(i), bdy(n_body)%y(i)
          END DO
          CLOSE(8)
       ELSE
          STOP 'foil file does not exist'
       END IF
    ! END DO

    ! look for actuators in input directory
    readinput = .TRUE.
    n_actuator = 0
    DO WHILE (readinput)
       WRITE(file_num,"(I3.3)") n_actuator+1
       INQUIRE(file="input/actuator."//file_num//".inp",exist=readinput)
       IF (readinput) THEN
          n_actuator=n_actuator+1
          OPEN(unit=8,file="input/actuator."//file_num//".inp",form='formatted',status='old')
          READ(8,*) act(n_actuator)%npts
          ALLOCATE( act(n_actuator)%x(act(n_actuator)%npts), act(n_actuator)%y(act(n_actuator)%npts) )
          ALLOCATE( act(n_actuator)%s(act(n_actuator)%npts) )
          DO i=1,act(n_actuator)%npts
             READ(8,*) act(n_actuator)%x(i), act(n_actuator)%y(i)
          END DO
          CLOSE(8)
       END IF
    END DO

    IF ( (n_body.gt.1).and.(Ncr.gt.1) ) THEN
      STOP 'Current code does not supoort multiple bodies with multiple center of rotations'
    END IF

    WRITE(*,*) '=> Read all bodies and actuators'
    ! accumulate all bodies and actuators into global vector xb
    nb = 0
    alloc_ibd = .TRUE.
    !stationary = .TRUE.
    DO i=1,n_body
       WRITE(*,*) '=> Body no.',i,'has',bdy(i)%npts,'points.'
       nb = nb + bdy(i)%npts
    END DO
    DO i=1,n_actuator
       WRITE(*,*) '=> Act. no.',i,'has',act(i)%npts,'points.' 
       nb = nb + act(i)%npts
    END DO
    WRITE(*,*) '=> There are',nb,'lagrangian points'
    nf = 2*nb     ! number of immersed body forces

    IF ( nb .eq. 0 ) THEN
      WRITE(*,*) '===> Note: there is no body or actuator <==='
      alloc_ibd = .FALSE.
    ELSE
      !indexing for vectors on the immersed boundary (force/position/velocity)
      ALLOCATE( f(1:nb,1:2)  )
      next = 0
      DO i=1,nb
        next = next + 1
        f(i,1) = next
      END DO
      DO i = 1,nb
        next = next + 1
        f(i,2) = next
      END DO
      IF (next.ne.nf) STOP "ERROR: error in setup_parms - f <==="

    END IF
   !
 END SUBROUTINE grid_setup_eulerian_grid_setup_geometry

! *****************************************************************************************
  SUBROUTINE grid_collect_bodies(xb,code)

    IMPLICIT NONE
    INTEGER :: code(:)
    REAL(KIND(0.D0)) :: xb(:)
    INTEGER :: i,i_bdy,i_act,next

    next = 1
    DO i_bdy=1,n_body
       bdy(i_bdy)%ipos = next
       DO i=1,bdy(i_bdy)%npts
          xb(f(next,1)) = bdy(i_bdy)%x(i)
          xb(f(next,2)) = bdy(i_bdy)%y(i)
          code(next) = i_bdy
          next = next+1
       END DO
    END DO
    DO i_act=1,n_actuator
       act(i_act)%ipos = next
       DO i=1,act(i_act)%npts
          xb(f(next,1)) = act(i_act)%x(i)
          xb(f(next,2)) = act(i_act)%y(i)
          code(next) = -i_act
          next = next+1
       END DO
    END DO

  END SUBROUTINE grid_collect_bodies

  ! *****************************************************************************************
    SUBROUTINE grid_destroy_grid
      IMPLICIT NONE
      INTEGER :: i

      DEALLOCATE( f, u, v)

      IF(ALLOCATED(indexx)) THEN
        DEALLOCATE(indexx)
      END IF

      IF(ALLOCATED(weight)) THEN
        DEALLOCATE(weight)
      END IF

      DO i=1,n_body
        DEALLOCATE( bdy(i)%x, bdy(i)%y )
      ENDDO

      DO i=1,n_actuator
        DEALLOCATE( act(i)%x, act(i)%y, act(i)%s )
      ENDDO

    END SUBROUTINE grid_destroy_grid

  ! *****************************************************************************************
  SUBROUTINE grid_setup_reg(xb)

    IMPLICIT NONE
    INTEGER          :: i, k, l, next
    REAL(KIND(0.D0)) :: x, y, d2
    REAL(KIND(0.D0)) :: xb(:)

    !    support2= INT(support/2) + MOD(INT(support),2)
    d2 = 0.5d0*delta

    IF(ALLOCATED(indexx)) THEN
      DEALLOCATE(indexx)
    END IF

    IF(ALLOCATED(weight)) THEN
      DEALLOCATE(weight)
    END IF

    ALLOCATE(weight((2*support+1)*(2*support+1),nf))
    ALLOCATE(indexx(nf))

    ! get index of body position relative to grid
    DO i = 1, nb
       indexx(i)    = INT((xb(i)+offsetx)/delta) ! body-x index
       indexx(i+nb) = INT((xb(i+nb)+offsety)/delta) ! body-y index
    END DO

      ! get regularized weight near ib points (u-vel points)
    DO i = 1,nb
       next = 0
       DO l = -support, support
          DO k = -support, support
            x = delta*(indexx(i)-1+k)-offsetx ! grid location x
            y = delta*(indexx(i+nb)-1+l)-offsety+d2 ! grid location y
            next = next+1
            weight(next,i)= delta * delta * &
            deltafnc(x,xb(i),delta) * &
            deltafnc(y,xb(i+nb),delta)
          END DO
       END DO
    END DO
    ! get regularized weight near ib points (v-vel points)
    DO i = 1, nb
       next = 0
       DO l = -support,support
         DO k = -support,support
            x = delta*(indexx(i)-1+k)-offsetx+d2 ! grid location x
            y = delta*(indexx(i+nb)-1+l)-offsety ! grid location y
            next = next+1
            weight(next,i+nb)= delta * delta * &
            deltafnc(x,xb(i),delta) * &
            deltafnc(y,xb(i+nb),delta)
          END DO
       END DO
    END DO

  END SUBROUTINE grid_setup_reg

  ! *****************************************************************************************
    FUNCTION deltafnc( r,r0, dr )
        IMPLICIT NONE

      REAL(KIND(0.D0)) :: r,r0,dr,deltafnc

      ! Roma, Peskin, & Berger (JCP 1999)
      IF (ABS(r-r0)<=(0.5D0*dr)) THEN
         deltafnc = (1.D0+SQRT(-3.D0*((r-r0)/dr)**2+1.D0))/(3.0D0*dr)
      ELSEIF (ABS(r-r0)<=(1.5D0*dr)) THEN
         deltafnc = (5.D0-3.D0*ABS((r-r0)/dr)-SQRT(-3.D0*(1.D0-ABS((r-r0)/dr))**2+1.D0))/(6.D0*dr)
      ELSE
         deltafnc = 0.D0
         !       write(*,*) " out of support - (deltafnc)"
      END IF

    END FUNCTION deltafnc

! *****************************************************************************************
  FUNCTION getbody( b, dir, bdyno, npts )
    IMPLICIT NONE
    INTEGER :: bdyno, dir, npts  ! body id, direction, number of points
    REAL(KIND(0.D0)) :: b(:)     ! ib data
    REAL(KIND(0.D0)) :: getbody(npts)

    getbody = b(f ( bdy(bdyno)%ipos:bdy(bdyno)%ipos+npts-1, dir ))

  END FUNCTION getbody

! *****************************************************************************************
  FUNCTION getact( b, dir, actno, npts )
    IMPLICIT NONE
    INTEGER :: actno, dir, npts  ! actuator id, direction, number of points
    REAL(KIND(0.D0)) :: b(:)     ! ib data
    REAL(KIND(0.D0)) :: getact(npts)

    getact = b(f ( act(actno)%ipos:act(actno)%ipos+npts-1, dir ))

  END FUNCTION getact

! *****************************************************************************************

END MODULE grid
