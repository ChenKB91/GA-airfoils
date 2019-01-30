MODULE mmisc_bc

  ! This module contains subroutines and functions relate to obtaining
  ! and applying boundary conditions on each grid level

  IMPLICIT NONE

CONTAINS

! *********************************************************************!
  SUBROUTINE mmisc_bc_getbc( r, rbc, fac)
    !***************************************************************!
    !*   given vorticity on a larger, coarser mesh, interpolate    *!
    !*   its values to the edge of a smaller, finer mesh           *!
    !*   Note: it is not used in current code. Instead, we use     *!
    !*         mmisc_bc_getbc_spline                               *!
    !***************************************************************!
    USE parameters
    USE grid
   IMPLICIT NONE

    REAL(KIND(0.D0)), DIMENSION(:,:) :: r
    REAL(KIND(0.D0)), DIMENSION(:) :: rbc
    REAL(KIND(0.D0)) :: fac
    INTEGER :: i,j

    ! get interpolated boundary conditions on finer grid

    DO i=0,m,2
       rbc(bottom + i+1) = r(m/4+i/2,n/4)
       rbc(top + i+1) = r(m/4+i/2,3*n/4)
    END DO
    DO i=1,m-1,2
       rbc(bottom +i+1)  = 0.5d0*( r(m/4+(i+1)/2,n/4) + r(m/4-1+(i+1)/2,n/4) )
       rbc(top + i+1) = 0.5d0*( r(m/4+(i+1)/2,3*n/4) + r(m/4-1+(i+1)/2,3*n/4) )
    END DO

    DO j=0,n,2
       rbc(left + j+1) = r(m/4, n/4+j/2)
       rbc(right + j+1) = r(3*m/4, n/4+j/2)
    END DO
    DO j=1,n-1,2
       rbc(left + j+1) = 0.5d0*( r(m/4, n/4+(j+1)/2) + r(m/4, n/4-1+(j+1)/2) )
       rbc(right + j+1) = 0.5d0*( r(3*m/4, n/4+(j+1)/2) + r(3*m/4, n/4-1+(j+1)/2) )
    END DO

    rbc = rbc*fac

  END SUBROUTINE mmisc_bc_getbc

  ! *********************************************************************!

    SUBROUTINE mmisc_bc_getbc_flux( r, rbc, fac)
      !***************************************************************!
      !*   given vorticity on a larger, coarser mesh, interpolate    *!
      !*   its values to the edge of a smaller, finer mesh           *!
      !*   Note: it is not used in current code. Instead, we use     *!
      !*         mmisc_bc_getbc_spline                               *!
      !***************************************************************!
      USE parameters
      USE grid
     IMPLICIT NONE

      REAL(KIND(0.D0)), DIMENSION(:) :: r
      REAL(KIND(0.D0)), DIMENSION(:) :: rbc
      REAL(KIND(0.D0)) :: fac
      INTEGER :: i,j

      ! get interpolated boundary conditions on finer grid

      DO i=0,m,2
         rbc(bottom + i+1) = 0.25D0*(        r(u(m/4+i/2,n/4  )) &
                                      + 3.D0*r(u(m/4+i/2,n/4-1))   )
         rbc(top + i+1) = 0.25D0*(   3.D0*r(u(m/4+i/2,3*n/4+1)) &
                                   +      r(u(m/4+i/2,3*n/4  )) )
      END DO
      DO i=1,m-1,2
        rbc(bottom + i+1) = 0.125D0*(        r(u(m/4+(i-1)/2,n/4  )) &
                                      + 3.D0*r(u(m/4+(i-1)/2,n/4-1)) &
                                      +      r(u(m/4+(i+1)/2,n/4  )) &
                                      + 3.D0*r(u(m/4+(i+1)/2,n/4-1))   )
        rbc(top + i+1) = 0.125D0*(   3.D0*r(u(m/4+(i-1)/2,3*n/4+1)) &
                                   +      r(u(m/4+(i-1)/2,3*n/4  )) &
                                   + 3.D0*r(u(m/4+(i+1)/2,3*n/4+1)) &
                                   +      r(u(m/4+(i+1)/2,3*n/4  ))  )
      END DO

      DO j=0,n,2
         rbc(left + j+1) = 0.25D0*(  3.D0*r(v(m/4-1, n/4+j/2)) &
                                   +      r(v(m/4  , n/4+j/2))  )
         rbc(right + j+1) = 0.25D0*(       r(v(3*m/4  , n/4+j/2)) &
                                    + 3.D0*r(v(3*m/4+1, n/4+j/2))  )
      END DO
      DO j=1,n-1,2
        rbc(left + j+1) = 0.125D0*(  3.D0*r(v(m/4-1, n/4+(j-1)/2)) &
                                   +      r(v(m/4  , n/4+(j-1)/2)) &
                                   + 3.D0*r(v(m/4-1, n/4+(j+1)/2)) &
                                   +      r(v(m/4  , n/4+(j+1)/2)) )
        rbc(right + j+1) = 0.125D0*(       r(v(3*m/4  , n/4+(j+1)/2)) &
                                    + 3.D0*r(v(3*m/4+1, n/4+(j+1)/2)) &
                                    +      r(v(3*m/4  , n/4+(j-1)/2)) &
                                    + 3.D0*r(v(3*m/4+1, n/4+(j-1)/2)) )
      END DO

      rbc = rbc*fac

    END SUBROUTINE mmisc_bc_getbc_flux

! ***********************************************************************!
  SUBROUTINE mmisc_bc_getbc_spline( r, rbc, fac)
    !***************************************************************!
    !*   given vorticity on a larger, coarser mesh, interpolate    *!
    !*   its values to the edge of a smaller, finer mesh           *!
    !***************************************************************!
    USE parameters
    USE grid
   IMPLICIT NONE

    REAL(KIND(0.D0)), DIMENSION(:,:) :: r
    REAL(KIND(0.D0)), DIMENSION(:) :: rbc
    REAL(KIND(0.D0)) :: fac
    INTEGER :: i,j

    REAL(KIND(0.D0)) :: xx(0:m/2),yy(0:n/2),rx(0:m/2),rx2(0:m/2),ry(0:n/2),ry2(0:n/2),xfine,yfine

   xx=0.d0
   yy=0.d0
   rx=0.d0
   rx2=0.d0
   ry=0.d0
   ry2=0.d0

   !location of the coarse grid points relative to the first point (counting in finer level grid spacing units)
   DO i=0,m/2
     xx(i)=2.d0*real(i)
   ENDDO

!bottom
   !corresponding values of the coarse grid array at fine grid boundary
   DO i=0,m/2
     rx(i)=r(m/4+i,n/4)
   ENDDO
   !get spline coefficients
   CALL mmisc_bc_spline(x = xx, y = rx, n = m/2+1, yp1 = 1.d31, ypn = 1.d31, y2 = rx2)
   !get interpolated value
   DO i=0,m
     xfine=real(i)
     CALL mmisc_bc_splint(xa = xx, ya = rx, y2a = rx2, n = m/2+1, x = xfine, y = rbc(bottom +i+1))
   ENDDO

!top
   !corresponding values of the coarse grid array at fine grid boundary
   DO i=0,m/2
     rx(i)=r(m/4+i,3*n/4)
   ENDDO
   !get spline coefficients
   CALL mmisc_bc_spline(x = xx, y = rx, n = m/2+1, yp1 = 1.d31, ypn = 1.d31, y2 = rx2)
   !get interpolated value
   DO i=0,m
     xfine=real(i)
     CALL mmisc_bc_splint(xa = xx, ya = rx, y2a = rx2, n = m/2+1, x = xfine, y = rbc(top +i+1))
   ENDDO

   !location of the coarse grid points relative to the first point (counting in finer level grid spacing units)
   DO j=0,n/2
     yy(j)=2.d0*real(j)
   ENDDO

!left
   !corresponding values of the coarse grid array at fine grid boundary
   DO j=0,n/2
     ry(j)=r(m/4,n/4+j)
   ENDDO
   !get spline coefficients
   CALL mmisc_bc_spline(x = yy, y = ry, n = n/2+1, yp1 = 1.d31, ypn = 1.d31, y2 = ry2)
   !get interpolated value
   DO j=0,n
     yfine=real(j)
     CALL mmisc_bc_splint(xa = yy, ya = ry, y2a = ry2, n = n/2+1, x = yfine, y = rbc(left +j+1))
   ENDDO

!right
   !corresponding values of the coarse grid array at fine grid boundary
   DO j=0,n/2
     ry(j)=r(3*m/4,n/4+j)
   ENDDO
   !get spline coefficients
   CALL mmisc_bc_spline(x = yy, y = ry, n = n/2+1, yp1 = 1.d31, ypn = 1.d31, y2 = ry2)
   !get interpolated value
   DO j=0,n
     yfine=real(j)
     CALL mmisc_bc_splint(xa = yy, ya = ry, y2a = ry2, n = n/2+1, x = yfine, y = rbc(right +j+1))
   ENDDO

    rbc = rbc*fac

  END SUBROUTINE mmisc_bc_getbc_spline

! **********************************************************************!
  SUBROUTINE mmisc_bc_applybc( r, rbc, fac)
    !***************************************************************!
    !*   given vorticity at edges of domain, rbc, (from larger,    *!
    !*   coraser mesh), add values to correct laplacian of         *!
    !*   vorticity  on the (smaller, finer) domain, r.             *!
    !***************************************************************!
    USE grid
    USE parameters
   IMPLICIT NONE

    REAL(KIND(0.D0)), DIMENSION(:,:) :: r
    REAL(KIND(0.D0)), DIMENSION(:) :: rbc
    REAL(KIND(0.D0)) :: fac
    INTEGER :: i,j

    ! add bc's from coarser grid

    DO i=1,m-1
       r(i,1) = r(i,1) + fac* rbc( bottom + i+1 )
       r(i,n-1) = r(i,n-1) + fac* rbc( top + i+1 )
    END DO
    DO j=1,n-1
       r(1,j) = r(1,j) + fac* rbc( left + j+1 )
       r(m-1,j) = r(m-1,j) + fac* rbc( right + j+1 )
    END DO

  END SUBROUTINE mmisc_bc_applybc

! *****************************************************************************************
  SUBROUTINE mmisc_bc_spline(x,y,n,yp1,ypn,y2)
  ! From Numerical Recipes in Fortran, Press et al.
  !Given arrays x(1:n) and y(1:n) containing a tabulated function, i.e., yi = f(xi), with x1< x2<...< xN,
  !and given values yp1 and ypn for the function at points 1 and n, respectively, this routine returns an array y2(1:n) of
  !length n which contains the second derivatives of the interpolating function at the tabulated
  !points xi. If yp1 and/or ypn are equal to 1x10^30 or larger, the routine is signaled to set
  !the corresponding boundary condition for a natural spline, with zero second derivative on that boundary

  ! Called once for a given set of data, need to call splint once for each desired
  ! spline interpolated point

    IMPLICIT NONE
    INTEGER n,NMAX
    REAL(KIND(0.D0)) :: yp1,ypn,x(n),y(n),y2(n)
    PARAMETER(NMAX=500) !Parameter: NMAX is the largest anticipated value of n.

    INTEGER i,k
    REAL p,qn,sig,un,u(NMAX)

    if (yp1.gt..99e30)then              !The lower boundary condition is set either to be
      y2(1)=0.                        !"natural"
      u(1)=0.
    else                            !or else to have a specified first derivative
      y2(1)=-0.5
      u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
    endif

    do i=2,n-1 !This is the decomposition loop of the tridiagonal algorithm. y2 and u are used for temporary storage of the decomposed factors.
      sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
      p=sig*y2(i-1)+2.
      y2(i)=(sig-1.)/p
      u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1)) &
       /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
    enddo

    if (ypn.gt..99e30) then !The upper boundary condition is set either to be
      qn=0.            !"natural"
      un=0.
    else               !or else to have a specified first derivative
      qn=0.5
      un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
    endif

    y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)

    do k=n-1,1,-1 !This is the backsubstitution loop of the tridiagonal algorithm
      y2(k)=y2(k)*y2(k+1)+u(k)
    enddo

    return
  END SUBROUTINE mmisc_bc_spline

  ! **********************************************************************!
  SUBROUTINE mmisc_bc_splint(xa,ya,y2a,n,x,y)
    ! From Numerical Recipes in Fortran, Press et al.
    !Given the arrays xa(1:n) and ya(1:n) of length n, which tabulate a function (with the
    !xai's in order), and given the array y2a(1:n), which is the output from spline above,
    !and given a value of x, this routine returns a cubic-spline interpolated value y.
    IMPLICIT NONE

    INTEGER :: n
    REAL(KIND(0.D0)) :: x,y,xa(n),y2a(n),ya(n)
    INTEGER :: k,khi,klo
    REAL(KIND(0.D0)) :: a,b,h

    klo=1 !We will find the right place in the table by means of bisection
    khi=n

    !This is optimal if sequential calls to this routine are at random
    !values of x. If sequential calls are in order, and closely
    !spaced, one would do better to store previous values of
    !klo and khi and test if they remain appropriate on the
    !next call.


    1 if (khi-klo.gt.1) then
      k=(khi+klo)/2
      if(xa(k).gt.x)then
        khi=k
      else
        klo=k
      endif
      goto 1
    endif !klo and khi now bracket the input value of x.

    h=xa(khi)-xa(klo)

    if (h.eq.0.) STOP 'bad xa input in splint' !The xa's must be distinct.
      a=(xa(khi)-x)/h !Cubic spline polynomial is now evaluated.
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+&
       ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.
    return
  END SUBROUTINE mmisc_bc_splint

 ! **********************************************************************!

! ********************************************************************!
! subroutines and functions below are pressure calculation related    !
! ********************************************************************!

!***************************************************************************************!
SUBROUTINE mmisc_bc_apply_bc_phi( r, rbc)

    !*****************************************************************!
    !*   given phi right outside of domain, phi_bc, (from larger,    *!
    !*   coraser mesh), add values to correct laplacian(D*D^T)of phi *!
    !*   , mass_rhs   on the (smaller, finer) domain, r.             *!
    !*****************************************************************!

  USE parameters
  USE grid
  IMPLICIT NONE
  REAL(KIND(0.D0)), DIMENSION(:,:) :: r
  REAL(KIND(0.D0)), DIMENSION(:)   :: rbc
  INTEGER                          :: i,j

  ! add bc's from coarser grid
  ! TOP and BOTTOM
  DO i=1,m
     r(i,1) = r(i,1) - rbc(bottom_phi + i)
     r(i,n) = r(i,n) - rbc(   top_phi + i)
  END DO

    ! LEFT and RIGHT
    DO j=1,n
       r(1,j) = r(1,j) - rbc( left_phi + j)
       r(m,j) = r(m,j) - rbc(right_phi + j)
    END DO

  END SUBROUTINE mmisc_bc_apply_bc_phi

!************************************************************************************!

SUBROUTINE mmisc_bc_get_bc_phi( r, rbc )

    !******************************************************************!
    !*   given potential,phi on a larger, coarser mesh,               *!
    !*   interpolate it's values to the edge of a smaller, finer mesh *!
    !******************************************************************!
  USE parameters
  USE grid
  IMPLICIT NONE
  REAL(KIND(0.D0)), DIMENSION(:,:) :: r
  REAL(KIND(0.D0)), DIMENSION(:)   :: rbc
  INTEGER                          :: i,j

  ! get interpolated boundary conditions on finer grid
  ! TOP and BOTTOM
  DO i=0,m-2,2
     rbc(bottom_phi + i+1) = .25d0*( r(m/4+i/2,n/4) + &
                               2.d0* r(m/4+1+i/2,n/4) + r(m/4+1+i/2,n/4+1) )
     rbc(top_phi    + i+1) = .25d0*( r(m/4+i/2,n/4+n/2+1) + &
                               2.d0* r(m/4+1+i/2,n/4+n/2+1) + r(m/4+1+i/2,n/4+n/2) )
  END DO

  DO i=1,m-1,2
     rbc(bottom_phi +i+1)  = .25d0*( r(m/4+2+(i-1)/2,n/4) + &
                               2.d0* r(m/4+1+(i-1)/2,n/4) + r(m/4+1+(i-1)/2,n/4+1) )
     rbc(top_phi    +i+1)  = .25d0*( r(m/4+2+(i-1)/2,n/4+n/2+1) + &
                               2.d0* r(m/4+1+(i-1)/2,n/4+n/2+1) + r(m/4+1+(i-1)/2,n/4+n/2) )
  END DO

  ! LEFT and RIGHT
  DO j=0,n-2,2
     rbc(left_phi  + j+1)  = .25d0*( r(m/4,n/4+j/2) + &
                               2.d0* r(m/4,n/4+1+j/2) + r(m/4+1,n/4+1+j/2) )
     rbc(right_phi + j+1)  = .25d0*( r(m/4+m/2+1,n/4+j/2) + &
                               2.d0* r(m/4+m/2+1,n/4+1+j/2) + r(m/4+m/2,n/4+1+j/2) )
  END DO
  DO j=1,n-1,2
     rbc(left_phi  + j+1)  = .25d0*( r(m/4,n/4+(j-1)/2) + &
                               2.d0* r(m/4,n/4+1+(j-1)/2) + r(m/4+1,n/4+1+(j-1)/2) )
     rbc(right_phi + j+1)  = .25d0*( r(m/4+m/2+1,n/4+2+(j-1)/2) + &
                               2.d0* r(m/4+m/2+1,n/4+1+(j-1)/2) + r(m/4+m/2,n/4+1+(j-1)/2) )
  END DO

END SUBROUTINE mmisc_bc_get_bc_phi

 ! **********************************************************************!
! ********************************************************************!
! subroutines and functions above are pressure calculation related    !
! ********************************************************************!
 ! **********************************************************************!

END MODULE mmisc_bc
