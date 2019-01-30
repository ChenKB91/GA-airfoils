MODULE mmisc

  ! This module contains subroutines and functions of numerical methods

  USE mmisc_bc
  USE mmisc_opr
  USE mmisc_motion

  IMPLICIT NONE

CONTAINS

  !*****************************************************************!

   SUBROUTINE mmisc_choldc(dtype)

    !***************************************************************!
    !*   Cholesky factorization of A                               *!
    !***************************************************************!
    USE dts
    USE grid
    IMPLICIT NONE
    TYPE(chol_data_t), INTENT(inout) :: dtype
    REAL(KIND(0.D0)) :: chsum
    INTEGER :: i,j,k
    DO i = 1, Nf
       DO j = i, Nf
          chsum = dtype%cholmat(i,j)
          DO k = i-1, 1, -1
             chsum = chsum - dtype%cholmat(i,k) * dtype%cholmat(j,k)
          END DO
          IF (i.eq.j) THEN
            IF (chsum.le.0.) STOP 'choldc failed'
            dtype%cholvec(i) = sqrt(chsum)
          ELSE
            dtype%cholmat(j,i) = chsum / dtype%cholvec(i)
          END IF
       END DO
    END DO
  END SUBROUTINE mmisc_choldc

  ! *******************************************************************!
    Function cholsl(b, cholmat, cholvec) Result(x)

      !***************************************************************!
      !*   Solve A x = b given it's Cholesky decomposition           *!
      !***************************************************************!
      USE grid
      IMPLICIT NONE

      REAL(KIND(0.D0)), DIMENSION(Nf,Nf)   :: cholmat
      REAL(KIND(0.D0)), DIMENSION(Nf)      :: cholvec
      REAL(KIND(0.D0)), DIMENSION(:)       :: b
      REAL(KIND(0.D0)), DIMENSION(SIZE(b)) :: x

      INTEGER :: i,k
      REAL(KIND(0.D0)) ::  chsum

      DO i=1,Nf
         chsum=b(i)
         DO  k=i-1,1,-1
            chsum=chsum-cholmat(i,k)*x(k)
         END DO
         x(i)=chsum/cholvec(i)
      END DO
      DO i=Nf,1,-1
         chsum=x(i)
         DO k=i+1,Nf
            chsum=chsum-cholmat(k,i)*x(k)
         END DO
         x(i)=chsum/cholvec(i)
      END DO

    END Function cholsl

! *******************************************************************!

 FUNCTION gaussianforce(x, y, xc, yc, wf)

    !************************************************************************!
    !*    creates a 2D Gaussian distribution of amplitude 1, centred at     *!
    !*    xc, yc, and with width defined by parameter a. It is used as a    *!
    !*    waveform for the control body force which multiplies it by an     *!
    !*    amplitude determined by the array control                         *!
    !************************************************************************!

    IMPLICIT NONE
    REAL(KIND(0.D0)) :: x,y, gaussianforce
    REAL(KIND(0.D0)) :: xc, yc,wf

   gaussianforce = 0.d0
   gaussianforce =               EXP(  -(x-xc)**2.D0 / (2.d0/9.d0*wf**2.d0)  ) ! x-dir function
   gaussianforce = gaussianforce*EXP(  -(y-yc)**2.D0 / (2.d0/9.d0*wf**2.d0)  ) ! y-dir function

 END FUNCTION gaussianforce

 !*****************************************************************!

 FUNCTION taylor_vortex(x, y, r, xc, yc, t )

    !***************************************************************!
    !* Given the radius, r, of the (Taylor) vortex at              *!
    !* nondimensionalized time t, the function returns the value   *!
    !* of vorticity at coordinate (x,y) with the center of vortex  *!
    !* at xc and y-yc                                              *!
    !***************************************************************!

   USE parameters
   IMPLICIT NONE

   REAL(KIND(0.D0)) :: x,y,xc,yc,t0,t,r,r2
   REAL(KIND(0.D0)) :: taylor_vortex

   t0=real(Re)*r**2/2.D0
   t0=t0/t
   r2 = (x-xc)**2.D0 + (y-yc)**2.D0
   taylor_vortex = 2.D0*t0**2/r*(1.D0-r2/2.D0*t0/r**2)*EXP( 5.D-1*(1.D0-r2*t0/r**2) )


 END FUNCTION taylor_vortex
! ****************************************************************!

END MODULE mmisc
