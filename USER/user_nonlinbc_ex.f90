
      MODULE User
      IMPLICIT NONE

      PUBLIC:: NextStep2layer, INIT, DONE, READ_DATA, Quantity, Norm, LSQ, EarlyQuit, IndependentVar, Diagnost
      PUBLIC:: N_params, N_int_params, params, int_params
      PUBLIC:: n_unknowns, NOP
      PUBLIC:: threshold, ex_data, norming_method, space, dt               !X: this line must not be changed!

      DOUBLE PRECISION:: threshold = -1., dt

      INTEGER, PARAMETER:: N_params = 2, N_int_params = 0
      INTEGER, PARAMETER:: N_fixed_params = 0
      INTEGER:: n_unknowns
      DOUBLE PRECISION, DIMENSION(N_params):: params
      INTEGER, DIMENSION(N_int_params):: int_params
      DOUBLE PRECISION, DIMENSION(:,:), POINTER:: ex_data   !must be Mx2 with M = as many as you need

      INTEGER:: NOP = 100

      CHARACTER(LEN=2), PARAMETER:: norming_method = 'L2' !L2 (sqrt int ex^2), C0 (max abs ex), 01 (no norming), L1 (int abs ex), my (user-provided function norm)
      CHARACTER(LEN=2), PARAMETER:: space = 'L2' !L2 (sqrt int ex^2, LSQ), C0 (max abs ex), L1 (int abs ex), my (user-provided function LSQ)

      PRIVATE

      DOUBLE PRECISION, PARAMETER:: pi = 3.1415926
      DOUBLE PRECISION, PARAMETER:: L= pi/2.
      DOUBLE PRECISION:: D, bb
      INTEGER, PARAMETER:: sigma = 0
      DOUBLE PRECISION:: time_limit
      DOUBLE PRECISION:: was, is,  gone
      INTEGER, PARAMETER:: notC = 0, left = 1+notC                             !notC variables beside distribution of c
      INTEGER:: right, xNOP
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: r0, U0

      CONTAINS

              PURE FUNCTION Z(t)
                      DOUBLE PRECISION:: Z
                      DOUBLE PRECISION, INTENT(IN):: t
                      Z = sqrt(D/bb*exp(-D*t))
              END FUNCTION Z

              PURE FUNCTION dZ(t)
                      DOUBLE PRECISION:: dZ
                      DOUBLE PRECISION, INTENT(IN):: t
                      dZ = -0.5*D*sqrt(D/bb*exp(-D*t))
              END FUNCTION dZ

      SUBROUTINE READ_DATA
      IMPLICIT NONE
         INTEGER:: i
         INTEGER, PARAMETER:: leng=100
         DOUBLE PRECISION:: t
         DOUBLE PRECISION, PARAMETER:: trueD = 31.07d-3, trueb = trueD, dt = 0.1/dble(leng-1)
         threshold = 0.5
         allocate(ex_data(leng,2))
         do i = 1,leng
               t = dble(i-1)*dt
               ex_data(i,1) = t
               ex_data(i,2) = trueD*exp(-trueD*t) 
         end do
         CONTINUE
      END SUBROUTINE

      SUBROUTINE INIT
      IMPLICIT NONE
         DOUBLE PRECISION,PARAMETER:: minD  = 1.d-3,  maxD = 100.d-3 !min, max for diffusivity
         DOUBLE PRECISION,PARAMETER:: minb  = 1.d-3,  maxb = 100.d-3 !min, max for desorption
         INTEGER:: i
	 DOUBLE PRECISION:: h
         xNOP = 101
         n_unknowns = xNOP+notC
         right = n_unknowns
         h = 1. / dble(xNOP-1)
         D  = linear(minD, maxD,  params(1))
 !D = 31.07d-3
         bb  = linear(minb, maxb,  params(2))
 !bb = D
       OPEN(42, file='diagn.out', status='replace')
       write(42, '("Diffusivity is equal to ", 5(ES16.6,1X), 3X, I3)') D
         time_limit = 0.1
         dt = 1.d-4
         NOP = ceiling(time_limit / dt)
         gone = 0.
         allocate(r0(xNOP), U0(xNOP))
         r0 = -L + [(i, i=0,xNOP-1)]*h*(2.*L)
         U0 = Steady(r0)
         was = Amount(r0, U0);     !how much hydrogen has left the particle at t=0
         CONTINUE
      END SUBROUTINE

      SUBROUTINE DONE
      IMPLICIT NONE
      DOUBLE PRECISION:: disb
         disb = (was-is-gone)/was*100.
         if (abs(disb) > 5.) write(42, *) 'BEWARE!!!', was, is, gone
         write(42, '("Matter disbalance is ", F6.2, "%")') disb
         CONTINUE
      END SUBROUTINE

      PURE FUNCTION IndependentVar(t)
      IMPLICIT NONE
         DOUBLE PRECISION, INTENT(IN):: t
         DOUBLE PRECISION:: IndependentVar
         IndependentVar = t
      END FUNCTION

      ELEMENTAL FUNCTION Steady(r)
      IMPLICIT NONE
         DOUBLE PRECISION, INTENT(IN):: r
         DOUBLE PRECISION:: Steady
         Steady = Z(0.d0) + cos(r)
      END FUNCTION

      PURE FUNCTION Quantity(t,U)
      IMPLICIT NONE
         DOUBLE PRECISION, INTENT(IN):: t
         DOUBLE PRECISION, INTENT(IN), DIMENSION(:):: U
         DOUBLE PRECISION:: Quantity
         Quantity = bb * U(right)**2
      END FUNCTION

      FUNCTION EarlyQuit(t,U)
      IMPLICIT NONE
         DOUBLE PRECISION, INTENT(IN):: t
         DOUBLE PRECISION, INTENT(IN), DIMENSION(:):: U
         LOGICAL:: EarlyQuit
         EarlyQuit = .false.
      END FUNCTION

      ELEMENTAL FUNCTION Vol(r)
      IMPLICIT NONE
        DOUBLE PRECISION, INTENT(IN):: r
        DOUBLE PRECISION:: Vol
         Vol  = Surf(r)*r/(sigma+1.)
      END FUNCTION

      ELEMENTAL FUNCTION Surf(r)
      IMPLICIT NONE
        DOUBLE PRECISION, INTENT(IN):: r
        DOUBLE PRECISION:: Surf
         select case(sigma)
         case(2)
                 Surf = 4.*pi*r**2
         case(1)
                 Surf = 2.*pi*r
         case(0)
                 Surf = 1.
         end select 
      END FUNCTION

      PURE FUNCTION Amount(r, U)
      USE Integration, ONLY: trapz
      IMPLICIT NONE
        DOUBLE PRECISION, INTENT(IN), DIMENSION(:):: r, U
        DOUBLE PRECISION:: Amount
        Amount = trapz(r, U*Surf(r))
      END FUNCTION

 SUBROUTINE NextStep2layer(t,t1,U)                             ! User provided 2-layer proceeder: takes current time, previous time, and previos state, return current state
   USE Parabolic1D, ONLY: derx1, LRnonlin, DIRICHLET_I
   use Integration, ONLY: trapz
   IMPLICIT NONE
     DOUBLE PRECISION, INTENT(IN):: t,t1                       !current time and previous time
     DOUBLE PRECISION, INTENT(INOUT), DIMENSION(:), TARGET:: U
     DOUBLE PRECISION, DIMENSION(xNOP):: r
     INTEGER:: i, ind
     U(left:right) = U(left:right) + U0
     call LRnonlin(dt, U(left:right), A,B,B1,E,F, LBC, RBC, L, -L, 0.0d0, 0.0d0)
     gone = gone + Quantity(t,U(left:right)) * Surf(L)*2 * dt
     r = -L + [(i, i=0,xNOP-1)]/dble(xNOP-1)*(2*L)
     is = Amount(r,U(left:right))
     U(left:right) = U(left:right) - U0
   RETURN

   CONTAINS

      PURE FUNCTION A(i)
       IMPLICIT NONE
         DOUBLE PRECISION:: A
         INTEGER, INTENT(IN):: i
         A = D
      END FUNCTION

      PURE FUNCTION B(i)
       IMPLICIT NONE
         DOUBLE PRECISION:: B
         INTEGER, INTENT(IN):: i
         B = 0.0d0
      END FUNCTION

      PURE FUNCTION B1(i)
       IMPLICIT NONE
         DOUBLE PRECISION:: B1
         INTEGER, INTENT(IN):: i
         B1 = 0.0d0
      END FUNCTION

      PURE FUNCTION E(i)
       IMPLICIT NONE
         DOUBLE PRECISION:: E
         INTEGER, INTENT(IN):: i
         E = 0.0d0
      END FUNCTION

      PURE FUNCTION F(i)
       IMPLICIT NONE
         DOUBLE PRECISION:: F
         INTEGER, INTENT(IN):: i
         F = dZ(t1)
      END FUNCTION

      PURE FUNCTION LBC(c)
       IMPLICIT NONE
         DOUBLE PRECISION:: LBC
         DOUBLE PRECISION, INTENT(IN):: c
         LBC = bb * c**2
      END FUNCTION

      PURE FUNCTION RBC(c)
       IMPLICIT NONE
         DOUBLE PRECISION:: RBC
         DOUBLE PRECISION, INTENT(IN):: c
         RBC = -bb * c**2
      END FUNCTION

 END SUBROUTINE

      FUNCTION norm(t,x)
      IMPLICIT NONE
         DOUBLE PRECISION:: norm
         DOUBLE PRECISION, INTENT(IN), DIMENSION(:):: t,x  !grid and data
         norm = 1.0
         ! do whatever you want to produce a positive real number. The L2 norm of the difference between the model curve and your experimental one will be divided on this quantity
         CONTINUE
      END FUNCTION

      SUBROUTINE Diagnost()
        implicit none
        return
      END SUBROUTINE

      FUNCTION LSQ(t,x)
      IMPLICIT NONE
         DOUBLE PRECISION:: LSQ
         DOUBLE PRECISION, INTENT(IN), DIMENSION(:):: t,x
         LSQ = 0.0
         ! do whatever you want to produce a positive real number. The L2 norm of the difference between the model curve and your experimental one will be divided on this quantity
         CONTINUE
      END FUNCTION

      FUNCTION linear(min_value, max_value, value01)
      IMPLICIT NONE
         DOUBLE PRECISION:: linear
         DOUBLE PRECISION, INTENT(IN):: min_value, max_value, value01
         linear = min_value + value01 * (max_value - min_value)
         return
      END FUNCTION

      END MODULE User

      !!perl -E 'for(1..10) {say $s if ($s=`./solver_example`)=~/^1/}'
