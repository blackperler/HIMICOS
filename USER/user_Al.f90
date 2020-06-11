! Solves the problem of dehydriding. Spherical particles.
! -bc**2 = Ddc/dr on the surface (P=0), equilibrium on the phaze boundary.
! No c in metal. Desorption flux density is given for both phazes.
! Quick diffusion. The skin grows both in tangent and radial directions, 
! Desorption is from the phaze boundary.

      MODULE User
      IMPLICIT NONE

      PUBLIC:: INIT, DONE, READ_DATA, Quantity, Norm, LSQ, EarlyQuit, IndependentVar, NextStep2layer
      PUBLIC:: N_params, N_int_params, params, int_params
      PUBLIC:: n_unknowns, NOP, dt
      PUBLIC:: threshold, ex_data, norming_method, space               !X: this line must not be changed!

      DOUBLE PRECISION:: threshold = -1.

      INTEGER, PARAMETER:: N_params = 5, N_int_params = 1
      INTEGER, PARAMETER:: N_fixed_params = 0
      INTEGER, PARAMETER:: n_unknowns = 2
      DOUBLE PRECISION, DIMENSION(N_params):: params
      INTEGER, DIMENSION(N_int_params):: int_params
      DOUBLE PRECISION, DIMENSION(:,:), POINTER:: ex_data   !must be Mx2 with M = as many as you need

      DOUBLE PRECISION:: dt
      INTEGER:: NOP = 100

      CHARACTER(LEN=2), PARAMETER:: norming_method = 'L2' !L2 (sqrt int ex^2), C0 (max abs ex), 01 (no norming), L1 (int abs ex), my (user-provided function norm)
      CHARACTER(LEN=2), PARAMETER:: space = 'L2' !L2 (sqrt int ex^2, LSQ), C0 (max abs ex), L1 (int abs ex), my (user-provided function LSQ)

      PRIVATE

      INTEGER, PARAMETER:: curves = 2
      DOUBLE PRECISION, DIMENSION(curves), PARAMETER:: time_limits = [1.4d4, 6392.d0]
      DOUBLE PRECISION, PARAMETER:: tol = 1e-6, cb = 1.,pi = 3.1415926 !tolerance constant, unit concentration in hydride
      DOUBLE PRECISION, PARAMETER:: c0 = 1d-3, ca = 1.0d0, L=8.d-4
      DOUBLE PRECISION:: Jb,Ja, R0, S0, nu
      INTEGER:: sigma=2, leng, curve !shape factor, Number Of spatial Points, aux
      DOUBLE PRECISION:: sc, time_limit, Surf, Vol

      CONTAINS

      SUBROUTINE READ_DATA
         INTEGER, DIMENSION(curves), PARAMETER:: lengs = [6862, 3148]   !the first is no UV activation, the second is 1hr UV exposued
         CHARACTER(len=64):: filename = "./curve#"
         INTEGER:: i
         threshold = 0.8
         leng = lengs(curve)
         write(filename(8:8),'(I1.1)') curve
         open(unit=2,file=filename,status='old',access='sequential',form='formatted')    !newunit?
         allocate(ex_data(leng,2))
         do i = 1,leng
               read(2,*) ex_data(i,:)
         end do
         close(2)
         ex_data(:,1) = ex_data(:,1) / ex_data(leng,1) !norming to [0,1] time
         CONTINUE
      END SUBROUTINE

      SUBROUTINE INIT
         DOUBLE PRECISION,PARAMETER:: bb = 0.1024711, ba = 3.593d-09     !an initial set of values
         DOUBLE PRECISION,PARAMETER:: maxbb = bb*2.0, maxba = ba*2.0, maxR0 = 0.9*L, maxS0 = 1.0/6.0, maxnu = 1.0     !maximal values
         DOUBLE PRECISION,PARAMETER:: minbb = bb*0.2, minba = ba*0.2, minR0 = 0.5*L, minS0 = 0.0, minnu = 0.1     !minimal values
         DOUBLE PRECISION:: bb_, ba_
         NOP = 1000
         curve = floor(real(int_params(1)) / 5.0) + 1
         bb_ = linear(minbb,maxbb, params(1))
         Jb = -bb_*c0**2
         ba_ = linear(minba,maxba, params(2))
         Ja = -ba_*ca**2
         R0 = linear(minR0,maxR0, params(3))
         S0 = linear(minS0,maxS0, params(4))
         nu = linear(minnu,maxnu, params(5))
         time_limit = time_limits(curve)
         dt = time_limit / NOP
         select case(sigma)
         case(2)
                 Surf = 4.*pi*L**2
         case(1)
                 Surf = 2.*pi*L
         case(0)
                 Surf = 1.
         case default
                 print*,'ERROR: sigma must be 2, 1, or 0!'
                 stop
         end select 
         Vol  = Surf*L/(sigma+1.)
         CONTINUE
      END SUBROUTINE

      SUBROUTINE DONE
         if(associated(ex_data)) deallocate(ex_data)
         CONTINUE
      END SUBROUTINE

      PURE FUNCTION EarlyQuit(t,U)
         IMPLICIT NONE
         DOUBLE PRECISION, INTENT(IN):: t
         DOUBLE PRECISION, INTENT(IN), DIMENSION(:):: U
         LOGICAL:: EarlyQuit
         EarlyQuit = .false.
      END FUNCTION

      PURE FUNCTION Quantity(t,U)
         IMPLICIT NONE
         DOUBLE PRECISION, INTENT(IN):: t
         DOUBLE PRECISION, INTENT(IN), DIMENSION(:):: U
         DOUBLE PRECISION:: Quantity
         DOUBLE PRECISION:: rho, Sb
         rho = U(1) + R0; Sb = U(2) + S0
         Quantity = (Sb*Surf/L**sigma*(L**(sigma+1)-rho**(sigma+1))/(sigma+1)) / Vol !reacted fraction
      END FUNCTION

      PURE FUNCTION IndependentVar(t)
         IMPLICIT NONE
         DOUBLE PRECISION, INTENT(IN):: t
         DOUBLE PRECISION:: IndependentVar
         IndependentVar = t
      END FUNCTION

     SUBROUTINE NextStep2layer(t,t1,U)        ! User provided 2-layer proceeder: takes current time, previous time, and previos state, return current state
        USE ODE, ONLY: PredCorr1step          ! use ODE solver(s)
        IMPLICIT NONE
        DOUBLE PRECISION, INTENT(IN):: t,t1
        DOUBLE PRECISION, INTENT(INOUT), DIMENSION(:):: U
        CALL  PredCorr1step(t, t-t1, U, RHS)  ! solving the ODEs under chosen values of parameters
     END SUBROUTINE

     SUBROUTINE RHS(t, U, dU)                !the right-hand side f(t,x) of the equations \dot x = f(t,x)
         IMPLICIT NONE
        DOUBLE PRECISION, INTENT(IN):: t                  !time t
        DOUBLE PRECISION, INTENT(IN), DIMENSION(:):: U    !the phase vector x(t)
        DOUBLE PRECISION, INTENT(OUT), DIMENSION(:):: dU  !dx/dt
        DOUBLE PRECISION:: rho, S, local_nu, W_           !variables to be used
        dU = 0.0d0
        rho = U(1) + R0; S = U(2) + S0                    !the core size, relative area of the skin
        if(rho <= 0.) rho = 0.                            !repair negative value
        S = max(min(S,1.),0.)                             !repair values outside [0,1]
        local_nu = nu                                     !shrinking core vs skin growth factor nu local to the sub 
        if(S.ge.1) local_nu = 0.                          !if skin is ready, pure shrinking core
        W_ = -cb*(L**(sigma+1)-rho**(sigma+1))/(sigma+1)  !geometrical factor
        if(rho>0) then                                    !two phases co-exist
                dU(1) = (1.-nu)*Jb / cb                   !d\rho/dt: note constnt fluxes Jb and Ja
                dU(2) = (Ja*(1.-S)*L**sigma + nu*S*Jb*rho**sigma) / W_ !dS/dt
        else                                              !single phase
                dU(1) = 0.                                !no core
                dU(2) = (Ja*(1.-S)*L**sigma + nu*S*Jb*rho**sigma) / W_ !S changes
        end if
        dU = dU * time_limit                              !
        RETURN
      END SUBROUTINE

      FUNCTION norm()
         DOUBLE PRECISION:: norm
         norm = 1.0
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



