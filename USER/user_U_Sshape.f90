      MODULE User                                             !user prepared module to get a new time step
      IMPLICIT NONE                                           !all names must be declared
                                                              !subs and variables that must be exported: used by the solver
      PUBLIC:: NextStep2layer                                 !take the current state and produce the next level
      PUBLIC:: INIT                                           !preparation routine
      PUBLIC:: DONE                                           !finalzing routine
      PUBLIC:: READ_DATA                                      !reading data from a file
      PUBLIC:: Quantity                                       !evaluates the quantity to be compared with measured data
      PUBLIC:: Norm                                           !norming user function (use the default unless youknow what you want)
      PUBLIC:: LSQ                                            !user norm to compare curves. Use the default unless you know what you want.
      PUBLIC:: EarlyQuit                                      !this sub allows early stop in case of events (like reaching zero).
      PUBLIC:: IndependentVar                                 !evaluate the independent var if this is not "time". 
      PUBLIC:: Diagnost                                       !do here whatever you want to see what is happening.
      PUBLIC:: N_params, N_int_params, params, int_params     !numbers of: real parameter, integer parameters; arrays of them
      PUBLIC:: n_unknowns                                     !number of unknowns (spatial grid points plus any additional variables)
      PUBLIC:: NOP                                            !number of time points
      PUBLIC:: threshold                                      !the threshold to pass/reject results (<0 means 'take all')
      PUBLIC:: ex_data, norming_method, space                 !array for measured data for comparing; methods of norming and of comparing
      PUBLIC:: dt                                             !the time step

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      DOUBLE PRECISION::   threshold = -1., dt                !negative threshold means "print the value", while positive: "say if ok"
      INTEGER, PARAMETER:: N_params = 10, N_int_params = 2    !parameters, real and int
      INTEGER, PARAMETER:: N_fixed_params = 0                 !non-varied paramerers, if any
      INTEGER:: n_unknowns
      DOUBLE PRECISION, DIMENSION(N_params):: params          !array of parameters
      INTEGER, DIMENSION(N_int_params):: int_params           !array of int parameters
      DOUBLE PRECISION, DIMENSION(:,:), POINTER:: ex_data     !data array: must be Mx2 with M = as many as you need; to be allocated later

      INTEGER:: NOP = 100                                     !number of time steps; to be changed

      !norms for compaing curves; space chooses the banach space, i.e. the norm to evaluate;
      !norming_method is also space, but to evaluate how 'big' is the experimental curve: to divide the norm of difference between model and experimental curves.
      CHARACTER(LEN=2), PARAMETER:: norming_method = 'L2' !L2 (sqrt int ex^2), C0 (max abs ex), 01 (no norming), L1 (int abs ex), my (user-provided function norm)
      CHARACTER(LEN=2), PARAMETER:: space = 'L2' !L2 (sqrt int ex^2, LSQ), C0 (max abs ex), L1 (int abs ex), my (user-provided function LSQ)

      PRIVATE                                                 !from here on everything is for use inside the module only

      DOUBLE PRECISION, PARAMETER:: pi = 3.1415926            !just pi
      DOUBLE PRECISION, PARAMETER:: mu = 1.46d21              !the H kinetic constant
      DOUBLE PRECISION:: c2, c1                               !H concentrations in two phases near the surface
      DOUBLE PRECISION:: P = 470. !torr                       !H pressure   
      DOUBLE PRECISION:: D, s1, bb1, s2, bb2, eta, L, rho0, S0, gamm, tgamm, multip, shift !parameters
      INTEGER, PARAMETER:: n_curves = 7                       !number of experimental curves
      INTEGER:: leng, curve_index, sigma                      !Number of points in experimental data, chosen curve, geometry parameter
      DOUBLE PRECISION, DIMENSION(3):: curva = [0,0,0]        !to be adjusted by putting 1 somewhere, see sigma
      DOUBLE PRECISION:: time_limit, min_rho=0.               !time limit, minimal hydride core
      DOUBLE PRECISION:: was, is,  gone                       !to control conservation
      INTEGER, PARAMETER:: notC = 3                           !variables beside distribution of H concentration
      INTEGER, PARAMETER:: left = 1+notC                      !all variables are in one array; left is the first concentration.
      INTEGER:: right, xNOP                                   !right is the last concentration, xNOP is amount of concentrations
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: r0, U0    !concentration distribution U0 against nodes r0
      DOUBLE PRECISION, DIMENSION(n_curves), PARAMETER:: time_limits = [6965., 2000., 967., 1954., 146., 486., 1953.] !durations of experiments
      DOUBLE PRECISION, DIMENSION(n_curves), PARAMETER:: pressures = [470., 470., 470., 470., 470., 470., 470.]       !pressures in different experiments

      CONTAINS                                                !subroutines

      SUBROUTINE READ_DATA                                    !sub for reading exper. data from a file
      IMPLICIT NONE
         INTEGER, DIMENSION(n_curves), PARAMETER:: lengs = [1401, 380, 169, 393, 50, 94, 389]   !lengths of files
         CHARACTER(len=64):: filename = "./curve#.dat"        !file name to be constructed
         INTEGER:: i
         threshold = +0.02                                    !threshold can be adjuced here, for each curve its own
         threshold = +0.12                                    !threshold can be adjuced here, for each curve its own
         threshold = +0.20                                    !threshold can be adjuced here, for each curve its own
         if(curve_index.eq.1) threshold=+0.02
         leng = lengs(curve_index)                            !length of the chosen (by curve_index) curve
         write(filename(8:8),'(I1.1)') curve_index            !construct the filename by putting number instead of #
         open(unit=2,file=filename,status='old',access='sequential',form='formatted')    !open the file
         if(.not.associated(ex_data)) allocate(ex_data(leng,2))                            !create an array for the data
         do i = 1,leng                                        !read data pairs
               read(2,*) ex_data(i,:)                         !let Fortran decide the format of the data
         end do
         close(2)                                             !close the file
         CONTINUE                                             !bye
      END SUBROUTINE

      SUBROUTINE INIT                                         !sub for preparing everything before calculation begins
      IMPLICIT NONE
      !initial guess of parameters
         DOUBLE PRECISION, PARAMETER:: D_ig = 7.2643d-12      !diffusivity, initial guess
         DOUBLE PRECISION, PARAMETER:: s2_ig = 2.5221d-28     !adhesion for phase2, initial guess
         DOUBLE PRECISION, PARAMETER:: s1_ig = 4.1994d-31     !adhesion for phase1, initial guess
         DOUBLE PRECISION, PARAMETER:: b2_ig = 1.2790d-04     !desorption for phase2, initial guess
         DOUBLE PRECISION, PARAMETER:: b1_ig = 0.0d0          !desorption for phase1, initial guess
         DOUBLE PRECISION, PARAMETER:: L_ig = 8.6913d-05      !radius of the particle, initial guess
         DOUBLE PRECISION, PARAMETER:: rho0_ig =  0.63322     !initial radius of the core, initial guess
         DOUBLE PRECISION, PARAMETER:: gamm_ig =  1.1525e-17  !gamma parameter for the nuceation phase, initial guess
         DOUBLE PRECISION, PARAMETER:: tgamm_ig =  29.785     !duration of the nucleation phase, initial guess
         DOUBLE PRECISION, PARAMETER:: m_ig =  0.81509        !linear transformation factor, initial guess
         DOUBLE PRECISION, PARAMETER:: shift_ig = 0.0         !linear transforation shift, initial guess
      !parameters, span
         DOUBLE PRECISION, PARAMETER:: minD = D_ig*0.1,        maxD = D_ig * 10.
         DOUBLE PRECISION, PARAMETER:: mins2 = s2_ig,          maxs2 = s2_ig
         DOUBLE PRECISION, PARAMETER:: mins1 = s1_ig,          maxs1 = s1_ig
         DOUBLE PRECISION, PARAMETER:: minb2 = b2_ig,          maxb2 = b2_ig
         DOUBLE PRECISION, PARAMETER:: minb1 = b1_ig,          maxb1 = b1_ig
         DOUBLE PRECISION, PARAMETER:: minL = L_ig,            maxL = L_ig
         DOUBLE PRECISION, PARAMETER:: minV0 = 0.5,            maxV0 = 0.99
         DOUBLE PRECISION, PARAMETER:: mingamm = 1.d-2*1.d-13*gamm_ig,maxgamm = 5.*1.d-13*gamm_ig
         DOUBLE PRECISION, PARAMETER:: mintgamm = 0.,          maxtgamm = 2.*tgamm_ig
         DOUBLE PRECISION, PARAMETER:: minm = m_ig*0.8,        maxm = m_ig*1.2
         DOUBLE PRECISION, PARAMETER:: minshift = shift_ig*2., maxshift = -shift_ig*2.
         DOUBLE PRECISION, PARAMETER:: minR0 = 0.20,           maxR0 = 0.90 !maximal values
     !the code
         INTEGER:: i
         xNOP = 50                                                            !spatial grid size
         n_unknowns = xNOP+notC                                               !grid + other variables
         right = n_unknowns                                                   !now grid is from U(left) to U(right)
         curve_index = mod(int_params(1), n_curves) + 1    ! ;curve_index=2   !choose a curve using the first int param
         time_limit = time_limits(curve_index)             !                  !time limit equal to that of the exper.curve 
         dt = 5.0d0                                        !                  !the time step
         NOP = ceiling(time_limit / dt)                    !                  !number of time steps
         P = pressures(curve_index)                        !                  !constant pressure
         sigma = mod(int_params(2), 3)                     ! ;sigma=2         !geometry parameter: 2 for spherical particles, 1 for cylinders, 0 for plates
         curva(3-sigma) = 1                                ! ;curva=[1,0,0]    !describes the shape function: area of the shape's surface
         D  = linear(minD, maxD,  params(1))               ! ;D =D_ig         !diffusivity
         s2 = linear(mins2, maxs2,  params(2))             ! ;s2=s2_ig        !adgesion for phase2
         s1 = linear(mins1, maxs1,  params(3))             ! ;s1=s1_ig        !adhesion for phase1
         bb2 = linear(minb2, maxb2,  params(4))            ! ;bb2=b2_ig       !desorption for phase2
!         bb1 = linear(minb1, maxb1,  params(5))           !  !;bb1=b1_ig      !desorption for phase1
         bb1=b1_ig !constant param                         ! 
         L  = linear(minL,maxL, params(5))                 ! ;L =L_ig         !radius of the particle
         rho0 = linear(minR0,maxR0, params(6))             ! ;rho0=rho0_ig    !radius of the core of phase2 relative to L
         rho0 = rho0 * L                                   !                  !radius of the core
         eta= 0.0d0                                        !                  !how much H shifts the boundary, else broadens the skin
         gamm=   linear(mingamm,maxgamm,params(9))*Surf(L) !;gamm= gamm_ig
         tgamm=  linear(mintgamm,maxtgamm,params(10))      !;tgamm= tgamm_ig
         multip= linear(minm,maxm,params(7))               !;multip= m_ig    !linear transformation
         shift=  linear(minshift,maxshift,params(8))       !;shift= shift_ig
         S0 = 0.0d0                                                           !initial skin area
         min_rho = 0.1 * (L-rho0)                                                   !size of the core to be ignored
       !OPEN(42, file='diagn.out', status='replace')                           !open the diagnostic output
!         c2 = density / molar_mass * Avo * stoichiometric !H/cm^3
         c2 = 1.                                                              !concentration of H in phase2 near the surface 
!         c1 = c1_div_c2*c2
         c1 = 0.                                                              !concentration of H in phase2 near the surface 
         gone = 0.                                                            !no hydride left/entered initially
         if(.not. allocated(r0)) allocate(r0(xNOP), U0(xNOP))                                         !prepare arrays for concentration distribution
         r0 = rho0 + dble([(i, i=0,xNOP-1)])/dble(xNOP-1)*(L-rho0)                  !fill the grid nodes
         U0 = steady(r0)                                                      !fill the steady distribution using the sub
         was = Amount(r0, U0, S0)                                             !how much hydrogen was in the particle at t=0
         CONTINUE
      END SUBROUTINE

 SUBROUTINE NextStep2layer(t,t1,U)              ! 2-layer proceeder: takes current time, previous time, and previos state; return current state
   USE Parabolic1D, ONLY: derx, DLlinRnonlin, DIRICHLET_I, NEUMANN_II !use some stuff for solving parabolic 1D PDE
   use Integration, ONLY: trapz                                      !use the trapz numerical integration routine
   IMPLICIT NONE
     DOUBLE PRECISION, INTENT(IN):: t,t1                             !the current time and the previous time
     DOUBLE PRECISION, INTENT(INOUT), DIMENSION(:), TARGET:: U       !the phase vector, including concentration distribution and other variables
     DOUBLE PRECISION:: rho, J, JS, JSdt, drho, S, dS, cL, J1,J2, W  !different physial quantities
     DOUBLE PRECISION, DIMENSION(xNOP):: r                           !spatial grid
     DOUBLE PRECISION, DIMENSION(:), POINTER:: V                     !concentration distribution
     INTEGER:: i, ind                                                !counters
     V => U(left:right)                                              !assign the concentration distribution
     JSdt = 0.                                                       !flux: zero initially
     rho = U(1) + rho0; rho = max(rho,0.)                            !current radius of the shrinking core (U(1) stores rho-rho0 because of zero initial conditions)
     S = U(2) + S0; S = max(min(S,1.),0.)                            !current area under the uncompleted skin
     V = V + U0                                                      !distribution: U-U0 was stored because of zero initial conditions
     if(t1 < tgamm) then                                             !the nucleation stage
        U(3) = gamm*t**3                                             !simple formula for the amount of hydrogen stored in U(3)
        gone = gone - gamm*(t**3.-(t-dt)**3.)                        !change the amount gone/come
        was = was - gamm*(t**3.-(t-dt)**3.)                          !a trick: this amount is not in 'is', so we substract it from was, as if it was not there
     else                                                            !skin formation or shrinking core
        WHERE(V < 0.) V = 0.0                                        !repair negative values, if any
        cL = U(right)                                                !concentration near the surface
        J1 = mu*s1*P - bb1*c1**2                                     !total flux density for phase1
        J2 = mu*s2*P - bb2*cL**2                                     !total flux density for phase2
        J = J2*S + J1*(1.-S)                                         !get the total flux density from/to the surface
        if(S<1.) then                                                !if skin is not complete yet
           r = rho + [(i, i=0,xNOP-1)]/dble(xNOP-1)*(L-rho)          !build the grid points
           W  = trapz(r,V*Surf(r)) / Surf(L)                         !calculate geometry factor using an integrator
           dS = ( J2*S + (1.-eta)*J1*(1.-S) ) / W                    !evaluate the derivative of surface of the skin
        else                                                         !in case of complete skin
           dS = 0.0d0                                                !do nothing, no further growth
           eta = 1.                                                  !all hydrogen now moves the boundary
        end if
        S = S + dS*dt                                                !the Euler time step for surface of the skin
        S = max(min(S,1.),0.)
        if(rho>min_rho) then                                         !if the core is not too small
           drho = eta * D*derx(V,1,L,rho)/(c1-V(1))               !evaluate its derivative
           rho  = rho + drho*dt                                      !the Euler time step for the core radius
           rho  = max(rho,0.)                                        !repair negative value
           call DLlinRnonlin(dt, V, D, curva, DIRICHLET_I, LBC, RBC, L, rho, 0.0d0, drho) !call the solver for a PDE with linear left BC and nonlinear right BC
        else                                                         !if there is no more core 
           drho = 0.0d0                                              !no further change
           rho  = 0.0d0                                              !no core means no core
           call DLlinRnonlin(dt, V, D, curva, NEUMANN_II, LBC2, RBC, L, rho, 0.0d0, drho) !call the same solver but with different left BC
        end if
        WHERE(V < 0) V = 0.0                                         !repair negative concentrations
        cL = V(xNOP)                                                !concentration near the surface
        J1 = mu*s1*P - bb1*c1**2                                     !total flux density for phase1
        J2 = mu*s2*P - bb2*cL**2                                     !total flux density for phase2
        J = (   J + J2*S + J1*(1.-S)   ) / 2.                        !get the total flux density from/to the surface
        JS = J*Surf(L)                                               !get the total flux
        JSdt = JS*dt                                                 !get the amount that crossed the surface during the time step
        U(3) = U(3) + JSdt                                           !change the amount left
        gone = gone - JSdt                                              !change the amount gone/come
     end if
     r = rho + [(i, i=0,xNOP-1)]/dble(xNOP-1)*(L-rho)                !get the grid points
     is = Amount(r,V,S)                                              !evaluate the current amount
     V = V - U0                                                      !substract initial data because the problem solved has zero initial conditions
     U(1) = rho - rho0
     U(2) = S - S0
   RETURN

   CONTAINS                                                          !PDE coefficients to be given to PDE solvers

      PURE FUNCTION LBC()                                            !the Dirichlet boundary condition
       IMPLICIT NONE
         DOUBLE PRECISION:: LBC
         LBC = c2                                                    !simply the equilibrium concentration
      END FUNCTION

      PURE FUNCTION LBC2()                                           !the Neumann BC (single-phase stage)
       IMPLICIT NONE
         DOUBLE PRECISION:: LBC2
         LBC2 = 0.0d0                                                !u'=0
      END FUNCTION

      PURE FUNCTION RBC(c)                                           !the right-hand side of nonlinar right BC: Du'=mu s P - b c^2
       IMPLICIT NONE
       DOUBLE PRECISION:: RBC
       DOUBLE PRECISION, INTENT(IN):: c                              !surface concentration in
         RBC = (mu*s2*P - bb2*c**2) / D                              !u'= this
      END FUNCTION

 END SUBROUTINE

      SUBROUTINE DONE                                                !finalising sub
      IMPLICIT NONE
      DOUBLE PRECISION:: disb                                        !disbalance: violation of the conservation
      DOUBLE PRECISION:: clim
         if(associated(ex_data)) deallocate(ex_data)                 !free allocated data
         if(allocated(r0)) deallocate(r0, U0)                                         !prepare arrays for concentration distribution
         clim = sqrt(mu*s2*P/bb2)
       !  disb = (was-is-gone)/(Vol(L)*clim)*100.                       !evaluate disbalance
       if(abs(was-gone)<1.d-32) return
         disb = (was-is-gone)/(was-gone)*100.                       !evaluate disbalance
         if (abs(disb) > 5.) write(*, *) 'BEWARE!!!'                !write a warning if too high
         write(*, '("Matter disbalance is ", F6.2, "%")') disb      !report the disbalance
         write(*, '("Was ", F10.5, " is ", F10.5, " gone ", F10.5)') [was,is,gone]/vol(L)     !report the amounts of was, is now, and gone/arrived
         CONTINUE
      END SUBROUTINE

      PURE FUNCTION IndependentVar(t)                                !sub for evaluating the independent variable if it is not time
      IMPLICIT NONE
         DOUBLE PRECISION, INTENT(IN):: t
         DOUBLE PRECISION:: IndependentVar
         IndependentVar = t                                          !it is time in our case
      END FUNCTION

      ELEMENTAL FUNCTION Steady(r)                                   !sub for the steady distribution: in a point or in many points
      IMPLICIT NONE
         DOUBLE PRECISION, INTENT(IN):: r                            !point in real coordinates or array of points
         DOUBLE PRECISION:: Steady
         DOUBLE PRECISION:: a,b,c,Discr, AA,BB                       !auxiliary variables
         SELECT CASE(sigma)                                          !steady distribution depends on geometry
         CASE(2)                                                     !if a sphere,
                a = bb2*((L-rho0)/(L*rho0))**2                       !prepare...
                b = D/L**2 + 2.*c2*bb2*(L-rho0)/(L*rho0)             !the square...
                c = -(mu*s2*P - bb2*c2**2)                           !equation,
                Discr = b**2 - 4.*a*c                                !calculate its descriminant,
                BB = (sqrt(Discr) - b) / (2.*a)                      !and solve it
                AA = c2 + BB/rho0                                    !use the other BC to build the curve
                Steady = AA - BB / r                                 !evaluate the distribution (hyperbola)
         CASE(1)                                                     !if a cylinder,                                 
                a = bb2*log(L/rho0)**2                               !prepare...
                b = D/L + 2.*bb2*c2*log(L/rho0)                      !the square...
                c = -(mu*s2*P - b*c2**2)                             !equation,
                Discr = b**2 - 4.*a*c                                !calculate its descriminant,
                BB = (sqrt(Discr) - b) / (2.*a)                      !and solve it
                AA = c2 - BB*log(rho0)                               !use the other BC to build the curve
                Steady = AA + BB * log(r)                            !evaluate the distribution (logarithm)
         CASE(0)                                                     !if a plate,                                 
                a = bb2*(L-rho0)**2                                  !prepare...
                b = D + 2.*c2*bb2*(L-rho0)                           !the square...
                c = -(mu*s2*P - bb2*c2**2)                           !equation,
                Discr = b**2 - 4.*a*c                                !calculate its descriminant,
                BB = (sqrt(Discr) - b) / (2.*a)                      !and solve it
                AA = c2 - BB*rho0                                    !use the other BC to build the curve
                Steady = AA + BB * r                                 !evaluate the distribution (linear)
         END SELECT
         Steady = max(Steady,0.)                                     !repair negative values, if any
      END FUNCTION

      PURE FUNCTION Quantity(t,U)                                    !sub for evaluating the quantity to be compared with mesurements
      IMPLICIT NONE
         DOUBLE PRECISION, INTENT(IN):: t                            !time
         DOUBLE PRECISION, INTENT(IN), DIMENSION(:):: U              !variables, including the concentration distribution
         DOUBLE PRECISION:: Quantity
         DOUBLE PRECISION:: clim
         clim = sqrt(mu*s2*P/bb2)
         Quantity = U(3)/(Vol(L)*clim) * multip + shift                !the quantity is the hydrogen amount stored in U(3), only normed and after the linear transform
         return
      END FUNCTION

      FUNCTION EarlyQuit(t,U)                                        !some event can, possibly, stop the calculation
      IMPLICIT NONE
         DOUBLE PRECISION, INTENT(IN):: t                            !time
         DOUBLE PRECISION, INTENT(IN), DIMENSION(:):: U              !variables, including the concentration distribution
         LOGICAL:: EarlyQuit
         EarlyQuit = .false.                                         !not in our case
      END FUNCTION

      ELEMENTAL FUNCTION Vol(r)                                      !volume of the particle (or similar bodies)
      IMPLICIT NONE
        DOUBLE PRECISION, INTENT(IN):: r                             !radius in this sub!!!
        DOUBLE PRECISION:: Vol
         Vol  = Surf(r)*r/(sigma+1.)                                 !general formula for balls, cylinders, plates
      END FUNCTION

      ELEMENTAL FUNCTION Surf(r)                                     !surface of the particle (or similar bodies)
      IMPLICIT NONE
        DOUBLE PRECISION, INTENT(IN):: r                             !radius in this sub!!!
        DOUBLE PRECISION:: Surf
         select case(sigma)                                          !depending on what body it is
         case(2)
                 Surf = 4.*pi*r**2                                   !sphere
         case(1)
                 Surf = 2.*pi*r                                      !cylinder (of unit length)
         case(0)
                 Surf = 1.                                           !plate (of unit surface)
         end select 
      END FUNCTION

      PURE FUNCTION Amount(r, U, S)                                  !amount of hydrogen in the particle at any time
      USE Integration, ONLY: trapz                                   !need the integrator
      IMPLICIT NONE
        DOUBLE PRECISION, INTENT(IN), DIMENSION(:):: r, U            !grid points and values of the concentration
        DOUBLE PRECISION, INTENT(IN):: S                             !area of the skin, also necessary
        DOUBLE PRECISION:: Amount
        DOUBLE PRECISION:: rho, L                                    !core radius and particle's size are in r already
        INTEGER:: NOP                                                !grid size
        NOP = size(r)                                                !taken from the grid
        rho = r(1); L = r(NOP)
        Amount = c1 * Vol(rho) + &                                   !H inside the core: concentration their is c1, volume is Vol
                & S*trapz(r, U*Surf(r)) + &                          !H inside the new-phase skin: relative area times integral over volume
                & (1.-S)*trapz(r, c1*Surf(r))                        !H inside the not-yet-formed skin: relative area ties integral over surface
      END FUNCTION

      FUNCTION norm(t,x)                                   !function for the norm of the experimental curve to norm the difference between model and measured data
      IMPLICIT NONE
      DOUBLE PRECISION:: norm                              !use this only if none of the pre-designed norms suit you
         DOUBLE PRECISION, INTENT(IN), DIMENSION(:):: t,x  !grid and data
         norm = 1.0                                        !no norming by default
         CONTINUE
      END FUNCTION

      SUBROUTINE Diagnost()                                !diagnostic sub. Do whatever you want to learn what is happening
        implicit none
        return
      END SUBROUTINE

      FUNCTION LSQ(t,x)                                    !function for the norm of difference between measured and model curves                                   
      IMPLICIT NONE
         DOUBLE PRECISION:: LSQ                            !use this only if none of the pre-designed norms suit you
         DOUBLE PRECISION, INTENT(IN), DIMENSION(:):: t,x  !grid and data
         LSQ = 0.0                                         !replace this with real calculations if necessary
         CONTINUE
      END FUNCTION

      FUNCTION linear(min_value, max_value, value01)       !linear interpolator: to choose a point between M1 and M2 given an m from [0,1]
      IMPLICIT NONE
         DOUBLE PRECISION:: linear
         DOUBLE PRECISION, INTENT(IN):: min_value, max_value, value01
         linear = min_value + value01 * (max_value - min_value)
         return
      END FUNCTION

      END MODULE User


