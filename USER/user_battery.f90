! Solves the problem of dehydriding. Spherical particles.
! -bc**2 = Ddc/dr on the surface (P=0), equilibrium on the phaze boundary.
! No c in metal. Desorption flux density is given for both phazes.
! Quick diffusion. The skin grows both in tangent and radial directions, 
! Desorption is from the phaze boundary.

      MODULE User
      IMPLICIT NONE

      PUBLIC:: NextStep2layer                                !take the current state and produce the next level
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

      DOUBLE PRECISION:: threshold = -1., dt                  !negative threshold means "print the value", while positive: "say if ok"
      INTEGER, PARAMETER:: N_params = 5, N_int_params = 1     !parameters, real and int
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
      DOUBLE PRECISION, PARAMETER:: electron = 1.60217657d-19, e_per_C = 1. / electron   !elementary charge and electrons per a Coulomb 
      INTEGER, PARAMETER:: ncurrents = 9                      !number of experimentl curves
      DOUBLE PRECISION, DIMENSION(ncurrents),PARAMETER:: currents = [29, 58, 87, 116, 148, 268, 358, 448, 597] !currents used to produce the exper.curves
      DOUBLE PRECISION, DIMENSION(ncurrents),PARAMETER:: charges_ex = [346.7, 329.8, 319.6, 303.9, 292.9, 269.9, 252.3, 213.2, 151.3] !charges passed
      DOUBLE PRECISION:: cb, calph                            !concentrations in beta and alpha phases
      DOUBLE PRECISION:: L=2.0D-003                           !size of hydride particles
      DOUBLE PRECISION:: rho0, C2U, I_, D                     !physical constants
      DOUBLE PRECISION:: tafela = 0.940, tafelb = 0.004, save_tafel_a = 0.940 !the Tafel equation ! for ex5 0.01 of b is the same as 0.05 of a
      INTEGER:: leng, I_index                                 !number of exper. points
      INTEGER, PARAMETER:: sigma=2                            !shape factor: sphere
      DOUBLE PRECISION:: time_limit, J                        !time limit and flux density
      DOUBLE PRECISION:: c0_div_cb                            !ratio of concentrations: equilibrium in metal to stoichiometric in hydride
      DOUBLE PRECISION:: min_rho = 0., N, was, is,  gone      !minimal core size, number of particles in powder, balance quantities
      INTEGER, PARAMETER:: notC = 1, left = 1+notC            !notC variables beside distribution of c, c distribution between left and right
      INTEGER:: right, xNOP                                   !xNOP points in spatial grid
      DOUBLE PRECISION:: molar_mass = 830.357                 !g/mole, see http://ru.webqc.org/molecular-weight-of-La2MgNi9.html 
      DOUBLE PRECISION:: stoichiometric = 0.015               !(0.11) !H/molecule !stoich. concentration in hydride; real c can be less.
      DOUBLE PRECISION:: density = 6.22d0                     !g/cm^3 ! value from the IE, for hydride La2MgNi9H13
      DOUBLE PRECISION, PARAMETER:: Avo = 6.d23               !molecules/mole, the Avogadro number
      DOUBLE PRECISION, PARAMETER:: minV = 0.74               !voltage cut-off
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: r0, U0    !grid and distribution

      CONTAINS                                                !subroutines

      SUBROUTINE READ_DATA                                    !sub for reading exper. data from a file
      IMPLICIT NONE
         INTEGER, DIMENSION(ncurrents), PARAMETER:: lengs = [355,172,112,81,62,124,185,90,95]   !size of experimental data files
         CHARACTER(len=64):: filename = "./curve#"            !filename template
         INTEGER:: i
         threshold = 0.5                                      !threshold to separate good results from bad ones
         leng = lengs(I_index)                                !get size of the file
         write(filename(8:8),'(I1.1)') I_index                !construct the file name
         open(unit=2,file=filename,status='old',access='sequential',form='formatted')  !open the file for formatted sequential reading
         allocate(ex_data(leng,2))                            !allocate the array for data
         do i = 1,leng                                        !read data pairs
               read(2,*) ex_data(i,:)                         !one by one
         end do
         close(2)                                             !close the file
         CONTINUE
      END SUBROUTINE

      SUBROUTINE INIT                                         !sub for preparing everything before calculation begins
      IMPLICIT NONE
         DOUBLE PRECISION,PARAMETER:: D_ig = 4.d-8, tafela_ig = 1.3780319040402054, L_ig = 2.0D-3     !initial guess of values
         DOUBLE PRECISION,PARAMETER:: minR0 = 0.80, maxR0 = 0.99 !span
         DOUBLE PRECISION,PARAMETER:: minst = 10., maxst = 13.   !span
         DOUBLE PRECISION,PARAMETER:: minD  = D_ig*0.1,  maxD = D_ig*3. !min, max for diffusivity
         DOUBLE PRECISION,PARAMETER:: minL  = L_ig*0.8,  maxL = L_ig*1.2 !min, max for the radius
         DOUBLE PRECISION,PARAMETER:: minTafela  = tafela_ig*0.5,  maxTafela = tafela_ig*2. !min, max for the Tafel coef
         INTEGER:: i
         xNOP = 100                                              !number of spatial points
         n_unknowns = xNOP+notC                                  !number of variables: grid distribution plus other variables 
         right = n_unknowns                                      !rightmost index of distribution
         I_index = int_params(1); if(I_index.eq.0) I_index=5     !choose the curve
         I_ = currents(I_index) * 1.d-3                          !get the current that produced that curve
         D  = linear(minD, maxD,  params(1))                     !get the diffusivity
         C2U=  3.0E+021                                          !this is constant: Fermi equation constant, converts H concentration to voltage
         c0_div_cb=  1.E-2                                       !constant: H concentration in metal equilibrium wrt hydride
         save_tafel_a= linear(minTafela,maxTafela, params(2))    !choose tafel constant A (V = C/F + A - B ln(I))
         tafelb= 0.1                                             !tafel constant B is constant
         L= linear(minL,maxL, params(3))                         !particle's radius
         rho0 = linear(minR0,maxR0, params(4)) * L               !initial radius of the hydride core
         min_rho = 0.01 * L                                      !minimal core to be taken into account
         stoichiometric= linear(minst,maxst, params(5))          !stoichiometric concentration in hydride: real C can be less.
       OPEN(42, file='diagn.out', status='replace')              !open diagnostic file
       write(42, '("Parameters are ", 5(ES16.6,1X), 3X, I3)') D,SAVE_TAFEL_A, L, rho0/L, STOICHIOMETRIC, I_index !initial report
! I_index = 5
! I_= 0.14799999999999999
! D=  4.2248058186732596E-008
! C0_DIV_CB=  9.9999997764825821E-003
! C2U=  3.0000000000000000E+021
! SAVE_TAFEL_A=  1.3780319040402054
! TAFELB= 0.10000000149011612
! L=  2.1859480123086784E-003
! rho0= 0.97872818016507612*L
! STOICHIOMETRIC=  10.234377413367133
         N = 1.d0 / ( density * Vol(L) )                         !amount of powder particles per 1 gramm
         J = I_ / N / electron  / Surf(L)                        !flux density of charges from hydride surface to electrolyte: proportional to the current
         cb = density / molar_mass * Avo * stoichiometric !H/cm^3!H concentration in hydride
         time_limit = cb * Vol(L) / J / Surf(L)                  !time limit: how much time is needed to drive out all charges, counting 1 charge per one H atom
         dt = 10.0d0                                             !time step
         NOP = ceiling(time_limit / dt)                          !number of time steps
         calph = c0_div_cb*cb                                    !equilibrium concentration of H in metal
         gone = 0.                                               !at first no H has left the particle
         allocate(r0(xNOP), U0(xNOP))                            !allocate the grid
         r0 = rho0 + [(i, i=0,xNOP-1)]/dble(xNOP-1)*(L-rho0)     !create the grid nodes
         U0 = steady(r0)                                         !evaluate the steady distribution
         was = Amount(r0, U0);                                   !how much hydrogen was in the particle at t=0
         CONTINUE
      END SUBROUTINE

      SUBROUTINE DONE                                            !finalising sub
      IMPLICIT NONE
      DOUBLE PRECISION:: disb                                    !disbalance: violation of the conservation
         if(associated(ex_data)) deallocate(ex_data)             !free allocated data
         disb = (was-is-gone)/was*100.                           !evaluate disbalance
         if (abs(disb) > 5.) write(42, *) 'BEWARE!!!'            !write a warning if too high
         write(42, '("Matter disbalance is ", F6.2, "%")') disb  !report the disbalance
         close(42)                                               !report the amounts of was, is now, and gone/arrived
         CONTINUE                                                !close files
      END SUBROUTINE

      PURE FUNCTION IndependentVar(t)                            !sub for evaluating the independent variable if it is not time
      IMPLICIT NONE
         DOUBLE PRECISION, INTENT(IN):: t                        !time variable, current value
         DOUBLE PRECISION:: IndependentVar
         DOUBLE PRECISION, PARAMETER:: sec_per_hour = 3600., milli = 1d3
         IndependentVar = t * I_ / sec_per_hour * milli          !independent variable is not time but charge in mA hr
      END FUNCTION

      ELEMENTAL FUNCTION Steady(r)                               !sub for the steady distribution: in a point or in many points
      IMPLICIT NONE
         DOUBLE PRECISION, INTENT(IN):: r                        !point in real coord or many points
         DOUBLE PRECISION:: Steady
         SELECT CASE(sigma)                                      !shape-dependent:
         CASE(2)                                                 !sphere
                Steady = calph - J*L**2/D * (1./rho0 - 1./r)     !hydeprbola that satisfies the BC: c(rho)=calph, Dc'(L)=-J
         CASE(1)                                                 !cylinder
                Steady = calph - J*L/D * log(r / rho0)           !log curve that satisfies the BC
         CASE(0)                                                 !plate
                Steady = calph - J/D * (r - rho0)                !line that satisfies the BC
         END SELECT
         Steady = max(Steady,0.)                                 !repair negative values. if any
      END FUNCTION

      PURE FUNCTION Quantity(t,U)                                !sub for evaluating the quantity to be compared with mesurements
      IMPLICIT NONE
         DOUBLE PRECISION, INTENT(IN):: t                        !time
         DOUBLE PRECISION, INTENT(IN), DIMENSION(:):: U          !variables including the concentration distribution
         DOUBLE PRECISION:: Quantity
         Quantity = (U(right)+steady(L))/c2U + tafel()           !H concentration at the surface converted to voltage plus the Tafel shift
      END FUNCTION

      FUNCTION EarlyQuit(t,U)                                    !some event can, possibly, stop the calculation
      IMPLICIT NONE
         DOUBLE PRECISION, INTENT(IN):: t                        !time
         DOUBLE PRECISION, INTENT(IN), DIMENSION(:):: U          !variables
         LOGICAL:: EarlyQuit
         EarlyQuit = .false.                                     !generally, do not quit
         EarlyQuit = EarlyQuit .or. (Quantity(t,U) < minV)       !but stop if reached the cut-off value of the voltage
         EarlyQuit = EarlyQuit .or. (U(right)+steady(L) .le. 0.) !and quit if the battery is totally discharged
      END FUNCTION

      ELEMENTAL FUNCTION Vol(r)                                  !volume of the particle (or similar bodies)
      IMPLICIT NONE
        DOUBLE PRECISION, INTENT(IN):: r
        DOUBLE PRECISION:: Vol
         Vol  = Surf(r)*r/(sigma+1.)                             !common formula for balls, cylinders, and plates
      END FUNCTION

      ELEMENTAL FUNCTION Surf(r)                                 !surface of the particle (or similar bodies)
      IMPLICIT NONE                                                                                           
        DOUBLE PRECISION, INTENT(IN):: r                         !radius in this sub!!!
        DOUBLE PRECISION:: Surf                                                                               
         select case(sigma)                                      !depending on what body it is
         case(2)                                                                                              
                 Surf = 4.*pi*r**2                               !sphere
         case(1)                                                                                              
                 Surf = 2.*pi*r                                  !cylinder (of unit length)
         case(0)                                                                                              
                 Surf = 1.                                       !plate (of unit surface)
         end select 
      END FUNCTION

      PURE FUNCTION Amount(r, U)                                 !amount of hydrogen in the particle at any time
      USE Integration, ONLY: trapz                               !use the integrator
      IMPLICIT NONE
        DOUBLE PRECISION, INTENT(IN), DIMENSION(:):: r, U        !grid distribution 
        DOUBLE PRECISION:: Amount
        DOUBLE PRECISION:: rho, L                                !radii, taken from the grid
        INTEGER:: NOP                                            !size of the grid
        NOP = size(r)
        rho = r(1); L = r(NOP)
        Amount = cb * Vol(rho) + trapz(r, U*Surf(r))             !H in the hydride core plus integral over the metal phase
      END FUNCTION

 SUBROUTINE NextStep2layer(t,t1,U)           ! User provided 2-layer proceeder: takes current time, previous time, and previos state, returns current state
   USE Parabolic1D, ONLY: derx,derx1, LRlin, DIRICHLET_I, NEUMANN_II  !need some stuff from the PDE solver module
   use Integration, ONLY: trapz                                  !need an integrator
   IMPLICIT NONE
     DOUBLE PRECISION, INTENT(IN):: t,t1                         !the current time and previous time
     DOUBLE PRECISION, INTENT(INOUT), DIMENSION(:), TARGET:: U   !all variables inluding the concentration distribution
     DOUBLE PRECISION:: rho, realJ, drho                         !core radius = free boundary, flux density, velocity of the free boundary
     DOUBLE PRECISION, DIMENSION(xNOP):: r                       !spatial grid
     INTEGER:: i, ind
     rho = U(1) + rho0; rho = max(rho,0.)                        !get the core radius, repair if negative
     U(left:right) = U(left:right) + U0                          !get the concentration distribution; note that initial condition is zero, thus we add the steady one
     WHERE(U(left:right) < 0) U(left:right) = 0.0                !repair negative values, if any
     tafela = save_tafel_a - U(right) / C2U                      !get the Tafel constant A
     if(rho>min_rho) then                                        !if there still is a hydride core 
        drho = D*derx(U(left:right),1, L, rho) / (cb-U(left))    !evaluate velocity of the free boundary
        rho  = rho + dt*drho                                     !Euler step: new free boundary position
        rho  = max(rho,0.)                                       !repair neg
        call LRlin(dt, U(left:right), A,B,B1,E,F, DIRICHLET_I, LBC, NEUMANN_II, RBC, L, rho, 0.0d0, drho) !call the PDE solver
     else                                                        !if no more hydride
        drho = 0.0d0                                             !no free boundary
        rho = 0.0d0                                              !no hydride
        call LRlin(dt, U(left:right), A,B,B1,E,F, NEUMANN_II, LBC2, NEUMANN_II, RBC, L, rho, 0.0d0, drho) !same solver but other BC
     end if
     WHERE(U(left:right) < 0) U(left:right) = 0.0                !repair neg
     realJ = -D*derx1(U(left:right),xNOP,L,rho)                  !real flux can be less than the nominal one: if C becomes zero
     gone = gone + Surf(L)*realJ*dt                              !how much H has left the particle during the step
     r = rho + [(i, i=0,xNOP-1)]/dble(xNOP-1)*(L-rho)            !spatial grid
     is = Amount(r,U(left:right))                                !current amount
     U(left:right) = U(left:right) - U0                          !solution of the problem with zero initial state
     U(1) = rho - rho0
   RETURN

   CONTAINS                                                      !PDE coefficients

      PURE FUNCTION A(i)                                         !AC''
       IMPLICIT NONE
         DOUBLE PRECISION:: A
         INTEGER, INTENT(IN):: i
         A = D                                                  !simply diffusivity
      END FUNCTION

      PURE FUNCTION B(i)                    !+BC'
       IMPLICIT NONE
         DOUBLE PRECISION:: B
         INTEGER, INTENT(IN):: i
         DOUBLE PRECISION:: ri, h                               
         h = (L-rho)/dble(xNOP-1)          !spatial step
         ri = rho + dble(i-1)*h            !spatial point in original coordinates
         B = D * sigma / ri                !\frac{2D}{r} for sphere, e.g.
      END FUNCTION

      PURE FUNCTION B1(i)                  !-BC', to get negative B
       IMPLICIT NONE                                                
         DOUBLE PRECISION:: B1                                      
         INTEGER, INTENT(IN):: i           !index of the grid point
         B1 = 0.0d0                        !nothing in our case
      END FUNCTION

      PURE FUNCTION E(i)                   !+EC
       IMPLICIT NONE
         DOUBLE PRECISION:: E
         INTEGER, INTENT(IN):: i
         E = 0.0d0
      END FUNCTION

      PURE FUNCTION F(i)                  !+F - free term
       IMPLICIT NONE
         DOUBLE PRECISION:: F
         INTEGER, INTENT(IN):: i
         F = 0.0d0
      END FUNCTION

      PURE FUNCTION LBC()                !righthand side of the Dirichlet BC on the left side, i.e., on the free boundary
       IMPLICIT NONE
         DOUBLE PRECISION:: LBC
         LBC = calph                    !C=C_a
      END FUNCTION

      PURE FUNCTION LBC2()             !righthand side of the Neumann BC on the left side, i.e. in the centre of symmetry
       IMPLICIT NONE
         DOUBLE PRECISION:: LBC2
         LBC2 = 0.0d0                 !C'=0
      END FUNCTION

      PURE FUNCTION RBC()             !righthand side of the Neumann BC on the right side, i.e. on the surface
       IMPLICIT NONE
         DOUBLE PRECISION:: RBC
         RBC = -J/D*(L-rho)          !DC'=-J=const
      END FUNCTION

 END SUBROUTINE

      PURE FUNCTION tafel(units)                                        !sub for the Tafel formula: evaluate the voltage drop due to current
      IMPLICIT NONE
       character(len=3), intent(in), optional:: units                   !unit of current can be specified, Amper by default
       DOUBLE PRECISION:: tafel
       tafel = tafela - tafelb * log(I_ * 1.e+3)                        !default behaviour, A; coefs and the current I are global: available in the unit
       if(present(units)) then                                          !unit conversion necessary
               select case(units)
                case('A')                                               !current in Ampers
                        tafel = tafela - tafelb * log(I_ * 1.e+3)       !default behaviour, A
                case('mA')
                        tafel = tafela - tafelb * log(I_ * 1.)          !milliA
                end select
       end if
      END FUNCTION tafel

      FUNCTION norm(t,x)                   !function for the norm of the experimental curve to norm the difference between model and measured data
      IMPLICIT NONE                                                                                                                                  
         DOUBLE PRECISION:: norm           !use this only if none of the pre-designed norms suit you
         DOUBLE PRECISION, INTENT(IN), DIMENSION(:):: t,x  !grid and data
         norm = 1.0                                        !no norming by default
         CONTINUE                                                                                                                                    
      END FUNCTION                                                                                                                                   
                                                                                                                                                     
      SUBROUTINE Diagnost()                !diagnostic sub. Do whatever you want to learn what is happening
        implicit none                                                                                                                                
        return                                                                                                                                       
      END SUBROUTINE                                                                                                                                 
                                                                                                                                                     
      FUNCTION LSQ(t,x)                    !function for the norm of difference between measured and model curves                                   
      IMPLICIT NONE                                                                                                                                  
         DOUBLE PRECISION:: LSQ            !use this only if none of the pre-designed norms suit you
         DOUBLE PRECISION, INTENT(IN), DIMENSION(:):: t,x !grid and data                                                                                           
         LSQ = 0.0                         !replace this with real calculations if necessary
         CONTINUE                                                                                                                                    
      END FUNCTION                                                                                                                                   
         
      FUNCTION linear(min_value, max_value, value01)              !linear interpolator: to choose a point between M1 and M2 given an m from [0,1]
      IMPLICIT NONE
         DOUBLE PRECISION:: linear
         DOUBLE PRECISION, INTENT(IN):: min_value, max_value, value01
         linear = min_value + value01 * (max_value - min_value)
         return
      END FUNCTION

      END MODULE User


!1 Value is 0.3717 Parameters: 0.5025 0.4127 0.5005 0.2709 0.4360 2.0000

