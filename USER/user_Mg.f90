! Solves the problem of dehydriding. Spherical particles.
! -bc**2 = Ddc/dr on the surface (P=0), equilibrium on the phaze boundary.
! No c in metal. Desorption flux density is given for both phazes.
! Quick diffusion. The skin grows both in tangent and radial directions, 
! Desorption is from the phaze boundary.

      MODULE User
      IMPLICIT NONE

      PUBLIC:: INIT, DONE, READ_DATA, Quantity, Norm, LSQ, EarlyQuit, IndependentVar, NextStep2layer, Diagnost, savenml 
      PUBLIC:: N_params, N_int_params, params, int_params
      PUBLIC:: n_unknowns, NOP, dt
      PUBLIC:: threshold, ex_data, norming_method, space               !X: this line must not be changed!
      PUBLIC:: get_params_values, PERTURBE

      DOUBLE PRECISION:: threshold = -1.

      INTEGER, PARAMETER:: N_params = 9, N_int_params = 1
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

      INTEGER, PARAMETER:: curves = 6
      DOUBLE PRECISION, DIMENSION(curves), PARAMETER:: time_limits = [14644.345, 13965.357, 23045.968, 66224.273, 89372.458, 142663.] 

      DOUBLE PRECISION, PARAMETER:: tol = 1e-6, cb = 1.,pi = 3.1415926 !tolerance constant, unit concentration in hydride
      DOUBLE PRECISION, PARAMETER:: c0 = 1d-3, ca = 1.0d0
      DOUBLE PRECISION:: Jb,Ja, R0, S0, nu, L=8.d-4, slope = 0.
      INTEGER:: sigma=2, leng, curve !shape factor, Number Of spatial Points, aux
      DOUBLE PRECISION:: sc, time_limit, Surf, Vol, time_shift, tgamma=0., gamm=0., SNi=0., JNi=0.

      NAMELIST /parameters/ NOP, curve, L, Jb, Ja, R0, S0, nu, time_shift, slope, time_limit, tgamma, gamm, SNi, JNi

      CONTAINS

      SUBROUTINE READ_DATA
         INTEGER, DIMENSION(curves) :: lengs = [ 7207, 6872, 11339, 20828, 14325, 151112]
         CHARACTER(len=64):: datafile = "./Mg.#.dat"
         INTEGER:: i
         NAMELIST /metadata/ lengs, threshold, datafile
         threshold = 0.05
         threshold = -0.05
         open(44,file='metadata.nml')
         read(44,nml=metadata)
         close(44)
         leng = lengs(curve)
         write(datafile(6:6),'(I1.1)') curve
         open(unit=2,file=datafile,status='old',access='sequential',form='formatted')    !newunit?
         allocate(ex_data(leng,2))
         open(unit=3,file='exper.dat',status='replace',access='sequential',form='formatted')    !newunit?
         do i = 1,leng
               read(2,*,end=42) ex_data(i,:)
               ex_data(i,2) = min(1., ex_data(i,2) + slope * ex_data(i,1))
               write(3,*)  ex_data(i,:), (1.-ex_data(i,2))**(1./3.)
         end do
42       close(2)
         close(3)
         CONTINUE
      END SUBROUTINE

      SUBROUTINE INIT(nmlfile)
         CHARACTER(len=64), INTENT(IN), OPTIONAL:: nmlfile
         CHARACTER(len=64):: defnml = 'default.Mg.nml'
         DOUBLE PRECISION,PARAMETER:: maxR0 = 0.9, maxS0 = 1.0/6.0, maxnu = 1.0     !maximal values
         DOUBLE PRECISION,PARAMETER:: minR0 = 0.5, minS0 = 0.0,     minnu = 0.1     !minimal values
         DOUBLE PRECISION,PARAMETER:: mintime_shift=0., relative_max_timeshift = 0.15
         DOUBLE PRECISION:: maxtime_shift
         DOUBLE PRECISION:: bb_, ba_
         if(present(nmlfile)) then
                 open(33,file=nmlfile,status='old',access='sequential',form='formatted')
                 read(33,nml=parameters)
                 close(33)
         else
                 open(33,file=defnml,status='old',access='sequential',form='formatted')
                 read(33,nml=parameters)
                 close(33)
                 if(curve.eq.0) curve = mod(int_params(1),curves) + 1
                 Jb = Jb * 2. * params(1)                !from 0 to 2*initial_flux_dens
                 Ja = Ja * 2. * params(2)                !from 0 to 2*initial_flux_dens
                 R0 = linear(minR0,maxR0, params(3))*L
                 S0 = linear(minS0,maxS0, params(4))
                 nu = linear(minnu,maxnu, params(5))
                 time_limit = time_limits(curve)
                 maxtime_shift = time_limit * relative_max_timeshift 
                 time_shift = linear(mintime_shift,maxtime_shift, params(6))
                 slope = linear(0.0d0,0.5d0, params(7))
                 tgamma = linear(0.0d0,time_shift, params(8))
                 gamm = linear(0.0d0,time_shift**(-3), params(9))
                 SNi = 0.
                 JNi = 0.
         end if
         R0 = R0 * L
         time_limit = time_limit + time_shift
         slope = slope / time_limit
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

      SUBROUTINE savenml(nmlfile)
         CHARACTER(len=64), INTENT(IN):: nmlfile
         open(33,file=nmlfile,status='replace',access='sequential',form='formatted')
         write(33,nml=parameters)
         close(33)
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
         Sb = Sb + gamm*(tgamma-0.)**3*Vol/Surf*L**sigma/(L**(sigma+1)-R0**(sigma+1))*(sigma+1)
         if(rho <= 0.) rho = 0.                            !repair negative value
         Sb = max(min(Sb,1.),0.)                             !repair values outside [0,1]
         if(t.ge.tgamma) then 
                 Quantity = ((Sb+SNi)*Surf/L**sigma*(L**(sigma+1)-rho**(sigma+1))/(sigma+1)) / Vol !reacted fraction
         else
                 Quantity = gamm*(t-0.)**3 + ((S0+SNi)*Surf/L**sigma*(L**(sigma+1)-R0**(sigma+1))/(sigma+1)) / Vol !reacted fraction
         endif
      END FUNCTION

      PURE SUBROUTINE phase_vec(t,U,V)
         IMPLICIT NONE
         DOUBLE PRECISION, INTENT(IN):: t
         DOUBLE PRECISION, INTENT(IN), DIMENSION(:):: U
         DOUBLE PRECISION, INTENT(OUT), DIMENSION(:):: V
         V(1) = max(0., U(1) + R0) / L
         V(2) = U(2) + S0 + gamm*(tgamma-0.)**3*Vol/Surf*L**sigma/(L**(sigma+1)-R0**(sigma+1))*(sigma+1)
         V(2) = max(min(V(2),1.),0.)                             !repair values outside [0,1]
      END SUBROUTINE

      PURE FUNCTION IndependentVar(t)
         IMPLICIT NONE
         DOUBLE PRECISION, INTENT(IN):: t
         DOUBLE PRECISION:: IndependentVar
         IndependentVar = t - time_shift
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
        DOUBLE PRECISION:: rho, S, local_nu, W_, S_grows  !variables to be used
        dU = 0.0d0
        if(t<tgamma) return
        rho = U(1) + R0; S = U(2) + S0                    !the core size, relative area of the skin
        S = S + gamm*(tgamma-0.)**3*Vol/Surf*L**sigma/(L**(sigma+1)-R0**(sigma+1))*(sigma+1)
        if(rho <= 0.) rho = 0.                            !repair negative value
        S = max(min(S,1.),0.)                             !repair values outside [0,1]
        local_nu = nu                                     !shrinking core vs skin growth factor nu local to the sub 
        S_grows = 1.                                      
        if(S.ge.1-SNi) then
                local_nu = 0.                          !if skin is ready, pure shrinking core
                S_grows = 0.                                      
        endif
        W_ = -cb*(L**(sigma+1)-rho**(sigma+1))/(sigma+1)  !geometrical factor
        if(rho>0) then                                    !two phases co-exist
                dU(1) = (1.-local_nu)*Jb / cb                   !d\rho/dt: note constnt fluxes Jb and Ja
                dU(2) = ((Ja*(1.-SNi-S) + JNi*SNi)*S_grows*L**sigma + local_nu*S*Jb*rho**sigma) / W_ !dS/dt
        else                                              !single phase
                dU(1) = 0.                                !no core
                dU(2) = ((Ja*(1.-SNi-S) + JNi*SNi)*S_grows*L**sigma + local_nu*S*Jb*rho**sigma) / W_ !dS/dt
        end if
        RETURN
      END SUBROUTINE

      FUNCTION norm(x,y)
         DOUBLE PRECISION, DIMENSION(:), INTENT(IN):: x,y
         DOUBLE PRECISION:: norm
         norm = 1.0
         CONTINUE
      END FUNCTION

      SUBROUTINE Diagnost(t,U,realU)
        implicit none
         DOUBLE PRECISION, INTENT(IN):: t
         DOUBLE PRECISION, INTENT(IN), DIMENSION(:):: U
         DOUBLE PRECISION, INTENT(OUT), DIMENSION(:), OPTIONAL:: realU
         DOUBLE PRECISION:: rho, S                         !variables to be used
         LOGICAL, SAVE:: nocore = .false., skin = .false.
         rho = U(1) + R0; S = U(2) + S0                    !the core size, relative area of the skin
         S = S + gamm*(tgamma-0.)**3*Vol/Surf*L**sigma/(L**(sigma+1)-R0**(sigma+1))*(sigma+1)
         if(.not.nocore.and.rho.le.0.+1.e-19) then
                 print'("rho(",F12.4,")=0, S=",F12.4)',t,S
                 nocore = .true.
         endif
         if(.not.skin.and.S.ge.1.-SNi-1.e-19) then
                 print'("S(",F12.4,")=1, rho=",F12.4)',t,rho
                 skin = .true.
         endif
         if(t.ge.time_limit-dt/2.) print'("S(t^*)=",F12.2,", R(t^*)=",F12.2)', max(0.,min(1.,S+SNi)),max(0.,min(1.,rho/L))
         if(present(realU)) call phase_vec(t,U,realU)
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

      subroutine get_params_values(val,ival)
              double precision, dimension(:), intent(out):: val
              integer, dimension(:), intent(out):: ival
              val = [Jb, Ja, R0, S0, nu, time_shift, slope, tgamma, gamm]
              ival = [curve]
      end subroutine

      subroutine PERTURBE(dval)
              double precision, dimension(:), intent(in):: dval
              Jb = Jb + dval(1)
              Ja = Ja + dval(2)
              R0 = R0 + dval(3)
              S0 = S0 + dval(4)
              nu = nu + dval(5)
              time_shift = time_shift + dval(6)
              slope = slope + dval(7)
              tgamma = tgamma + dval(8)
              gamm = gamm + dval(9)
      end subroutine

      END MODULE User



