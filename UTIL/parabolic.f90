      MODULE Parabolic1D                                !solver for 1d parabolic PDE with different BC; IC always zero
        USE equations                                   !uses algebraic equation solvers

      IMPLICIT NONE

      PUBLIC:: derx, derx1, &
              &  LRlin, LLinRlin3, LLin3Rlin,      &
            &  LRlin3,  LLinRnonlin, LLin3Rnonlin, &
            &  LnonlinRlin, LnonlinRlin3, LRnonlin, &
            &  DLRlin,       DLLinRlin3,    DLLin3Rlin,    &
            &  DLRlin3,      DLLinRnonlin,  DLLin3Rnonlin, &
            &  DLnonlinRlin, DLnonlinRlin3, DLRnonlin,     &
            &  DIRICHLET_I, NEUMANN_II, NEUMANN_III, MAXITER

      INTEGER, PARAMETER:: DIRICHLET_I = 1, NEUMANN_II = 2, NEUMANN_III = 3  !codes of BC type
      INTEGER, PARAMETER:: MAXITER = 1000000

     !interface blocks of functions to be passed as arguments

      INTERFACE                                    !RHS of linear BC: c=f or c'=f or even f1*c+f2*c'=0
       PURE FUNCTION RightHandSide()
         DOUBLE PRECISION:: RightHandSide
       END FUNCTION RightHandSide
      END INTERFACE

      INTERFACE ! 
       PURE FUNCTION NonlinRightHandSide(c)     !RHS of nonlinear BC: c'=f(c)
         DOUBLE PRECISION, INTENT(IN):: c
         DOUBLE PRECISION:: NonlinRightHandSide
       END FUNCTION NonlinRightHandSide
      END INTERFACE

      INTERFACE ! 
       PURE FUNCTION PDE_Coefficient(i)       !coefficients of the PDE
         INTEGER, INTENT(IN):: i              !index of the spatial grid point
         DOUBLE PRECISION:: PDE_Coefficient
       END FUNCTION PDE_Coefficient
      END INTERFACE

      CONTAINS

      PURE FUNCTION derx(U,i, L, rho, dL, drho )           !evaluates right first spatial grid derivative
        IMPLICIT NONE
        DOUBLE PRECISION:: derx
        DOUBLE PRECISION:: h                     !spatial grid step
        DOUBLE PRECISION, DIMENSION(:), INTENT(IN):: U       !array of the grid function
        DOUBLE PRECISION, INTENT(IN):: L, rho
        DOUBLE PRECISION, INTENT(IN), OPTIONAL:: dL, drho 
        INTEGER, INTENT(IN):: i                              !index of the point where der is calculated
        INTEGER:: xNOP, ii
        DOUBLE PRECISION:: LR, LL, dLR, dLL
          LR = L; LL = rho; dLR = 0.; dLL = 0.;
          if(present(dL))  dLR = dL
          if(present(drho)) dLL = drho
          xNOP = size(U)
          h = (LR-LL)/dble(xNOP-1)
          ii = max(1,i)
          if(i<xNOP) then
                  derx = ( U(ii+1)-U(ii) )
          else
                  derx = ( U(xNOP)-U(xNOP-1) )
          end if
          derx = derx / h
      END FUNCTION

      PURE FUNCTION derx1(U,i, L, rho, dL, drho )           !evaluates left first spatial grid derivative
        IMPLICIT NONE
        DOUBLE PRECISION:: derx1
        DOUBLE PRECISION:: h                      !spatial grid step
        DOUBLE PRECISION, DIMENSION(:), INTENT(IN):: U        !array of the grid function
        DOUBLE PRECISION, INTENT(IN):: L, rho
        DOUBLE PRECISION, INTENT(IN), OPTIONAL:: dL, drho 
        INTEGER, INTENT(IN):: i                               !index of the point where der is calculated
        INTEGER:: xNOP, ii
        DOUBLE PRECISION:: LR, LL, dLR, dLL
          dLR = 0.; dLL = 0.;
          LR = L
          LL = rho
          if(present(dL))  dLR = dL
          if(present(drho)) dLL = drho
          xNOP = size(U)
          h = (LR-LL)/dble(xNOP-1)
          ii = min(xNOP,i)
          if(i>1) then
                  derx1 = ( U(ii)-U(ii-1) )
          else
                  derx1 = ( U(2)-U(1) )
          end if
          derx1 = derx1 / h
      END FUNCTION

      
      !Solves the parabolic 1D equation using implicit scheme; both boundary conditions are either of the I or the II kind
      SUBROUTINE LRlin(dt, U, A,B,B1,E,F, LBCtype, LBC, RBCtype, RBC, L, rho, dL, drho )
          DOUBLE PRECISION, INTENT(IN):: dt
          DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT):: U
          PROCEDURE(PDE_Coefficient):: A,B,B1,E,F
          PROCEDURE(RightHandSide):: LBC, RBC
          INTEGER, INTENT(IN):: LBCtype, RBCtype
          DOUBLE PRECISION:: h
          DOUBLE PRECISION, INTENT(IN):: L, rho
          DOUBLE PRECISION, INTENT(IN), OPTIONAL:: dL, drho 
          DOUBLE PRECISION, DIMENSION(size(U)):: sweepa,sweepb
          DOUBLE PRECISION:: th, thh, den
          DOUBLE PRECISION:: Ai, Bi, B1i, LR, LL, dLR, dLL, add, xi, factor
          INTEGER:: NOP, i
          NOP = size(U)
          h = 1./dble(NOP-1)
          th = dt/h; thh = th/h
          dLR = 0.; dLL = 0.;
          LR = L
          LL = rho
          if(present(dL))   dLR = dL
          if(present(drho)) dLL = drho
          factor = 1. / (LR-LL)
          SELECT CASE(LBCtype)
          CASE(DIRICHLET_I)    ! boundary condition is linear: c = RHS
              sweepa(1) = 0.; sweepb(1) = LBC()
          CASE(NEUMANN_II)     ! boundary condition is linear: Dc' = BC_RHS2(time,parameters)
              sweepa(1) = 1.; sweepb(1) = -LBC()*h/factor
          END SELECT
          do i = 2,NOP-1
             Ai = A(i) *factor**2; Bi = B(i)*factor; B1i = B1(i)*factor;
             xi = LL + dble(i-1)*h
             add = (dLL*(1.-xi)+dLR*xi) *factor
             if(add > 0.) Bi = Bi + add
             if(add < 0.) B1i = B1i - add 
             den =  1. - Ai*thh*(sweepa(i-1) - 2.) + Bi*th + B1i*th*(1. - sweepa(i-1)) - dt*E(i)
             sweepa(i) = (Ai*thh + Bi*th) / den
             sweepb(i) = (Ai*thh*sweepb(i-1) + B1i*th*sweepb(i-1) + F(i)*dt + U(i)) / den
          end do
          SELECT CASE(RBCtype)
          CASE(DIRICHLET_I)    ! boundary condition is linear: c = BC_RHS1(time,parameters)
              U(NOP) = RBC()
          CASE(NEUMANN_II)     ! boundary condition is linear: Dc' = BC_RHS2(time,parameters)
              U(NOP) = (sweepb(NOP-1) + h*RBC()/factor) / (1 - sweepa(NOP-1))
          END SELECT
          do i = NOP-1,1,-1
             U(i) = sweepa(i)*U(i+1) + sweepb(i)
          end do
       END SUBROUTINE
      
      SUBROUTINE LLinRlin3(dt, U, A,B,B1,E,F, LBCtype, LBC, RBC, RBC1, L, rho, dL, drho  ) !c'=rbc*c+rbc1
          DOUBLE PRECISION, INTENT(IN):: dt
          DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT):: U
          PROCEDURE(PDE_Coefficient):: A,B,B1,E,F
          PROCEDURE(RightHandSide):: LBC, RBC, RBC1
          INTEGER, INTENT(IN):: LBCtype
          DOUBLE PRECISION:: h
          DOUBLE PRECISION, INTENT(IN):: L, rho
          DOUBLE PRECISION, INTENT(IN), OPTIONAL:: dL, drho 
          DOUBLE PRECISION, DIMENSION(size(U)):: sweepa,sweepb
          DOUBLE PRECISION::  th, thh, den, U1, UI
          DOUBLE PRECISION:: Ai, Bi, B1i, LR, LL, dLR, dLL, add, xi, factor
          INTEGER:: NOP, i
          NOP = size(U)
          h = 1./dble(NOP-1)
          th = dt/h; thh = th/h
          dLR = 0.; dLL = 0.;
          LR = L
          LL = rho
          if(present(dL))  dLR = dL
          if(present(drho)) dLL = drho
          factor = 1. / (LR-LL)
          SELECT CASE(LBCtype)
          CASE(DIRICHLET_I)    ! boundary condition is linear: c = RHS
              sweepa(1) = 0.; sweepb(1) = LBC()
          CASE(NEUMANN_II)     ! boundary condition is linear: Dc' = BC_RHS2(time,parameters)
              sweepa(1) = 1.; sweepb(1) = -LBC()*h/factor
          END SELECT
          do i = 2,NOP-1
             Ai = A(i) *factor**2; Bi = B(i)*factor; B1i = B1(i)*factor;
             xi = LL + dble(i-1)*h
             add = (dLL*(1.-xi)+dLR*xi) *factor
             if(add > 0.) Bi = Bi + add
             if(add < 0.) B1i = B1i - add 
             den = 1. - Ai*thh*(sweepa(i-1) - 2.) + Bi*th + B1i*th*(1. - sweepa(i-1)) - dt*E(i)
             sweepa(i) = (Ai*thh + Bi*th) / den
             sweepb(i) = (Ai*thh*sweepb(i-1) + B1i*th*sweepb(i-1) + F(i)*dt + U(i)) / den
          end do
          U(NOP) = (sweepb(I-1) + h*RBC1()/factor) / (1 - sweepa(NOP-1) - h*RBC()/factor)
          do i = NOP-1,1,-1
             U(i) = sweepa(i)*U(i+1) + sweepb(i)
          end do
       END SUBROUTINE
      
      SUBROUTINE LLin3Rlin(dt, U, A,B,B1,E,F, LBC, LBC1, RBCtype, RBC, L, rho, dL, drho  ) !c'=rbc*c+rbc1
          DOUBLE PRECISION, INTENT(IN):: dt
          DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT):: U
          PROCEDURE(PDE_Coefficient):: A,B,B1,E,F
          PROCEDURE(RightHandSide):: LBC, LBC1, RBC
          INTEGER, INTENT(IN):: RBCtype
          DOUBLE PRECISION:: h
          DOUBLE PRECISION, INTENT(IN):: L, rho
          DOUBLE PRECISION, INTENT(IN), OPTIONAL:: dL, drho 
          DOUBLE PRECISION, DIMENSION(size(U)):: sweepa,sweepb
          DOUBLE PRECISION::  th, thh, den, U1, UI
          DOUBLE PRECISION:: Ai, Bi, B1i, LR, LL, dLR, dLL, add, xi, factor
          INTEGER:: NOP, i
          NOP = size(U)
          h = 1./dble(NOP-1)
          th = dt/h; thh = th/h
          dLR = 0.; dLL = 0.;
          LR = L
          LL = rho
          if(present(dL))  dLR = dL
          if(present(drho)) dLL = drho
          factor = 1. / (LR-LL)
          den = 1. + h*LBC()/factor
          sweepa(1) = 1. / den; sweepb(1) = -LBC1()*h/factor / den
          do i = 2,NOP-1
             Ai = A(i) *factor**2; Bi = B(i)*factor; B1i = B1(i)*factor;
             xi = LL + dble(i-1)*h
             add = (dLL*(1.-xi)+dLR*xi) *factor
             if(add > 0.) Bi = Bi + add
             if(add < 0.) B1i = B1i - add 
             den = 1. - Ai*thh*(sweepa(i-1) - 2.) + Bi*th + B1i*th*(1. - sweepa(i-1)) - dt*E(i)
             sweepa(i) = (Ai*thh + Bi*th) / den
             sweepb(i) = (Ai*thh*sweepb(i-1) + B1i*th*sweepb(i-1) + F(i)*dt + U(i)) / den
          end do
          SELECT CASE(RBCtype)
          CASE(DIRICHLET_I)    ! boundary condition is linear: c = BC_RHS1(time,parameters)
              U(NOP) = RBC()
          CASE(NEUMANN_II)     ! boundary condition is linear: Dc' = BC_RHS2(time,parameters)
              U(NOP) = (sweepb(I-1) + h*RBC()/factor) / (1 - sweepa(NOP-1))
          END SELECT
          do i = NOP-1,1,-1
             U(i) = sweepa(i)*U(i+1) + sweepb(i)
          end do
       END SUBROUTINE

      SUBROUTINE LRlin3(dt, U, A,B,B1,E,F, LBC, LBC1, RBC, RBC1, L, rho, dL, drho  )
          DOUBLE PRECISION, INTENT(IN):: dt
          DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT):: U
          PROCEDURE(PDE_Coefficient):: A,B,B1,E,F
          PROCEDURE(RightHandSide):: LBC, LBC1, RBC, RBC1
          DOUBLE PRECISION:: h
          DOUBLE PRECISION, INTENT(IN):: L, rho
          DOUBLE PRECISION, INTENT(IN), OPTIONAL:: dL, drho 
          DOUBLE PRECISION, DIMENSION(size(U)):: sweepa,sweepb
          DOUBLE PRECISION:: th, thh, den, U1, UI
          DOUBLE PRECISION:: Ai, Bi, B1i, LR, LL, dLR, dLL, add, xi, factor
          INTEGER:: NOP, i
          NOP = size(U)
          h = 1./dble(NOP-1)
          dLR = 0.; dLL = 0.;
          LR = L
          LL = rho
          if(present(dL))   dLR = dL
          if(present(drho)) dLL = drho
          factor = 1. / (LR-LL)
          th = dt/h; thh = th/h
          den = 1. + h*LBC()/factor
          sweepa(1) = 1. / den; sweepb(1) = -LBC1()*h/factor / den
          do i = 2,NOP-1
             Ai = A(i) *factor**2; Bi = B(i)*factor; B1i = B1(i)*factor;
             xi = LL + dble(i-1)*h
             add = (dLL*(1.-xi)+dLR*xi) *factor
             if(add > 0.) Bi = Bi + add
             if(add < 0.) B1i = B1i - add 
             den = 1. - Ai*thh*(sweepa(i-1) - 2.) + Bi*th + B1i*th*(1. - sweepa(i-1)) - dt*E(i)
             sweepa(i) = (Ai*thh + Bi*th) / den
             sweepb(i) = (Ai*thh*sweepb(i-1) + B1i*th*sweepb(i-1) + F(i)*dt + U(i)) / den
          end do
          U(NOP) = (sweepb(I-1) + h*RBC1()/factor) / (1 - sweepa(NOP-1) - h*RBC()/factor)
          do i = NOP-1,1,-1
             U(i) = sweepa(i)*U(i+1) + sweepb(i)
          end do
       END SUBROUTINE
      
      SUBROUTINE LLinRnonlin(dt, U, A,B,B1,E,F, LBCtype, LBC, RBC, L, rho, dL, drho  ) !c'=rbc(c)
      USE Equations
      IMPLICIT NONE
          DOUBLE PRECISION, INTENT(IN):: dt
          DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT):: U
          PROCEDURE(PDE_Coefficient):: A,B,B1,E,F
          PROCEDURE(RightHandSide):: LBC
          PROCEDURE(NonlinRightHandSide):: RBC
          INTEGER, INTENT(IN):: LBCtype
          DOUBLE PRECISION:: h
          DOUBLE PRECISION, INTENT(IN):: L, rho
          DOUBLE PRECISION, INTENT(IN), OPTIONAL:: dL, drho 
          DOUBLE PRECISION, DIMENSION(size(U)):: sweepa,sweepb
          DOUBLE PRECISION::  th, thh, den, U1, UI
          DOUBLE PRECISION:: Ai, Bi, B1i, LR, LL, dLR, dLL, add, xi, factor
          INTEGER:: NOP, i, ierr
          NOP = size(U)
          h = 1./dble(NOP-1)
          dLR = 0.; dLL = 0.;
          LR = L
          LL = rho
          if(present(dL))  dLR = dL
          if(present(drho)) dLL = drho
          factor = 1. / (LR-LL)
          th = dt/h; thh = th/h
          SELECT CASE(LBCtype)
          CASE(DIRICHLET_I)    ! boundary condition is linear: c = RHS
              sweepa(1) = 0.; sweepb(1) = LBC()
          CASE(NEUMANN_II)     ! boundary condition is linear: c' = BC_RHS2(time,parameters)
              sweepa(1) = 1.; sweepb(1) = -LBC()*h/factor
          END SELECT
          do i = 2,NOP-1
             Ai = A(i) *factor**2; Bi = B(i)*factor; B1i = B1(i)*factor;
             xi = dble(i-1)*h
             add = (dLL*(1.-xi)+dLR*xi) *factor
             if(add > 0.) Bi = Bi + add
             if(add < 0.) B1i = B1i - add 
             den = 1. - Ai*thh*(sweepa(i-1) - 2.) + Bi*th + B1i*th*(1. - sweepa(i-1)) - dt*E(i)
             sweepa(i) = (Ai*thh + Bi*th) / den
             sweepb(i) = (Ai*thh*sweepb(i-1) + B1i*th*sweepb(i-1) + F(i)*dt + U(i)) / den
          end do
          caLL dichotomy(RP, U(NOP), ierr, 0.D0)
          if(ierr .ne. SOLUTION_FOUND) then
                  print*, 'Something wrong in solving RP(x)=0, ierr = ', ierr
                  stop
          end if
          do i = NOP-1,1,-1
             U(i) = sweepa(i)*U(i+1) + sweepb(i)
          end do
      CONTAINS
                  PURE FUNCTION RP(X)
                          DOUBLE PRECISION:: RP
                          DOUBLE PRECISION, INTENT(IN):: X
                          RP = (1.-sweepa(NOP-1))*X - sweepb(NOP-1) - h*RBC(X)/factor
                  END FUNCTION RP
       END SUBROUTINE
      
      SUBROUTINE LLin3Rnonlin(dt, U, A,B,B1,E,F, LBC, LBC1, RBC, L, rho, dL, drho ) 
          DOUBLE PRECISION, INTENT(IN):: dt
          DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT):: U
          PROCEDURE(PDE_Coefficient):: A,B,B1,E,F
          PROCEDURE(RightHandSide):: LBC, LBC1
          PROCEDURE(NonlinRightHandSide):: RBC
          DOUBLE PRECISION:: h
          DOUBLE PRECISION, INTENT(IN):: L, rho
          DOUBLE PRECISION, INTENT(IN), OPTIONAL:: dL, drho 
          DOUBLE PRECISION, DIMENSION(size(U)):: sweepa,sweepb
          DOUBLE PRECISION:: th, thh, den, U1, UI
          DOUBLE PRECISION:: Ai, Bi, B1i, LR, LL, dLR, dLL, add, xi, factor
          INTEGER:: NOP, i, ierr
          NOP = size(U)
          h = 1./dble(NOP-1)
          dLR = 0.; dLL = 0.;
          LR = L
          LL = rho
          if(present(dL))  dLR = dL
          if(present(drho)) dLL = drho
          factor = 1. / (LR-LL)
          th = dt/h; thh = th/h
          den = 1. + h*LBC()/factor
          sweepa(1) = 1. / den; sweepb(1) = -LBC1()*h/factor / den
          do i = 2,NOP-1
             Ai = A(i) * factor**2; Bi = B(i)*factor; B1i = B1(i)*factor;
             xi = LL + dble(i-1)*h
             add = (dLL*(1.-xi)+dLR*xi) * factor
             if(add > 0.) Bi = Bi + add
             if(add < 0.) B1i = B1i - add 
             den = 1. - Ai*thh*(sweepa(i-1) - 2.) + Bi*th + B1i*th*(1. - sweepa(i-1)) - dt*E(i)
             sweepa(i) = (Ai*thh + Bi*th) / den
             sweepb(i) = (Ai*thh*sweepb(i-1) + B1i*th*sweepb(i-1) + F(i)*dt + U(i)) / den
          end do
          caLL dichotomy(RP, U(NOP), ierr, 0.D0)
          if(ierr .ne. SOLUTION_FOUND) then
                  print*, 'Something wrong in solving RP(x)=0, ierr = ', ierr
                  stop
          end if
          do i = NOP-1,1,-1
             U(i) = sweepa(i)*U(i+1) + sweepb(i)
          end do
       CONTAINS
                  PURE FUNCTION RP(X)
                          DOUBLE PRECISION:: RP
                          DOUBLE PRECISION, INTENT(IN):: X
                          RP = (1.-sweepa(NOP-1))*X - sweepb(NOP-1) - h*RBC(X) / factor
                  END FUNCTION RP
       END SUBROUTINE

      SUBROUTINE LNonlinRlin(dt, U, A,B,B1,E,F, LBC, RBCtype, RBC, L, rho, dL, drho  ) 
          DOUBLE PRECISION, INTENT(IN):: dt
          DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT):: U
          PROCEDURE(PDE_Coefficient):: A,B,B1,E,F
          PROCEDURE(RightHandSide):: RBC
          PROCEDURE(NonlinRightHandSide):: LBC
          INTEGER, INTENT(IN):: RBCtype
          DOUBLE PRECISION:: h
          DOUBLE PRECISION, INTENT(IN):: L, rho
          DOUBLE PRECISION, INTENT(IN), OPTIONAL:: dL, drho 
          DOUBLE PRECISION, DIMENSION(size(U)):: sweepa,sweepb
          DOUBLE PRECISION::  th, thh, den
          DOUBLE PRECISION:: Ai, Bi, B1i, LR, LL, dLR, dLL, add, xi, factor
          INTEGER:: NOP, i, ierr
          NOP = size(U)
          h = 1./dble(NOP-1)
          dLR = 0.; dLL = 0.;
          LR = L
          LL = rho
          if(present(dL))   dLR = dL
          if(present(drho)) dLL = drho
          factor = 1. / (LR-LL)
          th = dt/h; thh = th/h
          SELECT CASE(RBCtype)
          CASE(DIRICHLET_I)    ! boundary condition is linear: c = RHS
              sweepa(NOP) = 0.; sweepb(NOP) = RBC()
          CASE(NEUMANN_II)     ! boundary condition is linear: Dc' = BC_RHS2(time,parameters)
              sweepa(NOP) = 1.; sweepb(NOP) = -RBC()*h/factor
          END SELECT
          do i = NOP-1,2,-1
             Ai = A(i) *factor**2; Bi = B(i)*factor; B1i = B1(i)*factor;
             xi = LL + dble(i-1)*h
             add = (dLL*(1.-xi)+dLR*xi) *factor
             if(add > 0.) Bi = Bi + add
             if(add < 0.) B1i = B1i - add 
             den = 1. - Ai*thh*(sweepa(i+1) - 2.) + Bi*th*(1. - sweepa(i+1)) + B1i*th - dt*E(i)
             sweepa(i) = (Ai*thh + B1i*th) / den
             sweepb(i) = (Ai*thh*sweepb(i+1) + Bi*th*sweepb(i+1) + F(i)*dt + U(i)) / den
          end do
          caLL dichotomy(LP, U(1), ierr, 0.D0)
          if(ierr .ne. SOLUTION_FOUND) then
                  print*, 'Something wrong in solving LP(x)=0, ierr = ', ierr
                  stop
          end if
          do i = 2,NOP
             U(i) = sweepa(i)*U(i-1) + sweepb(i)
          end do
       CONTAINS
                  PURE FUNCTION LP(X)
                          DOUBLE PRECISION:: LP
                          DOUBLE PRECISION, INTENT(IN):: X
                          LP = (1.-sweepa(2))*X - sweepb(2) - h*LBC(X)/factor
                  END FUNCTION LP
       END SUBROUTINE
      
      SUBROUTINE LNonlinRlin3(dt, U, A,B,B1,E,F, LBC, RBC, RBC1, L, rho, dL, drho ) 
          DOUBLE PRECISION, INTENT(IN):: dt
          DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT):: U
          PROCEDURE(PDE_Coefficient):: A,B,B1,E,F
          PROCEDURE(RightHandSide):: RBC, RBC1
          PROCEDURE(NonlinRightHandSide):: LBC
          DOUBLE PRECISION:: h
          DOUBLE PRECISION, INTENT(IN):: L, rho
          DOUBLE PRECISION, INTENT(IN), OPTIONAL:: dL, drho 
          DOUBLE PRECISION, DIMENSION(size(U)):: sweepa,sweepb
          DOUBLE PRECISION:: th, thh, den, U1, UI
          DOUBLE PRECISION:: Ai, Bi, B1i, LR, LL, dLR, dLL, add, xi, factor
          INTEGER:: NOP, i, ierr
          NOP = size(U)
          h = 1./dble(NOP-1)
          dLR = 0.; dLL = 0.;
          LR = L
          LL = rho
          if(present(dL))  dLR = dL
          if(present(drho)) dLL = drho
          factor = 1. / (LR-LL)
          th = dt/h; thh = th/h
          den = 1. + h*RBC()/factor
          sweepa(NOP) = 1. / den; sweepb(NOP) = -RBC1()*h/factor / den
          do i = NOP-1,2,-1
             Ai = A(i) *factor**2; Bi = B(i)*factor; B1i = B1(i)*factor;
             xi = LL + dble(i-1)*h
             add = (dLL*(1.-xi)+dLR*xi) *factor
             if(add > 0.) Bi = Bi + add
             if(add < 0.) B1i = B1i - add 
             den = 1. - Ai*thh*(sweepa(i+1) - 2.) + Bi*th*(1. - sweepa(i+1)) + B1i*th - dt*E(i)
             sweepa(i) = (Ai*thh + B1i*th) / den
             sweepb(i) = (Ai*thh*sweepb(i+1) + Bi*th*sweepb(i+1) + F(i)*dt + U(i)) / den
          end do
          caLL dichotomy(LP, U(1), ierr, 0.D0)
          if(ierr .ne. SOLUTION_FOUND) then
                  print*, 'Something wrong in solving LP(x)=0, ierr = ', ierr
                  stop
          end if
          do i = 2,NOP
             U(i) = sweepa(i)*U(i-1) + sweepb(i)
          end do
       CONTAINS
                  PURE FUNCTION LP(X)
                          DOUBLE PRECISION:: LP
                          DOUBLE PRECISION, INTENT(IN):: X
                          LP = (1.-sweepa(2))*X - sweepb(2) - h*LBC(X)/factor
                  END FUNCTION LP
       END SUBROUTINE

      SUBROUTINE LRnonlin(dt, U, A,B,B1,E,F, LBC, RBC, L, rho, dL, drho ) 
          DOUBLE PRECISION, INTENT(IN):: dt
          DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT):: U
          PROCEDURE(PDE_Coefficient):: A,B,B1,E,F
          PROCEDURE(NonlinRightHandSide):: RBC, LBC
          DOUBLE PRECISION:: h
          DOUBLE PRECISION, INTENT(IN):: L, rho
          DOUBLE PRECISION, INTENT(IN), OPTIONAL:: dL, drho 
          DOUBLE PRECISION, DIMENSION(2,size(U)):: sweepa,sweepb,sweepg
          DOUBLE PRECISION:: th, thh, den, U1, UI
          DOUBLE PRECISION:: Ai, Bi, B1i, LR, LL, dLR, dLL, add, xi, factor
          INTEGER:: NOP, i,j, ierr, iter
          NOP = size(U)
          h = 1./dble(NOP-1)
          dLR = 0.; dLL = 0.;
          LR = L
          LL = rho
          if(present(dL))  dLR = dL
          if(present(drho)) dLL = drho
          factor = 1. / (LR-LL)
          th = dt/h; thh = th/h
          sweepa(1,1)   = 0.; sweepb(1,1)   = 0.; sweepg(1,1)   = 1.
          sweepa(2,NOP) = 0.; sweepb(2,NOP) = 0.; sweepg(2,NOP) = 1.
          do i = 2,NOP-1
             Ai = A(i) *factor**2; Bi = B(i)*factor; B1i = B1(i)*factor;
             xi = LL + dble(i-1)*h
             add = (dLL*(1.-xi)+dLR*xi) *factor
             if(add > 0.) Bi = Bi + add
             if(add < 0.) B1i = B1i - add 
             den = 1. - Ai*thh*(sweepa(1,i-1) - 2.) + Bi*th + B1i*th*(1. - sweepa(1,i-1)) - dt*E(i)
             sweepa(1,i) = (Ai*thh + Bi*th) / den
             sweepb(1,i) = (Ai*thh*sweepb(1,i-1) + B1i*th*sweepb(1,i-1) + F(i)*dt + U(i)) / den
             sweepg(1,i) = sweepg(1,i-1) / den
             j = NOP-i+1
             Ai = A(j) *factor**2; Bi = B(j)*factor; B1i = B1(j)*factor;
             xi = LL + dble(j-1)*h
             add = (dLL*(1.-xi)+dLR*xi) *factor
             den = 1. - Ai*thh*(sweepa(2,j+1) - 2.) + Bi*th*(1. - sweepa(2,j+1)) + B1i*th - dt*E(j)
             sweepa(2,j) = (Ai*thh + B1i*th) / den
             sweepb(2,j) = (Ai*thh*sweepb(2,j+1) + Bi*th*sweepb(2,j+1) + F(j)*dt + U(j)) / den
             sweepg(2,j) = sweepg(2,j+1) / den
          end do
          caLL dichotomy(GL, U(1), ierr, 0.D0)
          if(ierr .ne. SOLUTION_FOUND) then
                  print*, 'Something wrong in solving GL(x)=0, ierr = ', ierr
                  stop
          end if
          caLL dichotomy(GR, U(NOP), ierr, 0.D0)
          if(ierr .ne. SOLUTION_FOUND) then
                  print*, 'Something wrong in solving GR(Y)=0, ierr = ', ierr
                  stop
          end if
          do iter = 1,MAXITER
            U1     = (sweepb(2,2) + sweepg(2,2)*U(NOP) + h*LBC(U(1))/factor) / (1.-sweepa(2,2))
            UI = (sweepb(1,NOP-1) + sweepg(1,NOP-1)*U(1) + h*RBC(U(NOP))/factor) / (1.-sweepa(1,NOP-1))
            if((abs(UI-U(NOP))+abs(U1-U(1)))*0.5 < tol) then
                    U(1) = U1
                    U(NOP) = UI
                    exit
            end if
            U(1) = U1
            U(NOP) = UI
          end do
          do i = 2,NOP
             U(i) = sweepa(2,i)*U(i+1) + sweepb(2,i)
          end do
       CONTAINS
                  PURE FUNCTION GL(X)
                          DOUBLE PRECISION:: GL
                          DOUBLE PRECISION, INTENT(IN):: X
                          GL = (1.-sweepa(2,2))*X - sweepb(2,2) - sweepg(2,2)*0.0d0  - h*LBC(X)/factor
                  END FUNCTION GL
                  PURE FUNCTION GR(Y)
                          DOUBLE PRECISION:: GR
                          DOUBLE PRECISION, INTENT(IN):: Y
                          GR = (1.-sweepa(1,NOP-1))*Y - sweepb(1,NOP-1) - sweepg(1,NOP-1)*0.0d0  - h*RBC(Y)/factor
                  END FUNCTION GR
       END SUBROUTINE

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! SIMPLE FICK EQUATION IN CURVILINEAR COODINATES !!!!!!!!!!!!!!!!!!!!!!!

      !Solves the parabolic 1D equation using implicit scheme; both boundary conditions are either of the I or the II kind
      SUBROUTINE DLRlin(dt, U, D, curva, LBCtype, LBC, RBCtype, RBC, L, rho, dL, drho )
              DOUBLE PRECISION, INTENT(IN):: dt
              DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT):: U
              DOUBLE PRECISION, INTENT(IN):: D                   !diffusivity
              DOUBLE PRECISION, DIMENSION(3), INTENT(IN):: curva !three coefficients for second order polynomial of surface wrt distance r
              PROCEDURE(RightHandSide):: LBC, RBC
              INTEGER, INTENT(IN):: LBCtype, RBCtype
              INTEGER:: xNOP
              DOUBLE PRECISION, INTENT(IN):: L, rho
              DOUBLE PRECISION, INTENT(IN), OPTIONAL:: dL, drho 
              xNOP = size(U)
              call LRlin(dt, U, A,B,ZEROF,ZEROF,ZEROF, LBCtype, LBC, RBCtype, RBC, L, rho, dL, drho )

      CONTAINS                                                          !PDE coefficients to be given to PDE solvers

              PURE FUNCTION A(i)                                             !Au''
                      IMPLICIT NONE
                      DOUBLE PRECISION:: A
                      INTEGER, INTENT(IN):: i                                     !index of the grid point
                      A = D                                                       !simply the diffusivity in our case
              END FUNCTION

              PURE FUNCTION B(i)                                             !Bu'
                      IMPLICIT NONE
                      DOUBLE PRECISION:: B
                      INTEGER, INTENT(IN):: i                                     !index of the grid point 
                      DOUBLE PRECISION:: ri, nom, den                                       !distance from the centre to the current point
                      ri = rho + (dble(i-1)/dble(xNOP-1))*(L-rho)                            !spatial point in original coordinates
                      nom = 2.*curva(1)*ri + curva(2)
                      den = curva(1)*ri*ri + curva(2)*ri + curva(3) 
                      B = D * nom / den                                         !\frac{2D}{r} for sphere, e.g.
              END FUNCTION

      END SUBROUTINE
      
      SUBROUTINE DLLinRlin3(dt, U, D, curva, LBCtype, LBC, RBC, RBC1, L, rho, dL, drho  ) !c'=rbc*c+rbc1
          DOUBLE PRECISION, INTENT(IN):: dt
          DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT):: U
              DOUBLE PRECISION, INTENT(IN):: D                   !diffusivity
              DOUBLE PRECISION, DIMENSION(3), INTENT(IN):: curva !three coefficients for second order polynomial of surface wrt distance r
          PROCEDURE(RightHandSide):: LBC, RBC, RBC1
          INTEGER, INTENT(IN):: LBCtype
          DOUBLE PRECISION:: h
          DOUBLE PRECISION, INTENT(IN):: L, rho
          DOUBLE PRECISION, INTENT(IN), OPTIONAL:: dL, drho 
          INTEGER:: xNOP
          xNOP = size(U)
          call LLinRlin3(dt, U, A,B,ZEROF,ZEROF,ZEROF, LBCtype, LBC, RBC, RBC1, L, rho, dL, drho  ) !c'=rbc*c+rbc1

      CONTAINS                                                          !PDE coefficients to be given to PDE solvers

              PURE FUNCTION A(i)                                             !Au''
                      IMPLICIT NONE
                      DOUBLE PRECISION:: A
                      INTEGER, INTENT(IN):: i                                     !index of the grid point
                      A = D                                                       !simply the diffusivity in our case
              END FUNCTION

              PURE FUNCTION B(i)                                             !Bu'
                      IMPLICIT NONE
                      DOUBLE PRECISION:: B
                      INTEGER, INTENT(IN):: i                                     !index of the grid point 
                      DOUBLE PRECISION:: ri, nom, den                                       !distance from the centre to the current point
                      ri = rho + (dble(i-1)/dble(xNOP-1))*(L-rho)                            !spatial point in original coordinates
                      nom = 2.*curva(1)*ri + curva(2)
                      den = curva(1)*ri*ri + curva(2)*ri + curva(3) 
                      B = D * nom / den                                         !\frac{2D}{r} for sphere, e.g.
              END FUNCTION

       END SUBROUTINE
      
      SUBROUTINE DLLin3Rlin(dt, U, D, curva, LBC, LBC1, RBCtype, RBC, L, rho, dL, drho  ) !c'=rbc*c+rbc1
          DOUBLE PRECISION, INTENT(IN):: dt
          DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT):: U
              DOUBLE PRECISION, INTENT(IN):: D                   !diffusivity
              DOUBLE PRECISION, DIMENSION(3), INTENT(IN):: curva !three coefficients for second order polynomial of surface wrt distance r
          PROCEDURE(RightHandSide):: LBC, LBC1, RBC
          INTEGER, INTENT(IN):: RBCtype
          DOUBLE PRECISION:: h
          DOUBLE PRECISION, INTENT(IN):: L, rho
          DOUBLE PRECISION, INTENT(IN), OPTIONAL:: dL, drho 
          INTEGER:: xNOP
          xNOP = size(U)
          call LLin3Rlin(dt, U, A,B,ZEROF,ZEROF,ZEROF, LBC, LBC1, RBCtype, RBC, L, rho, dL, drho  ) !c'=rbc*c+rbc1

      CONTAINS                                                          !PDE coefficients to be given to PDE solvers

              PURE FUNCTION A(i)                                             !Au''
                      IMPLICIT NONE
                      DOUBLE PRECISION:: A
                      INTEGER, INTENT(IN):: i                                     !index of the grid point
                      A = D                                                       !simply the diffusivity in our case
              END FUNCTION

              PURE FUNCTION B(i)                                             !Bu'
                      IMPLICIT NONE
                      DOUBLE PRECISION:: B
                      INTEGER, INTENT(IN):: i                                     !index of the grid point 
                      DOUBLE PRECISION:: ri, nom, den                                       !distance from the centre to the current point
                      ri = rho + (dble(i-1)/dble(xNOP-1))*(L-rho)                            !spatial point in original coordinates
                      nom = 2.*curva(1)*ri + curva(2)
                      den = curva(1)*ri*ri + curva(2)*ri + curva(3) 
                      B = D * nom / den                                         !\frac{2D}{r} for sphere, e.g.
              END FUNCTION

       END SUBROUTINE

      SUBROUTINE DLRlin3(dt, U, D, curva, LBC, LBC1, RBC, RBC1, L, rho, dL, drho  )
          DOUBLE PRECISION, INTENT(IN):: dt
          DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT):: U
              DOUBLE PRECISION, INTENT(IN):: D                   !diffusivity
              DOUBLE PRECISION, DIMENSION(3), INTENT(IN):: curva !three coefficients for second order polynomial of surface wrt distance r
          PROCEDURE(RightHandSide):: LBC, LBC1, RBC, RBC1
          DOUBLE PRECISION:: h
          DOUBLE PRECISION, INTENT(IN):: L, rho
          DOUBLE PRECISION, INTENT(IN), OPTIONAL:: dL, drho 
          INTEGER:: xNOP
          xNOP = size(U)
          call LRlin3(dt, U, A,B,ZEROF,ZEROF,ZEROF, LBC, LBC1, RBC, RBC1, L, rho, dL, drho  )

      CONTAINS                                                          !PDE coefficients to be given to PDE solvers

              PURE FUNCTION A(i)                                             !Au''
                      IMPLICIT NONE
                      DOUBLE PRECISION:: A
                      INTEGER, INTENT(IN):: i                                     !index of the grid point
                      A = D                                                       !simply the diffusivity in our case
              END FUNCTION

              PURE FUNCTION B(i)                                             !Bu'
                      IMPLICIT NONE
                      DOUBLE PRECISION:: B
                      INTEGER, INTENT(IN):: i                                     !index of the grid point 
                      DOUBLE PRECISION:: ri, nom, den                                       !distance from the centre to the current point
                      ri = rho + (dble(i-1)/dble(xNOP-1))*(L-rho)                            !spatial point in original coordinates
                      nom = 2.*curva(1)*ri + curva(2)
                      den = curva(1)*ri*ri + curva(2)*ri + curva(3) 
                      B = D * nom / den                                         !\frac{2D}{r} for sphere, e.g.
              END FUNCTION

       END SUBROUTINE
      
      SUBROUTINE DLLinRnonlin(dt, U, D, curva, LBCtype, LBC, RBC, L, rho, dL, drho  ) !c'=rbc(c)
      USE Equations
      IMPLICIT NONE
          DOUBLE PRECISION, INTENT(IN):: dt
          DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT):: U
              DOUBLE PRECISION, INTENT(IN):: D                   !diffusivity
              DOUBLE PRECISION, DIMENSION(3), INTENT(IN):: curva !three coefficients for second order polynomial of surface wrt distance r
          PROCEDURE(RightHandSide):: LBC
          PROCEDURE(NonlinRightHandSide):: RBC
          INTEGER, INTENT(IN):: LBCtype
          DOUBLE PRECISION:: h
          DOUBLE PRECISION, INTENT(IN):: L, rho
          DOUBLE PRECISION, INTENT(IN), OPTIONAL:: dL, drho 
          INTEGER:: xNOP
          xNOP = size(U)
          call LLinRnonlin(dt, U, A,B,ZEROF,ZEROF,ZEROF, LBCtype, LBC, RBC, L, rho, dL, drho  ) !c'=rbc(c)

      CONTAINS                                                          !PDE coefficients to be given to PDE solvers

              PURE FUNCTION A(i)                                             !Au''
                      IMPLICIT NONE
                      DOUBLE PRECISION:: A
                      INTEGER, INTENT(IN):: i                                     !index of the grid point
                      A = D                                                       !simply the diffusivity in our case
              END FUNCTION

              PURE FUNCTION B(i)                                             !Bu'
                      IMPLICIT NONE
                      DOUBLE PRECISION:: B
                      INTEGER, INTENT(IN):: i                                     !index of the grid point 
                      DOUBLE PRECISION:: ri, nom, den                                       !distance from the centre to the current point
                      ri = rho + (dble(i-1)/dble(xNOP-1))*(L-rho)                            !spatial point in original coordinates
                      nom = 2.*curva(1)*ri + curva(2)
                      den = curva(1)*ri*ri + curva(2)*ri + curva(3) 
                      B = D * nom / den                                         !\frac{2D}{r} for sphere, e.g.
              END FUNCTION

       END SUBROUTINE
      
      SUBROUTINE DLLin3Rnonlin(dt, U, D, curva, LBC, LBC1, RBC, L, rho, dL, drho ) 
          DOUBLE PRECISION, INTENT(IN):: dt
          DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT):: U
              DOUBLE PRECISION, INTENT(IN):: D                   !diffusivity
              DOUBLE PRECISION, DIMENSION(3), INTENT(IN):: curva !three coefficients for second order polynomial of surface wrt distance r
          PROCEDURE(RightHandSide):: LBC, LBC1
          PROCEDURE(NonlinRightHandSide):: RBC
          DOUBLE PRECISION:: h
          DOUBLE PRECISION, INTENT(IN):: L, rho
          DOUBLE PRECISION, INTENT(IN), OPTIONAL:: dL, drho 
          INTEGER:: xNOP
          xNOP = size(U)
          call LLin3Rnonlin(dt, U, A,B,ZEROF,ZEROF,ZEROF, LBC, LBC1, RBC, L, rho, dL, drho ) 

      CONTAINS                                                          !PDE coefficients to be given to PDE solvers

              PURE FUNCTION A(i)                                             !Au''
                      IMPLICIT NONE
                      DOUBLE PRECISION:: A
                      INTEGER, INTENT(IN):: i                                     !index of the grid point
                      A = D                                                       !simply the diffusivity in our case
              END FUNCTION

              PURE FUNCTION B(i)                                             !Bu'
                      IMPLICIT NONE
                      DOUBLE PRECISION:: B
                      INTEGER, INTENT(IN):: i                                     !index of the grid point 
                      DOUBLE PRECISION:: ri, nom, den                                       !distance from the centre to the current point
                      ri = rho + (dble(i-1)/dble(xNOP-1))*(L-rho)                            !spatial point in original coordinates
                      nom = 2.*curva(1)*ri + curva(2)
                      den = curva(1)*ri*ri + curva(2)*ri + curva(3) 
                      B = D * nom / den                                         !\frac{2D}{r} for sphere, e.g.
              END FUNCTION

       END SUBROUTINE

      SUBROUTINE DLNonlinRlin(dt, U, D, curva, LBC, RBCtype, RBC, L, rho, dL, drho  ) 
          DOUBLE PRECISION, INTENT(IN):: dt
          DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT):: U
              DOUBLE PRECISION, INTENT(IN):: D                   !diffusivity
              DOUBLE PRECISION, DIMENSION(3), INTENT(IN):: curva !three coefficients for second order polynomial of surface wrt distance r
          PROCEDURE(RightHandSide):: RBC
          PROCEDURE(NonlinRightHandSide):: LBC
          INTEGER, INTENT(IN):: RBCtype
          DOUBLE PRECISION:: h
          DOUBLE PRECISION, INTENT(IN):: L, rho
          DOUBLE PRECISION, INTENT(IN), OPTIONAL:: dL, drho 
          INTEGER:: xNOP
          xNOP = size(U)
          call LNonlinRlin(dt, U, A,B,ZEROF,ZEROF,ZEROF, LBC, RBCtype, RBC, L, rho, dL, drho  ) 

      CONTAINS                                                          !PDE coefficients to be given to PDE solvers

              PURE FUNCTION A(i)                                             !Au''
                      IMPLICIT NONE
                      DOUBLE PRECISION:: A
                      INTEGER, INTENT(IN):: i                                     !index of the grid point
                      A = D                                                       !simply the diffusivity in our case
              END FUNCTION

              PURE FUNCTION B(i)                                             !Bu'
                      IMPLICIT NONE
                      DOUBLE PRECISION:: B
                      INTEGER, INTENT(IN):: i                                     !index of the grid point 
                      DOUBLE PRECISION:: ri, nom, den                                       !distance from the centre to the current point
                      ri = rho + (dble(i-1)/dble(xNOP-1))*(L-rho)                            !spatial point in original coordinates
                      nom = 2.*curva(1)*ri + curva(2)
                      den = curva(1)*ri*ri + curva(2)*ri + curva(3) 
                      B = D * nom / den                                         !\frac{2D}{r} for sphere, e.g.
              END FUNCTION

       END SUBROUTINE
      
      SUBROUTINE DLNonlinRlin3(dt, U, D, curva, LBC, RBC, RBC1, L, rho, dL, drho ) 
          DOUBLE PRECISION, INTENT(IN):: dt
          DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT):: U
              DOUBLE PRECISION, INTENT(IN):: D                   !diffusivity
              DOUBLE PRECISION, DIMENSION(3), INTENT(IN):: curva !three coefficients for second order polynomial of surface wrt distance r
          PROCEDURE(RightHandSide):: RBC, RBC1
          PROCEDURE(NonlinRightHandSide):: LBC
          DOUBLE PRECISION:: h
          DOUBLE PRECISION, INTENT(IN):: L, rho
          DOUBLE PRECISION, INTENT(IN), OPTIONAL:: dL, drho 
          INTEGER:: xNOP
          xNOP = size(U)
          call LNonlinRlin3(dt, U, A,B,ZEROF,ZEROF,ZEROF, LBC, RBC, RBC1, L, rho, dL, drho ) 

      CONTAINS                                                          !PDE coefficients to be given to PDE solvers

              PURE FUNCTION A(i)                                             !Au''
                      IMPLICIT NONE
                      DOUBLE PRECISION:: A
                      INTEGER, INTENT(IN):: i                                     !index of the grid point
                      A = D                                                       !simply the diffusivity in our case
              END FUNCTION

              PURE FUNCTION B(i)                                             !Bu'
                      IMPLICIT NONE
                      DOUBLE PRECISION:: B
                      INTEGER, INTENT(IN):: i                                     !index of the grid point 
                      DOUBLE PRECISION:: ri, nom, den                                       !distance from the centre to the current point
                      ri = rho + (dble(i-1)/dble(xNOP-1))*(L-rho)                            !spatial point in original coordinates
                      nom = 2.*curva(1)*ri + curva(2)
                      den = curva(1)*ri*ri + curva(2)*ri + curva(3) 
                      B = D * nom / den                                         !\frac{2D}{r} for sphere, e.g.
              END FUNCTION

       END SUBROUTINE

      SUBROUTINE DLRnonlin(dt, U, D, curva, LBC, RBC, L, rho, dL, drho ) 
          DOUBLE PRECISION, INTENT(IN):: dt
          DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT):: U
              DOUBLE PRECISION, INTENT(IN):: D                   !diffusivity
              DOUBLE PRECISION, DIMENSION(3), INTENT(IN):: curva !three coefficients for second order polynomial of surface wrt distance r
          PROCEDURE(NonlinRightHandSide):: RBC, LBC
          DOUBLE PRECISION:: h
          DOUBLE PRECISION, INTENT(IN):: L, rho
          DOUBLE PRECISION, INTENT(IN), OPTIONAL:: dL, drho 
          INTEGER:: xNOP
          xNOP = size(U)
          call LRnonlin(dt, U, A,B,ZEROF,ZEROF,ZEROF, LBC, RBC, L, rho, dL, drho ) 

      CONTAINS                                                          !PDE coefficients to be given to PDE solvers

              PURE FUNCTION A(i)                                             !Au''
                      IMPLICIT NONE
                      DOUBLE PRECISION:: A
                      INTEGER, INTENT(IN):: i                                     !index of the grid point
                      A = D                                                       !simply the diffusivity in our case
              END FUNCTION

              PURE FUNCTION B(i)                                             !Bu'
                      IMPLICIT NONE
                      DOUBLE PRECISION:: B
                      INTEGER, INTENT(IN):: i                                     !index of the grid point 
                      DOUBLE PRECISION:: ri, nom, den                                       !distance from the centre to the current point
                      ri = rho + (dble(i-1)/dble(xNOP-1))*(L-rho)                            !spatial point in original coordinates
                      nom = 2.*curva(1)*ri + curva(2)
                      den = curva(1)*ri*ri + curva(2)*ri + curva(3) 
                      B = D * nom / den                                         !\frac{2D}{r} for sphere, e.g.
              END FUNCTION

       END SUBROUTINE

              PURE FUNCTION ZEROF(i)                                             !free term
                      IMPLICIT NONE
                      DOUBLE PRECISION:: ZEROF
                      INTEGER, INTENT(IN):: i                                     !index of the grid point    
                      ZEROF = 0.0d0                                                   !nothing in our case
              END FUNCTION

      END MODULE Parabolic1D

