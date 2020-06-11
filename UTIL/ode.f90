
      MODULE ODE

      IMPLICIT NONE

      PUBLIC PredCorr, Euler, PredCorr1step, Euler1step

      INTERFACE ! 
       SUBROUTINE RightHandSide(t, U, dU)
         DOUBLE PRECISION, INTENT(IN):: t
         DOUBLE PRECISION, INTENT(IN), DIMENSION(:):: U
         DOUBLE PRECISION, INTENT(OUT), DIMENSION(:):: dU
       END SUBROUTINE RightHandSide
      END INTERFACE


      CONTAINS
      
      SUBROUTINE PredCorr(tv, Y, RHS)
          DOUBLE PRECISION, DIMENSION(0:), INTENT(OUT):: tv
          DOUBLE PRECISION, DIMENSION(0:,:), INTENT(INOUT):: Y
          PROCEDURE(RightHandSide):: RHS
          DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:):: dU, dU_
          DOUBLE PRECISION:: dt
          INTEGER:: n, NOP, n_unknowns
          NOP = size(tv)-1;
          n_unknowns = size(Y,dim=2)
          allocate(dU(1:n_unknowns),dU_(1:n_unknowns))
          tv = 0.d0
          dt = 1. / dble(NOP)
!          Y = 0.d0
          do n = 0,NOP-1
                  tv(n+1) = tv(n) + dt
                  call RHS(tv(n), Y(n,:), dU)
                  call RHS(tv(n+1), Y(n,:) + dt*dU, dU_)
                  Y(n+1,:) = Y(n,:) + 0.5 * dt * (dU+dU_)
          end do
          deallocate(dU,dU_)
       END SUBROUTINE

      SUBROUTINE PredCorr1step(t, dt, Y, RHS)
          DOUBLE PRECISION, INTENT(IN):: t, dt
          DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT):: Y
          PROCEDURE(RightHandSide):: RHS
          DOUBLE PRECISION, DIMENSION(size(Y)):: dU, dU_
          INTEGER:: n_unknowns
          n_unknowns = size(Y)
!          Y = 0.d0
          call RHS(t, Y, dU)
          call RHS(t+dt, Y + dt*dU, dU_)
          Y = Y + 0.5 * dt * (dU+dU_)
       END SUBROUTINE

      SUBROUTINE Euler(tv, Y, RHS)
          DOUBLE PRECISION, DIMENSION(0:), INTENT(OUT):: tv
          DOUBLE PRECISION, DIMENSION(0:,:), INTENT(INOUT):: Y
          PROCEDURE(RightHandSide):: RHS
          DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:):: dU
          DOUBLE PRECISION:: dt
          INTEGER:: n, NOP, n_unknowns
          NOP = size(tv)-1;
          n_unknowns = size(Y,dim=2)
          allocate(dU(1:n_unknowns))
          tv = 0.d0
          dt = 1. / dble(NOP)
!          Y = 0.d0
          do n = 0,NOP-1
                  tv(n+1) = tv(n) + dt
                  call RHS(tv(n), Y(n,:), dU)
                  Y(n+1,:) = Y(n,:) + dU*dt
          end do
          deallocate(dU)
       END SUBROUTINE

      SUBROUTINE Euler1step(t, dt, Y, RHS)
          DOUBLE PRECISION, INTENT(IN):: t, dt
          DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT):: Y
          PROCEDURE(RightHandSide):: RHS
          DOUBLE PRECISION, DIMENSION(1:size(Y)):: dU
          INTEGER:: n_unknowns
          n_unknowns = size(Y)
!          Y = 0.d0
          call RHS(t, Y(:), dU)
          Y(:) = Y(:) + dU*dt
       END SUBROUTINE

      END MODULE ODE

