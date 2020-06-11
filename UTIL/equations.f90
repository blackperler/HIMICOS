MODULE Equations
        IMPLICIT NONE

        PUBLIC:: dichotomy, SOLUTION_FOUND, WRONG_SEGMENT, tol, tolX

        INTEGER, PARAMETER:: SOLUTION_FOUND = 0,WRONG_SEGMENT = 1
        DOUBLE PRECISION:: tol = 1.d-9, tolX = 1.d-9

      INTERFACE ! 
       PURE FUNCTION RHS(X)
         DOUBLE PRECISION:: RHS
         DOUBLE PRECISION,INTENT(IN):: X
       END FUNCTION RHS
      END INTERFACE

        CONTAINS

                SUBROUTINE dichotomy(func,X,ierr,left,right)
                !PURE SUBROUTINE dichotomy(func,X,ierr,left,right)
                        double precision, intent(in), optional:: left, right
                        procedure(RHS):: func
                        integer, intent(out):: ierr
                        double precision, intent(out):: X
                        double precision:: L,R,M, FL, FR
                        ierr = SOLUTION_FOUND 
                        if(present(left) .and. present(right)) then
                                L = left; R = right;
                        end if
                        if(present(left) .and. .not. present(right)) then
                                L = left; R = max(L,1.);
                                FL = func(L)
                                do while(func(R)*FL.gt.0.)
                                  R = 10.*R
                                end do
                        end if
                        if(.not.present(left) .and. present(right)) then
                                R = right; L = min(R,-1.);
                                FR = func(R)
                                do while(func(L)*FR.gt.0.)
                                  L = 10.*L
                                end do
                        end if
                        if(.not.present(left) .and. .not.present(right)) then
                                R = 1.; L = -1.;
                                do while(func(L)*func(R).gt.0.)
                                  L = 10.*L; R = 10.*R;
                                end do
                        end if
                        FL = func(L); FR = func(R); M = (L+R)*0.5d0
                        do while (abs(FL+FR)*0.5d0 > tol .and. abs(R-L) > tolX)
                                FL = func(L); FR = func(R); 
                                if(FL*FR > 0) then
                                        ierr = WRONG_SEGMENT
                                        return
                                end if
                                M = (L+R) * 0.5d0
                                if(FL*func(M).le.0.) R = M
                                if(func(M)*FR.le.0.) L = M
                        end do
                        X = M
                END SUBROUTINE


END MODULE Equations
