! An example to show how 0-layer problems are solved.
! A problem is F(p) = J, where p is a parameter vector, F is a function,
! J is a measured quantity. No time, so 0 time layers are used.

      MODULE User
      IMPLICIT NONE

      PUBLIC:: INIT, DONE, READ_DATA, Quantity, Diagnost
      PUBLIC:: N_params, N_int_params, params, int_params
      PUBLIC:: threshold, ex_data                           !X: this line must not be changed!

      DOUBLE PRECISION:: threshold = -1.

      INTEGER, PARAMETER:: N_params = 2, N_int_params = 0
      INTEGER, PARAMETER:: N_fixed_params = 0
      DOUBLE PRECISION, DIMENSION(N_params):: params
      INTEGER, DIMENSION(N_int_params):: int_params
      DOUBLE PRECISION:: ex_data

      PRIVATE

      DOUBLE PRECISION, PARAMETER:: pi = 3.1415926
      DOUBLE PRECISION:: x

      CONTAINS

      SUBROUTINE READ_DATA
         threshold = 0.01
         ex_data = 1.d0/sqrt(2.d0)
         CONTINUE
      END SUBROUTINE

      SUBROUTINE INIT
         x = (params(1)-0.5)*pi
         CONTINUE
      END SUBROUTINE

      SUBROUTINE DONE
         CONTINUE
      END SUBROUTINE

      SUBROUTINE Quantity(Q)
         IMPLICIT NONE
         DOUBLE PRECISION, INTENT(OUT):: Q
         Q = sin(x)
      END SUBROUTINE

      SUBROUTINE Diagnost()
        implicit none
        return
      END SUBROUTINE

      FUNCTION linear(min_value, max_value, value01)
         IMPLICIT NONE
         DOUBLE PRECISION:: linear
         DOUBLE PRECISION, INTENT(IN):: min_value, max_value, value01
         linear = min_value + value01 * (max_value - min_value)
         return
      END FUNCTION

      END MODULE User



