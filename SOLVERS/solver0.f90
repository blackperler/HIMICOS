
 PROGRAM TASK                                                            ! A task for BOINC grid umbrella computing
 
 USE user                                                                ! The user-provided module with parameters and right-hand side of the system of ODEs
 use interpolation                                                       ! interpolators
 use integration                                                         ! numerical integrators

 IMPLICIT NONE

 INTEGER:: npar, par, n                                                  ! loop counters
 DOUBLE PRECISION:: Q                                                    ! the compared quantity
 DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: model_data              ! model quantity interpolated
 DOUBLE PRECISION:: base, sq, h

 CALL INIT_RANDOM_SEED()                                                 !random numbers generator initialization
 
  do par = 1,N_params                                                    !random choice of the real parameters; each is from [0,1)
        call random_number(params(par))
  end do

  do npar = 1,N_int_params                                               !random choice of the integer parameters; each is 0-9
        call random_number(h)
        int_params(npar) = int(h*10)
  end do

 CALL INIT()                                                             ! Initializer (user-provided)

 call Quantity(Q)                       ! User provided 0-layer proceeder: returns the quantity
 call diagnost                                !call the diagnostic routine to let the user learn what is happening

 CALL READ_DATA()                                                        ! call the user-provided data loader so that experimental value be available

 if(abs(Q/ex_data) .le. HUGE(1.0d0)) then
         sq = abs(1. - Q/ex_data)                                                    ! the dimensionless value of interest
 else
         sq = abs(Q)
 end if

 if(threshold < 0.) then
         print 60, sq, params, real(int_params)                     ! if no threshold is given, return the value itself. !!! parameters???
         open(31, file='model.res',status='replace')
         write(31,*) Q, ex_data
         close(31)
 end if
 if(threshold > 0.) then                                                 ! if the threshold is given, 
         if(sq .le. threshold) then                                      ! and the value is better,
                 print 66, 1, sq, params, real(int_params)                                               ! return the success flag; !!! params???
                 open(31, file='model.res',status='replace')
                 write(31,*) Q, ex_data
                 close(31)
         else                                                            ! otherwise
                 print 66, 0, sq, params, real(int_params)               ! the fail flag.
         end if
 end if
 
 CALL DONE()                                                             ! allow the user to clear up

! Note that the output to STDERR always begins with:
! 1) A number (1 or 0): this means pass/fail result with respect to the user-privided threshold; then come the parameters after the word 'Parameters: ', all (including integer ones) as real numbers
! 2) The words 'Value is ': this means that the result follows an a real number; then come the parameters after the word 'Parameters: ', all (including integer ones) as real numbers
! 3) The word 'ERROR: ': this means that an error has occured. Its description follows.
60   FORMAT ("Value is ",F12.4,1X,"Parameters: ", *(F6.4,1X))
66   FORMAT (I1.1,1X,"Value is ",F6.4,1X,"Parameters: ", *(F6.4,1X))

 CONTAINS

 SUBROUTINE init_random_seed()
 IMPLICIT NONE
 INTEGER, ALLOCATABLE :: seed(:)
 INTEGER :: i, n, un, istat, dt(8), pid, t(2), s
 INTEGER(8) :: count, tms

 call random_seed(size = n)
 allocate(seed(n))
 ! First try if the OS provides a random number generator
 open(newunit=un, file="/dev/urandom", access="stream", &
 form="unformatted", action="read", status="old", iostat=istat)
 if (istat == 0) then
         read(un) seed
         close(un)
 else
         ! Fallback to XOR:ing the current time and pid. The PID is
         ! useful in case one launches multiple instances of the same
         ! program in parallel.
         call system_clock(count)
         if (count /= 0) then
                 t = transfer(count, t)
         else
                 call date_and_time(values=dt)
                 tms = (dt(1) - 1970) * 365_8 * 24 * 60 * 60 * 1000 &
                 + dt(2) * 31_8 * 24 * 60 * 60 * 1000 &
                 + dt(3) * 24 * 60 * 60 * 60 * 1000 &
                 + dt(5) * 60 * 60 * 1000 &
                 + dt(6) * 60 * 1000 + dt(7) * 1000 &
                 + dt(8)
                 t = transfer(tms, t)
         end if
         s = ieor(t(1), t(2))
         pid = getpid() + 1099279 ! Add a prime
         s = ieor(s, pid)
         if (n >= 3) then
                 seed(1) = t(1) + 36269
                 seed(2) = t(2) + 72551
                 seed(3) = pid
                 if (n > 3) then
                         seed(4:) = s + 37 * (/ (i, i = 0, n - 4) /)
                 end if
         else
                 seed = s + 37 * (/ (i, i = 0, n - 1 ) /)
         end if
 end if
 call random_seed(put=seed)
 END SUBROUTINE init_random_seed

 END
