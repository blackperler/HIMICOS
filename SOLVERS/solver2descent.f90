
 PROGRAM TASK                                                            ! A task for BOINC grid umbrella computing
 
 USE user                                                                ! The user-provided module with parameters and right-hand side of the system of ODEs
 use interpolation                                                       ! interpolators
 use integration                                                         ! numerical integrators
 use random                                                              ! random number generator initialization

 IMPLICIT NONE

 INTEGER:: npar, par, i                                                  ! loop counters
 DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE::   tv, f, wa, x0           ! time array, model-measurements difference vector, stuff for the descent
 DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: Q                       ! the compared quantity array (x-y columns)
 DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, TARGET:: Y                 ! solution array
 DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: model_data                ! model quantity interpolated
 INTEGER:: arsize, m, n, info, lwa                                       ! size of the ex_data array
 INTEGER, DIMENSION(:), ALLOCATABLE::   iwa                              ! special stuff for descent routine 
 DOUBLE PRECISION, PARAMETER:: tol = 1e-9!tolerance constant
 DOUBLE PRECISION:: base, sq, h

 CALL INIT_RANDOM_SEED()                                                 !random numbers generator initialization
 
  do par = 1,N_params                                                    !random choice of the real parameters; each is from [0,1)
        call random_number(params(par))
  end do

  do npar = 1,N_int_params                                               !random choice of the integer parameters; each is 0-9
        call random_number(h)
        int_params(npar) = int(h*10)
  end do
  allocate(x0(N_params))                                                  !allocate the array for the init values of the parameters
  x0 = params
  CALL INIT()                                                             !Initializer (user-provided)
  n = N_params
  CALL READ_DATA()                                                        !call the user-provided data loader so that experimental data be available
  arsize = size(ex_data, dim=1)                                           !get the size of the experimental data array
  m = arsize
  allocate(model_data(arsize))                                          !allocate the array for the model data
  allocate(tv(0:NOP), Q(0:NOP,2), f(0:NOP))                               !allocate arrays for the time grid a nd 2-column model data (time-value or smth of the sort)
  allocate(Y(n_unknowns))                                                 !allocate the phase vector: all variables at any time will be here
  lwa = m*n+5*n+m
  if(allocated(wa)) deallocate(wa,iwa)
  ALLOCATE(wa(lwa), iwa(n))
  print'("initial x is ",*(F9.3,1X))', params(:)
  call lmdif1(fcn,m,n,x0,f,tol,info,iwa,wa,lwa)
  select case(info)
     case(0)
             print*, 'improper input parameters.'
     case(1)
             print*, 'the relative error in the sum of squares is at most tol.'
             print('("final x is  ",*(F9.3,1X))'), params
     case(2)
             print*, 'the relative error between x and the solution is at most tol'
             print('("final x is  ",*(F9.3,1X))'), params
     case(3)
             print*, 'the relative error in the sum of squares and between x and sol is at most tol.'
             print('("final x is  ",*(F9.3,1X))'), params
     case(4)
             print*, 'f is orthogonal to jacobian columns.'
             print('("final x is  ",*(F9.3,1X))'), params
     case(5)
             print*, 'number of calls > 200(n+1) = ', 200*(n+1)
     case(6)
             print*, 'tol is too small, no reduction of f is possible'
     case(7)
             print*, 'tol is too small, no improvement of x possible'
     case default
             print*, 'Something strange'
  end select
  CALL INIT()                                                             !Initializer (user-provided)
  CALL READ_DATA()                                                        !call the user-provided data loader so that experimental data be available
  select case(norming_method)                                             !now norm the obtained discrepancy
     case('L2')
             base = sqrt(trapz(ex_data(:,1), ex_data(:,2)**2))            !L2 norm of the experimental curve
     case('L1')
             base = trapz(ex_data(:,1), abs(ex_data(:,2)))                !L1 squared norm of the experimental curve
     case('C0')
             base = maxval(abs(ex_data(:,2)))                             !C norm of the experimental curve: maximal value
     case('01')
             base = 1.0d0                                                 !no norming
     case('my')
             base = norm(ex_data(:,1), ex_data(:,2)**2)                   !user defined norming
             if(base .le. 0.0) then
                     print*, 'ERROR: User-supplied function base returned negative or zero value:', base   
                     STOP
             end if
  end select
  f = f / base                                                    !the dimensionless value of interest
  sq = sqrt(sum(f**2))
 if(threshold < 0.) then
         print 60, sq, params, real(int_params)                          !if no threshold is given, return the value itself. !!! parameters???
         open(31, file='model.res',status='replace')
         do i = 1,NOP                                                    !evaluate the quantity given the solution and the time array
           write(31,*) Q(i,1), Q(i,2)
         end do
         close(31)
 end if
 if(threshold > 0.) then                                                 !if the threshold is given, 
         open(31, file='model.res',status='replace')
         do i = 1,NOP                                            !evaluate the quantity given the solution and the time array
                 write(31,*) Q(i,1), Q(i,2)
         end do
         close(31)
         if(sq .le. threshold) then                                      !and the value is better,
!                 print 66, 1, sq, params, real(int_params)               !return the success flag; !!! params???
                 print'("1 value is",F9.3," base is ",F9.3," Fb ",F9.3," final x is  ",*(F9.3,1X))', sq,base,sq*base,params, real(int_params)
                 open(32, file='res',status='old',position='append')
                 write(32,'("value is",F9.3,"final x is  ",*(F9.3,1X))'), sq,params, real(int_params)
                 close(32)
         else                                                            !otherwise
                 print'("0 value is",F9.3," base is ",F9.3," Fb ",F9.3," final x is  ",*(F9.3,1X))', sq,base,sq*base,params, real(int_params)
         end if
 end if
 CALL DONE()                                                             !allow the user to clear up
!	 print*,allocated(tv),allocated(Q),allocated(Y)
 if(allocated(tv)) deallocate(tv)
!	 print*,2
 if(allocated(Q)) deallocate(Q)
!	 print*,3
 if(allocated(Y)) deallocate(Y)
!	 print*,4

! Note that the output to STDERR always begins with:
! 1) A number (1 or 0): this means pass/fail result with respect to the user-privided threshold; then come the parameters after the word 'Parameters: ', all (including integer ones) as real numbers
! 2) The words 'Value is ': this means that the result follows an a real number; then come the parameters after the word 'Parameters: ', all (including integer ones) as real numbers
! 3) The word 'ERROR: ': this means that an error has occured. Its description follows.
60   FORMAT ("Value is ",F12.4,1X,"Parameters: ", *(F6.4,1X))
66   FORMAT (I1.1,1X,"Value is ",F6.4,1X,"Parameters: ", *(F6.4,1X))

 CONTAINS

 SUBROUTINE fcn(m,n,x,fvec,iflag)
  use user
  integer, intent(in):: m,n
  integer, intent(inout):: iflag
  integer:: step
  double precision, dimension(n), intent(in):: x
  double precision, dimension(m), intent(out):: fvec
  params = x
  where(params>1.) params=1.
  where(params<0.) params=0.
             print('(" x is  ",*(F9.3,1X))'), params, real(int_params)
  CALL INIT()                                    !Initializer (user-provided): we rely on not changing integer values provided that int params did not change
  tv(0) = 0.                                     !time starts with 0
  Y(:) = 0.                                      !initial values are zero; note: the solver must take real IC into account!
  Q(0,1) = IndependentVar(tv(0))                 !get the initial value of the independent variable (be it time or may be charge)
  Q(0,2) = Quantity(tv(0),Y)                     !get the initial value of the quantity to be compared (real IC to be added inside)
  do step = 1,NOP                                !loop over the time grid
    tv(step) = tv(step-1) + dt                   !the time step
    call NextStep2layer(tv(step),tv(step-1),Y(:))!User provided 2-layer proceeder: 
                                                 !takes current time, previous time, and previos state, returns current state
    if(EarlyQuit(tv(step),Y(:))) then            !check for early stop
           Q(step:NOP,1) = Q(step-1,1);          !if stop, then quantity is unchanged further on
           Q(step:NOP,2) = Q(step-1,2);
           exit                                  !stop the loop
    end if
    Q(step,1) = IndependentVar(tv(step))         !evaluate the data pair
    Q(step,2) = Quantity(tv(step),Y)
    call diagnost                                !call the diagnostic routine to let the user lwarn what is happening
  end do
  CALL READ_DATA()                                                        !call the user-provided data loader so that experimental data be available
  do i = 1,m
        model_data(i) = linspline(Q(:,1),Q(:,2), ex_data(i,1))          !interpolation of the obtained model data to the experimental time grid
  end do
  do i = 1,m
        fvec(i) = (model_data(i)-ex_data(i,2))*sqrt(ex_data(i,1)-ex_data(max(1,i-1),1)) !interpolation of the obtained model data to the experimental time grid
  end do
 CALL DONE()                                                             !allow the user to clear up
 END SUBROUTINE

 END

