
 PROGRAM TASK                                                            ! A task for BOINC grid umbrella computing
 
 USE user                                                                ! The user-provided module with parameters and right-hand side of the system of ODEs
 use interpolation                                                       ! interpolators
 use integration                                                         ! numerical integrators
 use random                                                              ! random number generator initialization

 IMPLICIT NONE

 INTEGER:: npar, par, n                                                  ! loop counters
 DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: tv                      ! time array
 DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: Q                       ! the compared quantity array (x-y columns)
 DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE, TARGET:: Y               ! solution array
 DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: model_data              ! model quantity interpolated
 DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: val, dval,grad          ! real values of the parameters and their perturbations ***
 INTEGER:: arsize                                                        ! size of the ex_data array
 DOUBLE PRECISION:: base, sq, h
 CHARACTER(len=64):: nmlfile = ''
 LOGICAL, PARAMETER:: compare = .false.
 INTEGER:: gr!***

 CALL INIT_RANDOM_SEED()                                                 !random numbers generator initialization
 
  do par = 1,N_params                                                    !random choice of the real parameters; each is from [0,1)
        call random_number(params(par))
  end do

  do npar = 1,N_int_params                                               !random choice of the integer parameters; each is 0-9
        call random_number(h)
        int_params(npar) = int(h*10)
  end do
  if(command_argument_count()>0) call get_command_argument(1,nmlfile)

  if(len(trim(nmlfile)) > 0) then 
          CALL INIT(nmlfile)                                             !Initializer (user-provided)
  else
          CALL INIT()                                                    !Initializer (user-provided)
          nmlfile = 'param.nml'
  end if

 allocate(tv(0:NOP), Q(0:NOP,2))                                         !allocate arrays for the time grid and 2-column model data (time-value or smth of the sort)
 allocate(Y(n_unknowns))                                                 !allocate the phase vector: all variables at any time will be here

 allocate(val(N_params),dval(N_params),grad(0:N_params)) !***
 do gr = 0,N_params      !***
         call get_params_values(val) !***
         dval = 0.!***
         if(gr>0) dval(gr) = val(gr) * (-0.01)!***
         call PERTURBE(dval)!***

 tv(0) = 0.                                                              !time starts with 0
 Y(:) = 0.                                                               !initial values are zero; note: the solver must take real IC into account!
 Q(0,1) = IndependentVar(tv(0))                                          !get the initial value of the independent variable (be it time or may be charge)
 Q(0,2) = Quantity(tv(0),Y)                                              !get the initial value of the quantity to be compared (real IC to be added inside)

 do n = 1,NOP                                                            !loop over the time grid
   tv(n) = tv(n-1) + dt                                                  !the time step
   call NextStep2layer(tv(n),tv(n-1),Y(:))                               !User provided 2-layer proceeder: 
                                                                         !takes current time, previous time, and previos state, returns current state
   if(EarlyQuit(tv(n),Y(:))) then                                        !check for early stop
          Q(n:NOP,1) = Q(n-1,1);                                         !if stop, then quantity is unchanged further on
          Q(n:NOP,2) = Q(n-1,2);
          exit                                                           !stop the loop
   end if
   Q(n,1) = IndependentVar(tv(n))                                        !evaluate the data pair
   Q(n,2) = Quantity(tv(n),Y)
   call diagnost                                                         !call the diagnostic routine to let the user lwarn what is happening
 end do

 CALL READ_DATA()                                                        !call the user-provided data loader so that experimental data be available

 arsize = size(ex_data, dim=1)                                           !get the size of the experimental data array
 if(.not.allocated(model_data)) allocate(model_data(arsize,2))                                          !allocate the array for the model data
 model_data(:,1) = ex_data(:,1)                                          !the same x-values
 do n = 1,arsize
       model_data(n,2) = linspline(Q(:,1),Q(:,2), model_data(n,1))       !interpolation of the obtained model data to the experimental time grid
 end do

 select case(space)                                                      !now norm the obtained discrepancy
    case('L2')
            sq = sqrt(trapz(ex_data(:,1), (model_data(:,2) - ex_data(:,2))**2)) !the L2 norm of the difference between model and experimental data
    case('L1')
            sq = trapz(ex_data(:,1), abs(model_data(:,2) - ex_data(:,2)))       !the L1 norm of the difference between model and experimental data
    case('C0')
            sq = maxval(abs(model_data(:,2) - ex_data(:,2)))                    !the C norm of the difference between model and experimental data
    case('my')
            sq = LSQ(ex_data(:,1), model_data(:,2) - ex_data(:,2))       !user defined norm
 end select
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
 sq = sq / base                                                    !the dimensionless value of interest
 if(gr.eq.0) then!***
         grad(0) = sq!***
 else!***
         grad(gr) = (sq - grad(0)) / dval(gr)!***
 endif!***
 enddo !***
 sq = grad(0) !***
 print*, 'grad=', grad(1:N_params)!***
 print*, 'normed grad=', grad(1:N_params)*val(:)!***

 if(threshold < 0.) then
         print 60, sq, params, real(int_params)                          !if no threshold is given, return the value itself. !!! parameters???
         open(31, file='model.res',status='replace')
         do n = 0,NOP                                                    !evaluate the quantity given the solution and the time array
           write(31,*) Q(n,1), Q(n,2)
         end do
         close(31)
         call savenml(nmlfile)
 end if
 if(threshold > 0.) then                                                 !if the threshold is given, 
         if(sq .le. threshold) then                                      !and the value is better,
                 print 66, 1, sq, params, real(int_params)               !return the success flag; !!! params???
                 open(31, file='model.res',status='replace')
                 do n = 0,NOP                                            !evaluate the quantity given the solution and the time array
                     write(31,*) Q(n,1), Q(n,2)
                 end do
                 close(31)
                 if(compare) then
                         open(33, file='cmp.res',status='replace')
                         do n = 1,arsize
                                   write(33,*) model_data(n,:), ex_data(n,2), ex_data(n,1)!-model_data(n,1)
                             enddo
                         close(33)
                 endif
                 call savenml(nmlfile)
         else                                                            !otherwise
                 print 66, 0, sq, params, real(int_params)               !the fail flag.
         end if
 end if
 
 CALL DONE()                                                             !allow the user to clear up
 deallocate(tv, Q, Y)

! Note that the output to STDERR always begins with:
! 1) A number (1 or 0): this means pass/fail result with respect to the user-privided threshold; then come the parameters after the word 'Parameters: ', all (including integer ones) as real numbers
! 2) The words 'Value is ': this means that the result follows an a real number; then come the parameters after the word 'Parameters: ', all (including integer ones) as real numbers
! 3) The word 'ERROR: ': this means that an error has occured. Its description follows.
60   FORMAT ("Value is ",F12.4,1X,"Parameters: ", *(F6.4,1X))
66   FORMAT (I1.1,1X,"Value is ",F8.6,1X,"Parameters: ", *(F6.4,1X))

 END

