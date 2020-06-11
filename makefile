
all: Al battery U example plate 

go: solver2
	./solver2

solver2: user.o integration.o interpolation.o ode.o diffusion.o equations.o solver2.f90
	gfortran -o solver2 $^

user.o: user.f90
	gfortran -c user.f90

ode.o: ode.f90
	gfortran -c ode.f90

integration.o: integration.f90
	gfortran -c integration.f90

interpolation.o: interpolation.f90
	gfortran -c interpolation.f90

equations.o: equations.f90
	gfortran -c equations.f90

diffusion.o: equations.o diffusion.f90
	gfortran -c diffusion.f90

Al: user_Al.f90 integration.o interpolation.o ode.o solver2.f90
	gfortran -o $@ -ffree-line-length-none $^

Mg: user_Mg.f90 integration.o interpolation.o ode.o solver2.f90
	gfortran -o $@ -ffree-line-length-none $^

battery:  user_battery.f90 integration.o interpolation.o diffusion.o equations.o solver2.f90
	gfortran -fbounds-check -o $@ -ffree-line-length-none $^

U:        user.U.f90 integration.o interpolation.o diffusion.o equations.o solver2.f90
	gfortran -fbounds-check -o $@ -ffree-line-length-none $^
	./run.pl

Us:        user_U_Sshape.f90 integration.o interpolation.o parabolic.o equations.o solver2.f90
	gfortran -fbounds-check -o U -ffree-line-length-none $^
	./run.pl
	#./U
	#gplot -outfile cmp.png -type png model.res curve2.dat
	cat diagn.out

example:  user_example.f90 integration.o interpolation.o diffusion.o equations.o solver2.f90
	gfortran -fbounds-check -o solver_example -ffree-line-length-none $^

ex0lay:  user.0.f90 solver0.f90
	gfortran -fbounds-check -o ex0 -ffree-line-length-none $^

plate:  user_nonlinbc_ex.f90 integration.o interpolation.o diffusion.o equations.o solver2.f90
	gfortran -fbounds-check -o solver_plate -ffree-line-length-none $^

%.o: %.f90
	gfortran -fbounds-check -o $*.o -ffree-line-length-none -c $<

clear:
	rm *.o *.mod

clean: clear

arc:
	tar -zcf backup`date +%F`.tar.gz *.f90






