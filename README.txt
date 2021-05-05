To compile, if one has gcc and intel fortran binded with an MPI, for example, do the following:

gcc -c WELL19937a_new.c
mpif90 -c parModules.f90
mpif90 -o edipic1d.out par*f90 WELL19937a_new.o

To run, copy the executable file and the input data files (ssc_*dat, inluding the cross sections for the selected ion species) into some directory, edit the input data files as necessary (note that the input is formatted, so one has to follow patterns provided in corresponding description lines), then do for example:

mpirun -np 8 ./edipic1d.out > output.txt &



The code has been developed with the Government sponsorship (contract DE-AC02-09CH11466). The developer grants to the Government, and others acting on its behalf, a nonexclusive, paid-up, irrevocable, world-wide license in such copyrighted data to reproduce, prepare derivative works, distribute copies to the public, and perform publicly and display publicly, by or on behalf of the Government.
