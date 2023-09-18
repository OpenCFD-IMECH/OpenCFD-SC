 Mach 5 flow over a half-cylinder
 1) run Grid_and_init_HC3d.f90 for mesh and initial data
 2) ln -s opencfd0.dat opencfd.dat
 3) mpirun -n 24 ./opencfd-x.x.out 
 4) run postproc3d-1.0.f90 (in ../../) for output tecplot-type data
