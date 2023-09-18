DNS of boundary layer transition of Mach 2.25 flate-plate 
----------------------------------------
Ref: 
  Li XL et al. Chinese Physics Letters 26(9), 094701, 2009
  Gao H. et al. Chinese Physics Letters 22(7), 1709, 2005 
  Rai MM, et al. AIAA 95-0583
  Pirozzoli S. et al. Phys. Fluids 16, 530, 2004
-----------------------------------------
1) Generate mesh and initial data
    Compile and run Grid-and-init.f90
 2) ln -s opencfd0.dat opencfd.dat
 3) mpirun -n 64 ./opencfd-beta-4.out (or new versions)
 4) rm opencfd.dat
     ln -s OCFDxxxxxxxx.dat opencfd.dat  
 5) compile and run read-flatplate.f90
