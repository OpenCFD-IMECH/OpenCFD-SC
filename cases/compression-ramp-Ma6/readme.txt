DNS of Mach 6, 34° Compression-Ramp flow
Ref:  童福林，李欣，于长平，李新亮，高超声速激波湍流边界层干扰直接数值模拟研究，力学学报，50（2），197-208，2018
-----------------------------------
1） Mesh generation:   run mesh-conner-4.f90 
2)  Initial data  
     1）运行： init3d-1.2.f90 (initial as zero flow) 
     2）或者以已有的结果为初值（可以使用Interploation-corner3d-1.0.f90 对已有的结果进行插值， 插值到现有网格）
    ln -s opencfd0.dat opencfd.dat
3) run opencfd-2 
   mpirun -n 192 ./opencfd-beta-4.out

