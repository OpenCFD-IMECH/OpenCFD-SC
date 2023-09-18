! Transform OCFD-Jacobian2d.dat to OCFD-grid.dat -------------------------------------
! Copyright by Li Xinliang, Email : lixl@imech.ac.cn
! LHD, Institute of Mechanics, CAS
     implicit none 
	 real*8,allocatable:: x1(:),y1(:),z1(:),x2(:,:),y2(:,:)
	 real*8:: Lz
	 integer:: Iflag,nx,ny,nz,i,j,k 
	 print*, "Transform OCFD2d-Jacobi.dat/ocfd-grid.dat to OCFD-grid.dat"
     print*, " 1 transform 1d (ocfd-grid.dat),  2 transform 2d (OCFD2d-Jacobi.dat)"
     read(*,*) Iflag	 
	 print*, "please input nx,ny,nz"
	 read(*,*) nx,ny,nz 
	 print*, "please input Lz (length of spanwise domain)"
     read(*,*) Lz 
	 allocate(x1(nx),y1(ny),z1(nz),x2(nx,ny),y2(nx,ny))
	 do k=1,nz 
	 z1(k)=Lz*(k-1.d0)/nz
	 enddo 
	 open(99,file="OCFD-grid.dat",form="unformatted") 
	 if(Iflag .eq. 1) then 
	 open(100,file="ocfd-grid.dat",form="unformatted")            ! opencfd-ver 1.x data
	 read(100) x1 
	 read(100) y1
	 close(100)
	 write(99) x1 
	 write(99) y1 
	 write(99) z1 
	 else if (Iflag .eq. 2) then 
	 open(100,file="OCFD2d-Jacobi.dat",form="unformatted")
	 read(100) x2 
	 read(100) y2 
	 close(100)
	 write(99) x2 
	 write(99) y2 
	 write(99) z1 
	 endif 
	 close(99)
	 deallocate(x1,y1,z1,x2,y2)
	 end 
	 
	 
	 
!--------------------------------------------
	subroutine read_mesh(nx,ny,nz,x,y,z)
     implicit none 
     integer:: nx,ny,nz,mesh_type,i,j,k
	 integer,parameter:: GRID1D=10, GRID2D_PLANE=20, GRID2D_AXIAL_SYMM=21, GRID3D=30, GRID_AND_JACOBIAN3D=31 
	 real*8:: x(nx,ny,nz),y(nx,ny,nz),z(nx,ny,nz)
	 real*8,allocatable:: x1d(:),y1d(:),z1d(:),x2d(:,:),y2d(:,:)
	 allocate(x1d(nx),y1d(ny),z1d(nz),x2d(nx,ny),y2d(nx,ny)) 
	 print*, "read  mesh ......"
	 print*, "pleas input mesh_type,   10: 1D mesh;  20: 2D-plane mesh;  21: 2D-AxialSymmetry mesh;  30: 3D mesh"
	 read(*,*)  mesh_type 
	 
     open(56,file='OCFD-grid.dat',form='unformatted')
     if(mesh_type ==GRID1D) then 
	   read(56) x1d 
       read(56) y1d
       read(56) z1d 
	  do k=1,nz
	  do j=1,ny 
	  do i=1,nx 
	   x(i,j,k)=x1d(i) 
	   y(i,j,k)=y1d(j) 
	   z(i,j,k)=z1d(k) 
	  enddo 
	  enddo 
	  enddo 
	 else if( mesh_type ==GRID2D_PLANE) then 
	  read(56) x2d 
	  read(56) y2d 
	  read(56) z1d 
	  do k=1,nz
	  do j=1,ny 
	  do i=1,nx 
	   x(i,j,k)=x2d(i,j) 
	   y(i,j,k)=y2d(i,j) 
	   z(i,j,k)=z1d(k) 
	  enddo 
	  enddo 
	  enddo	  
	 else if( mesh_type ==GRID2D_AXIAL_SYMM) then  
	  read(56) x2d 
	  read(56) y2d 
	  read(56) z1d 
	  do k=1,nz
	  do j=1,ny 
	  do i=1,nx 
	   x(i,j,k)=x2d(i,j) 
	   y(i,j,k)=y2d(i,j)*cos(z1d(k))
	   z(i,j,k)=y2d(i,j)*sin(z1d(k)) 
	  enddo 
	  enddo 
	  enddo	    
	 else  
      call read3d(56,nx,ny,nz,x)
      call read3d(56,nx,ny,nz,y)
      call read3d(56,nx,ny,nz,z)
	 endif 
	 close(56)
	 deallocate(x1d,y1d,z1d,x2d,y2d)
     end 
