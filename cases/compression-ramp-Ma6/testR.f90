! OpenCFD  Post-analysis code  -----------------------------------------
! Copyright by Li Xinliang, Email : lixl@imech.ac.cn
! LHD, Institute of Mechanics, CAS
     implicit none 
	 real*8,allocatable,dimension(:,:,:):: x,y,z,R 
	 integer:: nx,ny,nz,i,j,k
	 print*, "Post-processing code for OpenCFD2,  output tecplot-type data"
	 print*, "Please input nx, ny, nz"
	 read(*,*) nx,ny,nz 
	 
      allocate(x(nx,ny,nz),y(nx,ny,nz),z(nx,ny,nz))
      allocate(R(nx,ny,nz))
      call read_mesh(nx,ny,nz,x,y,z)

      open(77,file="Rhybrid.dat",form="unformatted"	)
	  call read3d(77,nx,ny,nz,R) 
      close(77)
	  k=nz/2
	  open(44,file='R2d.dat')
      write(44,*) 'variables= x,y,z,R'
      write(44,*) 'zone i=',nx, 'j=',ny
      do j=1,ny
      do i=1,nx
      write(44,'(10f16.8)') x(i,j,k),y(i,j,k),z(i,j,k), R(i,j,k)
      enddo
      enddo
     close(44)

    deallocate(x,y,z,R)
   end
!==================================
       subroutine read3d(no,nx,ny,nz,u)
       implicit none
	   integer:: no,nx,ny,nz,k
       real*8:: u(nx,ny,nz)
	   print*, "read 3d data ..."
       do k=1,nz
       read (no)  u(:,:,k)
       enddo
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
