! -----  Init 3D ----
    implicit none
	integer,parameter:: nx=2193,ny=72,nz=128
    real*8,allocatable,dimension(:,:,:):: u
	real*8,allocatable,dimension(:):: x1d, y1d, z1d
	real*8:: Lz=0.175
	integer:: i,j,k
	
	allocate(u(nx,ny,nz))
    allocate(x1d(nx),y1d(ny),z1d(nz))
	open(99,file="ocfd-grid.dat",form="unformatted")
	read(99) x1d 
	read(99) y1d 
	close(99)

	do k=1,nz 
	z1d(k)=(k-1.d0)*Lz/nz 
	enddo 
	

	open(99,file="OCFD-grid.dat",form="unformatted")
    write(99) x1d 
	write(99) y1d 
	write(99) z1d 
	close(99)
	
    deallocate(x1d,y1d,z1d)

!---------------------------------------------------------
	open(99,file="opencfd0.dat",form="unformatted")
    write(99) 0,0.d0
    u(:,:,:)=1.d0   ! density
	call write3d(99,nx,ny,nz,u)
    u(:,:,:)=1.d0
	call write3d(99,nx,ny,nz,u)
    u(:,:,:)=0.d0
	call write3d(99,nx,ny,nz,u)
    u(:,:,:)=0.d0
	call write3d(99,nx,ny,nz,u)
    u(:,:,:)=1.d0
	call write3d(99,nx,ny,nz,u)
    close(99)
	deallocate(u)
 	end
	
	subroutine write3d(no,nx,ny,nz,u)
    implicit none
 	integer:: no,nx,ny,nz,k
	real*8:: u(nx,ny,nz)
	do k=1,nz
	write(no) u(:,:,k)
	enddo
	end

