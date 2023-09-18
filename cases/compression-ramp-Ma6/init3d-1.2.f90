! -----  Init 3D ----
! Generate initial data for opencfd-sc, Ver 1.1     
	implicit none
	integer:: nx,ny,nz,Iflag
    real*8,allocatable,dimension(:,:,:):: u
	real*8:: Pi, AoA, ux,uy
    print*, "Please input nx,ny,nz ?"
	read(*,*) nx,ny,nz
	allocate(u(nx,ny,nz))
!--------------------
        Pi=3.14159265358979d0
        print*, "Please input Iflag (0/1): 0 init wiht 0 flow, 1 with free-stream flow"
        read(*,*) Iflag
        if(Iflag ==0) then
          print*, "Init with 0 flow"
          ux=0.d0
          uy=0.d0
        else
         print*, "Init with free-stream flow"
	 print*, "please input AoA (Angle of attack, in degree) "
	 read(*,*) AoA
	 AoA=AoA*Pi/180.d0
         ux=cos(AoA)
	 uy=sin(AoA)
       endif


	open(99,file="opencfd0.dat",form="unformatted")
       write(99) 0,0.d0
      u(:,:,:)=1.d0   ! density
	call write3d(99,nx,ny,nz,u)
      u(:,:,:)=ux     ! u
	call write3d(99,nx,ny,nz,u)
      u(:,:,:)=uy     ! v
	call write3d(99,nx,ny,nz,u)
      u(:,:,:)=0.d0   ! w
	call write3d(99,nx,ny,nz,u)
      u(:,:,:)=1.d0   ! T
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

