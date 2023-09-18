! -----  Init 3D ----
! Generate initial data for opencfd-sc, Ver 1.3     
	implicit none
	integer:: nx,ny,nz,j,Iflag
    real*8,allocatable,dimension(:,:):: d,u,v,w,T
	real*8:: Pi, AoA, tmp,d1,u1,v1,T1
    print*, "Please input nx,ny,nz ?"
	read(*,*) nx,ny,nz
	allocate(d(nx,ny),u(nx,ny),v(nx,ny),w(nx,ny),T(nx,ny))

!--------------------
        Pi=3.14159265358979d0
        print*, "Please input Iflag (0/1): 0 init wiht free-stream flow, 1 with flow1d-inlet.dat"
        read(*,*) Iflag
        if(Iflag ==0) then
         print*, "Init with free-stream flow"
	     print*, "please input AoA (Angle of attack, in degree) "
	     read(*,*) AoA
	     AoA=AoA*Pi/180.d0
		 d=1.d0 
         u=cos(AoA)
	     v=sin(AoA)
         w=0.d0 
		 T=1.d0
        else 
		 open(99,file="flow1d-inlet.dat")
		 read(99,*)
		 do j=1,ny 
		 read(99,*) tmp, d1,u1,v1,T1 
         d(:,j)=d1 
		 u(:,j)=u1
		 v(:,j)=v1
		 w(:,j)=0.d0
		 T(:,j)=T1
		 enddo 
         close(99)
		endif


	open(99,file="opencfd0.dat",form="unformatted")
    write(99) 0,0.d0
	call write3d1(99,nx,ny,nz,d)
	call write3d1(99,nx,ny,nz,u)
	call write3d1(99,nx,ny,nz,v)
	call write3d1(99,nx,ny,nz,w)
	call write3d1(99,nx,ny,nz,T)
    close(99)
	deallocate(d,u,v,w,T)
 	end
	
	subroutine write3d1(no,nx,ny,nz,u)
    implicit none
 	integer:: no,nx,ny,nz,k
	real*8:: u(nx,ny)
	do k=1,nz
	write(no) u
	enddo
	end

