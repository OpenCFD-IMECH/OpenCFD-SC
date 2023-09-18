    implicit none 
    real*8,dimension(:,:,:),allocatable:: x,y,z,d,u,v,w,T
	real*8,dimension(:,:),allocatable:: d2d,u2d,v2d,T2d
	real*8,dimension(:),allocatable:: x1,y1
    integer:: nx,ny,nz,i,j,k,it,i0,k0
    real*8::  Tw,T_inf,PI,Tsb,Amu_W,Re,Ama,gamma,time
    real*8:: h1,h2,us1,us2,uy,taow,cf,Ue,De,sseta,ut,uvd,yp,us 
!------------------------------------------------------
    nx=2193
	ny=72
	nz=128
	Re=635000.d0
	Ama=2.25d0
	gamma=1.4d0
    T_inf=169.44d0
	Tw=1.9d0

    allocate(x(nx,ny,nz),y(nx,ny,nz),z(nx,ny,nz), &
	 d(nx,ny,nz),u(nx,ny,nz),v(nx,ny,nz),w(nx,ny,nz), &
	 T(nx,ny,nz),d2d(nx,ny),u2d(nx,ny),v2d(nx,ny),T2d(nx,ny))
    allocate(x1(nx),y1(ny))

    Tsb=110.4/T_inf
    Amu_w=1.0/Re*(1.0+Tsb)*sqrt(Tw**3)/(Tsb+Tw)
    call  read_mesh3d(nx,ny,nz,x,y,z)
	 
     do i=1,nx 
	 x1(i)=x(i,1,1)
	 enddo 
	 do j=1,ny
	 y1(j)=y(1,j,1)
	 enddo 

       print*, "x=",x1
       print*, "y=",y1

     open(70,file='opencfd.dat',form='unformatted')
     read(70) it,time
     print*, it,time
     call read3d(70,d,nx,ny,nz)
     call read3d(70,u,nx,ny,nz)
     call read3d(70,v,nx,ny,nz)
     call read3d(70,w,nx,ny,nz)
     call read3d(70,T,nx,ny,nz)
     close(70)
	  print*, 'read data ok'
!-------- y/2 section --------------
 
    call aver2d(nx,ny,nz,d,u,v,w,T,d2d,u2d,v2d,T2d)



!-------Cf ---------- 
        
     open(66,file='cf.dat')
	  h1=y1(2) 
	  h2=y1(3)
	  do i=1,nx
	    us1=u2d(i,2) ;  us2=u2d(i,3)  
	    uy=(h2*h2*us1-h1*h1*us2)/(h2*h2*h1-h1*h1*h2)  ! 2nd order
	    taow=Amu_w*uy
        ! taow=Amu_w*u2d(i,2)/hs   ! 1st order
	    Cf=taow*2.d0
          write(66,*) x(i,1,1),cf
	  enddo
	  close(66)	 
         
     Print*, "Computation of momentum thickness ..."
     Ue=1.d0
  	 De=1.d0
      open(99,file="momentum-thickness.dat")
	  do i=1,nx
       sseta=0.d0
	   do j=1,ny-1
	     us=u2d(i,j)
!        sseta=sseta+(us/(De*Ue)*(1.d0-us/Ue))*(yy(j+1)-yy(j))          ! ???
         sseta=sseta+(d2d(i,j)*u2d(i,j)/(De*Ue)*(1.d0-u2d(i,j)/Ue))*(y1(j+1)-y1(j))
	   enddo
	  write(99,*) i,x1(i),sseta
	  enddo
	  close(99)


!c--- section normal to the wall 
       print*, "please input i0"
	   read(*,*) i0    
       print*, "x=", x1(i0)          
	   us1=u2d(i0,2) ;  us2=u2d(i0,3) 
	   h1=y1(2) ; h2=y1(3)
	   uy=(h2*h2*us1-h1*h1*us2)/(h2*h2*h1-h1*h1*h2)  ! 2nd order
	   taow=Amu_w*uy
	   ut=sqrt(taow/d2d(i0,1)) 
	   uvd=0.d0
       open(55,file='us.dat')
       do j=2,ny-1
        yp=y1(j)*d2d(i0,1)*ut/Amu_w
        uvd=uvd+(sqrt(d2d(i0,j)/d2d(i0,1)) +  sqrt(d2d(i0,j-1)/d2d(i0,1)))/2.d0 *(u2d(i0,j)-u2d(i0,j-1))
	   write(55,"(4f16.8)") yp, u2d(i0,j)/ut,uvd/ut,  1.d0/(0.41)*log(yp)+5.1
	   enddo
	   close(55) 
	   print*, "Ret=", d2d(i0,1)*ut/Amu_w 

       k0=nz/2
       call write_xy(k0,nx,ny,nz,x,y,z,d,u,v,w,T,gamma,Ama)  
       deallocate( x,y,z,d,u,v,w,T, d2d,u2d,v2d,T2d, x1,y1)
	  end
!c==================================
      subroutine read3d(no,u,nx,ny,nz)
      implicit none
      integer:: no,nx,ny,nz,k
      real*8:: u(nx,ny,nz)
	  do k=1,nz
	   read(no) u(:,:,k)
	  enddo
	  end
!c========================================


     subroutine write_xy(k0,nx,ny,nz,x,y,z,d,u,v,w,T,gamma,Ama)  
     implicit none
     integer::nx,ny,nz,i,j,k0
     real*8,dimension(nx,ny,nz):: d,u,v,w,T,x,y,z
     real*8::  gamma,Ama,p00,p
	 p00=1.d0/(gamma*Ama*Ama)
	 
     open(44,file='flow2d_xy.dat')
	 write(44,*) 'variables= x,y,d,u,v,w,p,T'
	 write(44,*) 'zone i=',nx, 'j=',ny
	  do j=1,ny
      do i=1,nx
	  p=d(i,j,k0)*T(i,j,k0)*p00
	  write(44,'(8f16.8)') x(i,j,k0),y(i,j,k0), &
         d(i,j,k0),u(i,j,k0),v(i,j,k0),w(i,j,k0),p, T(i,j,k0)
	  enddo
	  enddo
     close(44)
	 end

   subroutine aver2d(nx,ny,nz,d,u,v,w,T,d2d,u2d,v2d,T2d)
   implicit none  
   integer:: nx,ny,nz,i,j,k   
   real*8,dimension(nx,ny,nz):: d,u,v,w,T
   real*8:: d2d(nx,ny),T2d(nx,ny),u2d(nx,ny),v2d(nx,ny)
   do j=1,ny
   do i=1,nx
	 d2d(i,j)=0.d0
	 u2d(i,j)=0.d0
	 v2d(i,j)=0.d0
	 T2d(i,j)=0.d0
    do k=1,nz
	  d2d(i,j)=d2d(i,j)+d(i,j,k)
	  u2d(i,j)=u2d(i,j)+u(i,j,k)
	  v2d(i,j)=v2d(i,j)+v(i,j,k)
	  T2d(i,j)=T2d(i,j)+T(i,j,k)
    enddo
	 d2d(i,j)=d2d(i,j)/nz
	 u2d(i,j)=u2d(i,j)/nz
	 v2d(i,j)=v2d(i,j)/nz
	 T2d(i,j)=T2d(i,j)/nz
     enddo
	 enddo
    end
	
!---------------------------------------------------------	
	subroutine read_mesh3d(nx,ny,nz,x,y,z)
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
