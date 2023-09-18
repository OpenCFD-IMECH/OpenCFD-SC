  ! Mesh and Initial condition for  flow over a Half-Cylinder
 	implicit none 
	integer,parameter:: nx=181, ny=101 , nz=16
	real*8,parameter::  PI=3.1415926535897d0 
    integer:: i,j,k
	real*8:: r,seta,R1,R2,hy1,Lz 
	real*8:: x2d(nx,ny),y2d(nx,ny),z1d(nz)
	real*8 d(nx,ny,nz),u(nx,ny,nz),v(nx,ny,nz),w(nx,ny,nz),T(nx,ny,nz),y1d(ny)
!---------------------------------------	
    R1=1.d0
    R2=3.d0 
	Lz=1.d0 
    hy1=0.001d0    
	call getsy(ny,y1d,R2-R1,hy1)
 !-------------------------------------

   do j=1,ny
   do i=1,nx 
    r=R1+y1d(j) 
	seta= 3.d0/2.d0*PI -PI*(i-1.d0)/(nx-1.d0)
    x2d(i,j)=r*cos(seta)
    y2d(i,j)=r*sin(seta)
   enddo
   enddo
   do k=1,nz 
 	z1d(k)=(k-1.d0)/nz*Lz  
   enddo 
   
   d=1.d0 
   u=1.d0
   v=0.d0 
   w=0.d0
   T=1.d0 
   
   open(99,file="OCFD-grid.dat",form="unformatted")
   write(99) x2d
   write(99) y2d 
   write(99) z1d 
   close(99) 

   open(100,file="opencfd0.dat",form="unformatted")
   write(100) 0, 0.d0 
   call write3d(100,nx,ny,nz,d)
   call write3d(100,nx,ny,nz,u)
   call write3d(100,nx,ny,nz,v)
   call write3d(100,nx,ny,nz,w)
   call write3d(100,nx,ny,nz,T)
   close(100) 
   end 


!c=================================================================
         subroutine getsy(ny,yy,SL,hy1)
 	     implicit none 
	     real*8 yy(ny),SL,hy1,delta,fb,fbx,bnew,a,s,dy
         real*8,save:: b=3.5d0
		 integer:: ny,j
		 
	     dy=1.d0/(ny-1.d0)
!---------------------------------------
         delta=hy1/SL
 ! using Newton method get coefficient
 100     continue
        fb=(exp(b/(ny-1.d0))-1.d0)/(exp(b)-1.d0)-delta
        fbx=(exp(b/(ny-1.d0))/(ny-1.d0)*(exp(b)-1.d0) -      &     
         (exp(b/(ny-1.d0))-1.d0)*exp(b)  )/((exp(b)-1.d0))**2
         bnew=b-fb/fbx
         if(abs(b-bnew) .gt. 1.d-6) then
	      b=bnew
		  goto 100
	     endif
	    
         b=bnew 
         a=1.d0/(exp(b)-1.d0)*SL 
         do j=1,ny
         s=(j-1.d0)*dy
         yy(j)=a*(exp(s*b)-1.d0)
         enddo
       end
!-----------------------------------------------
     subroutine write3d(file_no,nx,ny,nz,U)
	 implicit none 
	 integer:: file_no, nx,ny,nz, k
	 real*8:: U(nx,ny,nz) 
     do k=1,nz 
      write(file_no) U(:,:,k) 
     enddo 
     end 
	 



