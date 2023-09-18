! Finite difference for viscous terms  (Centeral Schemes)

    subroutine OCFD_dx0(f,fx, Num_Scheme)
    use flow_para
	implicit none
    integer Num_Scheme,i,j,k
	real(kind=OCFD_REAL_KIND)::   f(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP),fx(nx,ny,nz)
	real(kind=OCFD_REAL_KIND)::   a1,a2,a3,c1,c2,c3,c4
       a1=1.d0/(60.d0*hx)             ! 6th centeral scheme
       a2=-3.d0/(20.d0*hx)
       a3=3.d0/(4.d0*hx)	
	   
       c1=0.8d0/hx                    ! 8th centreral scheme
	   c2=-0.2d0/hx
	   c3=3.80952380952380952d-2/hx
	   c4=-3.571428571428571428d-3/hx 	

	   
!------Scheme for inner points  ( CD6, CD8)
    if(Num_Scheme == OCFD_Scheme_CD6 ) then 
        do k=1,nz
		do j=1,ny
        do i=1,nx
         fx(i,j,k)=a1*(f(i+3,j,k)-f(i-3,j,k)) +a2*(f(i+2,j,k)-f(i-2,j,k))  +a3*(f(i+1,j,k)-f(i-1,j,k))
        enddo
        enddo	
	    enddo 
		
	else if (Num_Scheme == OCFD_Scheme_CD8) then 
	     do k=1,nz
		 do j=1,ny
         do i=1,nx
          fx(i,j,k)=c1*(f(i+1,j,k)-f(i-1,j,k)) +c2*(f(i+2,j,k)-f(i-2,j,k))  &
                   +c3*(f(i+3,j,k)-f(i-3,j,k)) +c4*(f(i+4,j,k)-f(i-4,j,k))
         enddo
         enddo
	     enddo 
	else 
	  print*, 'This Numerical Scheme is not supported in viscous terms !'
	  print*, 'Only CD6 or CD8 can be used in viscous terms'
	  stop    
	endif

!---------Boundary Scheme ------- (low-order scheme)----
!---------- i- boundary  ---------------------------     
	 if(npx .eq. 0 .and. Para%Iperiodic_X .eq. 0) then 
	  do k=1,nz
	  do j=1,ny
	    fx(1,j,k)=(-3.d0*f(1,j,k)+4.d0*f(2,j,k)-f(3,j,k))/(2.d0*hx)           ! 2nd one-side scheme
	    fx(2,j,k)=(f(3,j,k)-f(1,j,k))/(2.d0*hx)                             ! 2nd centeral scheme 
		fx(3,j,k)=(8.d0*(f(4,j,k)-f(2,j,k)) - (f(5,j,k)-f(1,j,k)))/(12.d0*hx)   ! 4th central scheme
	  enddo 
	  enddo 
       if(Num_Scheme == OCFD_Scheme_CD8) then 
	   do k=1,nz
	   do j=1,ny
         fx(4,j,k)=a1*(f(7,j,k)-f(1,j,k)) +a2*(f(6,j,k)-f(2,j,k))  +a3*(f(5,j,k)-f(3,j,k))	  ! 6th centeral scheme
	   enddo
	   enddo 
	   endif 
	 endif 
!--------- i+ boundary ------------------------
	 if(npx .eq. npx0-1 .and. Para%Iperiodic_X .eq. 0) then 
	   do k=1,nz
	   do j=1,ny
	    fx(nx,j,k)=(f(nx-2,j,k)-4.d0*f(nx-1,j,k)  +3.d0*f(nx,j,k))/(2.d0*hx)  ! 2nd one-side scheme
	    fx(nx-1,j,k)=(f(nx,j,k)-f(nx-2,j,k))/(2.d0*hx)                             ! 2nd centeral scheme 
		fx(nx-2,j,k)=(8.d0*(f(nx-1,j,k)-f(nx-3,j,k)) - (f(nx,j,k)-f(nx-4,j,k)))/(12.d0*hx)   ! 4th central scheme
	   enddo 
	   enddo
       if(Num_Scheme == OCFD_Scheme_CD8) then 
	   do k=1,nz
	   do j=1,ny
         fx(nx-3,j,k)=a1*(f(nx,j,k)-f(nx-6,j,k)) +a2*(f(nx-1,j,k)-f(nx-5,j,k))  +a3*(f(nx-2,j,k)-f(nx-4,j,k))	  ! 6th centeral scheme
	   enddo
	   enddo 
	   endif 
	 endif  


     end

!c----------------------------------------------------------

	subroutine OCFD_dy0(f,fy,Num_Scheme )
     use flow_para
	 implicit none
     integer Num_Scheme,i,j,k	  
	 real(kind=OCFD_REAL_KIND):: f(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP), fy(nx,ny,nz)
	 real(kind=OCFD_REAL_KIND):: a1,a2,a3,c1,c2,c3,c4
      a1=1.d0/(60.d0*hy)
      a2=-3.d0/(20.d0*hy)
      a3=3.d0/(4.d0*hy)	
      c1=0.8d0/hy
	  c2=-0.2d0/hy
	  c3=3.80952380952380952d-2/hy
	  c4=-3.571428571428571428d-3/hy 
	
!------Scheme for inner points  ( CD6, CD8)
    if(Num_Scheme == OCFD_Scheme_CD6 ) then 
     do k=1,nz 
	 do j=1,ny
     do i=1,nx
       fy(i,j,k)=a1*(f(i,j+3,k)-f(i,j-3,k))  +a2*(f(i,j+2,k)-f(i,j-2,k)) +a3*(f(i,j+1,k)-f(i,j-1,k))
     enddo
     enddo 
     enddo 
	 
	else if (Num_Scheme == OCFD_Scheme_CD8) then 
	  do k=1,nz
	  do j=1,ny
      do i=1,nx
       fy(i,j,k)=c1*(f(i,j+1,k)-f(i,j-1,k)) +c2*(f(i,j+2,k)-f(i,j-2,k))  &
                +c3*(f(i,j+3,k)-f(i,j-3,k)) +c4*(f(i,j+4,k)-f(i,j-4,k))
      enddo
      enddo
	  enddo 
	else 
	  
	  print*, 'This Numerical Scheme is not supported in viscous terms !'
	  print*, 'Only CD6 or CD8 can be used in viscous terms'
	  stop    
	endif	 
     
!---------Boundary Scheme ------- (low-order scheme)----
!---------- j- boundary  ---------------------------     
   if(npy .eq. 0 .and. Para%Iperiodic_Y .eq. 0) then 
	 do k=1,nz  
	 do i=1,nx
	    fy(i,1,k)=(-3.d0*f(i,1,k)+4.d0*f(i,2,k)-f(i,3,k))/(2.d0*hy)           ! 2nd one-side scheme
	    fy(i,2,k)=(f(i,3,k)-f(i,1,k))/(2.d0*hy)                             ! 2nd centeral scheme 
		fy(i,3,k)=(8.d0*(f(i,4,k)-f(i,2,k)) - (f(i,5,k)-f(i,1,k)))/(12.d0*hy)   ! 4th central scheme
     enddo 
     enddo 
  if(Num_Scheme == OCFD_Scheme_CD8) then 
     do k=1,nz	 
	 do i=1,nx
        fy(i,4,k)=a1*(f(i,7,k)-f(i,1,k)) +a2*(f(i,6,k)-f(i,2,k))  +a3*(f(i,5,k)-f(i,3,k))	  ! 6th centeral scheme
	 enddo
	 enddo 
	 endif 
   endif 
!--------- j+ boundary ------------------------
	 if(npy .eq. npy0-1 .and. Para%Iperiodic_Y .eq. 0) then 
      do k=1,nz	  
	  do i=1,nx
	    fy(i,ny,k)=(f(i,ny-2,k)-4.d0*f(i,ny-1,k)  +3.d0*f(i,ny,k))/(2.d0*hy)  ! 2nd one-side scheme
	    fy(i,ny-1,k)=(f(i,ny,k)-f(i,ny-2,k))/(2.d0*hy)                             ! 2nd centeral scheme 
		fy(i,ny-2,k)=(8.d0*(f(i,ny-1,k)-f(i,ny-3,k)) - (f(i,ny,k)-f(i,ny-4,k)))/(12.d0*hy)   ! 4th central scheme
	  enddo 
      enddo      
	 if(Num_Scheme == OCFD_Scheme_CD8) then 
     do k=1,nz	 
	 do i=1,nx
       fy(i,ny-3,k)=a1*(f(i,ny,k)-f(i,ny-6,k)) +a2*(f(i,ny-1,k)-f(i,ny-5,k))  +a3*(f(i,ny-2,k)-f(i,ny-4,k))	  ! 6th centeral scheme
	 enddo
	 enddo 
	 endif 
	 endif  	 

	 end

 
!c----------------------------------------------------------

	subroutine OCFD_dz0(f,fz,Num_Scheme )
     use flow_para
	 implicit none
     integer Num_Scheme,i,j,k	  
	 real(kind=OCFD_REAL_KIND):: f(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP), fz(nx,ny,nz)
	 real(kind=OCFD_REAL_KIND):: a1,a2,a3,c1,c2,c3,c4
      a1=1.d0/(60.d0*hz)
      a2=-3.d0/(20.d0*hz)
      a3=3.d0/(4.d0*hz)	
      c1=0.8d0/hz
	  c2=-0.2d0/hz
	  c3=3.80952380952380952d-2/hz
	  c4=-3.571428571428571428d-3/hz 
	
!------Scheme for inner points  ( CD6, CD8)
    if(Num_Scheme == OCFD_Scheme_CD6 ) then 
        do k=1,nz
		do j=1,ny
        do i=1,nx
          fz(i,j,k)=a1*(f(i,j,k+3)-f(i,j,k-3))  +a2*(f(i,j,k+2)-f(i,j,k-2)) +a3*(f(i,j,k+1)-f(i,j,k-1))
        enddo
        enddo 
		enddo 
 
	else if (Num_Scheme == OCFD_Scheme_CD8) then 
	    do k=1,nz
		do j=1,ny
        do i=1,nx
         fz(i,j,k)=c1*(f(i,j,k+1)-f(i,j,k-1)) +c2*(f(i,j,k+2)-f(i,j,k-2))  &
                  +c3*(f(i,j,k+3)-f(i,j,k-3)) +c4*(f(i,j,k+4)-f(i,j,k-4))
        enddo
        enddo
		enddo 
	else 
	  
	  print*, 'This Numerical Scheme is not supported in viscous terms !'
	  print*, 'Only CD6 or CD8 can be used in viscous terms'
	  stop    
	endif	 
     
!---------Boundary Scheme ------- (low-order scheme)----
!---------- k- boundary  ---------------------------     
	 if(npz .eq. 0 .and. Para%Iperiodic_Z .eq. 0) then 
	   do j=1,ny 
	   do i=1,nx
	    fz(i,j,1)=(-3.d0*f(i,j,1)+4.d0*f(i,j,2)-f(i,j,3))/(2.d0*hz)           ! 2nd one-side scheme
	    fz(i,j,2)=(f(i,j,3)-f(i,j,1))/(2.d0*hz)                             ! 2nd centeral scheme 
		fz(i,j,3)=(8.d0*(f(i,j,4)-f(i,j,2)) - (f(i,j,5)-f(i,j,1)))/(12.d0*hz)   ! 4th central scheme
	   enddo 
	   enddo 
       if(Num_Scheme == OCFD_Scheme_CD8) then 
	   do j=1,ny 
	   do i=1,nx
         fz(i,j,4)=a1*(f(i,j,7)-f(i,j,1)) +a2*(f(i,j,6)-f(i,j,2))  +a3*(f(i,j,5)-f(i,j,3))	  ! 6th centeral scheme
	   enddo
	   enddo 
	   endif 
	 endif 
!--------- k+ boundary ------------------------
	 if(npz .eq. npz0-1 .and. Para%Iperiodic_Z .eq. 0) then 
	   do j=1,ny 
	   do i=1,nx
	    fz(i,j,nz)=(f(i,j,nz-2)-4.d0*f(i,j,nz-1)  +3.d0*f(i,j,nz))/(2.d0*hz)  ! 2nd one-side scheme
	    fz(i,j,nz-1)=(f(i,j,nz)-f(i,j,nz-2))/(2.d0*hz)                             ! 2nd centeral scheme 
		fz(i,j,nz-2)=(8.d0*(f(i,j,nz-1)-f(i,j,nz-3)) - (f(i,j,nz)-f(i,j,nz-4)))/(12.d0*hz)   ! 4th central scheme
	   enddo 
	   enddo 
       if(Num_Scheme == OCFD_Scheme_CD8) then 
	   do j=1,ny 
	   do i=1,nx
         fz(i,j,nz-3)=a1*(f(i,j,nz)-f(i,j,nz-6)) +a2*(f(i,j,nz-1)-f(i,j,nz-5))  +a3*(f(i,j,nz-2)-f(i,j,nz-4))	  ! 6th centeral scheme
	   enddo
	   enddo 
	   endif 
	 endif  	 
 
	 end

!c----------------------------------------------------------	