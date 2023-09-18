! Inviscous term for Jacobian-transformed N-S equation (Non-Characterical flux) ---------------------

module inviscous_data                 ! data used by inviscous term
 use OCFD_precision
 implicit none 
 real(kind=OCFD_REAL_KIND),allocatable,dimension(:,:,:):: cc          ! speed of sound
 real(kind=OCFD_REAL_KIND),allocatable,dimension(:) ::   &
       c1x,d1x,u1x,v1x,w1x,A1x,A2x,A3x,hhx1,hhx2, & 
	   c1y,d1y,u1y,v1y,w1y,A1y,A2y,A3y,hhy1,hhy2, &
  	   c1z,d1z,u1z,v1z,w1z,A1z,A2z,A3z,hhz1,hhz2
  real(kind=OCFD_REAL_KIND),allocatable,dimension(:,:):: fpx,fmx,fpy,fmy,fpz,fmz, hhx,hhy,hhz
  integer,allocatable,dimension(:)::Scm_Hbx,Scm_Hby,Scm_Hbz        ! hybrid scheme index
end 

subroutine allocate_inviscous_data
 use flow_para
 use inviscous_data 
 implicit none 
   allocate(cc(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP))
   allocate( fpx(1-LAP:nx+LAP,5),fmx(1-LAP:nx+LAP,5))
   allocate( c1x(1-LAP:nx+LAP),d1x(1-LAP:nx+LAP), u1x(1-LAP:nx+LAP),  &
             v1x(1-LAP:nx+LAP),w1x(1-LAP:nx+LAP), A1x(1-LAP:nx+LAP),  &
			 A2x(1-LAP:nx+LAP), A3x(1-LAP:nx+LAP), &
			 hhx1(0:nx),hhx2(0:nx), hhx(0:nx,5))
  
   allocate( fpy(1-LAP:ny+LAP,5),fmy(1-LAP:ny+LAP,5))
   allocate( c1y(1-LAP:ny+LAP),d1y(1-LAP:ny+LAP), u1y(1-LAP:ny+LAP),  &
             v1y(1-LAP:ny+LAP),w1y(1-LAP:ny+LAP), A1y(1-LAP:ny+LAP),  &
			 A2y(1-LAP:ny+LAP), A3y(1-LAP:ny+LAP), &
			 hhy1(0:ny),hhy2(0:ny), hhy(0:ny,5))

   allocate( fpz(1-LAP:nz+LAP,5),fmz(1-LAP:nz+LAP,5))
   allocate( c1z(1-LAP:nz+LAP),d1z(1-LAP:nz+LAP), u1z(1-LAP:nz+LAP),  &
             v1z(1-LAP:nz+LAP),w1z(1-LAP:nz+LAP), A1z(1-LAP:nz+LAP),  &
			 A2z(1-LAP:nz+LAP), A3z(1-LAP:nz+LAP), &
			 hhz1(0:nz),hhz2(0:nz),hhz(0:nz,5))			 
   
   allocate(Scm_Hbx(0:nx),Scm_Hby(0:ny),Scm_Hbz(0:nz))
   Scm_Hbx=0 ; Scm_Hby=0;  Scm_Hbz=0    ! 2021-6-7     (initialized as zero, in case of non-hybrid scheme, character-type flux used)
 end 
 
 subroutine deallocate_inviscous_data
 use flow_para
 use inviscous_data 
 implicit none 
 deallocate( c1x,d1x,u1x,v1x,w1x,A1x,A2x,A3x,hhx1,hhx2,hhx, & 
	   c1y,d1y,u1y,v1y,w1y,A1y,A2y,A3y,hhy1,hhy2,hhy, &
  	   c1z,d1z,u1z,v1z,w1z,A1z,A2z,A3z,hhz1,hhz2,hhz, &
	   fpx,fmx,fpy,fmy,fpz,fmz,Scm_Hbx,Scm_Hby,Scm_Hbz)
 end 
 
 
 
!-----------Inviscous term  (non-character type)----------------------------
  subroutine du_inviscous 
   use flow_data 
   use inviscous_data
   implicit none
   integer:: i,j,k,m,set_hybrid_scheme
   real(kind=OCFD_REAL_KIND):: hx_1, hy_1, hz_1,Rhb0
!c-------------------------------------------------
    hx_1=1.d0/hx
	hy_1=1.d0/hy
    hz_1=1.d0/hz      
	  
	do k=1-LAP,nz+LAP 
	do j=1-LAP,ny+LAP
	do i=1-LAP,nx+LAP
     cc(i,j,k)=sqrt(T(i,j,k))/Para%Ma
	enddo
	enddo
    enddo 
	  
!c---------x direction-------------------------
   do k=1,nz
   do j=1,ny   
    
	if(Scheme%Scheme_Invis == 	OCFD_Scheme_Hybrid) then 
     do i=0,nx
	  Rhb0=max(Rhybrid(i,j,k),Rhybrid(i+1,j,k))
	  Scm_Hbx(i)=set_hybrid_scheme(Rhb0)
	 enddo
	endif 
   
	  do i=1-LAP,nx+LAP
	   c1x(i)=cc(i,j,k)
	   d1x(i)=d(i,j,k)
	   u1x(i)=u(i,j,k)
	   v1x(i)=v(i,j,k)
	   w1x(i)=w(i,j,k)
 	   A1x(i)=Akx1(i,j,k)
	   A2x(i)=Aky1(i,j,k)
	   A3x(i)=Akz1(i,j,k)
	  enddo 

 ! Steger-Warming Flux Vector Splitting (SW-FVS)
      if(Para%Flux_Splitting .eq. OCFD_Split_SW) then 
	   call split_Stager_Warming(nx,LAP,c1x,d1x,u1x,v1x,w1x,A1x,A2x,A3x,fpx,fmx) 
      else 
       call split_Local_LaxFriedrichs(nx,LAP,c1x,d1x,u1x,v1x,w1x,A1x,A2x,A3x,fpx,fmx)
      endif 
	  
	  do m=1,5
 		! upwind finite-difference: dx1 for Positive flux ;  dx2 for Negative flux  
		call OCFD2d_flux1(fpx(1-LAP,m),hhx1,nx,LAP,Scheme%Bound_index(1,1),Scm_Hbx, Scheme%Scheme_boundary(1:2))        
        call OCFD2d_flux2(fmx(1-LAP,m),hhx2,nx,LAP,Scheme%Bound_index(1,1),Scm_Hbx, Scheme%Scheme_boundary(1:2))    

	   do i=1,nx
	    du(i,j,k,m)= (hhx1(i)-hhx1(i-1) + hhx2(i)-hhx2(i-1))*hx_1
	   enddo
	  enddo
   enddo 
   enddo 
!c-------y direction ------------------------
   do k=1,nz    
   do i=1,nx 
     
	 if(Scheme%Scheme_Invis == 	OCFD_Scheme_Hybrid) then 
     do j=0,ny
	  Rhb0=max(Rhybrid(i,j,k),Rhybrid(i,j+1,k))
	  Scm_Hby(j)=set_hybrid_scheme(Rhb0)
	 enddo
	endif   
   
	do j=1-LAP,ny+LAP
      c1y(j)=cc(i,j,k)
	  d1y(j)=d(i,j,k)
	  u1y(j)=u(i,j,k)
	  v1y(j)=v(i,j,k)
	  w1y(j)=w(i,j,k)
	  A1y(j)=Aix1(i,j,k)
	  A2y(j)=Aiy1(i,j,k)
	  A3y(j)=Aiz1(i,j,k)
	enddo 
	
	if(Para%Flux_Splitting .eq. OCFD_Split_SW) then 
     call split_Stager_Warming(ny,LAP,c1y,d1y,u1y,v1y,w1y,A1y,A2y,A3y,fpy,fmy) 
    else 
	 call split_Local_LaxFriedrichs(ny,LAP,c1y,d1y,u1y,v1y,w1y,A1y,A2y,A3y,fpy,fmy) 
	endif 
	
    do m=1,5
 	 ! dy1 for positive flux;  dy2 for negative flux   
	 call OCFD2d_flux1(fpy(1-LAP,m),hhy1,ny,LAP,   Scheme%Bound_index(1,2),Scm_Hby, Scheme%Scheme_boundary(3:4) )      
     call OCFD2d_flux2(fmy(1-LAP,m),hhy2,ny,LAP,   Scheme%Bound_index(1,2),Scm_Hby, Scheme%Scheme_boundary(3:4) )
	 do j=1,ny
	  du(i,j,k,m)= du(i,j,k,m)+ (hhy1(j)-hhy1(j-1)+hhy2(j)-hhy2(j-1))*hy_1                 ! df/dt=-du
     enddo
    enddo 
   enddo 
   enddo 

!c-------z direction ------------------------
   do j=1,ny    
   do i=1,nx 
    if(Scheme%Scheme_Invis == 	OCFD_Scheme_Hybrid) then 
     do k=0,nz
	  Rhb0=max(Rhybrid(i,j,k),Rhybrid(i,j,k+1))
	  Scm_Hbz(k)=set_hybrid_scheme(Rhb0)	 
	 enddo
	endif 
      
	do k=1-LAP,nz+LAP
      c1z(k)=cc(i,j,k)
	  d1z(k)=d(i,j,k)
	  u1z(k)=u(i,j,k)
	  v1z(k)=v(i,j,k)
	  w1z(k)=w(i,j,k)
	  A1z(k)=Asx1(i,j,k)
	  A2z(k)=Asy1(i,j,k)
	  A3z(k)=Asz1(i,j,k)
	enddo 
	
	if(Para%Flux_Splitting .eq. OCFD_Split_SW) then	
     call split_Stager_Warming(nz,LAP,c1z,d1z,u1z,v1z,w1z,A1z,A2z,A3z,fpz,fmz) 
    else 
	 call split_Local_LaxFriedrichs(nz,LAP,c1z,d1z,u1z,v1z,w1z,A1z,A2z,A3z,fpz,fmz)
    endif 
	
	do m=1,5
 	 ! dz1 for positive flux;  dz2 for negative flux  
	 call OCFD2d_flux1(fpz(1-LAP,m),hhz1,nz,LAP,   Scheme%Bound_index(1,3),Scm_Hbz, Scheme%Scheme_boundary(5:6) )    
     call OCFD2d_flux2(fmz(1-LAP,m),hhz2,nz,LAP,   Scheme%Bound_index(1,3),Scm_Hbz, Scheme%Scheme_boundary(5:6) )
	 do k=1,nz
	  du(i,j,k,m)=-(du(i,j,k,m)+ (hhz1(k)-hhz1(k-1)+hhz2(k)-hhz2(k-1))*hz_1 )*Ajac(i,j,k)                ! df/dt=-du
     enddo
    enddo 
   enddo 
   enddo 


    if(Para%IF_Mass_Force .eq. 1) then 
	  do k=1,nz
	  do j=1,ny
	  do i=1,nx
	     du(i,j,k,2)=du(i,j,k,2)+d(i,j,k)*Para%Mass_Force(1)
		 du(i,j,k,3)=du(i,j,k,3)+d(i,j,k)*Para%Mass_Force(2)
		 du(i,j,k,4)=du(i,j,k,4)+d(i,j,k)*Para%Mass_Force(3)
		 du(i,j,k,5)=du(i,j,k,5)+d(i,j,k)*(u(i,j,k)*Para%Mass_Force(1) &
		          +v(i,j,k)*Para%Mass_Force(2) +w(i,j,k)*Para%Mass_Force(3)) 
	  enddo
	  enddo 
	  enddo 
   endif 
	 

  end
 
!c-----------------------------------
!Scheme_boundary(:)==0  ! WENO5-type boundary scheme (default);  ==-1 Ghost Cell

! Finite difference Numerical flux for inviscous terms  (Upwind Schemes)
! In this version, Only WENO5/WENO7/OMP6 schemes are supported! 
    subroutine OCFD2d_flux1(f1d,hh,nx1,LAP1,Bound_index,Scm_Hy,Scheme_boundary)
    use OCFD_constants
	use Scheme_Para
	implicit none
    integer nx1,LAP1,Bound_index(2),Ka,Kb,i,Scm_Hy(0:nx1),Scheme_boundary(2),ib,ie
	real(kind=OCFD_REAL_KIND)::  f1d(1-LAP1:nx1+LAP1),hh(0:nx1)  
      Ka=Scheme%Ka1 ;  Kb=Scheme%Kb1
	  
	  if(Bound_index(1) .eq. 0  .or. Scheme_boundary(1) .eq. -1) then      ! inner point [ib,ie] 
	   ib=0
	  else  
	   ib=-Ka+1
	  endif 

	  if(Bound_index(2) .eq. 0 .or. Scheme_boundary(2) .eq. -1  ) then 
       ie=nx1     
	  else  
	   ie=nx1-Kb	  
	  endif 
	  
	   select case(Scheme%Scheme_Invis)
	    case (OCFD_Scheme_OMP6)
	       call OCFD_OMP6P(f1d,hh,nx1,LAP1,ib,ie)	
        case(OCFD_Scheme_WENO7)
	       call OCFD_WENO7P(f1d,hh,nx1,LAP1,ib,ie)
        case(OCFD_Scheme_WENO5)
	       call OCFD_WENO5P(f1d,hh,nx1,LAP1,ib,ie)
        case(OCFD_Scheme_UD7L)
	       call OCFD_UD7L_P(f1d,hh,nx1,LAP1,ib,ie)
        case(OCFD_Scheme_Hybrid)
	       call OCFD_HybridP(f1d,hh,nx1,LAP1,Scm_Hy,ib,ie)		   
        case(OCFD_Scheme_USER)
	       call OCFD_Scheme_USER_P(f1d,hh,nx1,LAP1,ib,ie)		   
		   
	   end select  
	   
!------------boundary scheme --------------------------------	   
	 if(Bound_index(1) .eq. 1 )then         ! i- (j-, k-) boundary 
       if(Scheme_boundary(1) .eq. 0) then  ! Default (WENO5-type) 
		 do i=0, -Ka
         if(i==0 .or. i==1) then 
 		     hh(i)=(2.d0*f1d(1)+5.d0*f1d(2)-f1d(3))/6.d0 
		 else if(i==2) then 
		    call hh_weno5P_boundary(Ka,Kb,f1d(i+Ka:i+Kb),hh(i), DEL_LIFT )
	     else 
		    call hh_weno5P(Ka,Kb,f1d(i+Ka:i+Kb),hh(i) )
		 endif 
		 enddo 
! -1  : WENO5 + Ghost Cell;  
       else if(Scheme_boundary(1) .eq. 1  ) then  ! WENO5 (with Ghost Cell) 
 		do i=0, -Ka
		 call hh_weno5P(Ka,Kb,f1d(i+Ka:i+Kb),hh(i) )
	    enddo 	 
	   endif 
	  endif 	

	   if(Bound_index(2) .eq. 1 )then  ! i+ (j+, k+) boundary 
       if(Scheme_boundary (2) .eq. 0) then  ! Default (WENO5 Del-substencil type) 
	    do i=nx1-Kb+1, nx1
         if(i==nx1) then 
 		   hh(i)=(2.d0*f1d(nx1-2)-7.d0*f1d(nx1-1)+11.d0*f1d(nx1))/6.d0   
		 else if(i==nx1-1) then 
		   call hh_weno5P_boundary(Ka,Kb,f1d(i+Ka:i+Kb),hh(i), DEL_RIGHT)
	     else 
		    call hh_weno5P(Ka,Kb,f1d(i+Ka:i+Kb),hh(i) )
		 endif 
        enddo 
	   else if(Scheme_boundary(2) .eq. 1 ) then  ! WENO5 with Ghost Cell 	
  	    do i=nx1-Kb+1, nx1   
         call hh_weno5P(Ka,Kb,f1d(i+Ka:i+Kb),hh(i) )
		enddo 		
	  endif 
	  endif 	
      end 

!---------------------------------------------------------	
! Scheme for flux- 
    subroutine OCFD2d_flux2(f1d,hh,nx1,LAP1, Bound_index,Scm_Hy,Scheme_boundary)
    use OCFD_constants
	use Scheme_Para
	implicit none
    integer nx1,LAP1,Bound_index(2),Ka,Kb,i,Scm_Hy(0:nx1),Scheme_boundary(2),ib,ie
	real(kind=OCFD_REAL_KIND)::  f1d(1-LAP1:nx1+LAP1),hh(0:nx1)
      Ka=Scheme%Ka2 ;  Kb=Scheme%Kb2	
	  
	  if(Bound_index(1) .eq. 0  .or. Scheme_boundary(1) .eq. -1) then      ! inner point [ib,ie] 
	   ib=0
	  else  
	   ib=-Ka+1
	  endif 

	  if(Bound_index(2) .eq. 0 .or. Scheme_boundary(2) .eq. -1  ) then 
       ie=nx1     
	  else  
	   ie=nx1-Kb	  
	  endif 
	  

	  
	  select case(Scheme%Scheme_Invis)
	    case (OCFD_Scheme_OMP6)
	       call OCFD_OMP6M(f1d,hh,nx1,LAP1,ib,ie)	
        case(OCFD_Scheme_WENO7)
	       call OCFD_WENO7M(f1d,hh,nx1,LAP1,ib,ie)
        case(OCFD_Scheme_WENO5)
	       call OCFD_WENO5M(f1d,hh,nx1,LAP1,ib,ie)
        case(OCFD_Scheme_UD7L)
	       call OCFD_UD7L_M(f1d,hh,nx1,LAP1,ib,ie)
        case(OCFD_Scheme_Hybrid)
	       call OCFD_HybridM(f1d,hh,nx1,LAP1,Scm_Hy,ib,ie)		   
        case(OCFD_Scheme_USER)
	       call OCFD_Scheme_USER_M(f1d,hh,nx1,LAP1,ib,ie)		   
	   end select  

!------------boundary scheme --------------	   
	   if(Bound_index(1) .eq. 1 )then       ! i- boundary
       if(Scheme_boundary(1) .eq. 0) then  ! Default (WENO5-type)	    
		do i=0, -Ka
         if(i==0 ) then 
 	       hh(i)=(2.d0*f1d(3)-7.d0*f1d(2)+11.d0*f1d(1))/6.d0  
		 else if(i==1) then 
		    call hh_weno5M_boundary(Ka,Kb,f1d(i+Ka:i+Kb),hh(i), DEL_LIFT )
		 else 
		    call hh_weno5M(Ka,Kb,f1d(i+Ka:i+Kb),hh(i) )
		 endif
        enddo
		
		else if(Scheme_boundary(1) .eq. 1 ) then  ! WENO5 + Ghost Cell	    
   		do i=0, -Ka
         call hh_weno5M(Ka,Kb,f1d(i+Ka:i+Kb),hh(i) )
        enddo 
		
	   endif 
  	   endif 	

	    if(Bound_index(2) .eq. 1 )then  ! i+ boundary
        if(Scheme_boundary(2) .eq. 0) then  ! Default (WENO5-type)			
	    do i=nx1-Kb+1, nx1
         if(i==nx1 .or. i==nx1-1) then 
 		   hh(i)=(2.d0*f1d(nx1)+5.d0*f1d(nx1-1)-f1d(nx1-2))/6.d0 
		 else if(i==nx1-2) then 
		   call hh_weno5M_boundary(Ka,Kb,f1d(i+Ka:i+Kb),hh(i), DEL_RIGHT)
	     else 
		  call hh_weno5M(Ka,Kb,f1d(i+Ka:i+Kb),hh(i) )
		 endif 
		enddo 
	    else if(Scheme_boundary(2) .eq. 1 ) then  !  (WENO5 +Ghost Cell )		
 		do i=nx1-Kb+1, nx1	
		 call hh_weno5M(Ka,Kb,f1d(i+Ka:i+Kb),hh(i) )
		enddo 
		endif 
        endif 	

     end	 
	 
!c----------------------------------------------------------
  
