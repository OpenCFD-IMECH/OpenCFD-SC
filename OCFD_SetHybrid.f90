 ! Scheme%Hybrid_para (:)   1: Shock sensor,  2 seta1,  3 seta2 , 4 Num_patch,  5 ib, 6 ie, 7, jb, 8 je, 9 kb, 10 ke 
  function set_hybrid_scheme(R0)
   use OCFD_precision
   use  Scheme_para 
   implicit none 
   real(kind=OCFD_REAL_KIND):: R0 
   integer:: set_hybrid_scheme
   if(R0 <=Scheme%Hybrid_para(2) ) then 
    set_hybrid_scheme=1            ! linear scheme
   else if (R0 <=Scheme%Hybrid_para(3) ) then 
     set_hybrid_scheme=2           ! shock caputure scheme 1 (WENO7)
   else 
    set_hybrid_scheme=3           ! shock capture scheme 2 (WENO5)
   endif 
   end 
!======================================================   
   subroutine comput_Rhybrid
   use flow_data 
   implicit none
   real(kind=OCFD_REAL_KIND),allocatable,dimension(:,:,:)::p,pk,pi,ps
   integer:: i,j,k,ierr,NS1(4),NS0(4) 
   integer:: m,mpatch, ib, ie, jb, je, kb,ke,Shock_Sensor
   real(kind=OCFD_REAL_KIND):: p00,px,py,pz,dp0,dp_av,dp_av1,epsl,R0
!------------------------------------------------------------------------
   Shock_sensor=nint(Scheme%Hybrid_para(1))
   if(Shock_sensor==1) then 
    allocate(p(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP)) 
   else    
    allocate(p(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP),pk(nx,ny,nz),pi(nx,ny,nz),ps(nx,ny,nz))  
   endif    
    
	p00=1.d0/(Para%gamma*Para%Ma*Para%Ma)      
	epsl=1.d-30
	 
	  do k=1-LAP, nz+LAP
	  do j=1-LAP, ny+LAP
	  do i=1-LAP, nx+LAP 
	   p(i,j,k)=p00*d(i,j,k)*T(i,j,k) 
	  enddo
	  enddo
	  enddo
	  
   if( Shock_sensor==1 ) then 
     do k=1,nz
	 do j=1,ny 
	 do i=1,nx
     Rhybrid(i,j,k)=abs((p(i-1,j,k)-2.d0*p(i,j,k)+p(i+1,j,k))/(p(i-1,j,k)+2.d0*p(i,j,k)+p(i+1,j,k)) ) &
	   + abs((p(i,j-1,k)-2.d0*p(i,j,k)+p(i,j+1,k))/(p(i,j-1,k)+2.d0*p(i,j,k)+p(i,j+1,k)) ) &
 	   + abs((p(i,j,k-1)-2.d0*p(i,j,k)+p(i,j,k+1))/(p(i,j,k-1)+2.d0*p(i,j,k)+p(i,j,k+1)) )  
     enddo 
	 enddo 
	 enddo 

	else	 
	  
	  call OCFD_dx0(p,pk,Scheme%Scheme_Vis)
      call OCFD_dy0(p,pi,Scheme%Scheme_Vis)
      call OCFD_dz0(p,ps,Scheme%Scheme_Vis)

	  dp0=0.d0 
  	  do k=1,nz
      do j=1,ny
      do i=1,nx
 	   px=pk(i,j,k)*Akx(i,j,k)+pi(i,j,k)*Aix(i,j,k)+ps(i,j,k)*Asx(i,j,k)
       py=pk(i,j,k)*Aky(i,j,k)+pi(i,j,k)*Aiy(i,j,k)+ps(i,j,k)*Asy(i,j,k)
       pz=pk(i,j,k)*Akz(i,j,k)+pi(i,j,k)*Aiz(i,j,k)+ps(i,j,k)*Asz(i,j,k)
       Rhybrid(i,j,k)=sqrt(px*px+py*py+pz*pz) 
 	   dp0=dp0+Rhybrid(i,j,k)
      enddo
      enddo
      enddo

	 call MPI_ALLREDUCE(dp0,dp_av,1,OCFD_DATA_TYPE,MPI_SUM,MPI_COMM_WORLD,ierr)


	 dp_av=dp_av/(1.d0*nx_global*ny_global*nz_global)    
     dp_av1=1.d0/(dp_av+epsl)
  
	  do k=1,nz 
	  do j=1,ny
	  do i=1,nx
	   Rhybrid(i,j,k)=Rhybrid(i,j,k)*dp_av1 
	  enddo
	  enddo
	  enddo 
	 endif
	 
!-----------------------------------------------	  
	 mpatch=nint(Scheme%Hybrid_para(4))
	 do m=1,mpatch
	  ib=nint(Scheme%Hybrid_para(5+(m-1)*6))
	  ie=nint(Scheme%Hybrid_para(6+(m-1)*6))
	  jb=nint(Scheme%Hybrid_para(7+(m-1)*6))	  
	  ie=nint(Scheme%Hybrid_para(8+(m-1)*6))
	  kb=nint(Scheme%Hybrid_para(9+(m-1)*6))
	  ke=nint(Scheme%Hybrid_para(10+(m-1)*6))	  
	  ib=max(ib-i_offset(npx)+1,1)               ! transform global index to local index 
	  ie=min(ie-i_offset(npx)+1,nx)
      jb=max(jb-j_offset(npy)+1,1)
      je=min(je-j_offset(npy)+1,ny)
      kb=max(kb-k_offset(npz)+1,1)
      ke=min(ke-k_offset(npz)+1,nz)	  
	  
	  do k=kb,ke 
      do j=jb,je
      do i=ib,ie 
	  Rhybrid(i,j,k)=Rhybrid(i,j,k)+100.d0         ! set to be a large number (using highest rubust scheme)
	  enddo 
	  enddo 
	  enddo 
	enddo   
	  	  
      call exchange_boundary_xyz(Rhybrid)
 !------------------------------------------
	  
  if (mod(Istep,Para%Istep_show).eq.0) then  
     NS1=0
     do k=1,nz 
     do j=1,ny 
     do i=0,nx 
	  R0=max(Rhybrid(i,j,k),Rhybrid(i+1,j,k))
      if( R0 <= Scheme%Hybrid_para(2) ) then 
	    NS1(1)=NS1(1)+1
	  else if(R0 <=Scheme%Hybrid_para(3)) then 
	    NS1(2)=NS1(2)+1
	  else 
	    NS1(3)=NS1(3)+1
      endif 
     enddo 
     enddo 
     enddo 
      NS1(4)=(nx+1)*ny*nz
     call  MPI_Reduce(NS1,NS0,4,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)	 
	 if(my_id .eq. 0) then 
      print*, "---the percent of 3 schemes (linear ,WENO7,WENO5 are: ", NS0(1:3)/(1.d0*NS0(4))
	 endif 
	 
    endif 	 
	
	if(Shock_sensor==1) then
      deallocate(p) 
	else 
	  deallocate(p,pk,pi,ps) 
    endif 
	
!   if(my_id .eq. 0)   open(99,file="Rhybrid.dat",form="unformatted")
!	call write_3d1(99,Rhybrid)
!   if(my_id .eq. 0) then 
!    close(99)
!    stop
!   endif 
!-----------------------------------
    end

