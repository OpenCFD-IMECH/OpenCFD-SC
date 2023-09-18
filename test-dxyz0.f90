! 2-D Compressible N-S Solver
! CopyRight by Li Xinliang  Email: lixl@imech.ac.cn
! Ver 2.0.x  2021-1
!---------------------------------------------------------------------
   subroutine   NS_solver
   use flow_data
   implicit none
   real(kind=OCFD_REAL_KIND),allocatable:: uu(:,:,:), ux1(:,:,:), ux2(:,:,:), ux0(:,:,:)   
   real(kind=OCFD_REAL_KIND):: PI,x,y,z, err1,err2,err3 ,u1 
   real*8:: wt(4)  
   integer:: i,j,k,i1,j1,k1,ierr 
   character (len=30):: filename
!-----------------------------------------------------------------------
       allocate(uu(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP),ux0(nx,ny,nz)) 
        PI=3.14159265358979d0
        hx=2.d0*PI/nx_global 
        hy=2.d0*PI/ny_global 
        hz=2.d0*PI/nz_global
	   do k=1,nz 
       do j=1,ny
       do i=1,nx
        i1=i_offset(npx)+i-1
        j1=j_offset(npy)+j-1
        k1=k_offset(npz)+k-1
		
        x=(i1-1.d0)*hx
        y=(j1-1.d0)*hy
        z=(k1-1.d0)*hz
		uu(i,j,k)=cos(x)*cos(y)*cos(z)
       enddo        
	   enddo
       enddo

       call exchange_boundary_xyz(uu)
	   
!--------x- --------------------------------------	   
	   if(my_id .eq. 0) print*, "--------test dx -----------"
    
  	   wt(1)=MPI_wtime() 
       call  OCFD2d_dx0(uu,ux0,Scheme%Scheme_Vis)
       call MPI_Barrier(MPI_COMM_WORLD,ierr)
       wt(2)=MPI_wtime()


        if(my_id.eq.0) then 
        print*, 'CPU for x-direction is ...'
        print*, wt(2)-wt(1)
        endif
	   


 	   err1=0.d0
       do k=1,nz 
       do j=1,ny
       do i=1,nx
        i1=i_offset(npx)+i-1
        j1=j_offset(npy)+j-1
        k1=k_offset(npz)+k-1
        x=(i1-1.d0)*hx
        y=(j1-1.d0)*hy
        z=(k1-1.d0)*hz 
		
        u1=-sin(x)*cos(y)*cos(z)

        if(abs(u1-ux0(i,j,k)).gt.err1) err1=abs(u1-ux0(i,j,k))
     
       enddo
       enddo
	   enddo 
	   
       print*, 'test x ... my_id,err1=',my_id,err1
        
      if( npy.eq.0  .and. npz .eq. 0 ) then
       write(filename,"('err-x'I3.3'.dat')") npx
        open(33,file=filename)
        do i=1,nx
        i1=i_offset(npx)+i-1
        x=(i1-1.d0)*hx 
        u1=-sin(x)
        write(33,'(I5,3E25.15)') i1,ux0(i,1,1)-u1
      enddo
      close(33)
      endif
      call MPI_Barrier(MPI_COMM_WORLD,ierr)

!--------y- --------------------------------------	   
	   if(my_id .eq. 0) print*, "--------test dx -----------"
    
  	   wt(1)=MPI_wtime() 
       call  OCFD2d_dy0(uu,ux0,Scheme%Scheme_Vis)
       call MPI_Barrier(MPI_COMM_WORLD,ierr)
       wt(2)=MPI_wtime()


        if(my_id.eq.0) then 
        print*, 'CPU for y-direction is ...'
        print*, wt(2)-wt(1)
        endif
	   


 	   err1=0.d0
       do k=1,nz 
       do j=1,ny
       do i=1,nx
        i1=i_offset(npx)+i-1
        j1=j_offset(npy)+j-1
        k1=k_offset(npz)+k-1
        x=(i1-1.d0)*hx
        y=(j1-1.d0)*hy
        z=(k1-1.d0)*hz 
		
        u1=-cos(x)*sin(y)*cos(z)

        if(abs(u1-ux0(i,j,k)).gt.err1) err1=abs(u1-ux0(i,j,k))
     
       enddo
       enddo
	   enddo 
	   
       print*, 'test y ... my_id,err1=',my_id,err1
        
      if( npx.eq.0  .and. npz .eq. 0 ) then
       write(filename,"('err-y'I3.3'.dat')") npy
        open(33,file=filename)
        do j=1,ny
        j1=j_offset(npy)+j-1
        y=(j1-1.d0)*hy 
        u1=-sin(y)
        write(33,'(I5,3E25.15)') j1,ux0(1,j,1)-u1
      enddo
      close(33)
      endif
      call MPI_Barrier(MPI_COMM_WORLD,ierr)

 !--------k- ------------------------------------------------------	   
	   if(my_id .eq. 0) print*, "--------test dz -----------"
    
  	   wt(1)=MPI_wtime() 
       call  OCFD2d_dz0(uu,ux0,Scheme%Scheme_Vis)
       call MPI_Barrier(MPI_COMM_WORLD,ierr)
       wt(2)=MPI_wtime()


        if(my_id.eq.0) then 
        print*, 'CPU for y-direction is ...'
        print*, wt(2)-wt(1)
        endif
	   


 	   err1=0.d0
       do k=1,nz 
       do j=1,ny
       do i=1,nx
        i1=i_offset(npx)+i-1
        j1=j_offset(npy)+j-1
        k1=k_offset(npz)+k-1
        x=(i1-1.d0)*hx
        y=(j1-1.d0)*hy
        z=(k1-1.d0)*hz 
		
        u1=-cos(x)*cos(y)*sin(z)

        if(abs(u1-ux0(i,j,k)).gt.err1) err1=abs(u1-ux0(i,j,k))
     
       enddo
       enddo
	   enddo 
	   
       print*, 'test z ... my_id,err1=',my_id,err1
        
      if( npx.eq.0  .and. npy .eq. 0 ) then
       write(filename,"('err-z'I3.3'.dat')") npz
        open(33,file=filename)
        do k=1,nz
        k1=k_offset(npz)+k-1
        z=(k1-1.d0)*hz 
        u1=-sin(z)
        write(33,'(I5,3E25.15)') k1,ux0(1,1,k)-u1
      enddo
      close(33)
      endif
      call MPI_Barrier(MPI_COMM_WORLD,ierr) 
	  deallocate(uu,ux0)
	  end 

  